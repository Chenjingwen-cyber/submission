# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import csv
import json
import random
import time
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from candidate_nodes import (
    CandidateLibrary,
    Matrix,
    Value,
    build_seed_candidate_library,
    add_td_otd_completion,
    load_one_matrix,
)
from milp_model import MilpSolution, SelectedGate, solve_layered_milp


# ---------------------------------------------------------------------------
# GF(2) evaluation and selected-circuit utilities
# ---------------------------------------------------------------------------

def matmul_gf2(matrix: Matrix, x_bits: Sequence[int]) -> List[int]:
    out: List[int] = []
    for row in matrix:
        s = 0
        for bit, xb in zip(row, x_bits):
            if bit:
                s ^= xb
        out.append(s)
    return out


def bitset_eval(value: Value, x_bits: Sequence[int]) -> int:
    s = 0
    i = 0
    v = int(value)
    while v:
        if v & 1:
            s ^= x_bits[i]
        v >>= 1
        i += 1
    return s


def build_step_order(lib: CandidateLibrary, sol: MilpSolution) -> List[SelectedGate]:
    seen = set(lib.inputs)
    out: List[SelectedGate] = []
    for gate in sorted(sol.selected_gates, key=lambda g: (g.layer, g.value, g.arity, g.impl)):
        if gate.value in seen:
            continue
        out.append(gate)
        seen.add(gate.value)
    return out


def get_gate_sources(lib: CandidateLibrary, gate: SelectedGate) -> List[str]:
    sources = getattr(lib, "impl_sources", {}).get((gate.value, tuple(sorted(gate.impl))))
    if not sources:
        return ["untracked"]
    return sorted(str(s) for s in sources)


def source_family(sources: Sequence[str]) -> str:
    """Coarse provenance labels kept inside the circuit file only."""
    if any(s.startswith("td_otd_terminal") for s in sources):
        return "Strategy2_terminal"
    if any(s.startswith("td_otd_bridge") for s in sources):
        return "Strategy2_bridge"
    if any(s.startswith("srtc") for s in sources):
        return "srtc_core"
    if any(s.startswith("input_pair_layer1") or s.startswith("input_triple_layer1") for s in sources):
        return "Strategy1_input"
    if any("shortest" in s or "triple_hint" in s or s.startswith("trace_pair") for s in sources):
        return "Strategy1_hint"
    return "other"


# ---------------------------------------------------------------------------
# Exports: one circuit file and one result file per Strategy
# ---------------------------------------------------------------------------

def export_circuit(lib: CandidateLibrary, sol: MilpSolution, path: Path) -> Dict[Value, str]:
    """Export only the selected circuit equations and their depths.

    Format example:
        Step 1: b32 = b1 + b9    depth=1
        Step 61: b92 = b33 + b27 + b67    depth=3    // y19
    """
    steps = build_step_order(lib, sol)
    value_to_wire: Dict[Value, str] = {v: f"b{i}" for i, v in enumerate(lib.inputs)}
    next_idx = len(lib.inputs)
    lines: List[str] = []

    for step_no, gate in enumerate(steps, start=1):
        for parent in gate.impl:
            if parent not in value_to_wire:
                value_to_wire[parent] = f"UNKNOWN({parent})"
        lhs = f"b{next_idx}"
        next_idx += 1
        value_to_wire[gate.value] = lhs
        rhs = " + ".join(value_to_wire[p] for p in gate.impl)
        target_suffix = ""
        if gate.value in lib.target_outputs:
            ys = ",".join(f"y{i}" for i in lib.target_outputs[gate.value])
            target_suffix = f"    // {ys}"
        lines.append(f"Step {step_no}: {lhs} = {rhs}    depth={gate.layer}{target_suffix}")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    return value_to_wire


def verify_symbolic(lib: CandidateLibrary, sol: MilpSolution, depth: int) -> Tuple[bool, List[str], Dict[Value, int]]:
    errors: List[str] = []
    available: Dict[Value, int] = {v: 0 for v in lib.inputs}
    for gate in build_step_order(lib, sol):
        if any(p not in available for p in gate.impl):
            errors.append(f"Parent not available for value {gate.value}: {gate.impl}")
            continue
        value = 0
        parent_depth = 0
        for p in gate.impl:
            value ^= p
            parent_depth = max(parent_depth, available[p])
        expected_depth = parent_depth + 1
        if value != gate.value:
            errors.append(f"Gate mismatch: value {gate.value} != xor{gate.impl}={value}")
        # The MILP layer may be a valid upper layer.  Only actual depth larger
        # than the MILP layer is invalid.
        if expected_depth > gate.layer:
            errors.append(f"Gate depth mismatch: value {gate.value}, actual {expected_depth} > MILP layer {gate.layer}")
        available[gate.value] = expected_depth

    for row_idx, target in enumerate(lib.row_targets):
        if target not in available:
            errors.append(f"Output y{row_idx} target value {target} is not produced by selected circuit")
        elif available[target] > depth:
            errors.append(f"Output y{row_idx} has depth {available[target]} > {depth}")
    return len(errors) == 0, errors, available


def verify_random(matrix: Matrix, lib: CandidateLibrary, sol: MilpSolution, trials: int, seed: int = 12345) -> Tuple[bool, List[str]]:
    errors: List[str] = []
    rng = random.Random(seed)
    steps = build_step_order(lib, sol)
    n = len(matrix[0])
    for trial in range(trials):
        x = [rng.randint(0, 1) for _ in range(n)]
        wire_values: Dict[Value, int] = {v: bitset_eval(v, x) for v in lib.inputs}
        for gate in steps:
            s = 0
            for p in gate.impl:
                if p not in wire_values:
                    errors.append(f"Trial {trial}: missing parent value {p}")
                    return False, errors
                s ^= wire_values[p]
            wire_values[gate.value] = s
        expected = matmul_gf2(matrix, x)
        actual: List[int] = []
        for target in lib.row_targets:
            if target not in wire_values:
                errors.append(f"Trial {trial}: output target {target} missing")
                return False, errors
            actual.append(wire_values[target])
        if actual != expected:
            errors.append(f"Trial {trial}: output mismatch. expected={expected}, actual={actual}")
            return False, errors
    return True, errors


def _meta_value(lib: CandidateLibrary, family: str, key: str, default: object = "") -> object:
    meta = lib.meta.get(family, {})
    if isinstance(meta, dict):
        return meta.get(key, default)
    return default


def write_result_file(
    path: Path,
    *,
    Strategy_tag: str,
    matrix_name: str,
    index: int,
    depth: int,
    num_runs: int,
    seed: int,
    lambda2: float,
    lambda3: float,
    terminal_width: int,
    terminal_parent_pool_limit: int,
    lib: CandidateLibrary,
    sol: MilpSolution,
    verified: bool,
    symbolic_ok: bool,
    symbolic_errors: List[str],
    random_ok: bool,
    random_errors: List[str],
    actual_depths: Dict[Value, int],
    circuit_path: Path,
) -> None:
    actual_max_depth = max((actual_depths.get(t, 0) for t in lib.row_targets), default=0)
    lines = [
        "[experiment]",
        f"Strategy={Strategy_tag}",
        f"matrix={matrix_name}",
        f"index={index}",
        f"depth={depth}",
        f"num_runs={num_runs}",
        f"seed={seed}",
        f"lambda2={lambda2}",
        f"lambda3={lambda3}",
        f"terminal_width={terminal_width}",
        f"terminal_parent_pool_limit={terminal_parent_pool_limit}",
        f"circuit_file={circuit_path}",
        "",
        "[solve]",
        f"status={sol.status_name}",
        f"objective={sol.objective}",
        f"xor2={sol.xor2_count}",
        f"xor3={sol.xor3_count}",
        f"runtime_sec={sol.runtime}",
        f"mip_gap={sol.mip_gap}",
        f"best_bound={sol.best_bound}",
        "",
        "[verification]",
        f"verified={verified}",
        f"symbolic_ok={symbolic_ok}",
        f"random_ok={random_ok}",
        f"actual_max_output_depth={actual_max_depth}",
        "",
        "[candidate]",
        f"candidate_values={lib.meta.get('num_values', '')}",
        f"pair_impls={lib.meta.get('num_pair_impls', '')}",
        f"triple_impls={lib.meta.get('num_triple_impls', '')}",
        f"feasible_targets={lib.meta.get('feasible_targets', '')}/{lib.meta.get('total_targets', '')}",
        f"infeasible_targets={lib.meta.get('infeasible_targets_under_depth', [])}",
    ]

    if Strategy_tag == "Strategy2":
        lines += [
            f"td_otd_terminal_width={_meta_value(lib, 'td_otd', 'terminal_width')}",
            f"td_otd_parent_pool={_meta_value(lib, 'td_otd', 'parent_pool_size')}",
            f"td_otd_bridge_candidates={_meta_value(lib, 'td_otd', 'bridge_candidate_values')}",
            f"td_otd_selected_bridge_values={_meta_value(lib, 'td_otd', 'selected_bridge_values')}",
            f"td_otd_added_generated_bridges={_meta_value(lib, 'td_otd', 'added_generated_bridges')}",
            f"td_otd_added_terminal_impls={_meta_value(lib, 'td_otd', 'added_terminal_impls')}",
            f"td_otd_added_direct_terminals={_meta_value(lib, 'td_otd', 'added_direct_terminals')}",
            f"td_otd_added_residual_terminals={_meta_value(lib, 'td_otd', 'added_residual_terminals')}",
            f"terminal_parent_scope={_meta_value(lib, 'td_otd', 'terminal_parent_completion_scope')}",
            f"terminal_parent_pool_limit={_meta_value(lib, 'td_otd', 'terminal_parent_pool_limit')}",
            f"terminal_parent_pool_size={_meta_value(lib, 'td_otd', 'terminal_parent_pool_size')}",
            f"terminal_parent_seed_existing={_meta_value(lib, 'td_otd', 'terminal_parent_seed_existing')}",
            f"terminal_parent_generated={_meta_value(lib, 'td_otd', 'terminal_parent_generated')}",
            f"terminal_parent_ranked_candidates={_meta_value(lib, 'td_otd', 'terminal_parent_ranked_candidates')}",
            f"topw_sort_rule={_meta_value(lib, 'td_otd', 'topw_sort_rule')}",
        ]

    if symbolic_errors:
        lines.append("\n[symbolic_errors]")
        lines.extend(f"  {e}" for e in symbolic_errors)
    if random_errors:
        lines.append("\n[random_errors]")
        lines.extend(f"  {e}" for e in random_errors)

    result_json = {
        "Strategy": Strategy_tag,
        "matrix": matrix_name,
        "index": index,
        "depth": depth,
        "num_runs": num_runs,
        "seed": seed,
        "lambda2": lambda2,
        "lambda3": lambda3,
        "terminal_width": terminal_width,
        "terminal_parent_pool_limit": terminal_parent_pool_limit,
        "status": sol.status_name,
        "objective": sol.objective,
        "xor2": sol.xor2_count,
        "xor3": sol.xor3_count,
        "runtime_sec": sol.runtime,
        "mip_gap": sol.mip_gap,
        "best_bound": sol.best_bound,
        "verified": verified,
        "symbolic_ok": symbolic_ok,
        "random_ok": random_ok,
        "actual_max_output_depth": actual_max_depth,
        "candidate_values": lib.meta.get("num_values", ""),
        "pair_impls": lib.meta.get("num_pair_impls", ""),
        "triple_impls": lib.meta.get("num_triple_impls", ""),
        "circuit_file": str(circuit_path),
    }
    lines += ["", "[json]", json.dumps(result_json, ensure_ascii=False, sort_keys=True)]

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def solve_export_one(
    Strategy_tag: str,
    matrix: Matrix,
    matrix_name: str,
    index: int,
    lib: CandidateLibrary,
    depth: int,
    num_runs: int,
    seed: int,
    lambda2: float,
    lambda3: float,
    terminal_width: int,
    terminal_parent_pool_limit: int,
    time_limit: float,
    output_flag: int,
    random_tests: int,
    prefix: Path,
) -> Tuple[MilpSolution, bool, Dict[str, Path]]:
    print(f"\n[solve] {Strategy_tag}: solving layered MILP")
    sol = solve_layered_milp(
        lib=lib,
        depth=depth,
        lambda2=lambda2,
        lambda3=lambda3,
        time_limit=time_limit,
        output_flag=output_flag,
    )
    print(f"[solve] {Strategy_tag}: status={sol.status_name}, objective={sol.objective}, xor2={sol.xor2_count}, xor3={sol.xor3_count}, gap={sol.mip_gap}")

    paths = {
        "circuit": prefix.with_name(prefix.name + "_circuit.txt"),
        "result": prefix.with_name(prefix.name + "_result.txt"),
    }

    verified = False
    symbolic_ok = False
    random_ok = False
    symbolic_errors: List[str] = []
    random_errors: List[str] = []
    actual_depths: Dict[Value, int] = {}

    if sol.selected_gates:
        export_circuit(lib, sol, paths["circuit"])
        symbolic_ok, symbolic_errors, actual_depths = verify_symbolic(lib, sol, depth)
        random_ok, random_errors = verify_random(matrix, lib, sol, random_tests)
        verified = symbolic_ok and random_ok
        print(f"[verify] {Strategy_tag}: verified={verified}")
        print(f"[write]  {Strategy_tag}: {paths['circuit']}")
    else:
        symbolic_errors = ["No selected gates"]
        random_errors = ["No selected gates"]

    write_result_file(
        paths["result"],
        Strategy_tag=Strategy_tag,
        matrix_name=matrix_name,
        index=index,
        depth=depth,
        num_runs=num_runs,
        seed=seed,
        lambda2=lambda2,
        lambda3=lambda3,
        terminal_width=terminal_width,
        terminal_parent_pool_limit=terminal_parent_pool_limit,
        lib=lib,
        sol=sol,
        verified=verified,
        symbolic_ok=symbolic_ok,
        symbolic_errors=symbolic_errors,
        random_ok=random_ok,
        random_errors=random_errors,
        actual_depths=actual_depths,
        circuit_path=paths["circuit"],
    )
    return sol, verified, paths


# ---------------------------------------------------------------------------
# Compact summary CSV across all experiment configurations
# ---------------------------------------------------------------------------

def append_summary(summary_path: Path, row: Dict[str, object]) -> None:
    fieldnames = [
        "matrix", "index", "depth", "num_runs", "seed", "Strategy", "status", "verified",
        "objective", "xor2", "xor3", "runtime_sec", "gap", "best_bound",
        "candidate_values", "pair_impls", "triple_impls",
        "terminal_width", "terminal_parent_pool_limit", "topw_sort_rule",
        "circuit_file", "result_file",
    ]
    write_header = not summary_path.exists()
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerow({k: row.get(k, "") for k in fieldnames})


def summary_row(
    *,
    matrix_name: str,
    index: int,
    depth: int,
    num_runs: int,
    seed: int,
    Strategy_tag: str,
    lib: CandidateLibrary,
    sol: MilpSolution,
    verified: bool,
    terminal_width: int,
    terminal_parent_pool_limit: int,
    paths: Dict[str, Path],
) -> Dict[str, object]:
    return {
        "matrix": matrix_name,
        "index": index,
        "depth": depth,
        "num_runs": num_runs,
        "seed": seed,
        "Strategy": Strategy_tag,
        "status": sol.status_name,
        "verified": verified,
        "objective": sol.objective if sol.objective is not None else "",
        "xor2": sol.xor2_count,
        "xor3": sol.xor3_count,
        "runtime_sec": sol.runtime,
        "gap": sol.mip_gap,
        "best_bound": sol.best_bound,
        "candidate_values": lib.meta.get("num_values", ""),
        "pair_impls": lib.meta.get("num_pair_impls", ""),
        "triple_impls": lib.meta.get("num_triple_impls", ""),
        "terminal_width": terminal_width,
        "terminal_parent_pool_limit": terminal_parent_pool_limit,
        "topw_sort_rule": _meta_value(lib, "td_otd", "topw_sort_rule", ""),
        "circuit_file": paths.get("circuit", ""),
        "result_file": paths.get("result", ""),
    }


# ---------------------------------------------------------------------------
# Main two-Strategy runner
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Final two-Strategy GDA-Seed / TD-OTD experiment runner.")
    parser.add_argument("--matrix", required=True, help="Path to BP-format matrix txt file")
    parser.add_argument("--index", type=int, default=0)
    parser.add_argument("--depth", type=int, default=3)
    parser.add_argument("--num-runs", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260424, help="Base seed. First run uses exactly this seed.")
    parser.add_argument("--lambda2", type=float, default=3.0)
    parser.add_argument("--lambda3", type=float, default=4.67)
    parser.add_argument("--time-limit-seed", type=float, default=1800)
    parser.add_argument("--time-limit-td-otd", type=float, default=1800)
    parser.add_argument("--terminal-width", type=int, default=512, help="Strategy-2 per-target terminal completion width W. <=0 keeps all plans.")
    parser.add_argument("--terminal-parent-pool-limit", type=int, default=3000, help="D=3 controlled layer-2 terminal-parent pool limit. Passed through to candidate_nodes.py.")
    parser.add_argument("--max-minlin-solutions", type=int, default=None)
    parser.add_argument("--output-dir", default="results")
    parser.add_argument("--gurobi-output", type=int, default=1)
    parser.add_argument("--random-tests", type=int, default=200)
    args = parser.parse_args()

    matrix_path = Path(args.matrix)
    matrix_name = matrix_path.stem
    matrix = load_one_matrix(matrix_path, args.index)
    root_out_dir = Path(args.output_dir)
    run_dir = root_out_dir / matrix_name / f"D{args.depth}" / f"runs{args.num_runs}_seed{args.seed}"
    run_dir.mkdir(parents=True, exist_ok=True)

    width_tag = "all" if args.terminal_width <= 0 else str(args.terminal_width)
    basename_common = f"{matrix_name}_idx{args.index}_D{args.depth}_runs{args.num_runs}_seed{args.seed}"

    print("[1/3] Building Strategy1 GDA-Seed candidate DAG once")
    t_build_seed = time.time()
    Strategy1_lib = build_seed_candidate_library(
        matrix=matrix,
        depth=args.depth,
        num_runs=args.num_runs,
        base_seed=args.seed,
        max_minlin_solutions=args.max_minlin_solutions,
        verbose=True,
    )
    print(f"[build] Strategy1 DAG built in {time.time() - t_build_seed:.2f}s")

    Strategy1_prefix = run_dir / f"{basename_common}_Strategy1"
    Strategy1_sol, Strategy1_verified, Strategy1_paths = solve_export_one(
        "Strategy1",
        matrix,
        matrix_name,
        args.index,
        Strategy1_lib,
        args.depth,
        args.num_runs,
        args.seed,
        args.lambda2,
        args.lambda3,
        args.terminal_width,
        args.terminal_parent_pool_limit,
        args.time_limit_seed,
        args.gurobi_output,
        args.random_tests,
        Strategy1_prefix,
    )

    print("\n[2/3] Extending the existing Strategy1 DAG by Strategy2 TD-OTD; no LCB-BP rerun")
    t_build_td = time.time()
    Strategy2_lib = add_td_otd_completion(
        seed_lib=Strategy1_lib,
        depth=args.depth,
        lambda2=args.lambda2,
        lambda3=args.lambda3,
        terminal_width=args.terminal_width,
        terminal_parent_pool_limit=args.terminal_parent_pool_limit,
        verbose=True,
    )
    print(f"[build] Strategy2 DAG built in {time.time() - t_build_td:.2f}s")

    Strategy2_prefix = run_dir / f"{basename_common}_Strategy2_W{width_tag}"
    Strategy2_sol, Strategy2_verified, Strategy2_paths = solve_export_one(
        "Strategy2",
        matrix,
        matrix_name,
        args.index,
        Strategy2_lib,
        args.depth,
        args.num_runs,
        args.seed,
        args.lambda2,
        args.lambda3,
        args.terminal_width,
        args.terminal_parent_pool_limit,
        args.time_limit_td_otd,
        args.gurobi_output,
        args.random_tests,
        Strategy2_prefix,
    )

    print("\n[3/3] Writing compact summary")
    summary_path = root_out_dir / "summary_all.csv"
    append_summary(summary_path, summary_row(
        matrix_name=matrix_name,
        index=args.index,
        depth=args.depth,
        num_runs=args.num_runs,
        seed=args.seed,
        Strategy_tag="Strategy1",
        lib=Strategy1_lib,
        sol=Strategy1_sol,
        verified=Strategy1_verified,
        terminal_width=args.terminal_width,
        terminal_parent_pool_limit=args.terminal_parent_pool_limit,
        paths=Strategy1_paths,
    ))
    append_summary(summary_path, summary_row(
        matrix_name=matrix_name,
        index=args.index,
        depth=args.depth,
        num_runs=args.num_runs,
        seed=args.seed,
        Strategy_tag="Strategy2",
        lib=Strategy2_lib,
        sol=Strategy2_sol,
        verified=Strategy2_verified,
        terminal_width=args.terminal_width,
        terminal_parent_pool_limit=args.terminal_parent_pool_limit,
        paths=Strategy2_paths,
    ))

    print(f"[summary] {summary_path}")
    print("\nDone.")


if __name__ == "__main__":
    main()
