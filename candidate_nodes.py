# -*- coding: utf-8 -*-
"""candidate_nodes.py

Paper-oriented candidate DAG construction for depth-constrained XOR2/XOR3
linear-layer synthesis.

This file implements two Strategys:

Strategy 1: GDA-Seed
    Global Depth-Aware LCB-BP Seed DAG.
    The LCB-BP traces are used only as candidate-discovery devices, not as
    fixed schedules.  Candidate producers are merged first, and depth feasibility
    is then recomputed globally over the merged DAG.

Strategy 2: SRTC
    Seed-DAG-Guided Shared Residual Terminal Completion.
    This Strategy extends a finished seed DAG without re-running LCB-BP.  It extracts
    a depth-aware parent pool from the seed DAG, constructs residual bridges from
    that pool, and admits only bridges that support at least two output targets.

The method is intentionally budget-free: no top-K terminal budget and no
per-target truncation are used in the SRTC Strategy.  The only structural control is
shared residual support.
"""

from __future__ import annotations

import random
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, List, Optional, Sequence, Set, Tuple

import gurobipy as gp
from gurobipy import GRB

Value = int
Pair = Tuple[Value, Value]
Triple = Tuple[Value, Value, Value]
Impl = Tuple[int, ...]
Matrix = List[List[int]]


@dataclass
class CandidateLibrary:
    bit_width: int
    inputs: List[Value]
    targets: List[Value]
    row_targets: List[Value]
    values: List[Value]
    pair_impls: Dict[Value, Set[Pair]]
    triple_impls: Dict[Value, Set[Triple]]
    names: Dict[Value, str]
    target_outputs: Dict[Value, List[int]]
    earliest_layer: Dict[Value, int]
    meta: Dict[str, object] = field(default_factory=dict)
    impl_sources: Dict[Tuple[Value, Impl], Set[str]] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# BP-format txt parser
# ---------------------------------------------------------------------------

def _clean_lines(path: str | Path) -> List[str]:
    lines: List[str] = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("#") or line.startswith("//"):
            continue
        lines.append(line)
    return lines


def load_bp_txt_matrices(path: str | Path) -> List[Matrix]:
    lines = _clean_lines(path)
    if not lines:
        raise ValueError(f"Empty matrix file: {path}")

    first = [int(x) for x in lines[0].split()]
    idx = 0
    if len(first) == 1 and len(lines) >= 2 and len(lines[1].split()) == 2:
        count = first[0]
        idx = 1
    elif len(first) == 2:
        count = 1
        idx = 0
    else:
        raise ValueError(f"Cannot parse first line of {path!s}: {lines[0]!r}")

    matrices: List[Matrix] = []
    for m_idx in range(count):
        if idx >= len(lines):
            raise ValueError(f"Unexpected EOF while reading dimensions of matrix {m_idx} in {path}")
        dims = [int(x) for x in lines[idx].split()]
        idx += 1
        if len(dims) != 2:
            raise ValueError(f"Expected two dimensions, got {dims} in {path}")
        rows, cols = dims
        mat: Matrix = []
        for r in range(rows):
            if idx >= len(lines):
                raise ValueError(f"Unexpected EOF while reading row {r} in {path}")
            row = [int(x) for x in lines[idx].split()]
            idx += 1
            if len(row) != cols:
                raise ValueError(f"Row length mismatch in {path}, row {r}: expected {cols}, got {len(row)}")
            if any(bit not in (0, 1) for bit in row):
                raise ValueError(f"Non-binary entry in {path}, row {r}: {row}")
            mat.append(row)
        matrices.append(mat)
    return matrices


def load_one_matrix(path: str | Path, index: int = 0) -> Matrix:
    mats = load_bp_txt_matrices(path)
    if index < 0 or index >= len(mats):
        raise IndexError(f"Matrix index {index} out of range for {path}; file contains {len(mats)} matrix/matrices")
    return mats[index]


# ---------------------------------------------------------------------------
# Matrix and LCB-BP seed helpers
# ---------------------------------------------------------------------------

def parse_matrix(matrix: Matrix) -> Tuple[List[Value], List[Value], List[Value], Set[Value], Dict[Value, str], Dict[Value, List[int]], int]:
    if not matrix or not matrix[0]:
        raise ValueError("matrix is empty")
    rows = len(matrix)
    cols = len(matrix[0])
    inputs = [1 << i for i in range(cols)]
    targets: List[Value] = []
    row_targets: List[Value] = []
    target_set: Set[Value] = set()
    names: Dict[Value, str] = {}
    target_outputs: Dict[Value, List[int]] = defaultdict(list)

    for r, row in enumerate(matrix):
        if len(row) != cols:
            raise ValueError(f"row {r} has length {len(row)}, expected {cols}")
        value = 0
        for c, bit in enumerate(row):
            if bit:
                value |= 1 << c
        row_targets.append(value)
        target_outputs[value].append(r)
        if value not in target_set:
            targets.append(value)
            target_set.add(value)
            names[value] = f"y{r}"
    return inputs, targets, row_targets, target_set, names, dict(target_outputs), cols


def xor_all(values: Iterable[Value]) -> Value:
    out = 0
    for v in values:
        out ^= int(v)
    return out


def build_minlin_model(current_basis: Sequence[Value], target: Value, bit_width: int, fixed_weight: Optional[int] = None):
    n = len(current_basis)
    model = gp.Model("MinLinComb")
    model.setParam("OutputFlag", 0)
    select = model.addVars(n, vtype=GRB.BINARY, name="c")

    for j in range(bit_width):
        idxs = [i for i, v in enumerate(current_basis) if (v >> j) & 1]
        target_bit = (target >> j) & 1
        parity_aux = model.addVar(vtype=GRB.INTEGER, lb=0, ub=n, name=f"N_{j}")
        model.addConstr(gp.quicksum(select[i] for i in idxs) + target_bit == 2 * parity_aux)

    weight_expr = gp.quicksum(select[i] for i in range(n))
    if fixed_weight is None:
        model.setObjective(weight_expr, GRB.MINIMIZE)
    else:
        model.addConstr(weight_expr == fixed_weight)
        model.setObjective(0, GRB.MINIMIZE)
    return model, select


def gen_minimal_lin_combs(
    current_basis: Sequence[Value],
    target: Value,
    bit_width: int,
    max_solutions: Optional[int] = None,
) -> List[Tuple[Value, ...]]:
    """Enumerate minimum-weight linear representations of target over current_basis."""
    model, select = build_minlin_model(current_basis, target, bit_width, fixed_weight=None)
    model.optimize()
    if model.status != GRB.OPTIMAL:
        return []
    optimal_weight = int(round(model.ObjVal))

    enum_model, enum_select = build_minlin_model(current_basis, target, bit_width, fixed_weight=optimal_weight)
    solutions: List[Tuple[Value, ...]] = []
    seen: Set[Tuple[Value, ...]] = set()
    n = len(current_basis)

    while True:
        enum_model.optimize()
        if enum_model.status != GRB.OPTIMAL:
            break
        combo = tuple(sorted(current_basis[i] for i in range(n) if enum_select[i].X > 0.5))
        if combo in seen:
            break
        seen.add(combo)
        solutions.append(combo)
        if max_solutions is not None and len(solutions) >= max_solutions:
            break
        ones = [i for i in range(n) if enum_select[i].X > 0.5]
        zeros = [i for i in range(n) if enum_select[i].X <= 0.5]
        enum_model.addConstr(
            gp.quicksum(1 - enum_select[i] for i in ones) + gp.quicksum(enum_select[i] for i in zeros) >= 1
        )
    return solutions


def find_best_pair(
    current_basis: Sequence[Value],
    remaining_targets: Sequence[Value],
    bit_width: int,
    rng: random.Random,
    max_minlin_solutions: Optional[int] = None,
) -> Optional[Pair]:
    pair_counts: Counter[Pair] = Counter()
    current_basis_set = set(current_basis)
    for target in remaining_targets:
        min_combs = gen_minimal_lin_combs(current_basis, target, bit_width, max_solutions=max_minlin_solutions)
        direct_pairs = [tuple(c) for c in min_combs if len(c) == 2 and xor_all(c) not in current_basis_set]
        if direct_pairs:
            return tuple(rng.choice(direct_pairs))  # type: ignore[return-value]
        current_target_support: Set[Pair] = set()
        for combo in min_combs:
            elems = tuple(sorted(combo))
            for pair in combinations(elems, 2):
                current_target_support.add(tuple(pair))  # type: ignore[arg-type]
        pair_counts.update(current_target_support)

    if not pair_counts:
        return None
    max_freq = max(pair_counts.values())
    candidates = [pair for pair, freq in pair_counts.items() if freq == max_freq and (pair[0] ^ pair[1]) not in current_basis_set]
    if not candidates:
        candidates = [pair for pair, freq in pair_counts.items() if freq == max_freq]
    return tuple(rng.choice(candidates)) if candidates else None


# ---------------------------------------------------------------------------
# Producer helpers + source recording
# ---------------------------------------------------------------------------

def _source_key(value: Value, impl: Impl) -> Tuple[Value, Impl]:
    return int(value), tuple(sorted(int(x) for x in impl))


def _record_source(impl_sources: DefaultDict[Tuple[Value, Impl], Set[str]], value: Value, impl: Impl, source: str) -> None:
    impl_sources[_source_key(value, impl)].add(source)


def add_pair_impl(
    pair_impls: DefaultDict[Value, Set[Pair]],
    value_pool: Set[Value],
    a: Value,
    b: Value,
    impl_sources: Optional[DefaultDict[Tuple[Value, Impl], Set[str]]] = None,
    source: str = "unknown_pair",
) -> Value:
    pair: Pair = tuple(sorted((int(a), int(b))))  # type: ignore[assignment]
    value = pair[0] ^ pair[1]
    if value != 0 and pair[0] != pair[1]:
        pair_impls[value].add(pair)
        value_pool.add(value)
        if impl_sources is not None:
            _record_source(impl_sources, value, pair, source)
    return value


def add_triple_impl(
    triple_impls: DefaultDict[Value, Set[Triple]],
    value_pool: Set[Value],
    a: Value,
    b: Value,
    c: Value,
    impl_sources: Optional[DefaultDict[Tuple[Value, Impl], Set[str]]] = None,
    source: str = "unknown_triple",
) -> Value:
    triple: Triple = tuple(sorted((int(a), int(b), int(c))))  # type: ignore[assignment]
    value = triple[0] ^ triple[1] ^ triple[2]
    if value != 0 and len(set(triple)) == 3:
        triple_impls[value].add(triple)
        value_pool.add(value)
        if impl_sources is not None:
            _record_source(impl_sources, value, triple, source)
    return value


def add_impl_by_arity(
    value: Value,
    impl: Impl,
    pair_impls: DefaultDict[Value, Set[Pair]],
    triple_impls: DefaultDict[Value, Set[Triple]],
    value_pool: Set[Value],
    impl_sources: Optional[DefaultDict[Tuple[Value, Impl], Set[str]]] = None,
    source: str = "unknown_impl",
) -> bool:
    value = int(value)
    impl = tuple(sorted(int(x) for x in impl))
    if value == 0 or len(set(impl)) != len(impl):
        return False
    added = False
    if len(impl) == 2:
        a, b = impl
        if a != b and (a ^ b) == value:
            old = len(pair_impls[value])
            pair_impls[value].add((a, b))
            added = len(pair_impls[value]) > old
            if impl_sources is not None:
                _record_source(impl_sources, value, (a, b), source)
    elif len(impl) == 3:
        a, b, c = impl
        if (a ^ b ^ c) == value:
            old = len(triple_impls[value])
            triple_impls[value].add((a, b, c))
            added = len(triple_impls[value]) > old
            if impl_sources is not None:
                _record_source(impl_sources, value, (a, b, c), source)
    else:
        raise ValueError(f"Unsupported implementation arity: {len(impl)}")
    if added:
        value_pool.add(value)
        value_pool.update(impl)
    return added


def collect_all_input_layer1_candidates(inputs: Sequence[Value], pair_impls, triple_impls, value_pool, value_tags, impl_sources) -> None:
    inputs_sorted = sorted(inputs)
    for a, b in combinations(inputs_sorted, 2):
        value = add_pair_impl(pair_impls, value_pool, a, b, impl_sources, "input_pair_layer1")
        value_tags[value].add("input_pair_layer1")
    for a, b, c in combinations(inputs_sorted, 3):
        value = add_triple_impl(triple_impls, value_pool, a, b, c, impl_sources, "input_triple_layer1")
        value_tags[value].add("input_triple_layer1")


def collect_shortest_hints(min_combs, target_set, value_pool, pair_impls, triple_impls, value_tags, impl_sources) -> None:
    for combo in min_combs:
        elems = tuple(sorted(combo))
        if len(elems) >= 2:
            for a, b in combinations(elems, 2):
                value = add_pair_impl(pair_impls, value_pool, a, b, impl_sources, "shortest_pair_hint")
                value_tags[value].add("shortest_pair_hint")
        if len(elems) >= 3:
            for a, b, c in combinations(elems, 3):
                value = a ^ b ^ c
                anchored = ((a ^ b) in value_pool) or ((a ^ c) in value_pool) or ((b ^ c) in value_pool)
                direct_target = value in target_set
                if anchored or direct_target:
                    parts: List[str] = []
                    if anchored:
                        parts.append("anchored_triple_hint")
                    if direct_target:
                        parts.append("direct_target_triple_hint")
                    src = "+".join(parts) if parts else "shortest_triple_hint"
                    add_triple_impl(triple_impls, value_pool, a, b, c, impl_sources, src)
                    if anchored:
                        value_tags[value].add("anchored_triple_hint")
                    if direct_target:
                        value_tags[value].add("direct_target_triple_hint")


def add_trace_pair(
    current_basis: Sequence[Value],
    unfinished_targets: Sequence[Value],
    bit_width: int,
    rng: random.Random,
    pair_impls,
    value_pool,
    value_tags,
    impl_sources,
    max_minlin_solutions: Optional[int] = None,
) -> Optional[Value]:
    best_pair = find_best_pair(current_basis, unfinished_targets, bit_width, rng, max_minlin_solutions=max_minlin_solutions)
    if best_pair is None:
        return None
    a, b = best_pair
    value = add_pair_impl(pair_impls, value_pool, a, b, impl_sources, "trace_pair")
    value_tags[value].add("trace_pair")
    return value


# ---------------------------------------------------------------------------
# Global depth-aware normalization
# ---------------------------------------------------------------------------

def compute_earliest_layers(inputs: Set[Value], pair_impls, triple_impls, max_depth: int) -> Dict[Value, int]:
    INF = max_depth + 1
    earliest: Dict[Value, int] = {value: 0 for value in inputs}
    changed = True
    while changed:
        changed = False
        for value, impls in pair_impls.items():
            best = earliest.get(value, INF)
            for a, b in impls:
                if a in earliest and b in earliest:
                    best = min(best, 1 + max(earliest[a], earliest[b]))
            if best < earliest.get(value, INF) and best <= INF:
                earliest[value] = best
                changed = True
        for value, impls in triple_impls.items():
            best = earliest.get(value, INF)
            for a, b, c in impls:
                if a in earliest and b in earliest and c in earliest:
                    best = min(best, 1 + max(earliest[a], earliest[b], earliest[c]))
            if best < earliest.get(value, INF) and best <= INF:
                earliest[value] = best
                changed = True
    return earliest


def impl_layer(impl: Impl, earliest: Dict[Value, int]) -> Optional[int]:
    if not all(parent in earliest for parent in impl):
        return None
    return 1 + max(earliest[parent] for parent in impl)


def build_parent_to_targets_map(targets: Sequence[Value], pair_impls, triple_impls, earliest, depth: int) -> Dict[Value, Set[Value]]:
    parent_to_targets: Dict[Value, Set[Value]] = defaultdict(set)
    for target in targets:
        for impl in pair_impls.get(target, set()):
            layer = impl_layer(impl, earliest)
            if layer is not None and layer <= depth:
                for parent in impl:
                    parent_to_targets[parent].add(target)
        for impl in triple_impls.get(target, set()):
            layer = impl_layer(impl, earliest)
            if layer is not None and layer <= depth:
                for parent in impl:
                    parent_to_targets[parent].add(target)
    return parent_to_targets


def backward_useful_closure(inputs: Set[Value], targets: Sequence[Value], pair_impls, triple_impls, earliest, depth: int) -> Set[Value]:
    useful: Set[Value] = set(inputs)
    worklist: List[Value] = []
    for target in targets:
        if earliest.get(target, depth + 1) <= depth:
            useful.add(target)
            worklist.append(target)
    while worklist:
        value = worklist.pop()
        impls_to_check: List[Impl] = []
        impls_to_check.extend(pair_impls.get(value, set()))
        impls_to_check.extend(triple_impls.get(value, set()))
        for impl in impls_to_check:
            layer = impl_layer(impl, earliest)
            if layer is None or layer > depth:
                continue
            for parent in impl:
                if parent not in useful:
                    useful.add(parent)
                    if parent not in inputs:
                        worklist.append(parent)
    return useful


def classify_value_role(value: Value, targets: Set[Value], parent_to_targets, earliest, depth: int) -> str:
    if value in targets:
        return "terminal_target"
    layer = earliest.get(value, depth + 1)
    if layer == 1:
        support = len(parent_to_targets.get(value, set()))
        return "layer1_skeleton_shared" if support >= 2 else "layer1_skeleton_local"
    if layer == 2:
        return "layer2_absorb" if len(parent_to_targets.get(value, set())) > 0 else "layer2_bridge"
    return "late_or_unknown"


def prune_to_depth_aware_library(inputs, targets, pair_impls, triple_impls, value_tags, depth: int):
    inputs_set = set(inputs)
    targets_set = set(targets)
    earliest = compute_earliest_layers(inputs_set, pair_impls, triple_impls, max_depth=depth)
    useful = backward_useful_closure(inputs_set, targets, pair_impls, triple_impls, earliest, depth)
    parent_to_targets = build_parent_to_targets_map(targets, pair_impls, triple_impls, earliest, depth)

    kept_pair_impls: DefaultDict[Value, Set[Pair]] = defaultdict(set)
    kept_triple_impls: DefaultDict[Value, Set[Triple]] = defaultdict(set)
    kept_values: Set[Value] = set(inputs_set) | set(targets_set)

    for value, impls in pair_impls.items():
        if value not in useful:
            continue
        for impl in impls:
            layer = impl_layer(impl, earliest)
            if layer is None or layer > depth:
                continue
            if value not in targets_set and layer == depth:
                # A non-target value generated at final depth cannot feed any later gate.
                continue
            if not all(parent in useful for parent in impl):
                continue
            kept_pair_impls[value].add(impl)
            kept_values.add(value)
            kept_values.update(impl)

    for value, impls in triple_impls.items():
        if value not in useful:
            continue
        for impl in impls:
            layer = impl_layer(impl, earliest)
            if layer is None or layer > depth:
                continue
            if value not in targets_set and layer == depth:
                continue
            if not all(parent in useful for parent in impl):
                continue
            kept_triple_impls[value].add(impl)
            kept_values.add(value)
            kept_values.update(impl)

    role_tags: Dict[Value, str] = {value: classify_value_role(value, targets_set, parent_to_targets, earliest, depth) for value in kept_values}
    return kept_values, kept_pair_impls, kept_triple_impls, earliest, useful, role_tags, parent_to_targets


# ---------------------------------------------------------------------------
# Strategy 1: GDA-Seed
# ---------------------------------------------------------------------------

def _finalize_library(
    inputs: Sequence[Value],
    targets: Sequence[Value],
    row_targets: Sequence[Value],
    names: Dict[Value, str],
    target_outputs: Dict[Value, List[int]],
    bit_width: int,
    kept_values: Set[Value],
    kept_pair_impls: Dict[Value, Set[Pair]],
    kept_triple_impls: Dict[Value, Set[Triple]],
    earliest: Dict[Value, int],
    impl_sources: DefaultDict[Tuple[Value, Impl], Set[str]] | Dict[Tuple[Value, Impl], Set[str]],
    meta: Dict[str, object],
) -> CandidateLibrary:
    final_values: Set[Value] = set(kept_values) | set(inputs) | set(targets)
    for value, impls in kept_pair_impls.items():
        final_values.add(value)
        for impl in impls:
            final_values.update(impl)
    for value, impls in kept_triple_impls.items():
        final_values.add(value)
        for impl in impls:
            final_values.update(impl)

    final_impl_sources: Dict[Tuple[Value, Impl], Set[str]] = {}
    for value, impls in kept_pair_impls.items():
        for impl in impls:
            key = _source_key(value, impl)
            final_impl_sources[key] = set(impl_sources.get(key, {"untracked_pair"}))
    for value, impls in kept_triple_impls.items():
        for impl in impls:
            key = _source_key(value, impl)
            final_impl_sources[key] = set(impl_sources.get(key, {"untracked_triple"}))

    return CandidateLibrary(
        bit_width=bit_width,
        inputs=list(sorted(inputs)),
        targets=list(targets),
        row_targets=list(row_targets),
        values=list(sorted(final_values)),
        pair_impls={v: set(impls) for v, impls in kept_pair_impls.items()},
        triple_impls={v: set(impls) for v, impls in kept_triple_impls.items()},
        names=dict(names),
        target_outputs={v: list(rows) for v, rows in target_outputs.items()},
        earliest_layer={v: int(layer) for v, layer in earliest.items() if v in final_values},
        meta=dict(meta),
        impl_sources=final_impl_sources,
    )


def build_seed_candidate_library(
    matrix: Matrix,
    depth: int = 3,
    num_runs: int = 1,
    base_seed: int = 20260424,
    max_minlin_solutions: Optional[int] = None,
    verbose: bool = True,
) -> CandidateLibrary:
    """Build the Strategy-1 Global Depth-Aware LCB-BP Seed DAG.

    The first LCB-BP run uses exactly ``base_seed``.  Run i uses
    ``base_seed + i`` for i=0..num_runs-1.
    """
    inputs, targets, row_targets, target_set, names, target_outputs, bit_width = parse_matrix(matrix)

    value_pool: Set[Value] = set(inputs) | set(targets)
    pair_impls: DefaultDict[Value, Set[Pair]] = defaultdict(set)
    triple_impls: DefaultDict[Value, Set[Triple]] = defaultdict(set)
    value_tags: DefaultDict[Value, Set[str]] = defaultdict(set)
    impl_sources: DefaultDict[Tuple[Value, Impl], Set[str]] = defaultdict(set)

    for value in inputs:
        value_tags[value].add("input")
    for target in targets:
        value_tags[target].add("target")

    if verbose:
        print(f"[candidate:seed] collect complete layer-1 input XOR2/XOR3 for depth={depth}")
    collect_all_input_layer1_candidates(inputs, pair_impls, triple_impls, value_pool, value_tags, impl_sources)

    run_stats: List[Dict[str, object]] = []
    total_trace_steps = 0

    for run_idx in range(num_runs):
        seed = base_seed + run_idx
        rng = random.Random(seed)
        current_basis = list(inputs)
        current_basis_set = set(current_basis)
        rounds = 0
        if verbose:
            print(f"[candidate:seed] LCB-BP hint collection run {run_idx + 1}/{num_runs}, seed={seed}")
        while not target_set.issubset(current_basis_set):
            rounds += 1
            unfinished = [t for t in targets if t not in current_basis_set]
            for target in unfinished:
                min_combs = gen_minimal_lin_combs(current_basis, target, bit_width, max_solutions=max_minlin_solutions)
                if min_combs:
                    collect_shortest_hints(min_combs, target_set, value_pool, pair_impls, triple_impls, value_tags, impl_sources)
            trace_value = add_trace_pair(
                current_basis, unfinished, bit_width, rng, pair_impls, value_pool, value_tags, impl_sources,
                max_minlin_solutions=max_minlin_solutions,
            )
            if trace_value is None:
                raise RuntimeError(f"Candidate generation stopped early at run {run_idx + 1}: no feasible pair")
            if trace_value in current_basis_set:
                raise RuntimeError(f"Candidate generation stopped early at run {run_idx + 1}: duplicate new value")
            current_basis.append(trace_value)
            current_basis_set.add(trace_value)
            value_tags[trace_value].add("bp_derived_skeleton")
        trace_steps = len(current_basis) - len(inputs)
        total_trace_steps += trace_steps
        run_stats.append({
            "run": run_idx + 1,
            "seed": seed,
            "rounds": rounds,
            "trace_steps": trace_steps,
            "trace_basis_size": len(current_basis),
        })

    kept_values, kept_pair_impls, kept_triple_impls, earliest, useful, role_tags, parent_to_targets = prune_to_depth_aware_library(
        inputs, targets, pair_impls, triple_impls, value_tags, depth
    )

    feasible_targets = [t for t in targets if kept_pair_impls.get(t) or kept_triple_impls.get(t)]
    infeasible_targets = [names[t] for t in targets if not (kept_pair_impls.get(t) or kept_triple_impls.get(t))]

    meta = {
        "method": "GDA-Seed",
        "depth": depth,
        "num_runs": num_runs,
        "base_seed": base_seed,
        "max_minlin_solutions": max_minlin_solutions,
        "num_values": len(kept_values),
        "num_pair_impls": sum(len(v) for v in kept_pair_impls.values()),
        "num_triple_impls": sum(len(v) for v in kept_triple_impls.values()),
        "feasible_targets": len(feasible_targets),
        "total_targets": len(targets),
        "infeasible_targets_under_depth": infeasible_targets,
        "layer1_skeleton_candidates": sum(1 for v in kept_values if role_tags.get(v, "").startswith("layer1_")),
        "layer2_candidates": sum(1 for v in kept_values if role_tags.get(v, "").startswith("layer2_")),
        "trace_steps_total": total_trace_steps,
        "runs": run_stats,
        "srtc": {"enabled": False},
    }

    if verbose:
        print("\n[candidate:seed] depth-aware seed summary")
        print(f"  values kept            : {meta['num_values']}")
        print(f"  pair implementations   : {meta['num_pair_impls']}")
        print(f"  triple implementations : {meta['num_triple_impls']}")
        print(f"  feasible targets       : {meta['feasible_targets']}/{meta['total_targets']}")
        print(f"  layer1 candidates      : {meta['layer1_skeleton_candidates']}")
        print(f"  layer2 candidates      : {meta['layer2_candidates']}")
        if infeasible_targets:
            print(f"  infeasible targets     : {infeasible_targets}")

    return _finalize_library(
        inputs, targets, row_targets, names, target_outputs, bit_width,
        kept_values, kept_pair_impls, kept_triple_impls, earliest, impl_sources, meta,
    )


# ---------------------------------------------------------------------------
# Strategy 2: Seed-DAG-Guided Shared Residual Terminal Completion (SRTC)
# ---------------------------------------------------------------------------

def impl_cost(impl: Impl, lambda2: float, lambda3: float) -> float:
    if len(impl) == 2:
        return lambda2
    if len(impl) == 3:
        return lambda3
    raise ValueError(f"Unsupported implementation arity: {len(impl)}")


def _canonical_better(cand: Dict[str, object], cur: Optional[Dict[str, object]]) -> bool:
    if cur is None:
        return True
    return (
        (float(cand["cost"]), int(cand["layer"]), tuple(cand["impl"]))
        <
        (float(cur["cost"]), int(cur["layer"]), tuple(cur["impl"]))
    )


def _as_defaultdict_pair(pair_impls: Dict[Value, Set[Pair]]) -> DefaultDict[Value, Set[Pair]]:
    out: DefaultDict[Value, Set[Pair]] = defaultdict(set)
    for v, impls in pair_impls.items():
        out[v].update(tuple(i) for i in impls)
    return out


def _as_defaultdict_triple(triple_impls: Dict[Value, Set[Triple]]) -> DefaultDict[Value, Set[Triple]]:
    out: DefaultDict[Value, Set[Triple]] = defaultdict(set)
    for v, impls in triple_impls.items():
        out[v].update(tuple(i) for i in impls)
    return out


def add_srtc_completion(
    seed_lib: CandidateLibrary,
    depth: int,
    lambda2: float,
    lambda3: float,
    min_bridge_support: int = 2,
    seed_parent_values: Optional[Set[Value]] = None,
    verbose: bool = True,
) -> CandidateLibrary:
    """Extend a seed DAG by Seed-Solution-Guided SRTC.

    No LCB-BP run is executed here.  The method uses only ``seed_lib`` and,
    when provided, the shallow values selected by the seed-DAG MILP solution.

    Paper-oriented depth-D rule:
        parent pool A  = inputs ∪ {selected seed-circuit values with layer ≤ D-2}
        bridge pool B  = seed/generated values available by layer D-1
        keep bridge r  iff it supports at least ``min_bridge_support`` outputs
        terminal forms = y=p^r and y=p^q^r, plus direct shallow terminals.

    If ``seed_parent_values`` is None, the function falls back to the broader
    seed-DAG-wide parent pool {v in seed DAG | earliest(v) ≤ D-2}.  The runner
    in this project passes the seed solution core by default.
    """
    if depth < 1:
        raise ValueError("depth must be at least 1")
    if min_bridge_support < 1:
        raise ValueError("min_bridge_support must be at least 1")

    inputs = list(seed_lib.inputs)
    targets = list(seed_lib.targets)
    target_set = set(targets)
    names = dict(seed_lib.names)
    row_targets = list(seed_lib.row_targets)
    target_outputs = {v: list(rows) for v, rows in seed_lib.target_outputs.items()}
    bit_width = int(seed_lib.bit_width)

    pair_impls = _as_defaultdict_pair(seed_lib.pair_impls)
    triple_impls = _as_defaultdict_triple(seed_lib.triple_impls)
    value_pool: Set[Value] = set(seed_lib.values) | set(inputs) | set(targets)
    value_tags: DefaultDict[Value, Set[str]] = defaultdict(set)
    impl_sources: DefaultDict[Tuple[Value, Impl], Set[str]] = defaultdict(set)
    for key, srcs in seed_lib.impl_sources.items():
        impl_sources[key].update(srcs)

    seed_earliest = dict(seed_lib.earliest_layer)
    for x in inputs:
        seed_earliest[x] = 0

    parent_cut = max(0, depth - 2)
    bridge_cut = max(0, depth - 1)

    if seed_parent_values is None:
        # Fallback mode: broad seed-DAG-wide shallow pool.  This is retained for
        # ablation only; run_compare.py uses the seed-solution-guided pool.
        parent_pool_source = "seed_dag_wide_earliest"
        parent_candidates = set(v for v in value_pool if seed_earliest.get(v, depth + 1) <= parent_cut)
    else:
        # Main paper mode: use the certified seed solution as the shallow
        # skeleton for terminal completion.  This controls the expansion without
        # top-K/budget truncation and ties Strategy 2 directly to Strategy 1.
        parent_pool_source = "seed_solution_core"
        parent_candidates = set(int(v) for v in seed_parent_values)

    parent_candidates.update(inputs)
    parent_pool = sorted(
        v for v in parent_candidates
        if v in value_pool and (v in inputs or seed_earliest.get(v, depth + 1) <= parent_cut)
    )
    parent_set = set(parent_pool)
    parent_layers = {v: seed_earliest.get(v, 0 if v in inputs else depth + 1) for v in parent_pool}

    # bridge_info maps bridge value -> canonical implementation info.  Existing
    # seed values can act as bridges without a new implementation.
    bridge_info: Dict[Value, Dict[str, object]] = {}

    def consider_bridge(value: Value, impl: Optional[Impl], layer: int, source: str) -> None:
        if value == 0 or layer > bridge_cut:
            return
        if impl is None:
            cand = {"impl": tuple(), "layer": int(layer), "cost": 0.0, "source": source, "existing": True}
        else:
            impl = tuple(sorted(int(x) for x in impl))
            if len(set(impl)) != len(impl):
                return
            cand = {"impl": impl, "layer": int(layer), "cost": impl_cost(impl, lambda2, lambda3), "source": source, "existing": False}
        cur = bridge_info.get(value)
        if _canonical_better(cand, cur):
            bridge_info[value] = cand

    # Existing seed values that are already available before the terminal layer.
    for v, layer in seed_earliest.items():
        if v in value_pool and layer <= bridge_cut:
            consider_bridge(v, None, layer, "srtc_existing_seed_bridge")

    # Generated bridge values from the parent pool, using one deterministic
    # canonical implementation per bridge value.  This is not a top-K budget;
    # it is a value-level normalization.
    for a, b in combinations(parent_pool, 2):
        layer = 1 + max(parent_layers[a], parent_layers[b])
        consider_bridge(a ^ b, (a, b), layer, "srtc_generated_pair_bridge")
    for a, b, c in combinations(parent_pool, 3):
        layer = 1 + max(parent_layers[a], parent_layers[b], parent_layers[c])
        consider_bridge(a ^ b ^ c, (a, b, c), layer, "srtc_generated_triple_bridge")

    bridge_values = set(bridge_info)

    # First collect all direct terminals and all residual terminal plans.  Bridge
    # plans are activated only if their bridge supports enough distinct targets.
    direct_terminal_impls: Dict[Value, Set[Impl]] = defaultdict(set)
    bridge_terminal_impls: Dict[Value, List[Tuple[Impl, Value, str]]] = defaultdict(list)
    bridge_support_targets: Dict[Value, Set[Value]] = defaultdict(set)

    def terminal_layer(impl: Impl, bridge_value: Optional[Value] = None) -> Optional[int]:
        layers = []
        for p in impl:
            if bridge_value is not None and p == bridge_value:
                if p not in bridge_info:
                    return None
                layers.append(int(bridge_info[p]["layer"]))
            else:
                if p not in seed_earliest:
                    return None
                layers.append(int(seed_earliest[p]))
        if not layers:
            return None
        return 1 + max(layers)

    for target in targets:
        # Direct shallow pair: y = p ^ q.
        for p in parent_pool:
            q = target ^ p
            if q in parent_set and p < q:
                impl: Impl = tuple(sorted((p, q)))
                layer = terminal_layer(impl)
                if layer is not None and layer <= depth:
                    direct_terminal_impls[target].add(impl)

        # Direct shallow triple: y = p ^ q ^ s.
        for p, q in combinations(parent_pool, 2):
            s = target ^ p ^ q
            if s in parent_set and q < s:
                impl = tuple(sorted((p, q, s)))
                layer = terminal_layer(impl)
                if layer is not None and layer <= depth:
                    direct_terminal_impls[target].add(impl)

        # Residual pair terminal: y = p ^ r.
        seen_bridge_for_target: Set[Value] = set()
        for p in parent_pool:
            r = target ^ p
            if r in bridge_values and r != p:
                impl = tuple(sorted((p, r)))
                layer = terminal_layer(impl, bridge_value=r)
                if layer is not None and layer <= depth:
                    bridge_terminal_impls[target].append((impl, r, "srtc_terminal_pair_bridge"))
                    seen_bridge_for_target.add(r)

        # Residual triple terminal: y = p ^ q ^ r.
        for p, q in combinations(parent_pool, 2):
            r = target ^ p ^ q
            if r in bridge_values and r not in (p, q):
                impl = tuple(sorted((p, q, r)))
                layer = terminal_layer(impl, bridge_value=r)
                if layer is not None and layer <= depth:
                    bridge_terminal_impls[target].append((impl, r, "srtc_terminal_triple_bridge"))
                    seen_bridge_for_target.add(r)

        for r in seen_bridge_for_target:
            bridge_support_targets[r].add(target)

    kept_bridges = {r for r, ys in bridge_support_targets.items() if len(ys) >= min_bridge_support}

    added_direct_terminals = 0
    added_bridge_terminals = 0
    added_generated_bridges = 0

    # Add direct shallow terminals.
    for target, impls in direct_terminal_impls.items():
        for impl in impls:
            if add_impl_by_arity(
                target, impl, pair_impls, triple_impls, value_pool,
                impl_sources=impl_sources,
                source=f"srtc_direct_terminal:arity{len(impl)}",
            ):
                added_direct_terminals += 1

    # Add kept shared bridge implementations and their terminal decompositions.
    for r in sorted(kept_bridges):
        info = bridge_info[r]
        impl = tuple(info["impl"])
        if impl:
            if add_impl_by_arity(
                r, impl, pair_impls, triple_impls, value_pool,
                impl_sources=impl_sources,
                source=f"srtc_shared_bridge:{info.get('source', 'unknown')}:support{len(bridge_support_targets[r])}",
            ):
                added_generated_bridges += 1

    for target, items in bridge_terminal_impls.items():
        seen_impls: Set[Impl] = set()
        for impl, r, kind in items:
            if r not in kept_bridges:
                continue
            if impl in seen_impls:
                continue
            seen_impls.add(impl)
            if add_impl_by_arity(
                target, impl, pair_impls, triple_impls, value_pool,
                impl_sources=impl_sources,
                source=f"{kind}:support{len(bridge_support_targets[r])}",
            ):
                added_bridge_terminals += 1

    kept_values, kept_pair_impls, kept_triple_impls, earliest, useful, role_tags, parent_to_targets = prune_to_depth_aware_library(
        inputs, targets, pair_impls, triple_impls, value_tags, depth
    )

    feasible_targets = [t for t in targets if kept_pair_impls.get(t) or kept_triple_impls.get(t)]
    infeasible_targets = [names[t] for t in targets if not (kept_pair_impls.get(t) or kept_triple_impls.get(t))]

    meta = dict(seed_lib.meta)
    meta.update({
        "method": "GDA-Seed+SRTC",
        "num_values": len(kept_values),
        "num_pair_impls": sum(len(v) for v in kept_pair_impls.values()),
        "num_triple_impls": sum(len(v) for v in kept_triple_impls.values()),
        "feasible_targets": len(feasible_targets),
        "total_targets": len(targets),
        "infeasible_targets_under_depth": infeasible_targets,
        "layer1_skeleton_candidates": sum(1 for v in kept_values if role_tags.get(v, "").startswith("layer1_")),
        "layer2_candidates": sum(1 for v in kept_values if role_tags.get(v, "").startswith("layer2_")),
        "srtc": {
            "enabled": True,
            "parent_pool_source": parent_pool_source,
            "parent_cut_layer": parent_cut,
            "bridge_cut_layer": bridge_cut,
            "parent_pool_size": len(parent_pool),
            "seed_parent_values_provided": seed_parent_values is not None,
            "bridge_candidate_values": len(bridge_values),
            "bridges_with_any_terminal_support": len(bridge_support_targets),
            "min_bridge_support": min_bridge_support,
            "kept_shared_bridges": len(kept_bridges),
            "added_generated_bridges": added_generated_bridges,
            "added_direct_terminals": added_direct_terminals,
            "added_bridge_terminals": added_bridge_terminals,
        },
    })

    if verbose:
        s = meta["srtc"]
        print("\n[candidate:srtc] Shared Residual Terminal Completion summary")
        print(f"  parent pool source          : {s['parent_pool_source']}")
        print(f"  parent pool earliest <=     : {s['parent_cut_layer']}")
        print(f"  bridge available <=         : {s['bridge_cut_layer']}")
        print(f"  parent pool size            : {s['parent_pool_size']}")
        print(f"  bridge candidate values     : {s['bridge_candidate_values']}")
        print(f"  bridges with support        : {s['bridges_with_any_terminal_support']}")
        print(f"  min bridge support          : {s['min_bridge_support']}")
        print(f"  kept shared bridges         : {s['kept_shared_bridges']}")
        print(f"  added generated bridges     : {s['added_generated_bridges']}")
        print(f"  added direct terminals      : {s['added_direct_terminals']}")
        print(f"  added bridge terminals      : {s['added_bridge_terminals']}")
        print("\n[candidate:srtc] final depth-aware summary")
        print(f"  values kept                 : {meta['num_values']}")
        print(f"  pair implementations        : {meta['num_pair_impls']}")
        print(f"  triple implementations      : {meta['num_triple_impls']}")
        print(f"  feasible targets            : {meta['feasible_targets']}/{meta['total_targets']}")
        if infeasible_targets:
            print(f"  infeasible targets          : {infeasible_targets}")

    return _finalize_library(
        inputs, targets, row_targets, names, target_outputs, bit_width,
        kept_values, kept_pair_impls, kept_triple_impls, earliest, impl_sources, meta,
    )


# ---------------------------------------------------------------------------
# Strategy 2b: Target-Driven Output-Terminal Decomposition (TD-OTD)
# ---------------------------------------------------------------------------

def _plan_sort_key(plan: Dict[str, object]) -> Tuple[object, ...]:
    """Deterministic support-first Top-W ranking for terminal plans.

    Smaller is better.  The ordering keeps the candidate-control rule simple:
    first prioritize cross-output support of terminal parents, then local cost,
    then terminal layer and lexical order for reproducibility.
    """
    return (
        -int(plan.get("support_score", 0)),
        float(plan["local_cost"]),
        int(plan["target_layer"]),
        tuple(plan["target_impl"]),
    )


def _better_plan_for_same_impl(new_plan: Dict[str, object], old_plan: Optional[Dict[str, object]]) -> bool:
    """Choose a canonical representative when several derivations give the same terminal impl."""
    if old_plan is None:
        return True
    return _plan_sort_key(new_plan) < _plan_sort_key(old_plan)


def add_td_otd_completion(
    seed_lib: CandidateLibrary,
    depth: int,
    lambda2: float,
    lambda3: float,
    terminal_width: int = 512,
    terminal_parent_pool_limit: int = 3000,
    verbose: bool = True,
) -> CandidateLibrary:
    """Target-Driven Output-Terminal Decomposition Completion.

    This is the paper-oriented aggressive Strategy-2 expansion.

    No LCB-BP run is executed here.  The method starts from a completed seed DAG
    and adds target-wise terminal decompositions.

    General depth-D rule:
        parent pool A = inputs ∪ {seed-DAG value v | earliest(v) ≤ D-2}
        bridge pool B = existing/generated values available by layer D-1
        terminal forms per target y:
            y = p ^ q
            y = p ^ q ^ s
            y = p ^ r
            y = p ^ q ^ r
        where p,q,s ∈ A and r ∈ B.

    For each target, all depth-feasible plans are ranked by a deterministic
    support-first rule: cross-output terminal-parent support is prioritized first,
    and local gate cost is used second.  The first ``terminal_width`` plans are
    retained.  If terminal_width <= 0, all plans are retained.

    The same terminal-parent completion rule is used for D=2 and D=3:
        y = t1 ^ t2 or y = t1 ^ t2 ^ t3, with layer(ti) <= D-1.
    D=2 uses all layer-1 terminal parents; D=3 uses a support-ranked layer-2
    parent pool capped by ``terminal_parent_pool_limit``.  Only auxiliary parents
    used by retained plans are added to the candidate DAG.

    """
    if depth < 1:
        raise ValueError("depth must be at least 1")

    inputs = list(seed_lib.inputs)
    targets = list(seed_lib.targets)
    names = dict(seed_lib.names)
    row_targets = list(seed_lib.row_targets)
    target_outputs = {v: list(rows) for v, rows in seed_lib.target_outputs.items()}
    bit_width = int(seed_lib.bit_width)

    pair_impls = _as_defaultdict_pair(seed_lib.pair_impls)
    triple_impls = _as_defaultdict_triple(seed_lib.triple_impls)
    value_pool: Set[Value] = set(seed_lib.values) | set(inputs) | set(targets)
    value_tags: DefaultDict[Value, Set[str]] = defaultdict(set)
    impl_sources: DefaultDict[Tuple[Value, Impl], Set[str]] = defaultdict(set)
    for key, srcs in seed_lib.impl_sources.items():
        impl_sources[key].update(srcs)

    seed_earliest = dict(seed_lib.earliest_layer)
    for x in inputs:
        seed_earliest[x] = 0

    parent_cut = max(0, depth - 2)
    bridge_cut = max(0, depth - 1)

    parent_pool = sorted(
        v for v in value_pool
        if v in inputs or seed_earliest.get(v, depth + 1) <= parent_cut
    )
    parent_set = set(parent_pool)
    parent_layers = {v: seed_earliest.get(v, 0 if v in inputs else depth + 1) for v in parent_pool}

    # bridge_info maps bridge value -> canonical implementation. Existing seed
    # values available before the terminal layer can be used as zero-cost bridges.
    bridge_info: Dict[Value, Dict[str, object]] = {}

    def consider_bridge(value: Value, impl: Optional[Impl], layer: int, source: str) -> None:
        value = int(value)
        if value == 0 or layer > bridge_cut:
            return
        if impl is None:
            cand = {
                "impl": tuple(),
                "layer": int(layer),
                "cost": 0.0,
                "source": source,
                "existing": True,
            }
        else:
            impl = tuple(sorted(int(x) for x in impl))
            if len(set(impl)) != len(impl):
                return
            cand = {
                "impl": impl,
                "layer": int(layer),
                "cost": impl_cost(impl, lambda2, lambda3),
                "source": source,
                "existing": False,
            }
        cur = bridge_info.get(value)
        if _canonical_better(cand, cur):
            bridge_info[value] = cand

    # Existing seed-DAG values usable before the terminal layer.
    for value, layer in seed_earliest.items():
        if value in value_pool and layer <= bridge_cut:
            consider_bridge(value, None, layer, "td_otd_existing_seed_bridge")

    # Generated one-gate bridges from the shallow parent pool.
    for a, b in combinations(parent_pool, 2):
        layer = 1 + max(parent_layers[a], parent_layers[b])
        consider_bridge(a ^ b, (a, b), layer, "td_otd_generated_pair_bridge")
    for a, b, c in combinations(parent_pool, 3):
        layer = 1 + max(parent_layers[a], parent_layers[b], parent_layers[c])
        consider_bridge(a ^ b ^ c, (a, b, c), layer, "td_otd_generated_triple_bridge")

    bridge_values = set(bridge_info)

    plans_by_target: Dict[Value, Dict[Impl, Dict[str, object]]] = defaultdict(dict)
    bridge_support_targets: Dict[Value, Set[Value]] = defaultdict(set)

    def plan_bridge_values(plan: Dict[str, object]) -> List[Value]:
        """Return all auxiliary bridge values that a terminal plan depends on."""
        if "bridge_values" in plan and plan["bridge_values"] is not None:
            return [int(x) for x in plan["bridge_values"]]
        bridge_value = plan.get("bridge_value")
        if bridge_value is None:
            return []
        return [int(bridge_value)]

    def prune_target_plan_map(target: Value, keep: int) -> None:
        """Keep only the locally best plans for one target during large D=2 rescue enumeration."""
        if keep <= 0:
            return
        plan_map = plans_by_target.get(target)
        if not plan_map or len(plan_map) <= keep:
            return
        best = sorted(plan_map.values(), key=_plan_sort_key)[:keep]
        plans_by_target[target] = {tuple(p["target_impl"]): p for p in best}

    def add_plan(target: Value, plan: Dict[str, object]) -> None:
        impl = tuple(sorted(int(x) for x in plan["target_impl"]))
        if len(set(impl)) != len(impl) or any(x == 0 for x in impl):
            return
        plan["target_impl"] = impl
        old = plans_by_target[target].get(impl)
        if _better_plan_for_same_impl(plan, old):
            plans_by_target[target][impl] = plan

    for target in targets:
        # Direct shallow pair: y = p ^ q.
        for p in parent_pool:
            q = target ^ p
            if q in parent_set and p < q:
                layer = 1 + max(parent_layers[p], parent_layers[q])
                if layer <= depth:
                    add_plan(target, {
                        "target_impl": (p, q),
                        "bridge_value": None,
                        "bridge_impl": None,
                        "bridge_layer": None,
                        "bridge_is_new": 0,
                        "terminal_cost": lambda2,
                        "bridge_cost": 0.0,
                        "local_cost": lambda2,
                        "target_layer": layer,
                        "kind": "direct_shallow_pair",
                    })

        # Direct shallow triple: y = p ^ q ^ s.
        for p, q in combinations(parent_pool, 2):
            s = target ^ p ^ q
            if s in parent_set and q < s:
                layer = 1 + max(parent_layers[p], parent_layers[q], parent_layers[s])
                if layer <= depth:
                    add_plan(target, {
                        "target_impl": (p, q, s),
                        "bridge_value": None,
                        "bridge_impl": None,
                        "bridge_layer": None,
                        "bridge_is_new": 0,
                        "terminal_cost": lambda3,
                        "bridge_cost": 0.0,
                        "local_cost": lambda3,
                        "target_layer": layer,
                        "kind": "direct_shallow_triple",
                    })

        # Residual pair bridge: y = p ^ r.
        for p in parent_pool:
            r = target ^ p
            if r in bridge_values and r != p:
                binfo = bridge_info[r]
                blayer = int(binfo["layer"])
                layer = 1 + max(parent_layers[p], blayer)
                if layer <= depth:
                    terminal_cost = lambda2
                    bridge_cost = float(binfo["cost"])
                    add_plan(target, {
                        "target_impl": (p, r),
                        "bridge_value": r,
                        "bridge_impl": tuple(binfo["impl"]),
                        "bridge_layer": blayer,
                        "bridge_is_new": 0 if bool(binfo.get("existing")) else 1,
                        "terminal_cost": terminal_cost,
                        "bridge_cost": bridge_cost,
                        "local_cost": terminal_cost + bridge_cost,
                        "target_layer": layer,
                        "kind": "residual_pair_bridge",
                    })

        # Residual triple bridge: y = p ^ q ^ r.
        for p, q in combinations(parent_pool, 2):
            r = target ^ p ^ q
            if r in bridge_values and r not in (p, q):
                binfo = bridge_info[r]
                blayer = int(binfo["layer"])
                layer = 1 + max(parent_layers[p], parent_layers[q], blayer)
                if layer <= depth:
                    terminal_cost = lambda3
                    bridge_cost = float(binfo["cost"])
                    add_plan(target, {
                        "target_impl": (p, q, r),
                        "bridge_value": r,
                        "bridge_impl": tuple(binfo["impl"]),
                        "bridge_layer": blayer,
                        "bridge_is_new": 0 if bool(binfo.get("existing")) else 1,
                        "terminal_cost": terminal_cost,
                        "bridge_cost": bridge_cost,
                        "local_cost": terminal_cost + bridge_cost,
                        "target_layer": layer,
                        "kind": "residual_triple_bridge",
                    })

    # Unified terminal-parent completion statistics and implementation.
    terminal_parent_completion_counts: Dict[str, int] = {}
    terminal_parent_completion_targets: List[str] = []
    terminal_parent_scope = "disabled"
    terminal_parent_max_layer = max(0, depth - 1)
    terminal_parent_pool_size = 0
    terminal_parent_seed_existing = 0
    terminal_parent_generated = 0
    terminal_parent_ranked_candidates = 0

    def current_terminal_parent_support() -> Dict[Value, Set[Value]]:
        """Support of terminal parents induced by already enumerated TD-OTD plans."""
        support: Dict[Value, Set[Value]] = defaultdict(set)
        for yy, plan_map in plans_by_target.items():
            for pp in plan_map.values():
                for v in tuple(pp.get("target_impl", tuple())):
                    if v in bridge_info and int(bridge_info[v].get("layer", depth + 1)) <= terminal_parent_max_layer:
                        support[int(v)].add(yy)
                for v in plan_bridge_values(pp):
                    if v in bridge_info and int(bridge_info[v].get("layer", depth + 1)) <= terminal_parent_max_layer:
                        support[int(v)].add(yy)
        return support

    if depth in (2, 3):
        parent_support_targets = current_terminal_parent_support()
        parent_support: Dict[Value, int] = {int(v): len(ts) for v, ts in parent_support_targets.items()}

        all_terminal_parents = sorted(
            v for v, info in bridge_info.items()
            if int(info.get("layer", depth + 1)) <= terminal_parent_max_layer
        )
        existing_terminal_parents = [
            int(v) for v in all_terminal_parents
            if bool(bridge_info[v].get("existing"))
        ]
        generated_terminal_parents = [
            int(v) for v in all_terminal_parents
            if not bool(bridge_info[v].get("existing")) and tuple(bridge_info[v].get("impl", tuple()))
        ]
        terminal_parent_ranked_candidates = len(generated_terminal_parents)

        if depth == 2:
            # D=2: the layer<=1 terminal-parent pool is still small enough to use fully.
            terminal_parent_pool = all_terminal_parents
            terminal_parent_scope = "all_targets_layer1_full"
        else:
            # D=3: layer<=2 terminal parents can be very numerous.  Keep all
            # existing seed parents, and fill the remaining capacity with
            # generated parents ranked by support first and local cost second.
            limit = int(terminal_parent_pool_limit) if terminal_parent_pool_limit and terminal_parent_pool_limit > 0 else 3000
            capacity = max(0, limit - len(existing_terminal_parents))
            generated_ranked = sorted(
                generated_terminal_parents,
                key=lambda v: (
                    -int(parent_support.get(v, 0)),
                    float(bridge_info[v].get("cost", 0.0)),
                    int(bridge_info[v].get("layer", depth + 1)),
                    int(v),
                ),
            )[:capacity]
            terminal_parent_pool = sorted(set(existing_terminal_parents) | set(generated_ranked))
            terminal_parent_scope = "all_targets_layer2_controlled"

        terminal_parent_pool_size = len(terminal_parent_pool)
        terminal_parent_seed_existing = sum(1 for v in terminal_parent_pool if bool(bridge_info[v].get("existing")))
        terminal_parent_generated = terminal_parent_pool_size - terminal_parent_seed_existing
        terminal_parent_set = set(terminal_parent_pool)
        terminal_parent_layers = {v: int(bridge_info[v]["layer"]) for v in terminal_parent_pool}

        def terminal_parent_extra_cost(v: Value) -> float:
            info = bridge_info[v]
            if bool(info.get("existing")):
                return 0.0
            return float(info.get("cost", 0.0))

        def generated_parent_bridges(impl: Impl) -> Tuple[Value, ...]:
            vals: List[Value] = []
            for v in impl:
                info = bridge_info.get(v)
                if info is None:
                    continue
                if bool(info.get("existing")):
                    continue
                if tuple(info.get("impl", tuple())):
                    vals.append(int(v))
            return tuple(sorted(set(vals)))

        def support_score_for_impl(impl: Impl) -> int:
            # Coarse, parameter-free reuse score: terminal parents already seen
            # across more output decompositions are preferred.  Existing seed
            # parents get at least support 1 so reusable existing structure is not
            # discarded merely because it did not appear as a residual bridge.
            score = 0
            for v in impl:
                if bool(bridge_info.get(v, {}).get("existing")):
                    score += max(1, int(parent_support.get(v, 0)))
                else:
                    score += int(parent_support.get(v, 0))
            return int(score)

        prune_keep = max(terminal_width if terminal_width > 0 else 4096, 1024)

        for target in targets:
            before = len(plans_by_target.get(target, {}))

            # y = t1 ^ t2, layer(ti) <= D-1.
            for a in terminal_parent_pool:
                b = target ^ a
                if b in terminal_parent_set and a < b:
                    layer = 1 + max(terminal_parent_layers[a], terminal_parent_layers[b])
                    if layer <= depth:
                        impl = (a, b)
                        bridges = generated_parent_bridges(impl)
                        extra = sum(terminal_parent_extra_cost(v) for v in impl)
                        add_plan(target, {
                            "target_impl": impl,
                            "bridge_value": None,
                            "bridge_values": bridges,
                            "bridge_impl": None,
                            "bridge_layer": None,
                            "bridge_is_new": 1 if bridges else 0,
                            "terminal_cost": lambda2,
                            "bridge_cost": extra,
                            "local_cost": lambda2 + extra,
                            "target_layer": layer,
                            "support_score": support_score_for_impl(impl),
                            "kind": f"depth{depth}_terminal_parent_pair",
                        })

            # y = t1 ^ t2 ^ t3, layer(ti) <= D-1.
            for i, a in enumerate(terminal_parent_pool):
                for b in terminal_parent_pool[i + 1:]:
                    c = target ^ a ^ b
                    if c in terminal_parent_set and b < c:
                        layer = 1 + max(terminal_parent_layers[a], terminal_parent_layers[b], terminal_parent_layers[c])
                        if layer <= depth:
                            impl = (a, b, c)
                            bridges = generated_parent_bridges(impl)
                            extra = sum(terminal_parent_extra_cost(v) for v in impl)
                            add_plan(target, {
                                "target_impl": impl,
                                "bridge_value": None,
                                "bridge_values": bridges,
                                "bridge_impl": None,
                                "bridge_layer": None,
                                "bridge_is_new": 1 if bridges else 0,
                                "terminal_cost": lambda3,
                                "bridge_cost": extra,
                                "local_cost": lambda3 + extra,
                                "target_layer": layer,
                                "support_score": support_score_for_impl(impl),
                                "kind": f"depth{depth}_terminal_parent_triple",
                            })
                    if len(plans_by_target.get(target, {})) > 4 * prune_keep:
                        prune_target_plan_map(target, prune_keep)
                if len(plans_by_target.get(target, {})) > 4 * prune_keep:
                    prune_target_plan_map(target, prune_keep)

            prune_target_plan_map(target, prune_keep)
            added = len(plans_by_target.get(target, {})) - before
            terminal_parent_completion_counts[names.get(target, str(target))] = max(0, added)
            if added > 0:
                terminal_parent_completion_targets.append(names.get(target, str(target)))

    # Compute bridge support after all canonical plans are known.
    for target, plan_map in plans_by_target.items():
        for plan in plan_map.values():
            for bridge_value in plan_bridge_values(plan):
                bridge_support_targets[int(bridge_value)].add(target)

    selected_plans_by_target: Dict[Value, List[Dict[str, object]]] = {}
    selected_bridge_values: Set[Value] = set()
    target_plan_counts: Dict[str, int] = {}
    enumerated_plan_counts: Dict[str, int] = {}

    for target in targets:
        plans = list(plans_by_target.get(target, {}).values())
        enumerated_plan_counts[names.get(target, str(target))] = len(plans)
        plans.sort(key=_plan_sort_key)
        if terminal_width > 0:
            plans = plans[:terminal_width]
        selected_plans_by_target[target] = plans
        target_plan_counts[names.get(target, str(target))] = len(plans)
        for plan in plans:
            for bridge_value in plan_bridge_values(plan):
                selected_bridge_values.add(int(bridge_value))

    added_generated_bridges = 0
    added_terminal_impls = 0
    added_direct_terminals = 0
    added_residual_terminals = 0

    # Add only bridges used by selected terminal plans.
    for r in sorted(selected_bridge_values):
        binfo = bridge_info[r]
        bimpl = tuple(binfo["impl"])
        if bimpl:
            source = "td_otd_bridge"
            if add_impl_by_arity(r, bimpl, pair_impls, triple_impls, value_pool, impl_sources, source):
                added_generated_bridges += 1

    # Add terminal decompositions.
    for target, plans in selected_plans_by_target.items():
        for rank, plan in enumerate(plans, start=1):
            impl = tuple(plan["target_impl"])
            kind = str(plan["kind"])
            source = "td_otd_terminal"
            added = add_impl_by_arity(target, impl, pair_impls, triple_impls, value_pool, impl_sources, source)
            if added:
                added_terminal_impls += 1
                if kind.startswith("direct"):
                    added_direct_terminals += 1
                else:
                    added_residual_terminals += 1

    kept_values, kept_pair_impls, kept_triple_impls, earliest, useful, role_tags, parent_to_targets = prune_to_depth_aware_library(
        inputs, targets, pair_impls, triple_impls, value_tags, depth
    )

    feasible_targets = [t for t in targets if kept_pair_impls.get(t) or kept_triple_impls.get(t)]
    infeasible_targets = [names[t] for t in targets if not (kept_pair_impls.get(t) or kept_triple_impls.get(t))]

    meta = dict(seed_lib.meta)
    meta.update({
        "method": "GDA-Seed+TD-OTD",
        "num_values": len(kept_values),
        "num_pair_impls": sum(len(v) for v in kept_pair_impls.values()),
        "num_triple_impls": sum(len(v) for v in kept_triple_impls.values()),
        "feasible_targets": len(feasible_targets),
        "total_targets": len(targets),
        "infeasible_targets_under_depth": infeasible_targets,
        "layer1_skeleton_candidates": sum(1 for v in kept_values if role_tags.get(v, "").startswith("layer1_")),
        "layer2_candidates": sum(1 for v in kept_values if role_tags.get(v, "").startswith("layer2_")),
        "td_otd": {
            "enabled": True,
            "terminal_width": terminal_width,
            "parent_cut_layer": parent_cut,
            "bridge_cut_layer": bridge_cut,
            "parent_pool_source": "seed_dag_shallow_earliest",
            "parent_pool_size": len(parent_pool),
            "bridge_candidate_values": len(bridge_values),
            "bridges_with_any_terminal_support": len(bridge_support_targets),
            "selected_bridge_values": len(selected_bridge_values),
            "added_generated_bridges": added_generated_bridges,
            "added_terminal_impls": added_terminal_impls,
            "added_direct_terminals": added_direct_terminals,
            "added_residual_terminals": added_residual_terminals,
            "enumerated_plan_counts": enumerated_plan_counts,
            "selected_plan_counts": target_plan_counts,
            "terminal_parent_completion_enabled": depth in (2, 3),
            "terminal_parent_completion_scope": terminal_parent_scope,
            "terminal_parent_max_layer": terminal_parent_max_layer if depth in (2, 3) else 0,
            "terminal_parent_pool_size": terminal_parent_pool_size,
            "terminal_parent_pool_limit": terminal_parent_pool_limit if depth == 3 else terminal_parent_pool_size,
            "terminal_parent_seed_existing": terminal_parent_seed_existing,
            "terminal_parent_generated": terminal_parent_generated,
            "terminal_parent_ranked_candidates": terminal_parent_ranked_candidates,
            "terminal_parent_completion_targets": terminal_parent_completion_targets,
            "terminal_parent_completion_counts": terminal_parent_completion_counts,
            "topw_sort_rule": "support_first_local_cost_second",
        },
    })

    if verbose:
        t = meta["td_otd"]
        print("\n[candidate:td-otd] Target-Driven Output-Terminal Decomposition summary")
        print(f"  terminal width W          : {t['terminal_width']}")
        print(f"  parent pool source        : {t['parent_pool_source']}")
        print(f"  parent pool earliest <=   : {t['parent_cut_layer']}")
        print(f"  bridge available <=       : {t['bridge_cut_layer']}")
        print(f"  parent pool size          : {t['parent_pool_size']}")
        print(f"  bridge candidate values   : {t['bridge_candidate_values']}")
        print(f"  bridges with support      : {t['bridges_with_any_terminal_support']}")
        print(f"  selected bridge values    : {t['selected_bridge_values']}")
        print(f"  added generated bridges   : {t['added_generated_bridges']}")
        print(f"  added terminal impls      : {t['added_terminal_impls']}")
        print(f"    direct terminals        : {t['added_direct_terminals']}")
        print(f"    residual terminals      : {t['added_residual_terminals']}")
        if t.get("terminal_parent_completion_enabled"):
            print(
                f"  terminal-parent completion: scope={t.get('terminal_parent_completion_scope', 'unknown')}, "
                f"max_layer={t.get('terminal_parent_max_layer', '')}, "
                f"parent_pool={t.get('terminal_parent_pool_size', '')}, "
                f"limit={t.get('terminal_parent_pool_limit', '')}"
            )
            print(f"    seed/existing parents   : {t.get('terminal_parent_seed_existing', 0)}")
            print(f"    generated parents       : {t.get('terminal_parent_generated', 0)}")
            print(f"    support-ranked values   : {t.get('terminal_parent_ranked_candidates', 0)}")
            print(f"  completion improved targets: {t.get('terminal_parent_completion_targets', [])}")
        print("\n[candidate:td-otd] final depth-aware summary")
        print(f"  values kept               : {meta['num_values']}")
        print(f"  pair implementations      : {meta['num_pair_impls']}")
        print(f"  triple implementations    : {meta['num_triple_impls']}")
        print(f"  feasible targets          : {meta['feasible_targets']}/{meta['total_targets']}")
        if infeasible_targets:
            print(f"  infeasible targets        : {infeasible_targets}")

    return _finalize_library(
        inputs, targets, row_targets, names, target_outputs, bit_width,
        kept_values, kept_pair_impls, kept_triple_impls, earliest, impl_sources, meta,
    )
