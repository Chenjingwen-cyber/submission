# -*- coding: utf-8 -*-
"""Layered MILP for depth-D XOR2/XOR3 h-XOR minimization over a candidate DAG."""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any
import gurobipy as gp
from gurobipy import GRB
from candidate_nodes import CandidateLibrary, Value, Impl


@dataclass
class SelectedGate:
    value: Value
    impl: Impl
    arity: int
    layer: int
    cost: float


@dataclass
class MilpSolution:
    status: int
    status_name: str
    objective: Optional[float]
    xor2_count: int
    xor3_count: int
    selected_gates: List[SelectedGate] = field(default_factory=list)
    runtime: Optional[float] = None
    mip_gap: Optional[float] = None
    best_bound: Optional[float] = None


def _require_gurobi() -> None:
    if gp is None or GRB is None:
        raise RuntimeError("gurobipy is required to solve the layered MILP.")


def _status_name(status: int) -> str:
    if GRB is None:
        return str(status)
    names = {
        GRB.OPTIMAL: "OPTIMAL",
        GRB.INFEASIBLE: "INFEASIBLE",
        GRB.TIME_LIMIT: "TIME_LIMIT",
        GRB.INTERRUPTED: "INTERRUPTED",
        GRB.SUBOPTIMAL: "SUBOPTIMAL",
        GRB.INF_OR_UNBD: "INF_OR_UNBD",
        GRB.UNBOUNDED: "UNBOUNDED",
    }
    return names.get(status, str(status))


def _normalize_earliest_layer(lib: CandidateLibrary, inputs: set[Value]) -> Dict[Value, int]:
    earliest: Dict[Value, int] = {v: 0 for v in inputs}
    raw: Optional[Any] = getattr(lib, "earliest_layer", None)
    if isinstance(raw, dict):
        for k, v in raw.items():
            try:
                earliest[int(k)] = int(v)
            except Exception:
                pass
    return earliest


def _min_possible_layer(impl: Impl, earliest: Dict[Value, int]) -> int:
    if all(parent in earliest for parent in impl):
        return 1 + max(earliest[parent] for parent in impl)
    return 1


def solve_layered_milp(
    lib: CandidateLibrary,
    depth: int,
    lambda2: float,
    lambda3: float,
    time_limit: Optional[float] = None,
    output_flag: int = 1,
    mip_focus: Optional[int] = None,
    heuristics: Optional[float] = None,
    cuts: Optional[int] = None,
    best_obj_stop: Optional[float] = None,
    use_layer_pruning: bool = True,
) -> MilpSolution:
    """Solve h-XOR minimization inside the generated candidate DAG.

    The model uses tight availability and first-generation variables:
        A[v,l] = value v is available by end of layer l
        G[v,l] = non-input value v is first generated at layer l
        X/Z    = selected XOR2/XOR3 producer at layer l

    For every non-input value and layer:
        A[v,l] = A[v,l-1] + G[v,l]
        sum(producers of v at l) = G[v,l]
    """
    _require_gurobi()

    inputs = set(lib.inputs)
    targets = list(lib.targets)
    values = sorted(set(lib.values) | inputs | set(targets) | set(lib.row_targets))
    earliest = _normalize_earliest_layer(lib, inputs)

    model = gp.Model("DepthD_hXOR_MILP")
    model.setParam("OutputFlag", output_flag)
    if time_limit is not None and time_limit > 0:
        model.setParam("TimeLimit", time_limit)
    if mip_focus is not None:
        model.setParam("MIPFocus", mip_focus)
    if heuristics is not None:
        model.setParam("Heuristics", heuristics)
    if cuts is not None:
        model.setParam("Cuts", cuts)
    if best_obj_stop is not None:
        model.setParam("BestObjStop", best_obj_stop)

    A: Dict[Tuple[Value, int], gp.Var] = {}
    for v in values:
        for l in range(depth + 1):
            A[(v, l)] = model.addVar(vtype=GRB.BINARY, name=f"A_{v}_{l}")

    G: Dict[Tuple[Value, int], gp.Var] = {}
    for v in values:
        if v in inputs:
            continue
        for l in range(1, depth + 1):
            G[(v, l)] = model.addVar(vtype=GRB.BINARY, name=f"G_{v}_{l}")

    X: Dict[Tuple[Value, Impl, int], gp.Var] = {}
    Z: Dict[Tuple[Value, Impl, int], gp.Var] = {}
    producers_by_value_layer: Dict[Tuple[Value, int], List[gp.Var]] = {}

    for v, impls in lib.pair_impls.items():
        if v in inputs or v not in values:
            continue
        for impl0 in sorted(impls):
            impl = tuple(impl0)
            if len(impl) != 2 or not all(p in values for p in impl):
                continue
            min_l = _min_possible_layer(impl, earliest) if use_layer_pruning else 1
            if min_l > depth:
                continue
            for l in range(max(1, min_l), depth + 1):
                var = model.addVar(vtype=GRB.BINARY, name=f"X_{v}_{'_'.join(map(str, impl))}_{l}")
                X[(v, impl, l)] = var
                producers_by_value_layer.setdefault((v, l), []).append(var)

    for v, impls in lib.triple_impls.items():
        if v in inputs or v not in values:
            continue
        for impl0 in sorted(impls):
            impl = tuple(impl0)
            if len(impl) != 3 or not all(p in values for p in impl):
                continue
            min_l = _min_possible_layer(impl, earliest) if use_layer_pruning else 1
            if min_l > depth:
                continue
            for l in range(max(1, min_l), depth + 1):
                var = model.addVar(vtype=GRB.BINARY, name=f"Z_{v}_{'_'.join(map(str, impl))}_{l}")
                Z[(v, impl, l)] = var
                producers_by_value_layer.setdefault((v, l), []).append(var)

    model.update()

    for v in values:
        if v in inputs:
            for l in range(depth + 1):
                model.addConstr(A[(v, l)] == 1, name=f"input_{v}_{l}")
        else:
            model.addConstr(A[(v, 0)] == 0, name=f"noninput0_{v}")

    for v in values:
        if v in inputs:
            continue
        for l in range(1, depth + 1):
            model.addConstr(A[(v, l)] == A[(v, l - 1)] + G[(v, l)], name=f"flow_{v}_{l}")

    for (v, impl, l), var in X.items():
        for parent in impl:
            model.addConstr(var <= A[(parent, l - 1)], name=f"X_parent_{v}_{parent}_{l}")
        model.addConstr(var <= G[(v, l)], name=f"X_gen_{v}_{l}")

    for (v, impl, l), var in Z.items():
        for parent in impl:
            model.addConstr(var <= A[(parent, l - 1)], name=f"Z_parent_{v}_{parent}_{l}")
        model.addConstr(var <= G[(v, l)], name=f"Z_gen_{v}_{l}")

    for v in values:
        if v in inputs:
            continue
        for l in range(1, depth + 1):
            producers = producers_by_value_layer.get((v, l), [])
            if producers:
                model.addConstr(gp.quicksum(producers) == G[(v, l)], name=f"one_producer_{v}_{l}")
            else:
                model.addConstr(G[(v, l)] == 0, name=f"no_producer_{v}_{l}")

    for target in targets:
        model.addConstr(A[(target, depth)] == 1, name=f"target_{target}")

    objective = gp.quicksum(lambda2 * var for var in X.values()) + gp.quicksum(lambda3 * var for var in Z.values())
    model.setObjective(objective, GRB.MINIMIZE)
    model.optimize()

    status = model.status
    selected: List[SelectedGate] = []
    objective_value: Optional[float] = None
    mip_gap: Optional[float] = None
    best_bound: Optional[float] = None
    try:
        best_bound = float(model.ObjBound)
    except Exception:
        best_bound = None

    if status in (GRB.OPTIMAL, GRB.TIME_LIMIT, GRB.SUBOPTIMAL, GRB.INTERRUPTED) and model.SolCount > 0:
        objective_value = float(model.ObjVal)
        try:
            mip_gap = float(model.MIPGap)
        except Exception:
            mip_gap = None
        for (v, impl, l), var in X.items():
            if var.X > 0.5:
                selected.append(SelectedGate(value=v, impl=impl, arity=2, layer=l, cost=lambda2))
        for (v, impl, l), var in Z.items():
            if var.X > 0.5:
                selected.append(SelectedGate(value=v, impl=impl, arity=3, layer=l, cost=lambda3))
        selected.sort(key=lambda g: (g.layer, g.value, g.arity, g.impl))

    return MilpSolution(
        status=status,
        status_name=_status_name(status),
        objective=objective_value,
        xor2_count=sum(1 for g in selected if g.arity == 2),
        xor3_count=sum(1 for g in selected if g.arity == 3),
        selected_gates=selected,
        runtime=float(model.Runtime) if hasattr(model, "Runtime") else None,
        mip_gap=mip_gap,
        best_bound=best_bound,
    )
