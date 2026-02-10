from __future__ import annotations

from typing import Dict, List, Optional, Tuple
import numpy as np

from ..flowsheet_tools import UnitOp


class MultiStageCentrifugalContactorStrip(UnitOp):
    def __init__(
        self,
        name: str,
        D: Dict[str, float],
        required_recovery_to_aq: Dict[str, float],
        nonextract_to_aq: Optional[List[str]] = None,
        nonextract_to_org: Optional[List[str]] = None,
        N_max: int = 50,
        AO_max: float = 1e3,
        N_balance_cap: int = 20,
    ):
        super().__init__(name)
        self.D = {k: float(v) for k, v in D.items()}
        self.required_recovery_to_aq = {k: float(v) for k, v in required_recovery_to_aq.items()}
        self.nonextract_to_aq = set(nonextract_to_aq or [])
        self.nonextract_to_org = set(nonextract_to_org or [])
        self.N_max = int(N_max)
        self.AO_max = float(AO_max)
        self.N_balance_cap = int(N_balance_cap)
        self.N_used: int = 1
        self.AO_used: float = 0.0

    @staticmethod
    def _yN_over_y0(S: float, N: int) -> float:
        if abs(S - 1.0) < 1e-12:
            return 1.0 / (N + 1.0)
        ssum = (S ** (N + 1) - 1.0) / (S - 1.0)
        return 1.0 / ssum

    @classmethod
    def _recovery_to_aq(cls, D: float, AO: float, N: int) -> float:
        S = AO / D
        return 1.0 - cls._yN_over_y0(S, N)

    @classmethod
    def _required_AO_for_recovery_to_aq(cls, D: float, N: int, recovery: float, AO_max: float) -> float:
        lo, hi = 0.0, 1.0
        while cls._recovery_to_aq(D, hi, N) < recovery:
            hi *= 2.0
            if hi > AO_max:
                raise ValueError("Required A/O exceeds AO_max.")
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if cls._recovery_to_aq(D, mid, N) >= recovery:
                hi = mid
            else:
                lo = mid
        return hi

    def _design_N_and_AO(self) -> Tuple[int, float]:
        N_cap = min(self.N_max, self.N_balance_cap)
        best_N, best_AO, best_J = None, None, float("inf")
        for N in range(1, N_cap + 1):
            AO_req = 0.0
            for sp, rec in self.required_recovery_to_aq.items():
                AO_req = max(AO_req, self._required_AO_for_recovery_to_aq(self.D[sp], N, rec, self.AO_max))
            J = AO_req * N
            if J < best_J:
                best_J, best_N, best_AO = J, N, AO_req
        if best_N is None or best_AO is None:
            raise ValueError("No feasible (N, A/O) found.")
        return best_N, best_AO

    def apply(self) -> None:
        org_in = self.inlets["org_in"]
        aq_in = self.inlets["aq_in"]
        org_out = self.outlets["org_out"]
        aq_out = self.outlets["aq_out"]

        reg = org_in.reg
        idx = reg.index()

        N, AO = self._design_N_and_AO()
        self.N_used = N
        self.AO_used = AO

        A = aq_in.to_dense()
        O = org_in.to_dense()

        Aq = np.zeros(reg.n(), dtype=float)
        Org = np.zeros(reg.n(), dtype=float)

        for sp in reg.species:
            j = idx[sp]
            nA, nO = A[j], O[j]
            nT = nA + nO
            if nT == 0.0:
                continue

            if sp in self.nonextract_to_aq:
                Aq[j] = nT
                continue
            if sp in self.nonextract_to_org:
                Org[j] = nT
                continue

            if sp in self.D and self.D[sp] > 0:
                S = AO / self.D[sp]
                yratio = self._yN_over_y0(S, N)
                nO_out = yratio * nO
                nA_out = (nO - nO_out) + nA
                Org[j] = max(nO_out, 0.0)
                Aq[j] = max(nA_out, 0.0)
            else:
                Aq[j] = nA
                Org[j] = nO

        aq_out.from_dense(Aq)
        org_out.from_dense(Org)
