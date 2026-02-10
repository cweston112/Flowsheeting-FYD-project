from __future__ import annotations

from typing import Dict, Iterable, Optional
import numpy as np

from ..flowsheet_tools import UnitOp


class SolventStripTransferSplitter(UnitOp):
    def __init__(
        self,
        name: str,
        *,
        frac_org_to_aq: Dict[str, float],
        force_org_species: Optional[Iterable[str]] = None,
        force_aq_species: Optional[Iterable[str]] = None,
    ):
        super().__init__(name)
        self.frac_org_to_aq = {k: float(v) for k, v in frac_org_to_aq.items()}
        for sp, f in self.frac_org_to_aq.items():
            if not (0.0 <= f <= 1.0):
                raise ValueError(f"{name}: frac_org_to_aq[{sp}] must be in [0,1]")
        self.force_org_species = set(force_org_species or [])
        self.force_aq_species = set(force_aq_species or [])

    def apply(self) -> None:
        org_in = self.inlets["org_in"]
        aq_in = self.inlets["aq_in"]
        org_out = self.outlets["org_out"]
        aq_out = self.outlets["aq_out"]

        reg = org_in.reg
        idx = reg.index()

        Oin = org_in.to_dense()
        Ain = aq_in.to_dense()

        Oout = Oin.copy()
        Aout = Ain.copy()

        for sp in reg.species:
            j = idx[sp]
            nO = Oin[j]
            if nO <= 0.0:
                continue

            if sp in self.force_org_species:
                continue

            if sp in self.force_aq_species:
                Aout[j] += Oout[j]
                Oout[j] = 0.0
                continue

            f = self.frac_org_to_aq.get(sp, 0.0)
            if f > 0.0:
                moved = f * nO
                Aout[j] += moved
                Oout[j] -= moved

        Oout[Oout < 0.0] = 0.0
        org_out.from_dense(Oout)
        aq_out.from_dense(Aout)
