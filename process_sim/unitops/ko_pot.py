from __future__ import annotations

from typing import Iterable, Optional
import numpy as np

from ..flowsheet_tools import UnitOp


class KOPot(UnitOp):
    def __init__(self, name: str, *, eta_KO: float = 0.999, droplet_set: Optional[Iterable[str]] = None):
        super().__init__(name)
        if not (0.0 <= eta_KO <= 1.0):
            raise ValueError("eta_KO must be in [0,1]")
        self.eta_KO = float(eta_KO)
        self.droplet_set = set(droplet_set or ["H2O", "HNO3", "I_aq"])

    def apply(self) -> None:
        Sin = self.inlets["in"]
        Sg = self.outlets["g"]
        Sl = self.outlets["l"]

        reg = Sin.reg
        idx = reg.index()
        Fin = Sin.to_dense()

        Fg = np.zeros(reg.n(), dtype=float)
        Fl = np.zeros(reg.n(), dtype=float)

        for sp in reg.species:
            j = idx[sp]
            n = Fin[j]
            if n == 0.0:
                continue
            if sp in self.droplet_set:
                Fl[j] = self.eta_KO * n
                Fg[j] = (1.0 - self.eta_KO) * n
            else:
                Fg[j] = n

        Sg.from_dense(Fg)
        Sl.from_dense(Fl)
