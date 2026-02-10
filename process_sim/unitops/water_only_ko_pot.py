from __future__ import annotations

import numpy as np

from ..flowsheet_tools import UnitOp


class WaterOnlyKOPot(UnitOp):
    """Removes all water to liquid; everything else stays in gas."""

    def __init__(self, name: str, *, water_name: str = "H2O"):
        super().__init__(name)
        self.water_name = water_name

    def apply(self) -> None:
        Sin = self.inlets["in"]
        Sg = self.outlets["g"]
        Sl = self.outlets["l"]

        reg = Sin.reg
        idx = reg.index()
        Fin = Sin.to_dense()

        Fg = Fin.copy()
        Fl = np.zeros(reg.n(), dtype=float)

        if self.water_name in idx:
            j = idx[self.water_name]
            Fl[j] = Fg[j]
            Fg[j] = 0.0

        Sg.from_dense(Fg)
        Sl.from_dense(Fl)
