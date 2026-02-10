from __future__ import annotations
import numpy as np
from ..flowsheet_tools import UnitOp

class Mixer(UnitOp):
    def apply(self) -> None:
        out = self.outlets["out"]
        fin = np.zeros(out.reg.n(), dtype=float)
        for s in self.inlets.values():
            fin += s.to_dense()
        out.from_dense(fin)
