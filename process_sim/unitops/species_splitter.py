# =============================================================================
# 1) unitops/species_splitter.py  (or just add into your existing unitops.py)
# =============================================================================
from __future__ import annotations

from typing import Dict
import numpy as np

from ..flowsheet_tools import UnitOp


class SpeciesSplitter(UnitOp):
    """
    Split selected species to outlet A by specified fractions.
    Unspecified species default to default_to_A (typically 0, i.e. stay in B).

    Ports:
      in  : inlet stream
      A   : outlet A
      B   : outlet B
    """

    def __init__(self, name: str, frac_to_A: Dict[str, float], default_to_A: float = 0.0):
        super().__init__(name)
        self.frac_to_A = {k: float(v) for k, v in frac_to_A.items()}
        self.default_to_A = float(default_to_A)
        for sp, f in self.frac_to_A.items():
            if not (0.0 <= f <= 1.0):
                raise ValueError(f"{name}: frac_to_A[{sp}] must be in [0,1]")
        if not (0.0 <= self.default_to_A <= 1.0):
            raise ValueError(f"{name}: default_to_A must be in [0,1]")

    def apply(self) -> None:
        Sin = self.inlets["in"]
        SA = self.outlets["A"]
        SB = self.outlets["B"]

        reg = Sin.reg
        idx = reg.index()
        Fin = Sin.to_dense()

        FA = np.zeros(reg.n(), dtype=float)
        FB = np.zeros(reg.n(), dtype=float)

        for sp in reg.species:
            j = idx[sp]
            n = Fin[j]
            if n == 0.0:
                continue
            fA = self.frac_to_A.get(sp, self.default_to_A)
            FA[j] = fA * n
            FB[j] = (1.0 - fA) * n

        SA.from_dense(FA)
        SB.from_dense(FB)


class UO2NitrateToUO3Reactor(UnitOp):
    """
    Simple "complete consumption" calciner-style reactor.

    User spec:
      - Feed: F24
      - All U in the form UO2(NO3)2 -> UO3 (same molar flow) exits solids stream F26
      - Offgas stream F29 contains NO2, O2, H2O with relevant stoichiometry
      - All other species assumed consumed (i.e. do not appear in outlets)

    Stoichiometry assumption (to include H2O in offgas):
      UO2(NO3)2·6H2O  ->  UO3  +  2 NO2  +  0.5 O2  +  6 H2O

    If you want the anhydrous nitrate decomposition instead (no H2O produced),
    set hydrate_n = 0.0.
    """

    def __init__(
        self,
        name: str,
        *,
        uo2n_name: str = "UO2(NO3)2",
        uo3_name: str = "UO3",
        no2_name: str = "NO2",
        o2_name: str = "O2",
        h2o_name: str = "H2O",
        hydrate_n: float = 6.0,
        print_diagnostics: bool = False,
    ):
        super().__init__(name)
        self.uo2n_name = str(uo2n_name)
        self.uo3_name = str(uo3_name)
        self.no2_name = str(no2_name)
        self.o2_name = str(o2_name)
        self.h2o_name = str(h2o_name)
        self.hydrate_n = float(hydrate_n)
        self.print_diagnostics = bool(print_diagnostics)

    def apply(self) -> None:
        Sin = self.inlets["in"]
        Ssol = self.outlets["solids"]   # F26
        Sgas = self.outlets["offgas"]   # F29

        n_uo2n = Sin.get(self.uo2n_name)

        # Clear outlets then write only what the spec allows
        Ssol.mol = {}
        Sgas.mol = {}

        if n_uo2n > 0.0:
            # solids: UO3 at same molar rate
            Ssol.set(self.uo3_name, n_uo2n)

            # offgas products
            Sgas.set(self.no2_name, 2.0 * n_uo2n)
            Sgas.set(self.o2_name, 0.5 * n_uo2n)
            if self.hydrate_n > 0.0:
                Sgas.set(self.h2o_name, self.hydrate_n * n_uo2n)

        n_h2o_in = Sin.get(self.h2o_name)
        if n_h2o_in > 0.0:
            Sgas.set(self.h2o_name, n_h2o_in)

        if self.print_diagnostics:
            print(f"[{self.name}] n(UO2(NO3)2)={n_uo2n:.6g} -> "
                  f"UO3={Ssol.get(self.uo3_name):.6g}, "
                  f"NO2={Sgas.get(self.no2_name):.6g}, "
                  f"O2={Sgas.get(self.o2_name):.6g}, "
                  f"H2O={Sgas.get(self.h2o_name):.6g}")
