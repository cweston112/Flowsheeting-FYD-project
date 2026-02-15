from __future__ import annotations

from dataclasses import dataclass

from ..flowsheet_tools import UnitOp, EPS

class NH3NOxReductionColumn(UnitOp):
    """
    Simple NH3-based NOx reduction bookkeeping.

    Removes a fixed fraction of NO and NO2 (separately) by reaction with NH3.
    Optionally limited by available O2 and NH3 (if feed undersupplied).

    Stoichiometries used:
      NO:  4 NH3 + 4 NO  + 1 O2 -> 4 N2 + 6 H2O
           per NO: NH3=1, O2=0.25, N2=1, H2O=1.5

      NO2: 4 NH3 + 2 NO2 + 1 O2 -> 3 N2 + 6 H2O
           per NO2: NH3=2, O2=0.5, N2=1.5, H2O=3
    """
    def __init__(
        self,
        name: str,
        *,
        removal_frac: float,
        no2_name: str = "NO2",
        no_name: str = "NO",
        nh3_name: str = "NH3",
        o2_name: str = "O2",
        n2_name: str = "N2",
        h2o_name: str = "H2O",
        print_diagnostics: bool = False,
    ):
        super().__init__(name)
        if not (0.0 <= removal_frac <= 1.0):
            raise ValueError("removal_frac must be in [0,1]")
        self.removal_frac = float(removal_frac)
        self.no2_name = no2_name
        self.no_name = no_name
        self.nh3_name = nh3_name
        self.o2_name = o2_name
        self.n2_name = n2_name
        self.h2o_name = h2o_name
        self.print_diagnostics = bool(print_diagnostics)

    def apply(self) -> None:
        gas = self.inlets["gas_in"]
        nh3 = self.inlets["nh3_in"]
        out = self.outlets["out"]

        # start from combined inlet contents
        out.mol = {}
        for s in (gas, nh3):
            for sp, n in s.mol.items():
                if n > EPS:
                    out.set(sp, out.get(sp) + float(n))

        NO2 = out.get(self.no2_name)
        NO  = out.get(self.no_name)
        NH3 = out.get(self.nh3_name)
        O2  = out.get(self.o2_name)

        # targets to remove
        dNO2_target = self.removal_frac * NO2
        dNO_target  = self.removal_frac * NO

        # Required reactants for those targets
        # NO2: NH3=2*dNO2, O2=0.5*dNO2
        # NO : NH3=1*dNO,  O2=0.25*dNO
        NH3_req = 2.0 * dNO2_target + 1.0 * dNO_target
        O2_req  = 0.5 * dNO2_target + 0.25 * dNO_target

        # Scale down if NH3 or O2 insufficient
        scale = 1.0
        if NH3_req > EPS:
            scale = min(scale, NH3 / NH3_req)
        if O2_req > EPS:
            scale = min(scale, O2 / O2_req)

        scale = max(min(scale, 1.0), 0.0)

        dNO2 = dNO2_target * scale
        dNO  = dNO_target  * scale

        # Consume reactants
        out.set(self.no2_name, NO2 - dNO2)
        out.set(self.no_name,  NO  - dNO)

        out.set(self.nh3_name, NH3 - (2.0 * dNO2 + 1.0 * dNO))
        out.set(self.o2_name,  O2  - (0.5 * dNO2 + 0.25 * dNO))

        # Produce products
        out.set(self.n2_name,  out.get(self.n2_name) + (1.5 * dNO2 + 1.0 * dNO))
        out.set(self.h2o_name, out.get(self.h2o_name) + (3.0 * dNO2 + 1.5 * dNO))

        if self.print_diagnostics:
            print(f"[{self.name}] removed NO2={dNO2:.6g} mol/s, NO={dNO:.6g} mol/s (scale={scale:.4g})")
