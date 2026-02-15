from __future__ import annotations
from ..flowsheet_tools import UnitOp, EPS

class AmmoniaAbsorptionColumn(UnitOp):
    """
    Simple NH3 absorber:
      - Transfers capture_frac of NH3 from gas_in to liq_in.
      - Everything else passes through unchanged (no chemistry).
    """
    def __init__(
        self,
        name: str,
        *,
        capture_frac: float = 0.99,
        nh3_name: str = "NH3",
        print_diagnostics: bool = False,
    ):
        super().__init__(name)
        if not (0.0 <= capture_frac <= 1.0):
            raise ValueError("capture_frac must be in [0,1]")
        self.capture_frac = float(capture_frac)
        self.nh3_name = str(nh3_name)
        self.print_diagnostics = bool(print_diagnostics)

    def apply(self) -> None:
        gas_in = self.inlets["gas_in"]
        liq_in = self.inlets["liq_in"]
        liq_out = self.outlets["liq_out"]
        gas_out = self.outlets["gas_out"]

        # start by passing both streams through
        liq_out.mol = {}
        gas_out.mol = {}

        # copy liquid inlet
        for sp, n in liq_in.mol.items():
            if n > EPS:
                liq_out.set(sp, float(n))

        # copy gas inlet
        for sp, n in gas_in.mol.items():
            if n > EPS:
                gas_out.set(sp, float(n))

        # transfer NH3
        n_nh3 = float(gas_out.get(self.nh3_name))
        d = self.capture_frac * n_nh3

        if d > EPS:
            gas_out.set(self.nh3_name, n_nh3 - d)
            liq_out.set(self.nh3_name, liq_out.get(self.nh3_name) + d)

        # phases
        liq_out.phase = "L"
        gas_out.phase = "G"

        if self.print_diagnostics:
            print(f"[{self.name}] captured NH3={d:.6g} mol/s (frac={self.capture_frac})")
