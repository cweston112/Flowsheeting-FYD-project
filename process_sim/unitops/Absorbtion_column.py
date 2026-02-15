from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from ..flowsheet_tools import UnitOp, EPS


class NO2AbsorptionColumn(UnitOp):
    """
    Simple NOx absorber "book-keeping" model.

    Inlets:
      - gas_in (bottom): NO2/NO/O2/N2/H2O etc.
      - liq_in (top): water (and possibly some acid)

    Outlets:
      - liq_out: dilute nitric acid
      - gas_out: remaining gas

    Chemistry handled:
      (A) NO2 absorption via overall: 3 NO2 + H2O -> 2 HNO3 + NO
          Implemented as: capture_frac_NO2 of NO2 is removed from gas.
          For each 1 mol NO2 captured: consumes (1/3) mol H2O, produces (2/3) mol HNO3,
          and produces (1/3) mol NO back to gas.

      (B) NO conversion (oxidation+absorption) overall:
          4 NO + 3 O2 + 2 H2O -> 4 HNO3
          Implemented as: capture_frac_NO of NO is converted to HNO3, limited by available O2 and H2O.
          For each 1 mol NO converted: consumes 0.75 mol O2 and 0.5 mol H2O, produces 1 mol HNO3.

    Notes:
      - This is not rate/HTU/NTU-based; it is a stoichiometric splitter for flowsheet closure.
      - If O2 or H2O is insufficient, NO conversion is reduced accordingly.
    """

    def __init__(
        self,
        name: str,
        *,
        no2_name: str = "NO2",
        no_name: str = "NO",
        o2_name: str = "O2",
        h2o_name: str = "H2O",
        hno3_name: str = "HNO3",
        capture_frac: float = 0.99,          # legacy: NO2 capture fraction
        capture_frac_no: float = 0.95,       # NEW: NO conversion fraction
        print_diagnostics: bool = False,
    ):
        super().__init__(name)
        self.no2_name = no2_name
        self.no_name = no_name
        self.o2_name = o2_name
        self.h2o_name = h2o_name
        self.hno3_name = hno3_name
        self.capture_frac = float(capture_frac)
        self.capture_frac_no = float(capture_frac_no)
        self.print_diagnostics = bool(print_diagnostics)

        if not (0.0 <= self.capture_frac <= 1.0):
            raise ValueError(f"{name}: capture_frac must be in [0,1]")
        if not (0.0 <= self.capture_frac_no <= 1.0):
            raise ValueError(f"{name}: capture_frac_no must be in [0,1]")

        # For reporting
        self.last_captured_no2 = 0.0
        self.last_converted_no = 0.0
        self.last_hno3_made = 0.0

    def apply(self) -> None:
        gas_in = self.inlets["gas_in"]
        liq_in = self.inlets["liq_in"]
        liq_out = self.outlets["liq_out"]
        gas_out = self.outlets["gas_out"]

        # Start as passthrough copies
        gas_out.mol = dict(gas_in.mol)
        liq_out.mol = dict(liq_in.mol)

        n_no2 = gas_in.get(self.no2_name)
        n_no  = gas_in.get(self.no_name)
        n_o2  = gas_in.get(self.o2_name)

        # Available water is in liquid inlet primarily, but allow any in gas too
        n_h2o_avail = liq_in.get(self.h2o_name) + gas_in.get(self.h2o_name)

        # ---------------------------------------------------------------------
        # (A) NO2 capture: remove fraction of NO2 from gas and turn into HNO3
        #     3 NO2 + H2O -> 2 HNO3 + NO
        # ---------------------------------------------------------------------
        cap_no2_target = self.capture_frac * max(n_no2, 0.0)

        # water required per 1 NO2 captured = 1/3
        h2o_req_no2 = (1.0 / 3.0) * cap_no2_target
        if h2o_req_no2 > n_h2o_avail + EPS:
            # limit capture by water availability
            cap_no2 = max(0.0, 3.0 * n_h2o_avail)
        else:
            cap_no2 = cap_no2_target

        # apply capture
        if cap_no2 > EPS:
            gas_out.set(self.no2_name, max(n_no2 - cap_no2, 0.0))

            # consume water, produce acid, produce NO (back to gas)
            h2o_used = (1.0 / 3.0) * cap_no2
            hno3_made = (2.0 / 3.0) * cap_no2
            no_made = (1.0 / 3.0) * cap_no2

            # allocate water consumption from liquid first, then gas
            n_h2o_liq = liq_out.get(self.h2o_name)
            take_liq = min(n_h2o_liq, h2o_used)
            liq_out.set(self.h2o_name, max(n_h2o_liq - take_liq, 0.0))
            rem = h2o_used - take_liq
            if rem > EPS:
                n_h2o_gas = gas_out.get(self.h2o_name)
                gas_out.set(self.h2o_name, max(n_h2o_gas - rem, 0.0))

            liq_out.set(self.hno3_name, liq_out.get(self.hno3_name) + hno3_made)
            gas_out.set(self.no_name, gas_out.get(self.no_name) + no_made)

            n_h2o_avail = n_h2o_avail - h2o_used  # update available pool
        else:
            hno3_made = 0.0

        # Update NO amount after NO2 step (includes NO produced)
        n_no_after = gas_out.get(self.no_name)
        n_o2_after = gas_out.get(self.o2_name)
        n_h2o_avail = (liq_out.get(self.h2o_name) + gas_out.get(self.h2o_name))

        # ---------------------------------------------------------------------
        # (B) NO conversion to HNO3:
        #     4 NO + 3 O2 + 2 H2O -> 4 HNO3
        #   per 1 NO: O2 0.75, H2O 0.5, HNO3 +1
        # ---------------------------------------------------------------------
        conv_no_target = self.capture_frac_no * max(n_no_after, 0.0)

        # limit by O2 and H2O availability
        if conv_no_target > EPS:
            max_by_o2 = (n_o2_after / 0.75) if n_o2_after > EPS else 0.0
            max_by_h2o = (n_h2o_avail / 0.5) if n_h2o_avail > EPS else 0.0
            conv_no = max(0.0, min(conv_no_target, max_by_o2, max_by_h2o))
        else:
            conv_no = 0.0

        if conv_no > EPS:
            # consume NO, O2, H2O; produce HNO3
            gas_out.set(self.no_name, max(n_no_after - conv_no, 0.0))

            o2_used = 0.75 * conv_no
            gas_out.set(self.o2_name, max(n_o2_after - o2_used, 0.0))

            h2o_used = 0.5 * conv_no
            # take from liquid first, then gas
            n_h2o_liq = liq_out.get(self.h2o_name)
            take_liq = min(n_h2o_liq, h2o_used)
            liq_out.set(self.h2o_name, max(n_h2o_liq - take_liq, 0.0))
            rem = h2o_used - take_liq
            if rem > EPS:
                n_h2o_gas = gas_out.get(self.h2o_name)
                gas_out.set(self.h2o_name, max(n_h2o_gas - rem, 0.0))

            liq_out.set(self.hno3_name, liq_out.get(self.hno3_name) + conv_no)

        # Phase/T/p stamping (keep whatever conditioning stamps later too)
        liq_out.phase = "L"
        gas_out.phase = "G"
        liq_out.T, liq_out.p = gas_in.T, gas_in.p
        gas_out.T, gas_out.p = gas_in.T, gas_in.p

        # Reporting
        self.last_captured_no2 = float(cap_no2)
        self.last_converted_no = float(conv_no)
        self.last_hno3_made = float((2.0 / 3.0) * cap_no2 + conv_no)

        if self.print_diagnostics:
            print(f"[{self.name}] captured NO2={self.last_captured_no2:.6g} mol/s, "
                  f"converted NO={self.last_converted_no:.6g} mol/s, "
                  f"HNO3 made={self.last_hno3_made:.6g} mol/s")
