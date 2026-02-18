from __future__ import annotations

from typing import Optional

from ..flowsheet_tools import EPS, Stream, UnitOp


class IdealTSAColumnEMM17(UnitOp):
    """
    Ideal temperature swing adsorption (TSA) column model for iodine on EMM-17.

    Two operating modes:
      - mode="adsorb":
          * removes iodine from process gas (perfect capture until working capacity is exceeded)
          * outputs cleaned gas (gas_out)

      - mode="regen":
          * desorbs iodine (time-averaged) into a regeneration gas (air)
          * outputs iodine-rich gas (regen_out)
    """

    def __init__(
        self,
        name: str,
        *,
        mode: str,  # "adsorb" or "regen"
        iodine_name: str = "I_g",
        bed_volume_m3: float = 1.0,
        cap_work_mol_per_m3: float = 0.0,
        t_ads_s: float = 8.0 * 3600.0,
        t_regen_s: float = 4.0 * 3600.0,
        regen_yI2_max: float = 0.02,
        print_diagnostics: bool = False,
    ):
        super().__init__(name)

        self.mode = str(mode).lower().strip()
        if self.mode not in ("adsorb", "regen"):
            raise ValueError(f"{name}: mode must be 'adsorb' or 'regen'")

        self.iodine_name = str(iodine_name)
        self.bed_volume_m3 = float(bed_volume_m3)
        self.cap_work_mol_per_m3 = float(cap_work_mol_per_m3)
        self.t_ads_s = float(t_ads_s)
        self.t_regen_s = float(t_regen_s)
        self.regen_yI2_max = float(regen_yI2_max)
        self.print_diagnostics = bool(print_diagnostics)

        # Optional coupling: in regen mode, pull captured inventory from this source column
        self.regen_source: Optional["IdealTSAColumnEMM17"] = None

        # --- calculated / reported (rates are mol/s, inventories are mol) ---
        self.I_in_mol_s: float = 0.0          # iodine into adsorber
        self.I_removed_mol_s: float = 0.0     # iodine captured during adsorption (average rate)
        self.I_slip_mol_s: float = 0.0        # iodine that breaks through during adsorption
        self.I_desorb_mol_s: float = 0.0      # iodine released during regen (average rate)

        self.capacity_cycle_mol: float = 0.0  # working capacity available this cycle (mol)
        self.loading_cycle_mol: float = 0.0   # iodine fed over adsorption half-cycle (mol)
        self.captured_cycle_mol: float = 0.0  # iodine actually captured this cycle (mol)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def apply(self) -> None:
        if self.mode == "adsorb":
            self._apply_adsorb()
        else:
            self._apply_regen()

    def set_regen_source(self, src: "IdealTSAColumnEMM17") -> None:
        """Convenience helper for A/B beds: regen bed reads captured iodine from adsorbing bed."""
        self.regen_source = src

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------
    def _copy_stream(self, src: Stream, dst: Stream) -> None:
        dst.mol = {}
        for sp, n in src.mol.items():
            dst.set(sp, n)

        # If your Stream class carries T/p/phase separately and you want to
        # preserve them, uncomment:
        # dst.T = src.T
        # dst.p = src.p
        # dst.phase = src.phase

    def _apply_adsorb(self) -> None:
        gas_in = self.inlets["gas_in"]
        gas_out = self.outlets["gas_out"]

        # Pass-through everything, then fix iodine
        self._copy_stream(gas_in, gas_out)

        I_in = float(gas_in.get(self.iodine_name))
        self.I_in_mol_s = I_in

        t_ads = max(self.t_ads_s, 0.0)
        t_reg = max(self.t_regen_s, EPS)

        # Working capacity available in this adsorption half-cycle (mol)
        cap_cycle = max(self.cap_work_mol_per_m3, 0.0) * max(self.bed_volume_m3, 0.0)
        load_cycle = I_in * t_ads

        self.capacity_cycle_mol = float(cap_cycle)
        self.loading_cycle_mol = float(load_cycle)

        if cap_cycle <= EPS or t_ads <= EPS or I_in <= EPS:
            captured_cycle = 0.0
        else:
            # Can only capture up to the available working capacity.
            captured_cycle = min(load_cycle, cap_cycle)

        # Convert captured per-cycle inventory back to an average capture rate (mol/s)
        I_removed = captured_cycle / max(t_ads, EPS)
        I_slip = max(I_in - I_removed, 0.0)

        # Update outlet iodine
        gas_out.set(self.iodine_name, I_slip)

        # Store results
        self.captured_cycle_mol = float(captured_cycle)
        self.I_removed_mol_s = float(I_removed)
        self.I_slip_mol_s = float(I_slip)

        # Time-averaged regen release rate (mol/s) for THIS cycle
        self.I_desorb_mol_s = float(captured_cycle / t_reg)

        if self.print_diagnostics:
            print(
                f"[{self.name}] ADSORB: I_in={I_in:.6g} mol/s, "
                f"I_removed={I_removed:.6g}, I_slip={I_slip:.6g}, "
                f"cap_cycle={cap_cycle:.6g} mol, load_cycle={load_cycle:.6g} mol, "
                f"captured_cycle={captured_cycle:.6g} mol, "
                f"I_desorb_avg={self.I_desorb_mol_s:.6g} mol/s"
            )

    def _apply_regen(self) -> None:
        regen_in = self.inlets["regen_in"]
        regen_out = self.outlets["regen_out"]

        # Regen inlet should be air-only; copy it to outlet then add desorbed iodine
        self._copy_stream(regen_in, regen_out)

        # Determine desorb rate:
        #  1) If coupled to an adsorbing bed, use its captured inventory for the cycle.
        #  2) Else fall back to this object's I_desorb_mol_s (may be set externally).
        I_rel = 0.0

        if self.regen_source is not None:
            # Prefer captured inventory if available, else fall back to its stored desorb rate
            captured = float(getattr(self.regen_source, "captured_cycle_mol", 0.0))
            if captured > EPS:
                I_rel = captured / max(self.t_regen_s, EPS)
            else:
                I_rel = float(getattr(self.regen_source, "I_desorb_mol_s", 0.0))
        else:
            I_rel = float(getattr(self, "I_desorb_mol_s", 0.0))

        I_rel = max(I_rel, 0.0)
        self.I_desorb_mol_s = float(I_rel)

        # Add iodine to regen outlet (air + iodine)
        if I_rel > EPS:
            regen_out.set(self.iodine_name, regen_out.get(self.iodine_name) + I_rel)

        if self.print_diagnostics:
            n_tot = regen_out.total_molar_flow()
            y = I_rel / max(n_tot, EPS)
            print(
                f"[{self.name}] REGEN: I_out={I_rel:.6g} mol/s, "
                f"total={n_tot:.6g} mol/s, y_I={y:.3e} "
                f"(limit={self.regen_yI2_max:.3e})"
            )
