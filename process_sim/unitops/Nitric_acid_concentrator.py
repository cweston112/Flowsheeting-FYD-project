from __future__ import annotations

from dataclasses import dataclass

from ..flowsheet_tools import UnitOp, EPS


class NitricAcidConcentrator(UnitOp):
    """
    Concentrate aqueous HNO3 by evaporating water only (no HNO3 loss).
    Target is implemented as mass fraction of HNO3 in liquid product, e.g. 0.45 ~ "45%".

    Inlet:  "in"
    Outlets:
      - "vap" : evaporated water (steam)
      - "liq" : concentrated acid

    Notes:
      - If you truly need *volume percent*, you’ll want a density model; this uses mass fraction.
      - All non-(H2O,HNO3) species are kept in the liquid.
    """

    def __init__(
        self,
        name: str,
        *,
        target_hno3_mass_frac: float = 0.45,
        hno3_name: str = "HNO3",
        h2o_name: str = "H2O",
        print_diagnostics: bool = False,
    ):
        super().__init__(name)
        self.target = float(target_hno3_mass_frac)
        if not (0.0 < self.target < 1.0):
            raise ValueError(f"{name}: target_hno3_mass_frac must be in (0,1)")
        self.hno3_name = hno3_name
        self.h2o_name = h2o_name
        self.print_diagnostics = bool(print_diagnostics)

        self.last_evap_h2o = 0.0

    def apply(self) -> None:
        sin = self.inlets["in"]
        svap = self.outlets["vap"]
        sliq = self.outlets["liq"]

        # Start: copy all to liquid; vapour empty
        sliq.mol = dict(sin.mol)
        svap.mol = {}

        reg = sin.reg
        mw = reg.mw  # expects dict-like: mw["H2O"], mw["HNO3"], ...

        n_hno3 = sin.get(self.hno3_name)
        n_h2o = sin.get(self.h2o_name)

        if n_hno3 <= EPS:
            # nothing to concentrate; just pass through, no evaporation
            svap.set(self.h2o_name, 0.0)
            sliq.set(self.h2o_name, n_h2o)
            sliq.phase = "L"
            svap.phase = "G"
            sliq.T, sliq.p = sin.T, sin.p
            svap.T, svap.p = sin.T, sin.p
            self.last_evap_h2o = 0.0
            return

        m_hno3 = n_hno3 * mw[self.hno3_name]
        # target mass fraction: m_hno3 / (m_hno3 + m_h2o_out) = target
        # => m_h2o_out = m_hno3*(1-target)/target
        m_h2o_req = m_hno3 * (1.0 - self.target) / self.target
        n_h2o_req = m_h2o_req / mw[self.h2o_name]

        # Evaporate excess water only
        n_h2o_out = min(n_h2o, max(n_h2o_req, 0.0))
        n_evap = max(n_h2o - n_h2o_out, 0.0)

        sliq.set(self.h2o_name, n_h2o_out)
        svap.set(self.h2o_name, n_evap)

        # Do NOT evaporate HNO3
        sliq.set(self.hno3_name, n_hno3)
        svap.set(self.hno3_name, 0.0)

        # Stamp conditions
        sliq.phase = "L"
        svap.phase = "G"
        sliq.T, sliq.p = sin.T, sin.p
        svap.T, svap.p = sin.T, sin.p

        self.last_evap_h2o = float(n_evap)

        if self.print_diagnostics:
            # compute resulting mass fraction
            m_h2o_out = n_h2o_out * mw[self.h2o_name]
            w = m_hno3 / max(m_hno3 + m_h2o_out, EPS)
            print(f"[{self.name}] evap_H2O={n_evap:.6g} mol/s -> w_HNO3={w:.3f}")
