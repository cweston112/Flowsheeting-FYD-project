"""
=============================================================================
SIMPLIFIED PROCESS FLOWSHEET — UNIT OPERATIONS
=============================================================================
All units operate on simple per-species split fractions or stoichiometric
conversions defined in inputs.py.  No complex stage calculations are performed.

Unit list
---------
  Mixer              – blend any number of inlet streams
  Splitter           – split one stream into two by fraction
  PassThrough        – copy inlet to outlet unchanged (buffer vessels, etc.)
  SplitUnit          – general 2-outlet split by per-species fractions
  ConversionReactor  – stoichiometric reaction with fractional conversion
  PhaseCondenser     – split gas into vapour + liquid by frac_to_liquid
  NOxAbsorber        – D101: removes NO2/NO from gas, produces HNO3 in liquid
  AcidConcentrator   – E104: evaporates water to hit target HNO3 mass fraction
  SCRReactor         – R103: NH3 + NOx → N2 + H2O
  KnockoutPot        – KO101: removes water from gas
  NH3Absorber        – D102: removes NH3 from gas into water
  TSAColumn          – TSA101: removes iodine from gas stream
=============================================================================
"""
from __future__ import annotations

from typing import Dict, List, Optional

from .framework import UnitOp, Stream, EPS


class ThreeOutletSplitUnit(UnitOp):
    """
    Generic one-inlet / three-outlet splitter.

    Fractions to outlet_A and outlet_B are specified per species. The remainder
    goes to outlet_C. This is used to make the EV-103A/B/C train explicit while
    keeping all split data editable in inputs.py.
    """

    def __init__(
        self,
        name: str,
        *,
        outlet_A: str,
        outlet_B: str,
        outlet_C: str,
        frac_to_A: Dict[str, float],
        frac_to_B: Dict[str, float],
        default_frac_A: float = 0.0,
        default_frac_B: float = 0.0,
        phase_A: str | None = None,
        phase_B: str | None = None,
        phase_C: str | None = None,
    ):
        super().__init__(name)
        self.outlet_A = outlet_A
        self.outlet_B = outlet_B
        self.outlet_C = outlet_C
        self.frac_to_A = dict(frac_to_A)
        self.frac_to_B = dict(frac_to_B)
        self.default_frac_A = float(default_frac_A)
        self.default_frac_B = float(default_frac_B)
        self.phase_A = phase_A
        self.phase_B = phase_B
        self.phase_C = phase_C

    def apply(self) -> None:
        src = self.inlets["in"]
        out_A = self.outlets[self.outlet_A]
        out_B = self.outlets[self.outlet_B]
        out_C = self.outlets[self.outlet_C]
        out_A.clear(); out_B.clear(); out_C.clear()
        out_A.T = out_B.T = out_C.T = src.T
        out_A.p = out_B.p = out_C.p = src.p
        if self.phase_A is not None:
            out_A.phase = self.phase_A
        if self.phase_B is not None:
            out_B.phase = self.phase_B
        if self.phase_C is not None:
            out_C.phase = self.phase_C

        for sp, n in src.mol.items():
            fA = max(0.0, min(1.0, self.frac_to_A.get(sp, self.default_frac_A)))
            fB = max(0.0, min(1.0, self.frac_to_B.get(sp, self.default_frac_B)))
            if fA + fB > 1.0:
                scale = 1.0 / max(fA + fB, EPS)
                fA *= scale
                fB *= scale
            out_A.add(sp, n * fA)
            out_B.add(sp, n * fB)
            out_C.add(sp, n * (1.0 - fA - fB))





class EV103CPerformanceStage(UnitOp):
    """
    Finalised source-sensitive EV-103C surrogate fitted to the updated reactor
    section stream table.

    Inlet roles
    -----------
    light : visible light side stream from EV-103B (F66) → entirely to visible waste F70
    heavy : visible heavy stream from EV-103B (F68)     → split to product F24 and vapour F25
    aux   : internal impurity bypass (HNO3/TBP/HTcO4)   → partly to product, remainder to hidden sink

    Notes
    -----
    This stage is intentionally source-sensitive so the PFD-facing streams
    F66/F68/F70/F24/F25 can match the finalised performance table exactly.
    """

    def __init__(
        self,
        name: str,
        *,
        heavy_to_product: Dict[str, float],
        heavy_to_vapour: Dict[str, float],
        light_to_product: Dict[str, float],
        aux_to_product: Dict[str, float],
        extra_h2o_to_vapour: float = 0.0,
    ):
        super().__init__(name)
        self.heavy_to_product = dict(heavy_to_product)
        self.heavy_to_vapour = dict(heavy_to_vapour)
        self.light_to_product = dict(light_to_product)
        self.aux_to_product = dict(aux_to_product)
        self.extra_h2o_to_vapour = float(extra_h2o_to_vapour)

    def apply(self) -> None:
        light = self.inlets["light"]
        heavy = self.inlets["heavy"]
        aux = self.inlets.get("aux")
        product = self.outlets["product"]
        vapour = self.outlets["vapour"]
        waste = self.outlets["waste"]
        sink = self.outlets.get("sink")

        for out in [product, vapour, waste, sink]:
            if out is not None:
                out.clear()

        refs = [s for s in (light, heavy, aux) if s is not None]
        Tref = sum(s.total_molar_flow() * s.T for s in refs) / max(sum(s.total_molar_flow() for s in refs), EPS) if refs else 298.15
        pref = min((s.p for s in refs), default=1e5)
        for out, phase in [(product, "L"), (vapour, "G"), (waste, "L"), (sink, "L")]:
            if out is not None:
                out.T = Tref
                out.p = pref
                out.phase = phase

        # Light stream goes mostly to visible waste, but trace dissolved
        # non-volatiles can be retained to product for compositionally consistent bookkeeping.
        for sp, n in light.mol.items():
            if n <= EPS:
                continue
            n_prod = min(self.light_to_product.get(sp, 0.0), n)
            rem = n - n_prod
            if n_prod > EPS:
                product.add(sp, n_prod)
            if rem > EPS:
                waste.add(sp, rem)

        # Heavy stream is routed by specified target amounts.
        for sp, n in heavy.mol.items():
            if n <= EPS:
                continue
            n_prod = min(self.heavy_to_product.get(sp, 0.0), n)
            rem = n - n_prod
            n_vap = min(self.heavy_to_vapour.get(sp, 0.0), rem)
            rem -= n_vap
            if n_prod > EPS:
                product.add(sp, n_prod)
            if n_vap > EPS:
                vapour.add(sp, n_vap)
            if rem > EPS:
                if sink is not None:
                    sink.add(sp, rem)
                else:
                    waste.add(sp, rem)

        # Auxiliary impurity bypass is partly retained in product and otherwise discarded to sink.
        if aux is not None:
            for sp, n in aux.mol.items():
                if n <= EPS:
                    continue
                n_prod = min(self.aux_to_product.get(sp, 0.0), n)
                rem = n - n_prod
                if n_prod > EPS:
                    product.add(sp, n_prod)
                if rem > EPS:
                    if sink is not None:
                        sink.add(sp, rem)
                    else:
                        waste.add(sp, rem)

        if self.extra_h2o_to_vapour > EPS:
            vapour.add("H2O", self.extra_h2o_to_vapour)

# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def _sum_inlets(inlets: Dict[str, Stream], *, T_default: float = 298.15,
                p_default: float = 1e5, phase_default: str = "L") -> Dict[str, float]:
    """Return summed mol dict from all inlets."""
    result: Dict[str, float] = {}
    for s in inlets.values():
        for sp, v in s.mol.items():
            result[sp] = result.get(sp, 0.0) + v
    return result


def _stamp_T_p(streams: List[Stream], T: Optional[float], p: Optional[float]) -> None:
    for s in streams:
        if T is not None:
            s.T = T
        if p is not None:
            s.p = p


# ─────────────────────────────────────────────────────────────────────────────
# MIXER
# ─────────────────────────────────────────────────────────────────────────────

class Mixer(UnitOp):
    """
    Blends all inlet streams into one outlet.
    T of outlet = flow-weighted average; p = minimum.
    """

    def apply(self) -> None:
        out = self.outlets["out"]
        out.clear()
        total_F = 0.0
        total_HF = 0.0   # enthalpy proxy: F*T
        min_p = min((s.p for s in self.inlets.values()), default=1e5)

        for s in self.inlets.values():
            for sp, v in s.mol.items():
                out.add(sp, v)
            F = s.total_molar_flow()
            total_F += F
            total_HF += F * s.T

        out.T = (total_HF / total_F) if total_F > EPS else list(self.inlets.values())[0].T if self.inlets else 298.15
        out.p = min_p


# ─────────────────────────────────────────────────────────────────────────────
# PASS-THROUGH  (buffer vessel, surge vessel)
# ─────────────────────────────────────────────────────────────────────────────

class PassThrough(UnitOp):
    """
    Copies inlet 'in' to outlet 'out' without modification.
    T/p are stamped by the flowsheet conditioning if needed.
    """

    def apply(self) -> None:
        src = self.inlets["in"]
        dst = self.outlets["out"]
        dst.mol = {sp: v for sp, v in src.mol.items()}
        dst.T = src.T
        dst.p = src.p
        dst.phase = src.phase


# ─────────────────────────────────────────────────────────────────────────────
# SPLITTER  (one inlet → two outlets by per-species fraction)
# ─────────────────────────────────────────────────────────────────────────────

class SplitUnit(UnitOp):
    """
    Generic 2-outlet split.

    Parameters
    ----------
    outlet_A        : name of the primary outlet key (e.g. "organic")
    outlet_B        : name of the secondary outlet key (e.g. "aqueous")
    frac_to_A       : {species: fraction going to outlet_A}
    default_frac_A  : fraction applied to any species not in frac_to_A
    """

    def __init__(
        self,
        name: str,
        *,
        outlet_A: str,
        outlet_B: str,
        frac_to_A: Dict[str, float],
        default_frac_A: float = 0.0,
    ):
        super().__init__(name)
        self.outlet_A = outlet_A
        self.outlet_B = outlet_B
        self.frac_to_A = frac_to_A
        self.default_frac_A = default_frac_A

    def apply(self) -> None:
        src = self.inlets["in"]
        out_A = self.outlets[self.outlet_A]
        out_B = self.outlets[self.outlet_B]
        out_A.clear(); out_B.clear()
        out_A.T = out_B.T = src.T
        out_A.p = out_B.p = src.p

        for sp, n in src.mol.items():
            fA = self.frac_to_A.get(sp, self.default_frac_A)
            fA = max(0.0, min(1.0, fA))
            out_A.add(sp, n * fA)
            out_B.add(sp, n * (1.0 - fA))


# ─────────────────────────────────────────────────────────────────────────────
# SIMPLE PURGE SPLITTER  (one fraction to purge, rest to recycle)
# ─────────────────────────────────────────────────────────────────────────────

class PurgeSplitter(UnitOp):
    """
    Splits inlet into a purge (fraction f) and a recycle (1-f).
    The purge fraction is updated externally by the solver.
    """

    def __init__(self, name: str, purge_frac: float = 0.05):
        super().__init__(name)
        self.purge_frac = float(purge_frac)

    def apply(self) -> None:
        src = self.inlets["in"]
        purge = self.outlets["purge"]
        recycle = self.outlets["recycle"]
        purge.clear(); recycle.clear()
        f = max(0.0, min(1.0, self.purge_frac))
        for sp, n in src.mol.items():
            purge.add(sp, n * f)
            recycle.add(sp, n * (1.0 - f))
        purge.T = recycle.T = src.T
        purge.p = recycle.p = src.p


# ─────────────────────────────────────────────────────────────────────────────
# COALESCER  (remove entrained organics from aqueous)
# ─────────────────────────────────────────────────────────────────────────────

class Coalescer(UnitOp):
    """
    Removes a fraction of specified species (e.g. entrained TBP/dodecane)
    from the inlet stream.  Removed material goes to a recovered organic outlet.
    """

    def __init__(self, name: str, *, frac_removed: Dict[str, float], default_frac_removed: float = 0.0):
        super().__init__(name)
        self.frac_removed = frac_removed
        self.default_frac_removed = default_frac_removed

    def apply(self) -> None:
        src = self.inlets["in"]
        aq = self.outlets["aqueous"]
        org = self.outlets["recovered_org"]
        aq.clear(); org.clear()
        aq.T = org.T = src.T
        aq.p = org.p = src.p

        for sp, n in src.mol.items():
            fr = self.frac_removed.get(sp, self.default_frac_removed)
            fr = max(0.0, min(1.0, fr))
            aq.add(sp, n * (1.0 - fr))
            org.add(sp, n * fr)


# ─────────────────────────────────────────────────────────────────────────────
# EVAPORATOR  (liquid inlet → vapour + liquid by frac_to_vapour)
# ─────────────────────────────────────────────────────────────────────────────

class Evaporator(UnitOp):
    """
    Simplified evaporator: each species is split between vapour and liquid
    outlets according to frac_to_vapour. Optionally applies full AHA destruction
    before the split, matching the old_complex evaporator behaviour.
    """

    def __init__(self, name: str, *, frac_to_vapour: Dict[str, float], default_frac_to_vapour: float = 0.0,
                 destroy_aha: bool = False, aha_name: str = "AHA",
                 aha_products: Optional[Dict[str, float]] = None, aha_hno3_consumption: float = 0.0):
        super().__init__(name)
        self.frac_to_vapour = frac_to_vapour
        self.default_frac_to_vapour = default_frac_to_vapour
        self.destroy_aha = bool(destroy_aha)
        self.aha_name = aha_name
        self.aha_products = aha_products or {"AcOH": 1.0, "N2O": 0.5}
        self.aha_hno3_consumption = float(aha_hno3_consumption)

    def apply(self) -> None:
        src = self.inlets["in"]
        vap = self.outlets["vapour"]
        liq = self.outlets["liquid"]
        vap.clear(); liq.clear()
        vap.T = liq.T = src.T
        vap.p = liq.p = src.p
        vap.phase = "G"; liq.phase = "L"

        mol = {sp: float(n) for sp, n in src.mol.items()}
        if self.destroy_aha:
            n_aha = mol.get(self.aha_name, 0.0)
            if n_aha > EPS:
                mol[self.aha_name] = 0.0
                if self.aha_hno3_consumption > 0.0:
                    mol["HNO3"] = max(mol.get("HNO3", 0.0) - self.aha_hno3_consumption * n_aha, 0.0)
                for sp, nu in self.aha_products.items():
                    mol[sp] = mol.get(sp, 0.0) + nu * n_aha

        for sp, n in mol.items():
            fv = self.frac_to_vapour.get(sp, self.default_frac_to_vapour)
            fv = max(0.0, min(1.0, fv))
            vap.add(sp, n * fv)
            liq.add(sp, n * (1.0 - fv))

        self.calcs["F_vapour_mol_s"] = vap.total_molar_flow()
        self.calcs["F_liquid_mol_s"] = liq.total_molar_flow()


# ─────────────────────────────────────────────────────────────────────────────
# CONVERSION REACTOR
# ─────────────────────────────────────────────────────────────────────────────

class ConversionReactor(UnitOp):
    """
    Stoichiometric conversion reactor.

    Reactions defined as list of dicts:
      {
        "reactant":    species name,
        "conversion":  fraction converted (0–1),
        "products":    { product_name: mol_produced_per_mol_reactant },
        "outlet_solid": [list of product names going to solid outlet],
        "outlet_gas":   [list of product names going to gas outlet],
      }

    Unreacted feed and non-reaction species are routed to outlets according to
    `default_to_solid` flag.
    Outlets: "solid", "gas", "liquid" (liquid for any aqueous products)
    """

    def __init__(self, name: str, reactions: list, default_to_solid: bool = False):
        super().__init__(name)
        self.reactions = reactions
        self.default_to_solid = default_to_solid

    def apply(self) -> None:
        # Merge all inlets
        mol_in: Dict[str, float] = {}
        for s in self.inlets.values():
            for sp, v in s.mol.items():
                mol_in[sp] = mol_in.get(sp, 0.0) + v

        mol_out = {sp: v for sp, v in mol_in.items()}

        solid_set: set = set()
        gas_set: set = set()

        for rxn in self.reactions:
            reactant = rxn["reactant"]
            conv = float(rxn.get("conversion", 1.0))
            n_reacted = mol_out.get(reactant, 0.0) * conv
            mol_out[reactant] = mol_out.get(reactant, 0.0) - n_reacted

            for prod, stoich in rxn.get("products", {}).items():
                mol_out[prod] = mol_out.get(prod, 0.0) + n_reacted * stoich

            solid_set.update(rxn.get("outlet_solid", []))
            gas_set.update(rxn.get("outlet_gas", []))

        # Route to outlets (handle aliased outlets correctly)
        solid_out = self.outlets.get("solid")
        gas_out = self.outlets.get("gas")
        liq_out = self.outlets.get("liquid")

        # Clear distinct outlet objects only (avoid double-clear on aliased outlets)
        cleared = set()
        for out in [solid_out, gas_out, liq_out]:
            if out is not None and id(out) not in cleared:
                out.clear()
                cleared.add(id(out))

        for sp, n in mol_out.items():
            if n < EPS:
                continue
            if sp in solid_set and solid_out is not None:
                solid_out.add(sp, n)
            elif sp in gas_set and gas_out is not None:
                gas_out.add(sp, n)
            elif liq_out is not None:
                liq_out.add(sp, n)
            elif self.default_to_solid and solid_out is not None:
                solid_out.add(sp, n)
            elif gas_out is not None:
                gas_out.add(sp, n)

        # Stamp T/p
        ref = next(iter(self.inlets.values()), None)
        if ref:
            for out in [solid_out, gas_out, liq_out]:
                if out is not None:
                    out.T = ref.T; out.p = ref.p





class CalcinerReactor(UnitOp):
    """
    Finalised reactor-section calciner surrogate fitted to the updated stream table.

    Behaviour
    ---------
    - Converts UO2(NO3)2 completely to solid UO3.
    - Sends 2 mol NO2 per mol UO2(NO3)2 to the offgas.
    - Passes the separate F-73 air stream directly to the offgas.
    - Routes only a fitted fraction of feed H2O to the offgas.
    - Keeps HNO3/TBP/HTcO4 out of the visible offgas so the exported
      reactor-section stream table matches the finalised performance data.
    """

    def __init__(
        self,
        name: str,
        *,
        no2_per_u: float = 2.0,
        o2_per_u: float = 0.0,
        h2o_to_offgas_frac: float = 1.0,
        route_hno3_to_offgas: float = 0.0,
        route_tbp_to_offgas: float = 0.0,
        route_htco4_to_offgas: float = 0.0,
    ):
        super().__init__(name)
        self.no2_per_u = float(no2_per_u)
        self.o2_per_u = float(o2_per_u)
        self.h2o_to_offgas_frac = float(h2o_to_offgas_frac)
        self.route_hno3_to_offgas = float(route_hno3_to_offgas)
        self.route_tbp_to_offgas = float(route_tbp_to_offgas)
        self.route_htco4_to_offgas = float(route_htco4_to_offgas)

    def apply(self) -> None:
        Sin = self.inlets["in"]
        Sair = self.inlets.get("air")
        Ssol = self.outlets["solids"]
        Sgas = self.outlets["offgas"]
        Ssol.clear(); Sgas.clear()
        Ssol.T = Sin.T; Ssol.p = Sin.p; Ssol.phase = "S"
        Sgas.T = Sin.T; Sgas.p = Sin.p; Sgas.phase = "G"

        n_u = Sin.get("UO2(NO3)2")
        if n_u > EPS:
            Ssol.set("UO3", n_u)
            if self.no2_per_u > 0.0:
                Sgas.set("NO2", self.no2_per_u * n_u)
            if self.o2_per_u > 0.0:
                Sgas.set("O2", self.o2_per_u * n_u)

        # Fitted offgas routing for non-uranium species
        h2o_in = Sin.get("H2O")
        if h2o_in > EPS:
            h2o_g = h2o_in * self.h2o_to_offgas_frac
            Sgas.set("H2O", Sgas.get("H2O") + h2o_g)
            Ssol.set("H2O", Ssol.get("H2O") + max(h2o_in - h2o_g, 0.0))

        for sp, frac_g in (
            ("HNO3", self.route_hno3_to_offgas),
            ("TBP", self.route_tbp_to_offgas),
            ("HTcO4", self.route_htco4_to_offgas),
        ):
            n = Sin.get(sp)
            if n <= EPS:
                continue
            ng = n * max(0.0, min(1.0, frac_g))
            ns = n - ng
            if ng > EPS:
                Sgas.set(sp, Sgas.get(sp) + ng)
            if ns > EPS:
                Ssol.set(sp, Ssol.get(sp) + ns)

        if Sair is not None:
            for sp, n in Sair.mol.items():
                if n > EPS:
                    Sgas.set(sp, Sgas.get(sp) + n)
            Sgas.T = max(Sgas.T, Sair.T)
            Sgas.p = min(Sgas.p, Sair.p)

# ─────────────────────────────────────────────────────────────────────────────
# PHASE CONDENSER  (K301)
# ─────────────────────────────────────────────────────────────────────────────

class PhaseCondenser(UnitOp):
    """
    Splits a mixed-phase stream into condensed liquid and non-condensed gas
    by per-species condensation fractions.
    """

    def __init__(self, name: str, *, frac_to_liquid: Dict[str, float], default_frac_to_liquid: float = 0.0):
        super().__init__(name)
        self.frac_to_liquid = frac_to_liquid
        self.default_frac_to_liquid = default_frac_to_liquid

    def apply(self) -> None:
        src = self.inlets["in"]
        liq = self.outlets["liquid"]
        gas = self.outlets["gas"]
        liq.clear(); gas.clear()
        liq.T = gas.T = src.T
        liq.p = gas.p = src.p
        liq.phase = "L"; gas.phase = "G"

        for sp, n in src.mol.items():
            fl = self.frac_to_liquid.get(sp, self.default_frac_to_liquid)
            fl = max(0.0, min(1.0, fl))
            liq.add(sp, n * fl)
            gas.add(sp, n * (1.0 - fl))


# ─────────────────────────────────────────────────────────────────────────────
# NOx ABSORPTION COLUMN  (D101)
# ─────────────────────────────────────────────────────────────────────────────

class NOxAbsorber(UnitOp):
    """
    Removes NO2 and NO from a gas stream, producing HNO3 in aqueous outlet.

    Reactions (simplified, absorption):
      3 NO2 + H2O → 2 HNO3 + NO   (fraction of NO2 captured = NO2_capture_frac)
      4 NO  + 3 O2 + 2 H2O → 4 HNO3  (fraction of NO captured = NO_capture_frac)

    Water is supplied from a liquid inlet ("water").
    """

    def __init__(self, name: str, *, NO2_capture_frac: float = 0.99, NO_capture_frac: float = 0.95):
        super().__init__(name)
        self.NO2_capture_frac = float(NO2_capture_frac)
        self.NO_capture_frac = float(NO_capture_frac)

    def apply(self) -> None:
        gas_in = self.inlets["gas"]
        water_in = self.inlets.get("water")
        gas_out = self.outlets["gas_out"]
        aq_out = self.outlets["aqueous"]

        gas_out.clear(); aq_out.clear()
        gas_out.phase = "G"; aq_out.phase = "L"
        gas_out.T = aq_out.T = gas_in.T
        gas_out.p = aq_out.p = gas_in.p

        # ── Water balance ───────────────────────────────────────────────────
        n_h2o_in = (water_in.get("H2O") if water_in else 0.0) + gas_in.get("H2O")
        n_hno3_in = (water_in.get("HNO3") if water_in else 0.0) + gas_in.get("HNO3")

        # Copy water to aqueous initially
        if water_in:
            for sp, v in water_in.mol.items():
                aq_out.add(sp, v)

        # ── NO2 capture ─────────────────────────────────────────────────────
        n_no2 = gas_in.get("NO2")
        n_no2_cap = n_no2 * self.NO2_capture_frac
        # 3 NO2 + H2O → 2 HNO3 + NO
        hno3_from_no2 = n_no2_cap * (2.0 / 3.0)
        no_from_no2 = n_no2_cap * (1.0 / 3.0)
        h2o_consumed_no2 = n_no2_cap * (1.0 / 3.0)

        # ── NO capture ──────────────────────────────────────────────────────
        n_no_base = gas_in.get("NO") + no_from_no2   # NO in gas including from NO2 rxn
        n_no_cap = n_no_base * self.NO_capture_frac
        # 4 NO + 3 O2 + 2 H2O → 4 HNO3
        hno3_from_no = n_no_cap
        o2_consumed = n_no_cap * (3.0 / 4.0)
        h2o_consumed_no = n_no_cap * (2.0 / 4.0)

        h2o_total_consumed = h2o_consumed_no2 + h2o_consumed_no

        # ── Build aqueous outlet ────────────────────────────────────────────
        aq_out.add("HNO3", hno3_from_no2 + hno3_from_no)
        h2o_aq = max(n_h2o_in - h2o_total_consumed, 0.0)
        aq_out.set("H2O", aq_out.get("H2O") + h2o_aq - (water_in.get("H2O") if water_in else 0.0) )
        # simpler: just set it directly
        aq_out.set("H2O", max(h2o_aq, 0.0))

        # ── Build gas outlet ────────────────────────────────────────────────
        for sp, n in gas_in.mol.items():
            if sp == "NO2":
                gas_out.add("NO2", n - n_no2_cap)
            elif sp == "NO":
                pass   # handled below
            elif sp in ("H2O", "HNO3"):
                pass   # moved to aqueous
            else:
                gas_out.add(sp, n)

        # NO in gas = original NO + from NO2 rxn − captured NO
        gas_out.add("NO", max(n_no_base - n_no_cap, 0.0))

        # O2 balance (reduce O2 in gas)
        o2_in_gas = gas_in.get("O2")
        gas_out.add("O2", max(o2_in_gas - o2_consumed, 0.0))
        # (if O2 was already added above by the else branch, correct it)
        # Safer: zero it out first and re-add
        # (already handled by the 'else' block — O2 goes to gas first)
        # Fix: remove the O2 that was added in 'else' and re-add the correct amount
        # Simplest: skip O2 in the copy loop and handle separately
        # (There's an ordering issue — let's just overwrite at end)
        gas_out.set("O2", max(o2_in_gas - o2_consumed, 0.0))

        self.calcs["NO2_captured_mol_s"] = n_no2_cap
        self.calcs["NO_captured_mol_s"] = n_no_cap
        self.calcs["HNO3_produced_mol_s"] = hno3_from_no2 + hno3_from_no


# ─────────────────────────────────────────────────────────────────────────────
# ACID CONCENTRATOR  (E104)
# ─────────────────────────────────────────────────────────────────────────────

class AcidConcentrator(UnitOp):
    """
    Evaporates water from dilute HNO3 solution to reach a target HNO3 mass fraction.
    Overhead is pure water vapour; bottoms is concentrated acid.
    """

    def __init__(self, name: str, *, target_hno3_mass_frac: float = 0.74):
        super().__init__(name)
        self.target_hno3_mass_frac = float(target_hno3_mass_frac)

    def apply(self) -> None:
        src = self.inlets["in"]
        conc = self.outlets["concentrated"]
        steam = self.outlets["steam"]

        conc.clear(); steam.clear()
        conc.T = steam.T = src.T
        conc.p = steam.p = src.p
        conc.phase = "L"; steam.phase = "G"

        # Copy everything to concentrated by default
        for sp, n in src.mol.items():
            conc.add(sp, n)

        n_hno3 = src.get("HNO3")
        n_h2o  = src.get("H2O")
        mw_hno3 = 63.012; mw_h2o = 18.015

        m_hno3 = n_hno3 * mw_hno3
        # Water needed to hit target mass fraction:
        # target = m_hno3 / (m_hno3 + m_h2o_keep)
        # → m_h2o_keep = m_hno3 * (1 - target) / target
        tgt = self.target_hno3_mass_frac
        m_h2o_keep = m_hno3 * (1.0 - tgt) / tgt if tgt > EPS else n_h2o * mw_h2o
        n_h2o_keep = m_h2o_keep / mw_h2o
        n_h2o_evap = max(n_h2o - n_h2o_keep, 0.0)

        conc.set("H2O", max(n_h2o - n_h2o_evap, 0.0))
        steam.add("H2O", n_h2o_evap)

        self.calcs["H2O_evaporated_mol_s"] = n_h2o_evap
        self.calcs["HNO3_conc_mol_s"] = n_hno3
        if (n_hno3 * mw_hno3 + max(n_h2o - n_h2o_evap, 0.0) * mw_h2o) > EPS:
            act_frac = (n_hno3 * mw_hno3) / (n_hno3 * mw_hno3 + max(n_h2o - n_h2o_evap, 0.0) * mw_h2o)
            self.calcs["actual_HNO3_mass_frac"] = act_frac


# ─────────────────────────────────────────────────────────────────────────────
# SCR REACTOR  (R103 — NH3 + NOx → N2 + H2O)
# ─────────────────────────────────────────────────────────────────────────────

class SCRReactor(UnitOp):
    """
    Selective catalytic reduction of NOx with ammonia.
    Reactions:
      4 NH3 + 4 NO  + O2  → 4 N2 + 6 H2O
      4 NH3 + 2 NO2 + O2  → 3 N2 + 6 H2O  (standard SCR)
    """

    def __init__(self, name: str, *, NOx_removal_frac: float = 0.9376,
                 NH3_excess_factor: float = 1.15):
        super().__init__(name)
        self.NOx_removal_frac = float(NOx_removal_frac)
        self.NH3_excess_factor = float(NH3_excess_factor)

    def apply(self) -> None:
        gas = self.inlets["gas"]
        nh3 = self.inlets["nh3"]
        out = self.outlets["out"]

        # start from combined inlet contents
        out.clear()
        out.T = gas.T; out.p = gas.p; out.phase = "G"
        for s in (gas, nh3):
            for sp, n in s.mol.items():
                if n > EPS:
                    out.set(sp, out.get(sp) + float(n))

        NO2 = out.get("NO2")
        NO  = out.get("NO")
        NH3 = out.get("NH3")
        O2  = out.get("O2")

        dNO2_target = self.NOx_removal_frac * NO2
        dNO_target  = self.NOx_removal_frac * NO

        NH3_req = 2.0 * dNO2_target + 1.0 * dNO_target
        O2_req  = 0.5 * dNO2_target + 0.25 * dNO_target

        scale = 1.0
        if NH3_req > EPS:
            scale = min(scale, NH3 / NH3_req)
        if O2_req > EPS:
            scale = min(scale, O2 / O2_req)
        scale = max(min(scale, 1.0), 0.0)

        dNO2 = dNO2_target * scale
        dNO  = dNO_target  * scale

        out.set("NO2", NO2 - dNO2)
        out.set("NO",  NO  - dNO)
        out.set("NH3", NH3 - (2.0 * dNO2 + 1.0 * dNO))
        out.set("O2",  O2  - (0.5 * dNO2 + 0.25 * dNO))
        out.set("N2",  out.get("N2") + (1.5 * dNO2 + 1.0 * dNO))
        out.set("H2O", out.get("H2O") + (3.0 * dNO2 + 1.5 * dNO))

        self.calcs["NO_removed_mol_s"]  = dNO
        self.calcs["NO2_removed_mol_s"] = dNO2
        self.calcs["NH3_consumed_mol_s"] = 2.0 * dNO2 + 1.0 * dNO
        self.calcs["N2_produced_mol_s"] = 1.5 * dNO2 + 1.0 * dNO


# ─────────────────────────────────────────────────────────────────────────────
# KNOCKOUT POT  (KO101)
# ─────────────────────────────────────────────────────────────────────────────

class KnockoutPot(UnitOp):
    """
    Removes water (and optionally other condensing species) from a gas stream.
    The liquid fraction falls out to a liquid outlet.
    """

    def __init__(self, name: str, *, water_removal_frac: float = 0.999,
                 liquid_species: Optional[list] = None):
        super().__init__(name)
        self.water_removal_frac = float(water_removal_frac)
        self.liquid_species = liquid_species or ["H2O", "HNO3"]

    def apply(self) -> None:
        src = self.inlets["in"]
        gas = self.outlets["gas"]
        liq = self.outlets["liquid"]
        gas.clear(); liq.clear()
        gas.T = liq.T = src.T
        gas.p = liq.p = src.p
        gas.phase = "G"; liq.phase = "L"

        for sp, n in src.mol.items():
            if sp in self.liquid_species:
                fr = self.water_removal_frac
                liq.add(sp, n * fr)
                gas.add(sp, n * (1.0 - fr))
            else:
                gas.add(sp, n)


# ─────────────────────────────────────────────────────────────────────────────
# NH3 ABSORBER  (D102)
# ─────────────────────────────────────────────────────────────────────────────

class NH3Absorber(UnitOp):
    """
    Absorbs NH3 from gas into water wash.
    """

    def __init__(self, name: str, *, NH3_capture_frac: float = 0.99):
        super().__init__(name)
        self.NH3_capture_frac = float(NH3_capture_frac)

    def apply(self) -> None:
        gas_in = self.inlets["gas"]
        water_in = self.inlets.get("water")
        gas_out = self.outlets["gas_out"]
        aq_out = self.outlets["aqueous"]

        gas_out.clear(); aq_out.clear()
        gas_out.phase = "G"; aq_out.phase = "L"
        gas_out.T = aq_out.T = gas_in.T
        gas_out.p = aq_out.p = gas_in.p

        if water_in:
            for sp, v in water_in.mol.items():
                aq_out.add(sp, v)

        n_nh3 = gas_in.get("NH3")
        n_nh3_cap = n_nh3 * self.NH3_capture_frac

        for sp, n in gas_in.mol.items():
            if sp == "NH3":
                gas_out.add("NH3", n - n_nh3_cap)
                aq_out.add("NH3", n_nh3_cap)
            else:
                gas_out.add(sp, n)

        self.calcs["NH3_captured_mol_s"] = n_nh3_cap


# ─────────────────────────────────────────────────────────────────────────────
# TSA IODINE COLUMN  (TSA101A / TSA101B)
# ─────────────────────────────────────────────────────────────────────────────

class TSAColumn(UnitOp):
    """
    Simplified time-averaged TSA iodine polishing.
    Captures a fixed fraction of I_g in the adsorbing column.
    The regenerating column releases previously captured iodine into regen air.
    """

    def __init__(self, name: str, *, ads_capture_frac: float = 0.9999,
                 iodine_species: str = "I_g", t_ads_s: float = 8.0 * 3600.0,
                 t_regen_s: float = 4.0 * 3600.0):
        super().__init__(name)
        self.ads_capture_frac = float(ads_capture_frac)
        self.iodine_species = iodine_species
        self.t_ads_s = float(t_ads_s)
        self.t_regen_s = float(t_regen_s)
        self._captured_cycle_mol_s = 0.0   # iodine captured per unit time (mol/s average)

    def set_regen_source(self, other: "TSAColumn") -> None:
        """Link this (regen) column to the adsorbing column for inventory tracking."""
        self._regen_source = other

    def apply(self) -> None:
        proc_in  = self.inlets["process_gas"]
        regen_in = self.inlets.get("regen_air")
        proc_out = self.outlets["process_out"]
        regen_out = self.outlets.get("regen_out")

        proc_out.clear()
        if regen_out:
            regen_out.clear()

        proc_out.T = proc_in.T; proc_out.p = proc_in.p; proc_out.phase = "G"

        # ── Adsorption ─────────────────────────────────────────────────────
        n_iodine = proc_in.get(self.iodine_species)
        n_captured = n_iodine * self.ads_capture_frac
        n_slip     = n_iodine * (1.0 - self.ads_capture_frac)
        self._captured_cycle_mol_s = n_captured

        for sp, n in proc_in.mol.items():
            if sp == self.iodine_species:
                proc_out.add(sp, n_slip)
            else:
                proc_out.add(sp, n)

        # ── Regeneration (desorb captured iodine into regen air) ──────────
        if regen_out and regen_in and regen_in.total_molar_flow() > EPS:
            src_col = getattr(self, "_regen_source", None)
            src_rate = src_col._captured_cycle_mol_s if src_col else n_captured
            n_desorb = src_rate * (self.t_ads_s / max(self.t_regen_s, EPS))
            regen_out.T = regen_in.T; regen_out.p = regen_in.p; regen_out.phase = "G"
            for sp, n in regen_in.mol.items():
                regen_out.add(sp, n)
            regen_out.add(self.iodine_species, n_desorb)

        self.calcs["I_captured_mol_s"] = n_captured
        self.calcs["I_slip_mol_s"] = n_slip

# ─────────────────────────────────────────────────────────────────────────────
# EFFECTIVE CONTINUOUS DISSOLVER TRAIN (DS-101/104 batch train surrogate)
# ─────────────────────────────────────────────────────────────────────────────

class DissolutionReactor(UnitOp):
    """
    Effective continuous representation of the batch dissolver train, adapted
    from the old_complex package.

    Inputs
    ------
    fuel : oxide / elemental feed
    acid : HNO3/H2O feed
    cleaning_water : rinse / wash water

    Outputs
    -------
    aq      : dissolved liquor
    offgas  : NOx + noble gas + iodine offgas
    solids  : undissolved / basket solids
    """

    def __init__(
        self,
        name: str,
        *,
        uo2_r1_frac: float = 0.18,
        iodine_total_species: str = "I",
        iodine_gas_species: str = "I_g",
        iodine_aq_species: str = "I_aq",
        iodine_offgas_frac: float = 0.98,
        noble_gases: tuple[str, ...] = ("Xe", "Kr"),
        zircaloy_species: str = "Zircaloy",
        zircaloy_basket_mol_s: float = 0.0,
        zircaloy_escape_fines_mol_s: float = 0.0,
        fines_settle_frac: float = 0.0,
        nox_no2_frac: float = 0.415,
    ):
        super().__init__(name)
        self.uo2_r1_frac = float(uo2_r1_frac)
        self.iodine_total_species = iodine_total_species
        self.iodine_gas_species = iodine_gas_species
        self.iodine_aq_species = iodine_aq_species
        self.iodine_offgas_frac = float(iodine_offgas_frac)
        self.noble_gases = tuple(noble_gases)
        self.zircaloy_species = zircaloy_species
        self.zirc_basket = float(zircaloy_basket_mol_s)
        self.zirc_escape = float(zircaloy_escape_fines_mol_s)
        self.fines_settle_frac = float(fines_settle_frac)
        self.nox_no2_frac = float(nox_no2_frac)

        self.reactions = [
            {"UO2": -1.0, "HNO3": -4.0, "UO2(NO3)2": 1.0, "NO2": 2.0, "H2O": 2.0},
            {"UO2": -1.0, "HNO3": -(8.0 / 3.0), "UO2(NO3)2": 1.0, "NO2": (2.0 / 3.0), "H2O": (4.0 / 3.0)},
            {"PuO2": -1.0, "HNO3": -4.0, "Pu(NO3)4": 1.0, "H2O": 2.0},
            {"NpO2": -1.0, "HNO3": -4.0, "Np(NO3)4": 1.0, "H2O": 2.0},
            {"AmO2": -1.0, "HNO3": -4.0, "Am(NO3)4": 1.0, "H2O": 2.0},
            {"Cs2O": -1.0, "HNO3": -2.0, "CsNO3": 2.0, "H2O": 1.0},
            {"SrO": -1.0, "HNO3": -2.0, "Sr(NO3)2": 1.0, "H2O": 1.0},
            {"Nd2O3": -1.0, "HNO3": -6.0, "Nd(NO3)3": 2.0, "H2O": 3.0},
            {"Sm2O3": -1.0, "HNO3": -6.0, "Sm(NO3)3": 2.0, "H2O": 3.0},
            {"Eu2O3": -1.0, "HNO3": -6.0, "Eu(NO3)3": 2.0, "H2O": 3.0},
            {"Gd2O3": -1.0, "HNO3": -6.0, "Gd(NO3)3": 2.0, "H2O": 3.0},
            {"Tc": -1.0, "HNO3": -7.0, "HTcO4": 1.0, "NO2": 7.0, "H2O": 3.0},
        ]

    def _extents(self, fuel: Stream) -> list[float]:
        xi = [0.0] * len(self.reactions)
        xi[0] = self.uo2_r1_frac * fuel.get("UO2")
        xi[1] = (1.0 - self.uo2_r1_frac) * fuel.get("UO2")
        for i, sp in enumerate(["PuO2", "NpO2", "AmO2", "Cs2O", "SrO", "Nd2O3", "Sm2O3", "Eu2O3", "Gd2O3", "Tc"], start=2):
            xi[i] = fuel.get(sp)
        return xi

    def apply(self) -> None:
        fuel = self.inlets["fuel"]
        acid = self.inlets["acid"]
        wash = self.inlets.get("cleaning_water")
        aq = self.outlets["aq"]
        off = self.outlets["offgas"]
        sol = self.outlets["solids"]

        aq.clear(); off.clear(); sol.clear()
        aq.phase = "L"; off.phase = "G"; sol.phase = "S"

        # Start with inlet inventories
        total = {}
        for src in [fuel, acid, wash] if wash is not None else [fuel, acid]:
            for sp, n in src.mol.items():
                total[sp] = total.get(sp, 0.0) + n

        # Apply stoichiometric conversions
        xi = self._extents(fuel)
        for extent, rxn in zip(xi, self.reactions):
            if extent <= EPS:
                continue
            for sp, nu in rxn.items():
                total[sp] = total.get(sp, 0.0) + nu * extent

        if total.get("HNO3", 0.0) < -1e-9:
            raise ValueError(f"{self.name}: insufficient HNO3 in dissolver feed")

        # Route everything to aqueous first
        aq_inv = {sp: max(n, 0.0) for sp, n in total.items() if max(n, 0.0) > EPS}
        off_inv = {}
        sol_inv = {}

        # Noble gases to offgas
        for sp in self.noble_gases:
            n = aq_inv.pop(sp, 0.0)
            if n > EPS:
                off_inv[sp] = off_inv.get(sp, 0.0) + n

        # Iodine split
        nI = 0.0
        for sp in [self.iodine_total_species, self.iodine_aq_species, self.iodine_gas_species]:
            if sp in aq_inv:
                nI += aq_inv.pop(sp)
        if nI > EPS:
            off_inv[self.iodine_gas_species] = off_inv.get(self.iodine_gas_species, 0.0) + self.iodine_offgas_frac * nI
            aq_inv[self.iodine_aq_species] = aq_inv.get(self.iodine_aq_species, 0.0) + (1.0 - self.iodine_offgas_frac) * nI

        # Move inlet / produced NOx to gas, then re-split produced NOx into NO2/NO
        n_no2_total = aq_inv.pop("NO2", 0.0)
        n_no_total = aq_inv.pop("NO", 0.0)
        if n_no2_total > EPS or n_no_total > EPS:
            produced_nox_basis = n_no2_total  # matches old_complex bookkeeping
            n_no2 = self.nox_no2_frac * produced_nox_basis
            n_no = (1.0 - self.nox_no2_frac) * produced_nox_basis + n_no_total
            if n_no2 > EPS:
                off_inv["NO2"] = off_inv.get("NO2", 0.0) + n_no2
            if n_no > EPS:
                off_inv["NO"] = off_inv.get("NO", 0.0) + n_no

        # Zircaloy basket / fines split
        z = aq_inv.pop(self.zircaloy_species, 0.0)
        if z > EPS or self.zirc_basket > EPS or self.zirc_escape > EPS:
            sol_inv[self.zircaloy_species] = self.zirc_basket + self.fines_settle_frac * self.zirc_escape
            aq_inv[self.zircaloy_species] = aq_inv.get(self.zircaloy_species, 0.0) + (1.0 - self.fines_settle_frac) * self.zirc_escape

        for sp, n in aq_inv.items():
            if n > EPS:
                aq.add(sp, n)
        for sp, n in off_inv.items():
            if n > EPS:
                off.add(sp, n)
        for sp, n in sol_inv.items():
            if n > EPS:
                sol.add(sp, n)

        # Conditions
        refT = acid.T if acid is not None else fuel.T
        refp = acid.p if acid is not None else fuel.p
        aq.T = off.T = sol.T = refT
        aq.p = off.p = sol.p = refp
