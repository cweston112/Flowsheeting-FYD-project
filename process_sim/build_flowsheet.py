from __future__ import annotations

from typing import Dict

from . import inputs as cfg
from .framework import ComponentRegistry, Flowsheet, Stream, Conditioner, solve_recycles, EPS
from .unitops import (
    Mixer, PassThrough, SplitUnit, ThreeOutletSplitUnit, EV103CPerformanceStage, PurgeSplitter, Coalescer, Evaporator,
    ConversionReactor, PhaseCondenser, NOxAbsorber, AcidConcentrator,
    SCRReactor, KnockoutPot, NH3Absorber, TSAColumn, DissolutionReactor, CalcinerReactor,
)


# ---------------------------------------------------------------------------
# REGISTRY
# ---------------------------------------------------------------------------
def _build_registry() -> ComponentRegistry:
    reg = ComponentRegistry(cfg.ATOMIC_WEIGHTS)
    for name, formula in cfg.SPECIES_FORMULAS.items():
        if name in cfg.SPECIES_MW_OVERRIDE:
            reg.register(name, mw=cfg.SPECIES_MW_OVERRIDE[name])
        else:
            try:
                reg.register(name, formula=formula)
            except Exception:
                reg.register(name, mw=100.0)
    for name, mw in cfg.SPECIES_MW_OVERRIDE.items():
        reg.register(name, mw=mw)
    return reg


def _stamp(unit) -> None:
    uc = cfg.UNIT_CONDITIONS.get(unit.name, {})
    T = uc.get("T")
    p = uc.get("p")
    for s in unit.outlets.values():
        if T is not None:
            s.T = T
        if p is not None:
            s.p = p


def _run(unit) -> None:
    unit.apply()
    _stamp(unit)


# ---------------------------------------------------------------------------
# SIZING HELPERS
# ---------------------------------------------------------------------------
def _set_aqueous_molarities(s: Stream, *, Vdot_L_s: float, molarity: Dict[str, float], water_molarity_pure: float = 55.5) -> None:
    s.clear()
    n_sol = 0.0
    for sp, M in molarity.items():
        n = float(M) * float(Vdot_L_s)
        if n > EPS:
            s.set(sp, n)
            n_sol += n
    n_w = max(float(water_molarity_pure) * float(Vdot_L_s) - n_sol, 0.0)
    if n_w > EPS:
        s.set("H2O", n_w)


def _size_acid_makeup(fs: Flowsheet) -> None:
    rec = fs.streams["F2R"]
    f2m = fs.streams["F2M"]
    f2w = fs.streams["F2W"]
    tgt_hno3 = cfg.ACID_TARGETS["hno3_mol_s"]
    tgt_h2o = cfg.ACID_TARGETS["h2o_mol_s"]
    r = cfg.ACID_STOCK["water_per_hno3_mol"]

    d_hno3 = max(tgt_hno3 - rec.get("HNO3"), 0.0)
    stock_h2o = d_hno3 * r
    d_h2o = max(tgt_h2o - rec.get("H2O") - stock_h2o, 0.0)

    f2m.clear(); f2w.clear()
    if d_hno3 > EPS:
        f2m.set("HNO3", d_hno3)
    if stock_h2o > EPS:
        f2m.set("H2O", stock_h2o)
    if d_h2o > EPS:
        f2w.set("H2O", d_h2o)


def _size_solvent_makeup(fs: Flowsheet, n_org_req_total: float) -> None:
    rec = fs.streams["F15"]
    f16 = fs.streams["F16"]
    reg = fs.reg
    tbp = cfg.SOLVENT["TBP_name"]
    dil = cfg.SOLVENT["diluent_name"]
    tbp_wt = cfg.SOLVENT["TBP_wt_frac"]

    M_tbp = reg.get_mw(tbp)
    M_dil = reg.get_mw(dil)
    M_mix = 1.0 / (tbp_wt / M_tbp + (1.0 - tbp_wt) / M_dil)
    mdot_req = n_org_req_total * M_mix
    mdot_tbp_req = tbp_wt * mdot_req
    mdot_dil_req = (1.0 - tbp_wt) * mdot_req
    mdot_tbp_rec = rec.get(tbp) * M_tbp
    mdot_dil_rec = rec.get(dil) * M_dil

    f16.clear()
    if mdot_tbp_req > mdot_tbp_rec + EPS:
        f16.set(tbp, (mdot_tbp_req - mdot_tbp_rec) / M_tbp)
    if mdot_dil_req > mdot_dil_rec + EPS:
        f16.set(dil, (mdot_dil_req - mdot_dil_rec) / M_dil)


def _compute_acid_purge_frac(fs: Flowsheet) -> float:
    s38 = fs.streams["F38"]
    tgt_hno3 = cfg.ACID_TARGETS["hno3_mol_s"]
    tgt_h2o = cfg.ACID_TARGETS["h2o_mol_s"]
    r = cfg.ACID_STOCK["water_per_hno3_mol"]
    H38 = s38.get("HNO3")
    W38 = s38.get("H2O")
    if H38 <= EPS and W38 <= EPS:
        return 0.0
    keep_max = 1.0
    if H38 > EPS:
        keep_max = min(keep_max, tgt_hno3 / H38)
    denom = W38 + r * H38
    rhs = tgt_h2o - r * tgt_hno3
    if denom > EPS:
        keep_max = min(keep_max, rhs / denom)
    keep = max(min(keep_max, 1.0), 0.0)
    return 1.0 - keep


def _size_f43(fs: Flowsheet) -> None:
    f42 = fs.streams["F42"]
    eta = float(cfg.R103["NOx_removal_frac"])
    nh3_stoich = eta * f42.get("NO") + 2.0 * eta * f42.get("NO2")
    nh3_feed = nh3_stoich * cfg.R103["NH3_excess_factor"]
    s = fs.streams["F43"]
    s.clear()
    if nh3_feed > EPS:
        s.set("NH3", nh3_feed)


def _size_f47(fs: Flowsheet) -> None:
    n_gas = fs.streams["F46"].total_molar_flow()
    s = fs.streams["F47"]
    s.clear()
    s.set("H2O", n_gas * cfg.D102["LG_water_per_gas"])


def _size_f39(fs: Flowsheet) -> None:
    n_gas = fs.streams["F41"].total_molar_flow()
    s = fs.streams["F39"]
    s.clear()
    s.set("H2O", n_gas * cfg.D101["LG_water_per_gas"])


def _size_tsa_regen_air(fs: Flowsheet) -> None:
    y_max = cfg.TSA["regen_max_y_I2"]
    active = cfg.TSA["active_column"].upper()
    ads = fs.units["TSA101A"] if active == "A" else fs.units["TSA101B"]
    regen_name = "F52B" if active == "A" else "F52A"
    t_ads = float(cfg.TSA.get("t_ads_s", 8.0 * 3600.0))
    t_reg = float(cfg.TSA.get("t_regen_s", 4.0 * 3600.0))
    regen_factor = t_ads / max(t_reg, EPS)
    n_captured_avg = getattr(ads, "_captured_cycle_mol_s", 0.0)
    n_desorb = n_captured_avg * regen_factor
    n_air = n_desorb * (1.0 - y_max) / max(y_max, EPS) if n_desorb > EPS else 0.0
    for name in ["F52A", "F52B"]:
        s = fs.streams[name]
        s.clear()
        if name == regen_name and n_air > EPS:
            s.set("O2", n_air * cfg.TSA["regen_O2_frac"])
            s.set("N2", n_air * cfg.TSA["regen_N2_frac"])


# ---------------------------------------------------------------------------
# BUILDER
# ---------------------------------------------------------------------------
def build_flowsheet(*, max_iter=None, tol=None, relax=None, verbose=None) -> Flowsheet:
    solver_cfg = cfg.SOLVER
    max_iter = solver_cfg["max_iter"] if max_iter is None else max_iter
    tol = solver_cfg["tol"] if tol is None else tol
    relax = solver_cfg["relax"] if relax is None else relax
    verbose = solver_cfg["verbose"] if verbose is None else verbose

    reg = _build_registry()
    fs = Flowsheet(reg, default_T=298.15, default_p=1e5)

    def S(name: str, phase="L", mol=None, T=None, p=None):
        return fs.add_stream(name, phase=phase, mol=mol or {}, T=T, p=p)

    # exogenous feeds
    for sname, spec in cfg.EXOGENOUS_STREAMS.items():
        S(sname, phase=spec["phase"], mol=spec["mol"], T=spec["T"], p=spec["p"])

    # tears / recycle seeds
    S("F15", "L", {"TBP": 0.5, "Dodecane": 2.0})
    S("F2R", "L", {"HNO3": 1.0, "H2O": 2.0})
    S("F45R", "L", {"H2O": 0.5})

    # process streams
    for name, phase in [
        ("F2WM", "L"), ("F2", "L"), ("F2S", "L"), ("F2SP", "L"), ("F2STP", "L"),
        ("F3", "G"), ("F4", "S"), ("F5", "L"), ("F6", "L"), ("F6P", "L"),
        ("F7", "L"), ("F8", "L"), ("F17", "L"), ("F53", "L"), ("F53PT", "L"), ("F54", "L"),
        ("F19", "L"), ("F20", "G"), ("F18", "L"), ("F10", "L"), ("F10T", "L"), ("F55", "L"), ("F55PT", "L"),
        ("F56", "L"), ("F21", "L"), ("F22", "G"), ("F23T", "G"), ("F23TP", "G"),
        ("F11T", "L"), ("F12", "L"), ("F12P", "L"), ("F13", "L"), ("F13P", "L"), ("F57", "L"), ("F58", "L"), ("F57T", "L"), ("F67", "L"),
        ("F64_raw", "L"), ("F64_aux", "L"), ("F64", "L"), ("F65", "G"), ("F66", "L"), ("F68", "L"), ("F69", "L"),
        ("F70", "L"), ("F70_sink", "L"), ("F71", "L"), ("F72", "L"), ("F24", "L"), ("F24P", "L"),
        ("F24PT", "L"), ("F25", "G"), ("F26", "S"), ("F26T", "S"), ("F29", "G"), ("F29P", "G"), ("F29TP", "G"), ("F27", "S"),
        ("F28", "G"), ("F30", "G"), ("F30P", "G"), ("F30TP", "G"), ("F31", "L"), ("F32", "G"), ("F63", "G"),
        ("F63P", "G"), ("F41", "G"), ("F41T", "G"), ("F39", "L"), ("F35", "L"),
        ("F42", "G"), ("F36", "L"), ("F36T", "L"), ("F37", "G"), ("F38", "L"),
        ("F59", "L"), ("F44", "G"), ("F45", "L"), ("F46", "G"), ("F47T", "L"),
        ("F48", "L"), ("F49", "G"), ("F60", "L"), ("F49A", "G"), ("F49B", "G"),
        ("F50A", "G"), ("F50B", "G"), ("F51A", "G"), ("F51B", "G"), ("F61", "L"),
        ("F62", "L"), ("F13T", "L"), ("F14", "L"), ("F_StackGas", "G"),
    ]:
        S(name, phase)

    # front-end / dissolver train
    m105 = Mixer("M105_WaterTrimMixer"); m105.connect_inlet("fresh", fs.streams["F2W"]); m105.connect_inlet("ko", fs.streams["F45R"]); m105.connect_outlet("out", fs.streams["F2WM"])
    m102 = Mixer("M102_AcidMakeupMixer"); m102.connect_inlet("acid", fs.streams["F2M"]); m102.connect_inlet("water", fs.streams["F2WM"]); m102.connect_inlet("recycle", fs.streams["F2R"]); m102.connect_outlet("out", fs.streams["F2"])
    v101 = PassThrough("V101_AcidTank"); v101.connect_inlet("in", fs.streams["F2"]); v101.connect_outlet("out", fs.streams["F2S"])
    p107 = Conditioner("P107_AcidTransfer", target_T=cfg.UNIT_CONDITIONS["P107_AcidTransfer"]["T"], target_p=cfg.UNIT_CONDITIONS["P107_AcidTransfer"]["p"]); p107.connect_inlet("in", fs.streams["F2S"]); p107.connect_outlet("out", fs.streams["F2SP"])
    e101acid = Conditioner("E101_AcidHeater", target_T=cfg.UNIT_CONDITIONS["E101_AcidHeater"]["T"]); e101acid.connect_inlet("in", fs.streams["F2SP"]); e101acid.connect_outlet("out", fs.streams["F2STP"])
    ds = DissolutionReactor("DS101_104_DissolverTrain", **cfg.DISSOLVER)
    ds.connect_inlet("fuel", fs.streams["F1"]); ds.connect_inlet("acid", fs.streams["F2STP"]); ds.connect_inlet("cleaning_water", fs.streams["F2B"])
    ds.connect_outlet("offgas", fs.streams["F3"]); ds.connect_outlet("aq", fs.streams["F5"]); ds.connect_outlet("solids", fs.streams["F4"])
    v102 = PassThrough("V102_DissolverHold"); v102.connect_inlet("in", fs.streams["F5"]); v102.connect_outlet("out", fs.streams["F6"])
    p108 = Conditioner("P108_ToX101", target_T=cfg.UNIT_CONDITIONS["P108_ToX101"]["T"], target_p=cfg.UNIT_CONDITIONS["P108_ToX101"]["p"]); p108.connect_inlet("in", fs.streams["F6"]); p108.connect_outlet("out", fs.streams["F6P"])
    v103 = PassThrough("V103_OffgasSurge"); v103.connect_inlet("in", fs.streams["F3"]); v103.connect_outlet("out", fs.streams["F63"])
    p104 = Conditioner("P104_OffgasBlower", target_T=cfg.UNIT_CONDITIONS["P104_OffgasBlower"]["T"], target_p=cfg.UNIT_CONDITIONS["P104_OffgasBlower"]["p"]); p104.connect_inlet("in", fs.streams["F63"]); p104.connect_outlet("out", fs.streams["F63P"])

    # solvent extraction
    m114 = Mixer("M114_SolventMakeupMixer"); m114.connect_inlet("makeup", fs.streams["F16"]); m114.connect_inlet("recycle", fs.streams["F15"]); m114.connect_outlet("out", fs.streams["F7"])
    S("F_X101", "L"); mx101 = Mixer("M_X101_feed"); mx101.connect_inlet("aq", fs.streams["F6P"]); mx101.connect_inlet("org", fs.streams["F7"]); mx101.connect_outlet("out", fs.streams["F_X101"])
    x101 = SplitUnit("X101_TBP_Extraction", outlet_A="organic", outlet_B="aqueous", frac_to_A=cfg.X101["frac_to_organic"], default_frac_A=cfg.X101["default_frac_to_organic"]); x101.connect_inlet("in", fs.streams["F_X101"]); x101.connect_outlet("organic", fs.streams["F8"]); x101.connect_outlet("aqueous", fs.streams["F17"])
    cf101 = Coalescer("CF101_Coalescer", frac_removed=cfg.CF101["frac_removed"], default_frac_removed=cfg.CF101["default_frac_removed"]); cf101.connect_inlet("in", fs.streams["F17"]); cf101.connect_outlet("aqueous", fs.streams["F53"]); cf101.connect_outlet("recovered_org", fs.streams["F54"])
    e103 = Conditioner("E103_RaffinateHeater", target_T=cfg.UNIT_CONDITIONS["E103_RaffinateHeater"]["T"], target_p=cfg.UNIT_CONDITIONS["E103_RaffinateHeater"].get("p")); e103.connect_inlet("in", fs.streams["F53"]); e103.connect_outlet("out", fs.streams["F53PT"])
    e102raff = Evaporator("E102_RaffinateEvaporator", frac_to_vapour=cfg.E101["frac_to_vapour"], default_frac_to_vapour=cfg.E101["default_frac_to_vapour"]); e102raff.connect_inlet("in", fs.streams["F53PT"]); e102raff.connect_outlet("vapour", fs.streams["F20"]); e102raff.connect_outlet("liquid", fs.streams["F19"])

    S("F_X102", "L"); mx102 = Mixer("M_X102_feed"); mx102.connect_inlet("org", fs.streams["F8"]); mx102.connect_inlet("aq", fs.streams["F9"]); mx102.connect_outlet("out", fs.streams["F_X102"])
    x102 = SplitUnit("X102_AHA_Strip", outlet_A="aqueous", outlet_B="organic", frac_to_A=cfg.X102["frac_to_aqueous"], default_frac_A=cfg.X102["default_frac_to_aqueous"]); x102.connect_inlet("in", fs.streams["F_X102"]); x102.connect_outlet("aqueous", fs.streams["F18"]); x102.connect_outlet("organic", fs.streams["F10"])
    cf102 = Coalescer("CF102_Coalescer", frac_removed=cfg.CF102["frac_removed"], default_frac_removed=cfg.CF102["default_frac_removed"]); cf102.connect_inlet("in", fs.streams["F18"]); cf102.connect_outlet("aqueous", fs.streams["F55"]); cf102.connect_outlet("recovered_org", fs.streams["F56"])
    e104s = Conditioner("E104_PuNpFeedHeater", target_T=cfg.UNIT_CONDITIONS["E104_PuNpFeedHeater"]["T"], target_p=cfg.UNIT_CONDITIONS["E104_PuNpFeedHeater"].get("p")); e104s.connect_inlet("in", fs.streams["F55"]); e104s.connect_outlet("out", fs.streams["F55PT"])
    ev101 = Evaporator("EV101_PuNpEvaporator", frac_to_vapour=cfg.E102["frac_to_vapour"], default_frac_to_vapour=cfg.E102["default_frac_to_vapour"], destroy_aha=True, aha_products=cfg.EVAP_AHA["products"], aha_hno3_consumption=cfg.EVAP_AHA["hno3_consumption"]); ev101.connect_inlet("in", fs.streams["F55PT"]); ev101.connect_outlet("vapour", fs.streams["F22"]); ev101.connect_outlet("liquid", fs.streams["F21"])
    m201 = Mixer("M201_EvapOverheadMixer"); m201.connect_inlet("e101", fs.streams["F20"]); m201.connect_inlet("e102", fs.streams["F22"]); m201.connect_outlet("out", fs.streams["F23T"])
    b101 = Conditioner("B101_OverheadBlower", target_T=cfg.UNIT_CONDITIONS["B101_OverheadBlower"]["T"], target_p=cfg.UNIT_CONDITIONS["B101_OverheadBlower"]["p"]); b101.connect_inlet("in", fs.streams["F23T"]); b101.connect_outlet("out", fs.streams["F23TP"])

    e106 = Conditioner("E106_OrgStripPreheater", target_T=cfg.UNIT_CONDITIONS["E106_OrgStripPreheater"]["T"], target_p=cfg.UNIT_CONDITIONS["E106_OrgStripPreheater"].get("p")); e106.connect_inlet("in", fs.streams["F10"]); e106.connect_outlet("out", fs.streams["F10T"])
    e105 = Conditioner("E105_AqStripPreheater", target_T=cfg.UNIT_CONDITIONS["E105_AqStripPreheater"]["T"], target_p=cfg.UNIT_CONDITIONS["E105_AqStripPreheater"].get("p")); e105.connect_inlet("in", fs.streams["F11"]); e105.connect_outlet("out", fs.streams["F11T"])
    S("F_X103", "L"); mx103 = Mixer("M_X103_feed"); mx103.connect_inlet("org", fs.streams["F10T"]); mx103.connect_inlet("aq", fs.streams["F11T"]); mx103.connect_outlet("out", fs.streams["F_X103"])
    x103 = SplitUnit("X103_Final_Strip", outlet_A="aqueous", outlet_B="organic", frac_to_A=cfg.X103["frac_to_aqueous"], default_frac_A=cfg.X103["default_frac_to_aqueous"]); x103.connect_inlet("in", fs.streams["F_X103"]); x103.connect_outlet("aqueous", fs.streams["F12"]); x103.connect_outlet("organic", fs.streams["F13"])
    p109 = Conditioner("P109_UFeedPump", target_T=cfg.UNIT_CONDITIONS["P109_UFeedPump"]["T"], target_p=cfg.UNIT_CONDITIONS["P109_UFeedPump"]["p"]); p109.connect_inlet("in", fs.streams["F12"]); p109.connect_outlet("out", fs.streams["F12P"])
    cf103 = Coalescer("CF103_Coalescer", frac_removed=cfg.CF103["frac_removed"], default_frac_removed=cfg.CF103["default_frac_removed"]); cf103.connect_inlet("in", fs.streams["F12P"]); cf103.connect_outlet("aqueous", fs.streams["F57"]); cf103.connect_outlet("recovered_org", fs.streams["F58"])
    e107 = Conditioner("E107_UFeedHeater", target_T=cfg.UNIT_CONDITIONS["E107_UFeedHeater"]["T"], target_p=cfg.UNIT_CONDITIONS["E107_UFeedHeater"]["p"]); e107.connect_inlet("in", fs.streams["F57"]); e107.connect_outlet("out", fs.streams["F57T"])
    e107_u = PassThrough("E107_SteamCondensate"); e107_u.connect_inlet("in", fs.streams["U67_src"]); e107_u.connect_outlet("out", fs.streams["F67"])
    ev103a = SplitUnit("EV103A_UConcentrator1", outlet_A="vapour", outlet_B="to_raw", frac_to_A=cfg.EV103A["frac_to_A"], default_frac_A=cfg.EV103A["default_frac_to_A"]); ev103a.connect_inlet("in", fs.streams["F57T"]); ev103a.connect_outlet("vapour", fs.streams["F65"]); ev103a.connect_outlet("to_raw", fs.streams["F64_raw"])
    ev103a_clean = SplitUnit("EV103A_InternalCleaner", outlet_A="aux", outlet_B="clean", frac_to_A={"HNO3": 1.0, "HTcO4": 1.0}, default_frac_A=0.0); ev103a_clean.connect_inlet("in", fs.streams["F64_raw"]); ev103a_clean.connect_outlet("aux", fs.streams["F64_aux"]); ev103a_clean.connect_outlet("clean", fs.streams["F64"])
    ev103a_aux = PassThrough("EV103A_CondensatePass"); ev103a_aux.connect_inlet("in", fs.streams["F67"]); ev103a_aux.connect_outlet("out", fs.streams["F69"])
    ev103b = SplitUnit("EV103B_UConcentrator2", outlet_A="to_C1", outlet_B="to_C2", frac_to_A=cfg.EV103B["frac_to_A"], default_frac_A=cfg.EV103B["default_frac_to_A"]); ev103b.connect_inlet("in", fs.streams["F64"]); ev103b.connect_outlet("to_C1", fs.streams["F66"]); ev103b.connect_outlet("to_C2", fs.streams["F68"])
    ev103b_aux = PassThrough("EV103B_OverheadPass"); ev103b_aux.connect_inlet("in", fs.streams["F65"]); ev103b_aux.connect_outlet("out", fs.streams["F72"])
    ev103c = EV103CPerformanceStage(
        "EV103C_UConcentrator3",
        heavy_to_product=cfg.EV103C_STAGE["heavy_to_product"],
        heavy_to_vapour=cfg.EV103C_STAGE["heavy_to_vapour"],
        light_to_product=cfg.EV103C_STAGE["light_to_product"],
        aux_to_product=cfg.EV103C_STAGE["aux_to_product"],
        extra_h2o_to_vapour=cfg.EV103C_STAGE["extra_h2o_to_vapour"],
    ); ev103c.connect_inlet("light", fs.streams["F66"]); ev103c.connect_inlet("heavy", fs.streams["F68"]); ev103c.connect_inlet("aux", fs.streams["F64_aux"]); ev103c.connect_outlet("product", fs.streams["F24"]); ev103c.connect_outlet("vapour", fs.streams["F25"]); ev103c.connect_outlet("waste", fs.streams["F70"]); ev103c.connect_outlet("sink", fs.streams["F70_sink"])
    m272 = Mixer("M272_EV103_WasteMixer"); m272.connect_inlet("a", fs.streams["F69"]); m272.connect_inlet("b", fs.streams["F72"]); m272.connect_inlet("c", fs.streams["F70"]); m272.connect_outlet("out", fs.streams["F71"])
    p105 = Conditioner("P105_UProductPump", target_T=cfg.UNIT_CONDITIONS["P105_UProductPump"]["T"], target_p=cfg.UNIT_CONDITIONS["P105_UProductPump"]["p"]); p105.connect_inlet("in", fs.streams["F24"]); p105.connect_outlet("out", fs.streams["F24P"])
    e109 = Conditioner("E109_UProductHeater", target_T=cfg.UNIT_CONDITIONS["E109_UProductHeater"]["T"], target_p=cfg.UNIT_CONDITIONS["E109_UProductHeater"]["p"]); e109.connect_inlet("in", fs.streams["F24P"]); e109.connect_outlet("out", fs.streams["F24PT"])
    r101 = CalcinerReactor(
        "R101_ThermalReactor",
        no2_per_u=cfg.R201["no2_per_u"],
        o2_per_u=cfg.R201["o2_per_u"],
        h2o_to_offgas_frac=cfg.R201["h2o_to_offgas_frac"],
        route_hno3_to_offgas=cfg.R201["route_hno3_to_offgas"],
        route_tbp_to_offgas=cfg.R201["route_tbp_to_offgas"],
        route_htco4_to_offgas=cfg.R201["route_htco4_to_offgas"],
    ); r101.connect_inlet("in", fs.streams["F24PT"]); r101.connect_inlet("air", fs.streams["F73"]); r101.connect_outlet("solids", fs.streams["F26"]); r101.connect_outlet("offgas", fs.streams["F29"])
    p106 = Conditioner("P106_ReactorOffgasBlower", target_T=cfg.UNIT_CONDITIONS["P106_ReactorOffgasBlower"]["T"], target_p=cfg.UNIT_CONDITIONS["P106_ReactorOffgasBlower"]["p"]); p106.connect_inlet("in", fs.streams["F29"]); p106.connect_outlet("out", fs.streams["F29P"])
    e111 = Conditioner("E111_ReactorOffgasCooler", target_T=cfg.UNIT_CONDITIONS["E111_ReactorOffgasCooler"]["T"], target_p=cfg.UNIT_CONDITIONS["E111_ReactorOffgasCooler"]["p"]); e111.connect_inlet("in", fs.streams["F29P"]); e111.connect_outlet("out", fs.streams["F29TP"])
    e110 = Conditioner("E110_ProductCooler", target_T=cfg.UNIT_CONDITIONS["E110_ProductCooler"]["T"], target_p=cfg.UNIT_CONDITIONS["E110_ProductCooler"]["p"]); e110.connect_inlet("in", fs.streams["F26"]); e110.connect_outlet("out", fs.streams["F26T"])
    k101 = ConversionReactor("K101_UProductSeparator", reactions=cfg.R202["reactions"], default_to_solid=cfg.R202["default_to_solid"]); k101.connect_inlet("in", fs.streams["F26T"]); k101.connect_outlet("solid", fs.streams["F27"]); k101.connect_outlet("gas", fs.streams["F28"]); k101.connect_outlet("liquid", fs.streams["F28"])
    m301 = Mixer("M301_HotVapourMixer"); m301.connect_inlet("evaps", fs.streams["F23TP"]); m301.connect_inlet("ucon", fs.streams["F25"]); m301.connect_inlet("calc", fs.streams["F29TP"]); m301.connect_outlet("out", fs.streams["F30"])
    p101 = Conditioner("P101_HotVapourBlower", target_T=cfg.UNIT_CONDITIONS["P101_HotVapourBlower"]["T"], target_p=cfg.UNIT_CONDITIONS["P101_HotVapourBlower"]["p"]); p101.connect_inlet("in", fs.streams["F30"]); p101.connect_outlet("out", fs.streams["F30P"])
    e112 = Conditioner("E112_HotVapourCooler", target_T=cfg.UNIT_CONDITIONS["E112_HotVapourCooler"]["T"], target_p=cfg.UNIT_CONDITIONS["E112_HotVapourCooler"].get("p")); e112.connect_inlet("in", fs.streams["F30P"]); e112.connect_outlet("out", fs.streams["F30TP"])
    k301 = PhaseCondenser("K301_Condenser", frac_to_liquid=cfg.K301["frac_to_liquid"], default_frac_to_liquid=cfg.K301["default_frac_to_liquid"]); k301.connect_inlet("in", fs.streams["F30TP"]); k301.connect_outlet("liquid", fs.streams["F31"]); k301.connect_outlet("gas", fs.streams["F32"])

    # offgas / acid recovery
    m302 = Mixer("M302_OffgasMixer"); m302.connect_inlet("cond", fs.streams["F32"]); m302.connect_inlet("diss", fs.streams["F63P"]); m302.connect_inlet("air", fs.streams["F40"]); m302.connect_outlet("out", fs.streams["F41"])
    e113 = Conditioner("E113_OffgasCooler", target_T=cfg.UNIT_CONDITIONS["E113_OffgasCooler"]["T"]); e113.connect_inlet("in", fs.streams["F41"]); e113.connect_outlet("out", fs.streams["F41T"])
    d101 = NOxAbsorber("D101_NO2_Absorber", NO2_capture_frac=cfg.D101["NO2_capture_frac"], NO_capture_frac=cfg.D101["NO_capture_frac"]); d101.connect_inlet("gas", fs.streams["F41T"]); d101.connect_inlet("water", fs.streams["F39"]); d101.connect_outlet("gas_out", fs.streams["F42"]); d101.connect_outlet("aqueous", fs.streams["F35"])
    m303 = Mixer("M303_AcidMixer"); m303.connect_inlet("abs", fs.streams["F35"]); m303.connect_inlet("cond", fs.streams["F31"]); m303.connect_outlet("out", fs.streams["F36"])
    e118 = Conditioner("E118_AcidTrimHeater", target_T=cfg.UNIT_CONDITIONS["E118_AcidTrimHeater"]["T"]); e118.connect_inlet("in", fs.streams["F36"]); e118.connect_outlet("out", fs.streams["F36T"])
    e104 = AcidConcentrator("E104_HNO3_Concentrator", target_hno3_mass_frac=cfg.E104["target_hno3_mass_frac"]); e104.connect_inlet("in", fs.streams["F36T"]); e104.connect_outlet("concentrated", fs.streams["F38"]); e104.connect_outlet("steam", fs.streams["F37"])
    v104 = PurgeSplitter("V104_AcidPurge", purge_frac=cfg.V104["purge_frac"]); v104.connect_inlet("in", fs.streams["F38"]); v104.connect_outlet("purge", fs.streams["F59"]); v104.connect_outlet("recycle", fs.streams["F2R"])
    r103 = SCRReactor("R103_NOxReduction", NOx_removal_frac=cfg.R103["NOx_removal_frac"], NH3_excess_factor=cfg.R103["NH3_excess_factor"]); r103.connect_inlet("gas", fs.streams["F42"]); r103.connect_inlet("nh3", fs.streams["F43"]); r103.connect_outlet("out", fs.streams["F44"])
    ko101 = KnockoutPot("KO101_KnockoutPot", water_removal_frac=cfg.KO101["water_removal_frac"]); ko101.connect_inlet("in", fs.streams["F44"]); ko101.connect_outlet("gas", fs.streams["F46"]); ko101.connect_outlet("liquid", fs.streams["F45"])
    v105 = PurgeSplitter("V105_KOWaterPurge", purge_frac=cfg.V105["purge_frac"]); v105.connect_inlet("in", fs.streams["F45"]); v105.connect_outlet("purge", fs.streams["F60"]); v105.connect_outlet("recycle", fs.streams["F45R"])
    e116 = Conditioner("E116_D102WaterHeater", target_T=cfg.UNIT_CONDITIONS["E116_D102WaterHeater"]["T"]); e116.connect_inlet("in", fs.streams["F47"]); e116.connect_outlet("out", fs.streams["F47T"])
    d102 = NH3Absorber("D102_NH3_Absorber", NH3_capture_frac=cfg.D102["NH3_capture_frac"]); d102.connect_inlet("gas", fs.streams["F46"]); d102.connect_inlet("water", fs.streams["F47T"]); d102.connect_outlet("gas_out", fs.streams["F49"]); d102.connect_outlet("aqueous", fs.streams["F48"])

    # TSA and solvent recovery
    active_A = 1.0 if cfg.TSA["active_column"].upper() == "A" else 0.0
    v201 = SplitUnit("V201_TSA_Split", outlet_A="A", outlet_B="B", frac_to_A={sp: active_A for sp in cfg.SPECIES_FORMULAS}, default_frac_A=active_A); v201.connect_inlet("in", fs.streams["F49"]); v201.connect_outlet("A", fs.streams["F49A"]); v201.connect_outlet("B", fs.streams["F49B"])
    tsaA = TSAColumn("TSA101A", ads_capture_frac=cfg.TSA["ads_capture_frac"], iodine_species=cfg.TSA["iodine_species"], t_ads_s=cfg.TSA["t_ads_s"], t_regen_s=cfg.TSA["t_regen_s"])
    tsaB = TSAColumn("TSA101B", ads_capture_frac=cfg.TSA["ads_capture_frac"], iodine_species=cfg.TSA["iodine_species"], t_ads_s=cfg.TSA["t_ads_s"], t_regen_s=cfg.TSA["t_regen_s"])
    if cfg.TSA["active_column"].upper() == "A":
        tsaB.set_regen_source(tsaA)
    else:
        tsaA.set_regen_source(tsaB)
    tsaA.connect_inlet("process_gas", fs.streams["F49A"]); tsaA.connect_inlet("regen_air", fs.streams["F52A"]); tsaA.connect_outlet("process_out", fs.streams["F50A"]); tsaA.connect_outlet("regen_out", fs.streams["F51A"])
    tsaB.connect_inlet("process_gas", fs.streams["F49B"]); tsaB.connect_inlet("regen_air", fs.streams["F52B"]); tsaB.connect_outlet("process_out", fs.streams["F50B"]); tsaB.connect_outlet("regen_out", fs.streams["F51B"])
    mstack = Mixer("M_StackGasMixer"); mstack.connect_inlet("A", fs.streams["F50A"]); mstack.connect_inlet("B", fs.streams["F50B"]); mstack.connect_outlet("out", fs.streams["F_StackGas"])
    p110 = Conditioner("P110_SolventRecyclePump", target_T=cfg.UNIT_CONDITIONS["P110_SolventRecyclePump"]["T"], target_p=cfg.UNIT_CONDITIONS["P110_SolventRecyclePump"]["p"]); p110.connect_inlet("in", fs.streams["F13"]); p110.connect_outlet("out", fs.streams["F13P"])
    m132 = Mixer("M132_SolventRecoveryMixer"); m132.connect_inlet("c1", fs.streams["F54"]); m132.connect_inlet("c2", fs.streams["F56"]); m132.connect_inlet("c3", fs.streams["F58"]); m132.connect_outlet("out", fs.streams["F61"])
    m133 = Mixer("M133_SolventMixer"); m133.connect_inlet("main", fs.streams["F13P"]); m133.connect_inlet("rec", fs.streams["F61"]); m133.connect_outlet("out", fs.streams["F62"])
    hx133 = Conditioner("HX133_SolventCooler", target_T=298.15); hx133.connect_inlet("in", fs.streams["F62"]); hx133.connect_outlet("out", fs.streams["F13T"])
    v133 = PurgeSplitter("V133_SolventPurge", purge_frac=cfg.V133["purge_frac"]); v133.connect_inlet("in", fs.streams["F13T"]); v133.connect_outlet("purge", fs.streams["F14"]); v133.connect_outlet("recycle", fs.streams["F15"])

    units = [m105, m102, v101, p107, e101acid, ds, v102, p108, v103, p104, m114, mx101, x101, cf101, e103, e102raff,
             mx102, x102, cf102, e104s, ev101, m201, b101, e106, e105, mx103, x103, p109, cf103, e107, e107_u, ev103a, ev103a_clean, ev103a_aux, ev103b, ev103b_aux, ev103c, m272, p105, e109, r101, p106, e111, e110, k101, m301, p101, e112, k301,
             m302, e113, d101, m303, e118, e104, v104, r103, ko101, v105, e116, d102, v201, tsaA, tsaB, mstack, p110, m132,
             m133, hx133, v133]
    for u in units:
        fs.add_unit(u)

    def evaluate_once() -> None:
        # A) size acid make-up from acid recycle
        _size_acid_makeup(fs)
        _run(m105); _run(m102); _run(v101); _run(p107); _run(e101acid); _run(ds); _run(v102); _run(p108)

        # B) size solvent make-up from current aqueous throughput entering X101
        n_aq = fs.streams["F6P"].total_molar_flow()
        _size_solvent_makeup(fs, n_aq * max(cfg.X101_SIZING["OA_target"], cfg.X101_SIZING["OA_min"]))
        _run(m114); _run(mx101); _run(x101); _run(cf101); _run(e103); _run(e102raff)

        # C) size AHA strip feed from organic throughput to X102
        O2 = fs.streams["F8"].total_molar_flow()
        Vdot2 = cfg.X102_SIZING["AO_target"] * O2 / max(cfg.X102_SIZING["water_molarity"], EPS)
        _set_aqueous_molarities(fs.streams["F9"], Vdot_L_s=Vdot2, molarity=cfg.X102_SIZING["molarities"], water_molarity_pure=cfg.X102_SIZING["water_molarity"])
        _run(mx102); _run(x102); _run(cf102); _run(e104s); _run(ev101); _run(m201); _run(b101)

        # D) size final strip feed from organic throughput to X103
        O3 = fs.streams["F10"].total_molar_flow()
        Vdot3 = cfg.X103_SIZING["AO_target"] * O3 / max(cfg.X103_SIZING["water_molarity"], EPS)
        _set_aqueous_molarities(fs.streams["F11"], Vdot_L_s=Vdot3, molarity=cfg.X103_SIZING["molarities"], water_molarity_pure=cfg.X103_SIZING["water_molarity"])
        _run(e106); _run(e105); _run(mx103); _run(x103); _run(p109); _run(cf103); _run(e107); _run(e107_u); _run(ev103a); _run(ev103a_clean); _run(ev103a_aux); _run(ev103b); _run(ev103b_aux); _run(ev103c); _run(m272); _run(p105); _run(e109); _run(r101); _run(p106); _run(e111); _run(e110); _run(k101); _run(m301); _run(p101); _run(e112); _run(k301)

        # E) offgas side and acid recycle
        _run(v103); _run(p104); _run(m302)
        _size_f39(fs)
        _run(e113); _run(d101); _run(m303); _run(e118); _run(e104)
        v104.purge_frac = _compute_acid_purge_frac(fs)
        _run(v104)

        # F) SCR / KO / NH3 absorber / KO recycle
        _size_f43(fs)
        _run(r103); _run(ko101)
        _run(v105)
        _size_f47(fs)
        _run(e116); _run(d102)

        # G) TSA + solvent recycle
        _run(v201)
        if cfg.TSA["active_column"].upper() == "A":
            _run(tsaA); _size_tsa_regen_air(fs); _run(tsaB)
        else:
            _run(tsaB); _size_tsa_regen_air(fs); _run(tsaA)
        _run(mstack)
        _run(p110); _run(m132); _run(m133); _run(hx133); _run(v133)

    evaluate_once()
    result = solve_recycles(
        fs,
        tears=cfg.SOLVER["tear_streams"],
        evaluate=evaluate_once,
        max_iter=max_iter,
        tol=tol,
        relax=relax,
        use_anderson=cfg.SOLVER.get("use_anderson", True),
        anderson_m=cfg.SOLVER.get("anderson_m", 6),
        verbose=verbose,
        print_every=cfg.SOLVER.get("print_every", 10),
    )
    fs.converged = result.converged
    fs.recycle_iter = result.iterations
    fs.recycle_err = result.error
    return fs
