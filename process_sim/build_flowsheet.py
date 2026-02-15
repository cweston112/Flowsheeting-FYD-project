from __future__ import annotations

from typing import List, Iterable, Optional

from .inputs import InputParameters
from .flowsheet_tools import (
    EPS,
    ComponentRegistry,
    Flowsheet,
    Stream,
    stamp_outlets_to_unit_conditions,
    snapshot_streams,
    max_rel_change_all_streams,
    _needs_T,
    _needs_p,
    size_solvent_makeup_for_required_org_flow,
    set_aqueous_molarities,
    get_or_create_conditioned_stream,
    size_acid_stock_and_water_makeup_for_targets, compute_required_acid_purge_fraction,
    compute_required_purge_fraction_with_stock_acid,
)
from .unitops import (
    Mixer,
    DissolutionReactor,
    HeaterCooler,
    PressureChanger,
    MultiStageCentrifugalContactorTBP,
    MultiStageCentrifugalContactorAHA_Strip,
    BufferVessel,
    Splitter,
    ContinuousEvaporatorAHA,
    CoalescingFilterTBP,
    SpeciesSplitter,
    UO2NitrateToUO3Reactor,
    UO3ToU3O8Reactor,
    NitricAcidConcentrator,
    NO2AbsorptionColumn,
    NH3NOxReductionColumn,
    KnockOutPotWater,
    AmmoniaAbsorptionColumn,
    IdealTSAColumnEMM17,
)


# ------------------------------------------------------------------------------
# REGISTRY
# ------------------------------------------------------------------------------

def build_registry() -> ComponentRegistry:
    reg = ComponentRegistry()
    reg.atomic_weights.update(InputParameters.ATOMIC_WEIGHTS)

    formula_species = {
        "UO2": "UO2", "PuO2": "PuO2", "NpO2": "NpO2", "AmO2": "AmO2",
        "Cs2O": "Cs2O", "SrO": "SrO", "Nd2O3": "Nd2O3", "Sm2O3": "Sm2O3",
        "Eu2O3": "Eu2O3", "Gd2O3": "Gd2O3", "Tc": "Tc",
        "Xe": "Xe", "Kr": "Kr",
        "I": "I", "I_g": "I2", "I_aq": "I",
        "Zircaloy": "Zr",
        "HNO3": "HNO3", "NO2": "NO2", "H2O": "H2O",
        "NH3": "NH3",
        "NO": "NO", "N2": "N2", "O2": "O2",
        "TBP": "C12H27O4P", "Dodecane": "C12H26",
        "AHA": "C2H5NO2", "AcOH": "C2H4O2", "N2O": "N2O",
        "UO2(NO3)2": "UO2(NO3)2",
        "Pu(NO3)4": "Pu(NO3)4",
        "Np(NO3)4": "Np(NO3)4",
        "Am(NO3)4": "Am(NO3)4",
        "Nd(NO3)3": "Nd(NO3)3",
        "Sm(NO3)3": "Sm(NO3)3",
        "Eu(NO3)3": "Eu(NO3)3",
        "Gd(NO3)3": "Gd(NO3)3",
        "CsNO3": "CsNO3",
        "Sr(NO3)2": "Sr(NO3)2",
        "HTcO4": "HTcO4",
        "UO3": "UO3",
        "U3O8": "U3O8",
    }
    for name, formula in formula_species.items():
        reg.add_species(name, formula=formula)

    # Explicit aliases (ensure MW exists even if formula differs)
    reg.add_species("I_g", mw_g_mol=reg.mw["I_g"])
    reg.add_species("I_aq", mw_g_mol=reg.mw["I"])
    return reg



# ------------------------------------------------------------------------------
# FLOWSHEET
# ------------------------------------------------------------------------------

def build_flowsheet(*,
                    max_iter: int = 200,
                    tol: float = 1e-8,
                    relax: float = 0.5,
                    include_Tp: bool = True,
                    tear_stream_names: Optional[Iterable[str]] = None) -> Flowsheet:
    """
    Full dissolver + TBP extraction + AHA strip + evaporators + calcination + offgas + HNO3 recovery
    + KO water recycle + TSA iodine polish + solvent purge/recycle.

    Key fixes vs your current version:
      1) **Remove the premature V104 run at the start of each iteration.**
         Running V104 before E104 has produced F38 causes an inconsistent F2_R at iteration start.
         (This is what was biting you: you were “using F38 before it exists”.)

      2) **Move M105 (water-trim mixer) BEFORE M102 in the unit_sequence.**
         M102 consumes F2_WM, but in your code M105 was executed late (after KO/V105),
         so F2_WM was stale when M102 ran.

      3) Purge fraction for V104 is computed **after E104** in each iteration, then V104 is run.
         This updates F2_R for the next iteration’s acid makeup sizing (standard fixed-point).
    """
    params = InputParameters()
    reg = build_registry()
    fs = Flowsheet(reg, default_T=params.DEFAULT_T_K, default_p=params.DEFAULT_P_PA)
    unit_cond = params.UNIT_CONDITIONS

    unit_sequence: List[object] = []

    def _add_unit(u):
        fs.add_unit(u)
        unit_sequence.append(u)
        return u

    def _add_unit_once(u):
        """Add unit to fs + unit_sequence only if it's not already present."""
        if u.name in fs.units:
            u_existing = fs.units[u.name]
            if u_existing not in unit_sequence:
                unit_sequence.append(u_existing)
            return u_existing
        fs.add_unit(u)
        unit_sequence.append(u)
        return u

    def run_unit(u) -> None:
        u.apply()
        stamp_outlets_to_unit_conditions(u, unit_cond)

    def relax_stream(name: str, prev_snap, curr_snap, alpha: float) -> None:
        s = fs.streams[name]
        x0, _, _ = prev_snap[name]
        x1, _, _ = curr_snap[name]
        s.from_dense((1.0 - alpha) * x0 + alpha * x1)
        if include_Tp:
            _, T0, p0 = prev_snap[name]
            _, T1, p1 = curr_snap[name]
            s.T = (1.0 - alpha) * T0 + alpha * T1
            s.p = (1.0 - alpha) * p0 + alpha * p1

    # -------------------------------------------------------------------------
    # DETERMINISTIC T/P CONDITIONING
    # -------------------------------------------------------------------------
    def condition_to_unit(src: Stream, unit_name: str) -> Stream:
        """
        Deterministic conditioning with producer-aware bypass:
          - If src.producer stamps exactly the target conditions, skip conditioning.
          - Else create/reuse HX/PC (via get_or_create_conditioned_stream).
        """
        cond = unit_cond.get(unit_name, {})
        target_T = cond.get("T", None)
        target_p = cond.get("p", None)

        if target_T is None and target_p is None:
            return src

        prod = getattr(src, "producer", None)
        if prod and prod in unit_cond:
            prod_cond = unit_cond.get(prod, {})
            prod_T = prod_cond.get("T", None)
            prod_p = prod_cond.get("p", None)

            same_T = (target_T is None) or (prod_T is not None and abs(float(prod_T) - float(target_T)) <= 1e-9)
            same_p = (target_p is None) or (prod_p is not None and abs(float(prod_p) - float(target_p)) <= 1e-6)
            if same_T and same_p:
                return src

        need_T = (target_T is not None) and _needs_T(src, float(target_T))
        need_p = (target_p is not None) and _needs_p(src, float(target_p))
        if not need_T and not need_p:
            return src

        def make_HX(name: str, T: float):
            return HeaterCooler(name, set_T=T)

        def make_PC(name: str, p: float):
            return PressureChanger(name, set_p=p)

        out = get_or_create_conditioned_stream(
            fs,
            src=src,
            unit_name=unit_name,
            target_T=float(target_T) if need_T else None,
            target_p=float(target_p) if need_p else None,
            make_HX=make_HX,
            make_PC=make_PC,
        )

        # Deterministic execution ordering: ensure any conditioner units are in unit_sequence
        for uname in getattr(fs, "_conditioning_units_created", []):
            if uname in fs.units:
                _add_unit_once(fs.units[uname])

        return out

    # -------------------------------------------------------------------------
    # EXOGENOUS STREAMS
    # -------------------------------------------------------------------------
    F1 = fs.new_exogenous_stream("F1", mol=params.FEED_COMPOSITION)  # fuel feed (solid)

    F2_M = fs.new_exogenous_stream("F2_M", mol={})  # 70% acid stock, sized each iteration
    F2_W = fs.new_exogenous_stream("F2_W", mol={})  # pure water trim, sized each iteration

    F2_B = fs.new_exogenous_stream("F2_B", mol={})  # cleaning water
    F2_B.set("H2O", params.F2B_flow_rate)

    F16 = fs.new_exogenous_stream("F16", mol={})  # organic solvent makeup (TBP+diluent)

    F9 = fs.new_exogenous_stream("F9", mol={})  # AHA feed
    set_aqueous_molarities(
        F9,
        Vdot_L_s=params.X2_STRIP_VDOT_L_S,
        molarity=params.X2_STRIP_MOLARITIES,
        water_name="H2O",
        water_molarity_pure=params.X2_STRIP_WATER_MOLARITY,
    )

    F11 = fs.new_exogenous_stream("F11", mol={})  # dilute nitric acid feed
    set_aqueous_molarities(
        F11,
        Vdot_L_s=params.X3_STRIP_VDOT_L_S,
        molarity=params.X3_STRIP_MOLARITIES,
        water_name="H2O",
        water_molarity_pure=params.X3_STRIP_WATER_MOLARITY,
    )

    F40 = fs.new_exogenous_stream("F40", mol={})
    F40.set("O2", params.F40_AIR_TOTAL_MOL_S * params.F40_AIR_O2_FRAC)
    F40.set("N2", params.F40_AIR_TOTAL_MOL_S * params.F40_AIR_N2_FRAC)

    F43 = fs.new_exogenous_stream("F43", mol={})  # pure NH3 (set each iteration)
    F47 = fs.new_exogenous_stream("F47", mol={})  # water wash (set each iteration)

    # TSA regeneration air feeds (exogenous, sized to be minimal)
    F52A = fs.new_exogenous_stream("F52A", mol={})
    F52B = fs.new_exogenous_stream("F52B", mol={})

    # -------------------------------------------------------------------------
    # PROCESS STREAMS
    # -------------------------------------------------------------------------
    F2 = fs.new_process_stream("F2", phase="L", mol={})
    F3 = fs.new_process_stream("F3", phase="G", mol={})

    F4 = fs.new_process_stream("F4", phase="S", mol={})
    F5 = fs.new_process_stream("F5", phase="L", mol={})
    F6 = fs.new_process_stream("F6", phase="L", mol={})

    F7 = fs.new_process_stream("F7", phase="L", mol={})
    F8 = fs.new_process_stream("F8", phase="L", mol={})
    F17 = fs.new_process_stream("F17", phase="L", mol={})

    F10 = fs.new_process_stream("F10", phase="L", mol={})
    F18 = fs.new_process_stream("F18", phase="L", mol={})

    F12 = fs.new_process_stream("F12", phase="L", mol={})
    F13 = fs.new_process_stream("F13", phase="L", mol={})

    F14 = fs.new_process_stream("F14", phase="L", mol={})
    F15 = fs.new_process_stream("F15", phase="L", mol={})

    F19 = fs.new_process_stream("F19", phase="L", mol={})  # E101 liquor
    F20 = fs.new_process_stream("F20", phase="G", mol={})  # E101 vapour

    F21 = fs.new_process_stream("F21", phase="L", mol={})  # E102 liquor
    F22 = fs.new_process_stream("F22", phase="G", mol={})  # E102 vapour

    F23 = fs.new_process_stream("F23", phase="G", mol={})  # combined vapour

    F101 = fs.new_process_stream("F101", phase="L", mol={})  # CF101 aqueous -> E101
    F105 = fs.new_process_stream("F105", phase="L", mol={})  # CF101 recovered organics

    F102 = fs.new_process_stream("F102", phase="L", mol={})  # CF102 aqueous -> E102
    F106 = fs.new_process_stream("F106", phase="L", mol={})  # CF102 recovered organics

    F24 = fs.new_process_stream("F24", phase="L", mol={})  # E103 liquid
    F25 = fs.new_process_stream("F25", phase="G", mol={})  # E103 vapour

    F26 = fs.new_process_stream("F26", phase="S", mol={})  # UO3 solids
    F29 = fs.new_process_stream("F29", phase="G", mol={})  # calciner offgas

    F30 = fs.new_process_stream("F30", phase="G", mol={})  # hot vapours before condenser

    F31 = fs.new_process_stream("F31", phase="L", mol={})  # condenser liquid (H2O+HNO3)
    F32 = fs.new_process_stream("F32", phase="G", mol={})  # condenser gas (NO2 + others)

    F41 = fs.new_process_stream("F41", phase="G", mol={})  # mixed offgas (F32 + F3 + F40)

    F112 = fs.new_process_stream("F112", phase="L", mol={})  # post-coalescer aqueous to E103
    F113 = fs.new_process_stream("F113", phase="L", mol={})  # recovered organics from CF103

    F27 = fs.new_process_stream("F27", phase="S", mol={})  # U3O8 product
    F28 = fs.new_process_stream("F28", phase="G", mol={})  # O2 offgas from R202

    F2_R = fs.new_process_stream("F2_R", phase="L", mol={})     # recycle acid back to M102
    F138 = fs.new_process_stream("F138", phase="L", mol={})     # purge from F38

    F44 = fs.new_process_stream("F44", phase="G", mol={})       # treated gas

    F45 = fs.new_process_stream("F45", phase="L", mol={})       # KO water
    F46 = fs.new_process_stream("F46", phase="G", mol={})       # dry gas

    F145 = fs.new_process_stream("F145", phase="L", mol={})     # purge from KO water
    F45_R = fs.new_process_stream("F45_R", phase="L", mol={})   # recycle KO water
    F2_WM = fs.new_process_stream("F2_WM", phase="L", mol={})   # mixed water trim (F2_W + F45_R)

    F2_S = fs.new_process_stream("F2_S", phase="L", mol={})     # surge vessel outlet to dissolver

    F48 = fs.new_process_stream("F48", phase="L", mol={})       # absorber liquid (NH3 in water)
    F49 = fs.new_process_stream("F49", phase="G", mol={})       # treated gas

    # -------------------------------------------------------------------------
    # TSA system streams (two columns A/B)
    # -------------------------------------------------------------------------
    F49A = fs.new_process_stream("F49A", phase="G", mol={})
    F49B = fs.new_process_stream("F49B", phase="G", mol={})

    F50A = fs.new_process_stream("F50A", phase="G", mol={})
    F50B = fs.new_process_stream("F50B", phase="G", mol={})

    F51A = fs.new_process_stream("F51A", phase="G", mol={})
    F51B = fs.new_process_stream("F51B", phase="G", mol={})

    # -------------------------------------------------------------------------
    # (FIX) M105 MUST COME BEFORE M102 IN THE UNIT SEQUENCE
    # -------------------------------------------------------------------------
    M105 = Mixer("M105_WaterTrimMixer")
    M105.add_inlet("fresh_water", F2_W)
    M105.add_inlet("recycle_ko_water", F45_R)
    M105.add_outlet("out", F2_WM)
    _add_unit(M105)

    # -------------------------------------------------------------------------
    # UNITS
    # -------------------------------------------------------------------------
    M102 = Mixer("M102_AcidMixer")
    M102.add_inlet("recycle", F2_R)
    M102.add_inlet("acid_stock", F2_M)
    M102.add_inlet("water_trim", F2_WM)
    M102.add_outlet("out", F2)
    _add_unit(M102)

    V102 = BufferVessel("V102_AcidSurge")
    V102.add_inlet("in", condition_to_unit(F2, "V102_AcidSurge"))
    V102.add_outlet("out", F2_S)
    _add_unit(V102)

    R101 = DissolutionReactor(
        "R101_Dissolver",
        uo2_r1_frac=params.R1_UO2_R1_FRACTION,
        iodine_total_species=getattr(params, "R1_IODINE_TOTAL_SPECIES", "I"),
        iodine_gas_species=getattr(params, "R1_IODINE_GAS_SPECIES", "I_g"),
        iodine_aq_species=getattr(params, "R1_IODINE_AQ_SPECIES", "I_aq"),
        iodine_offgas_frac=params.R1_IODINE_OFFGAS_FRACTION,
        noble_gases=getattr(params, "R1_NOBLE_GASES", ("Xe", "Kr")),
        zircaloy_species="Zircaloy",
        zircaloy_basket_mol_s=params.R1_ZIRCALOY_BASKET_MOL_S,
        zircaloy_escape_fines_mol_s=params.R1_ZIRCALOY_ESCAPE_FINES_MOL_S,
        fines_settle_frac=params.R1_FINES_SETTLE_FRACTION,
        cleaning_water_port="cleaning_water",
        nox_no2_frac=params.R1_NOX_NO2_FRACTION,
    )
    R101.add_inlet("fuel", F1)
    R101.add_inlet("acid", condition_to_unit(F2_S, "R101_Dissolver"))
    R101.add_inlet("cleaning_water", F2_B)
    R101.add_outlet("offgas", F3)
    R101.add_outlet("aq", F5)
    R101.add_outlet("solids", F4)
    _add_unit(R101)

    V101 = BufferVessel("V101_Liquid_Buffer")
    V101.add_inlet("in", condition_to_unit(F5, "V101_Liquid_Buffer"))
    V101.add_outlet("out", F6)
    _add_unit(V101)

    M114 = Mixer("M114_OrgFeedMixer")
    M114.add_inlet("recycle", F15)
    M114.add_inlet("makeup", F16)
    M114.add_outlet("out", F7)
    _add_unit(M114)

    X101 = MultiStageCentrifugalContactorTBP(
        "X101_TBP_Extraction",
        D=params.X1_DISTRIBUTION_COEFFICIENTS,
        required_recovery=params.X1_REQUIRED_RECOVERY,
        nonextract_to_aq=getattr(params, "X1_NONEXTRACT_TO_AQ", [
            "H2O", "HNO3", "CsNO3", "Sr(NO3)2",
            "Nd(NO3)3", "Sm(NO3)3", "Eu(NO3)3", "Gd(NO3)3", "I_aq",
        ]),
        nonextract_to_org=getattr(params, "X1_NONEXTRACT_TO_ORG", ["Dodecane"]),
        N_max=params.X1_N_MAX,
        OA_max=params.X1_OA_MAX,
        N_balance_cap=getattr(params, "X1_N_BALANCE_CAP", 20),
    )
    X101.add_inlet("aq_in", condition_to_unit(F6, "X101_TBP_Extraction"))
    X101.add_inlet("org_in", condition_to_unit(F7, "X101_TBP_Extraction"))
    X101.add_outlet("org_out", F8)
    X101.add_outlet("aq_out", F17)
    _add_unit(X101)

    CF101 = CoalescingFilterTBP(
        "CF101_Coalescer_X101_to_E101",
        tbp_name=params.COAL_TBP_NAME,
        diluent_name=params.COAL_DILUENT_NAME,
        water_name=params.COAL_WATER_NAME,
        hno3_name=params.COAL_HNO3_NAME,
        water_molarity_pure=params.COAL_WATER_MOLARITY_PURE,
        removal_eff=params.COAL_REMOVAL_EFF,
        tbp_solubility_weak_mol_L=params.COAL_TBP_SOLUBILITY_WEAK_MOL_L,
        tbp_solubility_strong_mol_L=params.COAL_TBP_SOLUBILITY_STRONG_MOL_L,
        dil_solubility_weak_mol_L=params.COAL_DIL_SOLUBILITY_WEAK_MOL_L,
        dil_solubility_strong_mol_L=params.COAL_DIL_SOLUBILITY_STRONG_MOL_L,
        hno3_strong_M=params.COAL_HNO3_STRONG_M,
        assume_pure_water_for_weak_stream=False,
        print_diagnostics=params.COAL_PRINT_DIAGNOSTICS,
    )
    CF101.add_inlet("in", condition_to_unit(F17, "CF101_Coalescer_X101_to_E101"))
    CF101.add_outlet("aq_out", F101)
    CF101.add_outlet("org_recovered", F105)
    _add_unit(CF101)

    E101 = ContinuousEvaporatorAHA(
        "E101_Evaporator",
        T_set_K=params.EVAP_T_K,
        P_set_torr=params.EVAP_P_TORR,
        use_equilibrium=params.EVAP_USE_EQUILIBRIUM,
        frac_vap=params.E101_FRAC_VAP,
        min_liq_mol_s=params.E101_MIN_LIQ_MOL_S,
        latent_kJ_per_mol=params.EVAP_LATENT_KJ_PER_MOL,
        cp_liq_J_per_molK=params.EVAP_CP_LIQ_J_PER_MOLK,
        aha_name="AHA",
        aha_products=params.EVAP_AHA_PRODUCTS,
        aha_hno3_consumption=params.EVAP_AHA_HNO3_CONSUMPTION,
        coupled_by_water=params.EVAP_COUPLED_BY_WATER,
        water_name=params.EVAP_WATER_NAME,
        print_diagnostics=params.EVAP_PRINT_DIAGNOSTICS,
        print_every=params.EVAP_PRINT_EVERY,
    )
    E101.add_inlet("in", condition_to_unit(F101, "E101_Evaporator"))
    E101.add_outlet("liq", F19)
    E101.add_outlet("vap", F20)
    _add_unit(E101)

    X102 = MultiStageCentrifugalContactorAHA_Strip(
        "X102_AHA_Strip",
        D_org_over_aq=params.X2_DISTRIBUTION_COEFFICIENTS,
        required_recovery_to_aq=params.X2_REQUIRED_RECOVERY_TO_AQ,
        nontransfer_keep_in_org=getattr(params, "X2_NONTRANSFER_KEEP_IN_ORG", None),
        nontransfer_keep_in_aq=getattr(params, "X2_NONTRANSFER_KEEP_IN_AQ", None),
        N_max=params.X2_N_MAX,
        AO_max=params.X2_AO_MAX,
        N_balance_cap=getattr(params, "X2_N_BALANCE_CAP", 20),
    )
    X102.add_inlet("org_in", condition_to_unit(F8, "X102_AHA_Strip"))
    X102.add_inlet("aq_in", condition_to_unit(F9, "X102_AHA_Strip"))
    X102.add_outlet("org_out", F10)
    X102.add_outlet("aq_out", F18)
    _add_unit(X102)

    CF102 = CoalescingFilterTBP(
        "CF102_Coalescer_X102_to_E102",
        tbp_name=params.COAL_TBP_NAME,
        diluent_name=params.COAL_DILUENT_NAME,
        water_name=params.COAL_WATER_NAME,
        hno3_name=params.COAL_HNO3_NAME,
        water_molarity_pure=params.COAL_WATER_MOLARITY_PURE,
        removal_eff=params.COAL_REMOVAL_EFF,
        tbp_solubility_weak_mol_L=params.COAL_TBP_SOLUBILITY_WEAK_MOL_L,
        tbp_solubility_strong_mol_L=params.COAL_TBP_SOLUBILITY_STRONG_MOL_L,
        dil_solubility_weak_mol_L=params.COAL_DIL_SOLUBILITY_WEAK_MOL_L,
        dil_solubility_strong_mol_L=params.COAL_DIL_SOLUBILITY_STRONG_MOL_L,
        hno3_strong_M=params.COAL_HNO3_STRONG_M,
        assume_pure_water_for_weak_stream=params.COAL_ASSUME_PURE_WATER_FOR_X102_AQ,
        print_diagnostics=params.COAL_PRINT_DIAGNOSTICS,
    )
    CF102.add_inlet("in", condition_to_unit(F18, "CF102_Coalescer_X102_to_E102"))
    CF102.add_outlet("aq_out", F102)
    CF102.add_outlet("org_recovered", F106)
    _add_unit(CF102)

    E102 = ContinuousEvaporatorAHA(
        "E102_Evaporator",
        T_set_K=params.EVAP_T_K,
        P_set_torr=params.EVAP_P_TORR,
        use_equilibrium=params.EVAP_USE_EQUILIBRIUM,
        frac_vap=params.E102_FRAC_VAP,
        min_liq_mol_s=params.E102_MIN_LIQ_MOL_S,
        latent_kJ_per_mol=params.EVAP_LATENT_KJ_PER_MOL,
        cp_liq_J_per_molK=params.EVAP_CP_LIQ_J_PER_MOLK,
        aha_name="AHA",
        aha_products=params.EVAP_AHA_PRODUCTS,
        aha_hno3_consumption=params.EVAP_AHA_HNO3_CONSUMPTION,
        coupled_by_water=params.EVAP_COUPLED_BY_WATER,
        water_name=params.EVAP_WATER_NAME,
        print_diagnostics=params.EVAP_PRINT_DIAGNOSTICS,
        print_every=params.EVAP_PRINT_EVERY,
    )
    E102.add_inlet("in", condition_to_unit(F102, "E102_Evaporator"))
    E102.add_outlet("liq", F21)
    E102.add_outlet("vap", F22)
    _add_unit(E102)

    M201 = Mixer("M201_EvapOverheadMixer")
    M201.add_inlet("E101_vap", condition_to_unit(F20, "M201_EvapOverheadMixer"))
    M201.add_inlet("E102_vap", condition_to_unit(F22, "M201_EvapOverheadMixer"))
    M201.add_outlet("out", F23)
    _add_unit(M201)

    X103 = MultiStageCentrifugalContactorAHA_Strip(
        "X103_Final_Strip",
        D_org_over_aq=params.X3_DISTRIBUTION_COEFFICIENTS,
        required_recovery_to_aq=params.X3_REQUIRED_RECOVERY_TO_AQ,
        nontransfer_keep_in_org=getattr(params, "X3_NONTRANSFER_KEEP_IN_ORG", None),
        nontransfer_keep_in_aq=getattr(params, "X3_NONTRANSFER_KEEP_IN_AQ", None),
        N_max=params.X3_N_MAX,
        AO_max=params.X3_AO_MAX,
        N_balance_cap=getattr(params, "X3_N_BALANCE_CAP", 20),
    )
    X103.add_inlet("org_in", condition_to_unit(F10, "X103_Final_Strip"))
    X103.add_inlet("aq_in", condition_to_unit(F11, "X103_Final_Strip"))
    X103.add_outlet("org_out", F13)
    X103.add_outlet("aq_out", F12)
    _add_unit(X103)

    CF103 = CoalescingFilterTBP(
        "CF103_Coalescer_X103_to_E103",
        tbp_name=params.COAL_TBP_NAME,
        diluent_name=params.COAL_DILUENT_NAME,
        water_name=params.COAL_WATER_NAME,
        hno3_name=params.COAL_HNO3_NAME,
        water_molarity_pure=params.COAL_WATER_MOLARITY_PURE,
        removal_eff=params.COAL_X103_TO_E103_REMOVAL_EFF,
        tbp_solubility_weak_mol_L=params.COAL_X103_TO_E103_TBP_SOLUBILITY_MOL_L,
        tbp_solubility_strong_mol_L=params.COAL_X103_TO_E103_TBP_SOLUBILITY_MOL_L,
        dil_solubility_weak_mol_L=params.COAL_DIL_SOLUBILITY_WEAK_MOL_L,
        dil_solubility_strong_mol_L=params.COAL_DIL_SOLUBILITY_STRONG_MOL_L,
        hno3_strong_M=params.COAL_HNO3_STRONG_M,
        assume_pure_water_for_weak_stream=False,
        print_diagnostics=params.COAL_PRINT_DIAGNOSTICS,
    )
    CF103.add_inlet("in", condition_to_unit(F12, "CF103_Coalescer_X103_to_E103"))
    CF103.add_outlet("aq_out", F112)
    CF103.add_outlet("org_recovered", F113)
    _add_unit(CF103)

    E103 = SpeciesSplitter(
        "E103_Evaporator100C",
        frac_to_A=params.E103_FRAC_TO_VAP,
        default_to_A=0.0,
    )
    E103.add_inlet("in", condition_to_unit(F112, "E103_Evaporator100C"))
    E103.add_outlet("A", F25)
    E103.add_outlet("B", F24)
    _add_unit(E103)

    R201 = UO2NitrateToUO3Reactor(
        "R201_UO3_Calciner",
        uo2n_name="UO2(NO3)2",
        uo3_name="UO3",
        no2_name="NO2",
        o2_name="O2",
        h2o_name="H2O",
        hydrate_n=params.R201_HYDRATE_N,
        print_diagnostics=True,
    )
    R201.add_inlet("in", condition_to_unit(F24, "R201_UO3_Calciner"))
    R201.add_outlet("solids", F26)
    R201.add_outlet("offgas", F29)
    _add_unit(R201)

    R202 = UO3ToU3O8Reactor(
        "R202_U3O8_Converter",
        uo3_name="UO3",
        u3o8_name="U3O8",
        o2_name="O2",
        print_diagnostics=False,
    )
    R202.add_inlet("in", condition_to_unit(F26, "R202_U3O8_Converter"))
    R202.add_outlet("solids", F27)
    R202.add_outlet("offgas", F28)
    _add_unit(R202)

    M301 = Mixer("M301_HotVapMixer")
    M301.add_inlet("F23", condition_to_unit(F23, "M301_HotVapMixer"))
    M301.add_inlet("F25", condition_to_unit(F25, "M301_HotVapMixer"))
    M301.add_inlet("F29", condition_to_unit(F29, "M301_HotVapMixer"))
    M301.add_outlet("out", F30)
    _add_unit(M301)

    K301 = SpeciesSplitter(
        "K301_Condenser",
        frac_to_A={"H2O": 1.0, "HNO3": 1.0},  # A=liquid
        default_to_A=0.0,
    )
    K301.add_inlet("in", condition_to_unit(F30, "K301_Condenser"))
    K301.add_outlet("A", F31)
    K301.add_outlet("B", F32)
    _add_unit(K301)

    M302 = Mixer("M302_OffgasMixer")
    M302.add_inlet("cond_gas", condition_to_unit(F32, "M302_OffgasMixer"))
    M302.add_inlet("dissolver_offgas", condition_to_unit(F3, "M302_OffgasMixer"))
    M302.add_inlet("air", condition_to_unit(F40, "M302_OffgasMixer"))
    M302.add_outlet("out", F41)
    _add_unit(M302)

    # Acid recovery section streams
    F35 = fs.new_process_stream("F35", phase="L", mol={})
    F39 = fs.new_process_stream("F39", phase="L", mol={})
    F36 = fs.new_process_stream("F36", phase="L", mol={})
    F37 = fs.new_process_stream("F37", phase="G", mol={})
    F38 = fs.new_process_stream("F38", phase="L", mol={})
    F42 = fs.new_process_stream("F42", phase="G", mol={})

    F39.T, F39.p = params.D101_ABS_T_K, params.D101_ABS_P_PA
    F35.T, F35.p = params.D101_ABS_T_K, params.D101_ABS_P_PA
    F42.T, F42.p = params.D101_ABS_T_K, params.D101_ABS_P_PA

    D101 = NO2AbsorptionColumn(
        "D101_NO2_Absorber",
        capture_frac=params.D101_NO2_CAPTURE_FRAC,
        capture_frac_no=params.D101_NO_CAPTURE_FRAC,
        print_diagnostics=getattr(params, "D101_PRINT_DIAGNOSTICS", False),
    )
    D101.add_inlet("gas_in", condition_to_unit(F41, "D101_NO2_Absorber"))
    D101.add_inlet("liq_in", condition_to_unit(F39, "D101_NO2_Absorber"))
    D101.add_outlet("liq_out", F35)
    D101.add_outlet("gas_out", F42)
    _add_unit(D101)

    M303 = Mixer("M303_AcidToConcentratorMixer")
    M303.add_inlet("absorber_bottoms", condition_to_unit(F35, "M303_AcidToConcentratorMixer"))
    M303.add_inlet("condenser_liquid", condition_to_unit(F31, "M303_AcidToConcentratorMixer"))
    M303.add_outlet("out", F36)
    _add_unit(M303)

    E104 = NitricAcidConcentrator(
        "E104_HNO3_Concentrator",
        target_hno3_mass_frac=params.E104_TARGET_HNO3_MASS_FRAC,
        print_diagnostics=getattr(params, "E104_PRINT_DIAGNOSTICS", False),
    )
    E104.add_inlet("in", condition_to_unit(F36, "E104_HNO3_Concentrator"))
    E104.add_outlet("vap", F37)
    E104.add_outlet("liq", F38)
    _add_unit(E104)

    V104 = Splitter("V104_AcidPurge", frac_to_A=0.0)  # frac updated each iteration (after E104)
    V104.add_inlet("in", condition_to_unit(F38, "V104_AcidPurge"))
    V104.add_outlet("A", F138)  # purge
    V104.add_outlet("B", F2_R)  # recycle
    _add_unit(V104)

    R103 = NH3NOxReductionColumn(
        "R103_NOxReduction",
        removal_frac=params.R103_NOX_REMOVAL_FRAC,
        print_diagnostics=False,
    )
    R103.add_inlet("gas_in", condition_to_unit(F42, "R103_NOxReduction"))
    R103.add_inlet("nh3_in", condition_to_unit(F43, "R103_NOxReduction"))
    R103.add_outlet("out", F44)
    _add_unit(R103)

    KO101 = KnockOutPotWater(
        "KO101_KnockOutPot",
        water_name="H2O",
        print_diagnostics=False,
    )
    KO101.add_inlet("in", condition_to_unit(F44, "KO101_KnockOutPot"))
    KO101.add_outlet("liq", F45)
    KO101.add_outlet("gas", F46)
    _add_unit(KO101)

    D102 = AmmoniaAbsorptionColumn(
        "D102_NH3_Absorber",
        capture_frac=params.D102_NH3_CAPTURE_FRAC,
        nh3_name="NH3",
        print_diagnostics=False,
    )
    D102.add_inlet("gas_in", condition_to_unit(F46, "D102_NH3_Absorber"))
    D102.add_inlet("liq_in", condition_to_unit(F47, "D102_NH3_Absorber"))
    D102.add_outlet("liq_out", F48)
    D102.add_outlet("gas_out", F49)
    _add_unit(D102)

    V105 = Splitter("V105_KOWaterPurge", frac_to_A=float(params.KO_WATER_PURGE_FRAC))
    V105.add_inlet("in", condition_to_unit(F45, "V105_KOWaterPurge"))
    V105.add_outlet("A", F145)   # purge
    V105.add_outlet("B", F45_R)  # recycle
    _add_unit(V105)

    # TSA split + columns
    V201 = Splitter("V201_TSA_F49_Split", frac_to_A=1.0)
    V201.add_inlet("in", condition_to_unit(F49, "V201_TSA_F49_Split"))
    V201.add_outlet("A", F49A)
    V201.add_outlet("B", F49B)
    _add_unit(V201)

    active = str(getattr(params, "TSA_ACTIVE_COLUMN", "A")).upper().strip()
    V201.frac_to_A = 1.0 if active == "A" else 0.0

    def _emm17_cap_work_mol_m3() -> float:
        uptake_g_g = float(params.TSA_EMM17_UPTAKE_G_PER_G)
        rho_kg_m3 = float(params.TSA_EMM17_PACKING_DENSITY_G_CM3) * 1000.0
        mtz = float(params.TSA_MTZ_UTILIZATION)
        MW_I2 = float(fs.reg.mw[params.TSA_IODINE_NAME])
        return (uptake_g_g * 1000.0 / max(MW_I2, EPS)) * rho_kg_m3 * mtz

    TSA_cap = _emm17_cap_work_mol_m3()
    t_ads_s = float(params.TSA_T_ADS_H) * 3600.0
    t_reg_s = float(params.TSA_T_REGEN_H) * 3600.0

    TSA101A = IdealTSAColumnEMM17(
        "TSA101A_ColA",
        mode="adsorb" if active == "A" else "regen",
        iodine_name=params.TSA_IODINE_NAME,
        bed_volume_m3=1.0,
        cap_work_mol_per_m3=TSA_cap,
        t_ads_s=t_ads_s,
        t_regen_s=t_reg_s,
        regen_yI2_max=float(params.TSA_REGEN_MAX_Y_I2),
        print_diagnostics=False,
    )
    TSA101B = IdealTSAColumnEMM17(
        "TSA101B_ColB",
        mode="regen" if active == "A" else "adsorb",
        iodine_name=params.TSA_IODINE_NAME,
        bed_volume_m3=1.0,
        cap_work_mol_per_m3=TSA_cap,
        t_ads_s=t_ads_s,
        t_regen_s=t_reg_s,
        regen_yI2_max=float(params.TSA_REGEN_MAX_Y_I2),
        print_diagnostics=False,
    )

    TSA101A.add_inlet("gas_in", condition_to_unit(F49A, "TSA101A_ColA"))
    TSA101A.add_outlet("gas_out", F50A)
    TSA101B.add_inlet("gas_in", condition_to_unit(F49B, "TSA101B_ColB"))
    TSA101B.add_outlet("gas_out", F50B)

    TSA101A.add_inlet("regen_in", condition_to_unit(F52A, "TSA101A_ColA"))
    TSA101A.add_outlet("regen_out", F51A)
    TSA101B.add_inlet("regen_in", condition_to_unit(F52B, "TSA101B_ColB"))
    TSA101B.add_outlet("regen_out", F51B)

    _add_unit(TSA101A)
    _add_unit(TSA101B)

    # Solvent recovery + purge
    F107 = fs.new_process_stream("F107", phase="L", mol={})

    M132 = Mixer("M132_CoalescerOrgMixer")
    M132.add_inlet("cf101_rec", F105)
    M132.add_inlet("cf102_rec", F106)
    M132.add_inlet("cf103_rec", F113)
    M132.add_outlet("out", F107)
    _add_unit(M132)

    F107_to_F13 = condition_to_unit(F107, "M133_OrgRecoveryMixer")

    F13_mix = fs.new_process_stream("F13_mix", phase="L", mol={})
    M133 = Mixer("M133_OrgRecoveryMixer")
    M133.add_inlet("main_org", F13)
    M133.add_inlet("recovered_org", F107_to_F13)
    M133.add_outlet("out", F13_mix)
    _add_unit(M133)

    F13_T = fs.new_process_stream("F13_T", phase="L", mol={})
    HX133 = HeaterCooler("HX133_SolventCooldown", set_T=unit_cond["V133_SolventPurge"]["T"])
    HX133.add_inlet("in", F13_mix)
    HX133.add_outlet("out", F13_T)
    _add_unit(HX133)

    V133 = Splitter("V133_SolventPurge", frac_to_A=params.PURGE_FRACTION_SOLVENT_LOOP)
    V133.add_inlet("in", condition_to_unit(F13_T, "V133_SolventPurge"))
    V133.add_outlet("A", F14)
    V133.add_outlet("B", F15)
    _add_unit(V133)

    # Freeze after wiring
    reg.freeze()

    # -------------------------------------------------------------------------
    # Seed recycle solvent (F15) if needed
    # -------------------------------------------------------------------------
    if F15.total_molar_flow() <= EPS:
        M_tbp = fs.reg.mw["TBP"]
        M_dil = fs.reg.mw["Dodecane"]
        tbp_wt = params.X1_TBP_WT_FRACTION
        M_mix = 1.0 / (tbp_wt / M_tbp + (1.0 - tbp_wt) / M_dil)
        mdot_seed = 1.0 * M_mix
        F15.set("TBP", (tbp_wt * mdot_seed) / M_tbp)
        F15.set("Dodecane", ((1.0 - tbp_wt) * mdot_seed) / M_dil)

    # -------------------------------------------------------------------------
    # RECYCLE LOOP
    # -------------------------------------------------------------------------
    DEFAULT_TEARS = ["F15", "F2_R", "F45_R"]
    tears = list(tear_stream_names) if tear_stream_names is not None else DEFAULT_TEARS

    prev_snap = snapshot_streams(fs)
    converged = False

    last_OA1 = last_AO2 = last_AO3 = None
    last_Vdot2 = last_Vdot3 = None
    last_norg_req = None
    last_LG_water = None
    last_f38_purge = None
    last_f2r = None

    def run_until(target_unit_name: str) -> None:
        found = False
        for u in unit_sequence:
            run_unit(u)
            if u.name == target_unit_name:
                found = True
                break
        if not found:
            raise KeyError(f"Unit '{target_unit_name}' not found in unit_sequence")

    def run_range(start_unit_name: str, stop_before_unit_name: str) -> None:
        started = False
        for u in unit_sequence:
            if u.name == start_unit_name:
                started = True
            if started:
                if u.name == stop_before_unit_name:
                    break
                run_unit(u)

    def size_F43_for_R103():
        NO2 = fs.streams["F42"].get("NO2")
        NO = fs.streams["F42"].get("NO")
        eta = float(params.R103_NOX_REMOVAL_FRAC)
        dNO2 = eta * float(NO2)
        dNO = eta * float(NO)
        nh3_req = 2.0 * dNO2 + 1.0 * dNO
        nh3_feed = float(params.R103_NH3_EXCESS_FACTOR) * nh3_req

        F43.mol = {}
        if nh3_feed > EPS:
            F43.set("NH3", nh3_feed)
        F43.phase = "G"
        F43.T, F43.p = 298.15, 1e5

    def size_F47_for_D102():
        n_gas = fs.streams["F46"].total_molar_flow()
        n_water = float(params.D102_LG_WATER_PER_GAS) * float(n_gas)
        F47.mol = {}
        if n_water > EPS:
            F47.set("H2O", n_water)
        F47.phase = "L"
        F47.T, F47.p = params.D102_T_K, params.D102_P_PA

    def size_TSA_beds_and_regen_air():
        iodine = str(params.TSA_IODINE_NAME)
        y_max = float(params.TSA_REGEN_MAX_Y_I2)

        col_ads = TSA101A if TSA101A.mode == "adsorb" else TSA101B
        col_reg = TSA101A if TSA101A.mode == "regen" else TSA101B

        I_in = fs.streams[col_ads.inlets["gas_in"].name].get(iodine)

        cap = max(float(col_ads.cap_work_mol_per_m3), EPS)
        n_I_cycle = float(I_in) * max(float(col_ads.t_ads_s), 0.0)

        V_bed = float(params.TSA_BED_VOL_SF) * (n_I_cycle / cap) if n_I_cycle > EPS else 1e-6
        col_ads.bed_volume_m3 = V_bed
        col_reg.bed_volume_m3 = V_bed

        I_desorb = (float(I_in) * float(col_ads.t_ads_s) / max(float(col_ads.t_regen_s), EPS)) if I_in > EPS else 0.0
        n_air_min = (I_desorb * (1.0 - y_max) / max(y_max, EPS)) if I_desorb > EPS else 0.0

        def set_air(stream: Stream, n_air: float, add_iodine: float):
            stream.mol = {}
            if n_air > EPS:
                stream.set("O2", n_air * float(params.TSA_REGEN_O2_FRAC))
                stream.set("N2", n_air * float(params.TSA_REGEN_N2_FRAC))
            if add_iodine > EPS:
                stream.set(iodine, add_iodine)
            stream.phase = "G"
            stream.T, stream.p = float(params.TSA_REGEN_T_K), float(params.TSA_REGEN_P_PA)

        if col_reg is TSA101A:
            set_air(F52A, n_air_min, I_desorb)
            set_air(F52B, 0.0, 0.0)
        else:
            set_air(F52B, n_air_min, I_desorb)
            set_air(F52A, 0.0, 0.0)

    for it in range(1, max_iter + 1):
        # ---------------------------------------------------------------------
        # (A) SIZE ACID MAKEUPS USING CURRENT RECYCLE (F2_R) STATE
        #     (V104 will be updated later in this iteration after E104 runs)
        # ---------------------------------------------------------------------
        size_acid_stock_and_water_makeup_for_targets(
            recycle=fs.streams["F2_R"],
            acid_stock=F2_M,
            water_makeup=F2_W,
            target_hno3_mol_s=params.R1_F2_TARGET_HNO3_MOL_S,
            target_h2o_mol_s=params.R1_F2_TARGET_H2O_MOL_S,
            stock_water_per_hno3_mol=params.ACID_STOCK_WATER_PER_HNO3_MOL,
        )

        # ---------------------------------------------------------------------
        # (B) RUN FRONT-END TO GET AQUEOUS FLOW FOR X101
        # ---------------------------------------------------------------------
        run_until("V101_Liquid_Buffer")

        # ---------------------------------------------------------------------
        # (C) SIZE X101 ORGANIC MAKEUP FOR CURRENT AQUEOUS FEED
        # ---------------------------------------------------------------------
        _, OA_design = X101.design_N_and_OA()
        OA1 = max(float(OA_design), float(getattr(params, "X1_OA_MIN", 0.0)))
        A_tot = X101.inlets["aq_in"].total_molar_flow()
        n_org_req_total = OA1 * A_tot

        last_OA1 = float(OA1)
        last_norg_req = float(n_org_req_total)

        size_solvent_makeup_for_required_org_flow(
            reg=fs.reg,
            recycle=F15,
            makeup=F16,
            n_org_req_total=n_org_req_total,
            tbp_wt=params.X1_TBP_WT_FRACTION,
            tbp_name="TBP",
            dil_name="Dodecane",
        )

        # ---------------------------------------------------------------------
        # (D) RUN THROUGH E101 SO F8/F17 ARE CURRENT
        # ---------------------------------------------------------------------
        run_until("E101_Evaporator")

        # ---------------------------------------------------------------------
        # (E) SIZE X102 STRIP FEED (F9) BASED ON CURRENT ORG_IN (F8)
        # ---------------------------------------------------------------------
        _, AO2_design = X102.design_N_and_AO()
        AO2 = max(float(AO2_design), float(getattr(params, "X2_AO_MIN", 0.0)))

        O2_tot = X102.inlets["org_in"].total_molar_flow()
        A2_target = AO2 * O2_tot
        Vdot2 = A2_target / max(params.X2_STRIP_WATER_MOLARITY, EPS)

        last_AO2 = float(AO2)
        last_Vdot2 = float(Vdot2)

        set_aqueous_molarities(
            F9,
            Vdot_L_s=Vdot2,
            molarity=params.X2_STRIP_MOLARITIES,
            water_name="H2O",
            water_molarity_pure=params.X2_STRIP_WATER_MOLARITY,
        )

        # ---------------------------------------------------------------------
        # (F) RUN THROUGH M201 SO X102/E102/M201 ARE CURRENT
        # ---------------------------------------------------------------------
        run_until("M201_EvapOverheadMixer")

        # ---------------------------------------------------------------------
        # (G) ENSURE CONDITIONING BETWEEN M201 AND X103 EXECUTED, THEN SIZE F11 FOR AO3
        # ---------------------------------------------------------------------
        run_range("M201_EvapOverheadMixer", "X103_Final_Strip")

        _, AO3_design = X103.design_N_and_AO()
        AO3 = max(float(AO3_design), float(getattr(params, "X3_AO_MIN", 0.0)))

        O3_tot = X103.inlets["org_in"].total_molar_flow()
        A3_target = AO3 * O3_tot
        Vdot3 = A3_target / max(params.X3_STRIP_WATER_MOLARITY, EPS)

        last_AO3 = float(AO3)
        last_Vdot3 = float(Vdot3)

        set_aqueous_molarities(
            F11,
            Vdot_L_s=Vdot3,
            molarity=params.X3_STRIP_MOLARITIES,
            water_name="H2O",
            water_molarity_pure=params.X3_STRIP_WATER_MOLARITY,
        )

        # ---------------------------------------------------------------------
        # (H) RUN THROUGH OFFGAS MIXING SO F41 IS CURRENT
        # ---------------------------------------------------------------------
        run_until("M302_OffgasMixer")

        # (H1) SIZE ABSORBER WATER FEED F39 BASED ON CURRENT F41 TOTAL FLOW (molar L/G)
        n_gas = fs.streams["F41"].total_molar_flow()
        n_water = float(params.D101_LG_WATER_PER_GAS) * float(n_gas)
        last_LG_water = float(n_water)

        F39.mol = {}
        if n_water > EPS:
            F39.set("H2O", n_water)
        F39.phase = "L"
        F39.T, F39.p = params.D101_ABS_T_K, params.D101_ABS_P_PA

        # ---------------------------------------------------------------------
        # (I) RUN ABSORBER -> CONCENTRATOR (E104) TO GET F38, THEN UPDATE V104 PURGE
        # ---------------------------------------------------------------------
        run_until("E104_HNO3_Concentrator")

        f_purge = compute_required_purge_fraction_with_stock_acid(
            fs.streams["F38"],
            target_hno3_mol_s=params.R1_F2_TARGET_HNO3_MOL_S,
            target_h2o_mol_s=params.R1_F2_TARGET_H2O_MOL_S,
            stock_water_per_hno3_mol=params.ACID_STOCK_WATER_PER_HNO3_MOL,
            f_min=0.0,
            f_max=0.999,
        )
        V104.frac_to_A = float(f_purge)
        run_unit(V104)  # update F2_R for next iteration sizing

        last_f38_purge = float(f_purge)
        last_f2r = float(fs.streams["F2_R"].total_molar_flow())

        # ---------------------------------------------------------------------
        # (J) NOW RUN NOx REDUCTION / KO / NH3 ABS / KO WATER PURGE (updates F45_R)
        # ---------------------------------------------------------------------
        run_until("D101_NO2_Absorber")
        size_F43_for_R103()
        run_unit(R103)
        run_until("KO101_KnockOutPot")
        size_F47_for_D102()
        run_unit(D102)
        run_unit(V105)  # update F45_R for next iteration’s M105->M102

        # ---------------------------------------------------------------------
        # (K) CONTINUE THROUGH SOLVENT LOOP TO UPDATE F15 (tear)
        # ---------------------------------------------------------------------
        run_until("V133_SolventPurge")

        # ---------------------------------------------------------------------
        # (L) TSA sizing + execution (downstream polish)
        # ---------------------------------------------------------------------
        run_unit(V201)
        size_TSA_beds_and_regen_air()
        run_unit(TSA101A)
        run_unit(TSA101B)

        # ---------------------------------------------------------------------
        # (M) Convergence on ALL streams (as you originally did)
        # ---------------------------------------------------------------------
        curr_snap = snapshot_streams(fs)
        err = max_rel_change_all_streams(prev_snap, curr_snap, include_Tp=include_Tp)
        if err < tol:
            converged = True
            break

        # Relax tears
        if tears:
            for sname in tears:
                if sname not in fs.streams:
                    raise KeyError(f"Tear stream '{sname}' not found in flowsheet.streams")
                relax_stream(sname, prev_snap, curr_snap, alpha=relax)
            prev_snap = snapshot_streams(fs)
        else:
            prev_snap = curr_snap

    if not converged:
        raise RuntimeError(
            f"Recycle did not converge after {max_iter} iterations "
            f"(tol={tol}, relax={relax}, tears={tears})."
        )

    # -------------------------------------------------------------------------
    # Post-convergence: recompute required N at ACTUAL OA/AO used
    # -------------------------------------------------------------------------
    def _xN_over_x0(E: float, N: int) -> float:
        if abs(E - 1.0) < 1e-12:
            return 1.0 / (N + 1.0)
        s = (E ** (N + 1) - 1.0) / (E - 1.0)
        return 1.0 / s

    def _min_N_for_target_recovery(E: float, recovery: float, *, N_max: int = 200) -> int:
        if recovery <= 0.0:
            return 1
        if recovery >= 1.0:
            return N_max
        target_xratio = max(1.0 - recovery, 0.0)
        for N in range(1, N_max + 1):
            if _xN_over_x0(E, N) <= target_xratio:
                return N
        return N_max

    def _as_fraction(v: float) -> float:
        v = float(v)
        return v / 100.0 if v > 1.0 else v

    def recompute_N_for_X101_at_OA_used() -> tuple[int, dict[str, int]]:
        OA = float(getattr(X101, "OA_used", 0.0))
        D = getattr(X101, "D", {})
        req = getattr(X101, "required_recovery", {})
        per_species: dict[str, int] = {}
        N_req = 1
        for sp, r in req.items():
            r = _as_fraction(r)
            if sp not in D:
                continue
            Di = float(D[sp])
            if Di <= 0.0 or OA <= 0.0:
                n_sp = 200
            else:
                E = Di * OA
                n_sp = _min_N_for_target_recovery(E, r, N_max=max(getattr(X101, "N_max", 50), 200))
            per_species[sp] = int(n_sp)
            N_req = max(N_req, int(n_sp))
        return int(N_req), per_species

    def recompute_N_for_strip_at_AO_used(unit) -> tuple[int, dict[str, int]]:
        AO = float(getattr(unit, "AO_used", 0.0))
        D = getattr(unit, "D", {})
        req = getattr(unit, "required_recovery_to_aq", {})
        per_species: dict[str, int] = {}
        N_req = 1
        for sp, r in req.items():
            r = _as_fraction(r)
            if sp not in D:
                continue
            Di = float(D[sp])
            if Di <= 0.0:
                n_sp = 1
            elif AO <= 0.0:
                n_sp = 200
            else:
                E = AO / Di
                n_sp = _min_N_for_target_recovery(E, r, N_max=max(getattr(unit, "N_max", 50), 200))
            per_species[sp] = int(n_sp)
            N_req = max(N_req, int(n_sp))
        return int(N_req), per_species

    X101_N_req_at_OA_used, X101_N_req_by_species = recompute_N_for_X101_at_OA_used()
    X102_N_req_at_AO_used, X102_N_req_by_species = recompute_N_for_strip_at_AO_used(X102)
    X103_N_req_at_AO_used, X103_N_req_by_species = recompute_N_for_strip_at_AO_used(X103)

    print("\n--- Stage back-calculation at actual ratios ---")
    print(f"X101: OA_used={getattr(X101, 'OA_used', None)} -> N_required={X101_N_req_at_OA_used} (spec: {X101_N_req_by_species})")
    print(f"X102: AO_used={getattr(X102, 'AO_used', None)} -> N_required={X102_N_req_at_AO_used} (spec: {X102_N_req_by_species})")
    print(f"X103: AO_used={getattr(X103, 'AO_used', None)} -> N_required={X103_N_req_at_AO_used} (spec: {X103_N_req_by_species})")
    print("---------------------------------------------\n")

    fs.design = {
        "converged": bool(converged),
        "iterations": int(it),
        "tol": float(tol),
        "relax": float(relax),
        "tear_streams": list(tears),

        "X101_OA_target_for_sizing": float(last_OA1) if last_OA1 is not None else None,
        "X101_n_org_req_total": float(last_norg_req) if last_norg_req is not None else None,
        "X102_AO_target_for_sizing": float(last_AO2) if last_AO2 is not None else None,
        "X102_Vdot_L_s": float(last_Vdot2) if last_Vdot2 is not None else None,
        "X103_AO_target_for_sizing": float(last_AO3) if last_AO3 is not None else None,
        "X103_Vdot_L_s": float(last_Vdot3) if last_Vdot3 is not None else None,

        "D101_water_feed_mol_s": float(last_LG_water) if last_LG_water is not None else None,
        "D101_LG_water_per_gas": float(params.D101_LG_WATER_PER_GAS),

        "X101_N_used": int(X101.N_used),
        "X101_OA_design": float(X101.OA_design),
        "X101_OA_used": float(X101.OA_used),

        "X102_N_used": int(X102.N_used),
        "X102_AO_design": float(X102.AO_design),
        "X102_AO_used": float(X102.AO_used),

        "X103_N_used": int(X103.N_used),
        "X103_AO_design": float(X103.AO_design),
        "X103_AO_used": float(X103.AO_used),

        "X101_N_required_at_OA_used": int(X101_N_req_at_OA_used),
        "X101_N_required_at_OA_used_by_species": dict(X101_N_req_by_species),

        "X102_N_required_at_AO_used": int(X102_N_req_at_AO_used),
        "X102_N_required_at_AO_used_by_species": dict(X102_N_req_by_species),

        "X103_N_required_at_AO_used": int(X103_N_req_at_AO_used),
        "X103_N_required_at_AO_used_by_species": dict(X103_N_req_by_species),

        "V104_F38_purge_frac": float(last_f38_purge) if last_f38_purge is not None else None,
        "F2_R_total_mol_s": float(last_f2r) if last_f2r is not None else None,
    }

    print("\n--- Design / sizing decisions (final iteration) ---")
    print(f"converged: {fs.design['converged']}")
    print(f"iterations: {fs.design['iterations']}")
    print(f"tear_streams: {fs.design['tear_streams']}")
    print(f"X101: N_used={fs.design['X101_N_used']}, OA_design={fs.design['X101_OA_design']:.6g}, "
          f"OA_used={fs.design['X101_OA_used']:.6g}, OA_target_for_sizing={fs.design['X101_OA_target_for_sizing']:.6g}, "
          f"n_org_req_total={fs.design['X101_n_org_req_total']:.6g} mol/s")
    print(f"X102: N_used={fs.design['X102_N_used']}, AO_design={fs.design['X102_AO_design']:.6g}, "
          f"AO_used={fs.design['X102_AO_used']:.6g}, AO_target_for_sizing={fs.design['X102_AO_target_for_sizing']:.6g}, "
          f"Vdot={fs.design['X102_Vdot_L_s']:.6g} L/s")
    print(f"X103: N_used={fs.design['X103_N_used']}, AO_design={fs.design['X103_AO_design']:.6g}, "
          f"AO_used={fs.design['X103_AO_used']:.6g}, AO_target_for_sizing={fs.design['X103_AO_target_for_sizing']:.6g}, "
          f"Vdot={fs.design['X103_Vdot_L_s']:.6g} L/s")
    print(f"D101: water_feed={fs.design['D101_water_feed_mol_s']:.6g} mol/s "
          f"(L/G={fs.design['D101_LG_water_per_gas']:.6g} mol/mol)")
    print(f"V104: purge_frac={fs.design['V104_F38_purge_frac']:.6g}, F2_R={fs.design['F2_R_total_mol_s']:.6g} mol/s")
    print("----------------------------------------------\n")

    return fs

