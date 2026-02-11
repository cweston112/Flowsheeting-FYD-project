from __future__ import annotations

from typing import List, Iterable, Optional
import numpy as np

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
    }
    for name, formula in formula_species.items():
        reg.add_species(name, formula=formula)

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
    Dissolver + TBP extraction (X101) + AHA strip (X102) + final strip (X103) + purge/recycle (V133)
    - Recycle-capable execution loop.
    - Organic feed to X101 is mixed from recycle (F15) + makeup (F16).
      Makeup is sized each iteration to meet X101 recovery/OA design requirements.
    - Tear streams: pass e.g. tear_stream_names=["F15"] (default).
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

    def run_unit(u) -> None:
        u.apply()
        stamp_outlets_to_unit_conditions(u, unit_cond)

    def relax_stream(name: str, prev_snap, curr_snap, alpha: float) -> None:
        s = fs.streams[name]
        x0, _, _ = prev_snap[name]
        x1, _, _ = curr_snap[name]
        s.from_dense((1.0 - alpha) * x0 + alpha * x1)

    # -------------------------------------------------------------------------
    # T/P INTEGRATION
    # -------------------------------------------------------------------------
    def condition_to_unit(src: Stream, unit_name: str) -> Stream:
        cond = unit_cond.get(unit_name, {})

        # Predict src T/p from producing unit (if known)
        src_T = src.T
        src_p = src.p
        if src.producer is not None:
            prod_cond = unit_cond.get(src.producer, {})
            if "T" in prod_cond:
                src_T = float(prod_cond["T"])
            if "p" in prod_cond:
                src_p = float(prod_cond["p"])

        class _View:
            def __init__(self, T, p):
                self.T = T
                self.p = p

        view = _View(src_T, src_p)
        out = src

        need_T = _needs_T(view, cond.get("T", None))
        need_p = _needs_p(view, cond.get("p", None))

        if need_T:
            sT = fs.new_process_stream(f"{out.name}_T", phase=out.phase, mol={})
            sT.T = float(cond["T"])
            sT.p = src_p
            sT.phase = out.phase

            u = HeaterCooler(f"HX__{out.name}__to__{unit_name}", set_T=float(cond["T"]))
            u.add_inlet("in", out)
            u.add_outlet("out", sT)

            fs.add_unit(u)
            unit_sequence.append(u)

            out = sT
            view = _View(sT.T, sT.p)

        if need_p:
            sp = fs.new_process_stream(f"{out.name}_p", phase=out.phase, mol={})
            sp.p = float(cond["p"])
            sp.T = out.T
            sp.phase = out.phase

            u = PressureChanger(f"PC__{out.name}__to__{unit_name}", set_p=float(cond["p"]))
            u.add_inlet("in", out)
            u.add_outlet("out", sp)

            fs.add_unit(u)
            unit_sequence.append(u)

            out = sp

        return out

    # -------------------------------------------------------------------------
    # EXOGENOUS STREAMS
    # -------------------------------------------------------------------------
    F1 = fs.new_exogenous_stream("F1", mol=params.FEED_COMPOSITION)  # fuel feed (solid)

    F2_M = fs.new_exogenous_stream("F2_M", mol={})  # acid makeup
    F2_M.set("HNO3", params.R1_F2_TARGET_HNO3_MOL_S)
    F2_M.set("H2O", params.R1_F2_TARGET_H2O_MOL_S)

    F2_B = fs.new_exogenous_stream("F2_B", mol={})  # cleaning water
    F2_B.set("H2O", params.F2B_flow_rate)

    F16 = fs.new_exogenous_stream("F16", mol={})  # organic solvent makeup (TBP+diluent)

    F9 = fs.new_exogenous_stream("F9", mol={}) # AHA feed
    set_aqueous_molarities(
        F9,
        Vdot_L_s=params.X2_STRIP_VDOT_L_S,
        molarity=params.X2_STRIP_MOLARITIES,
        water_name="H2O",
        water_molarity_pure=params.X2_STRIP_WATER_MOLARITY,
    )

    F11 = fs.new_exogenous_stream("F11", mol={}) # Dilute nitric acid feed
    set_aqueous_molarities(
        F11,
        Vdot_L_s=params.X3_STRIP_VDOT_L_S,
        molarity=params.X3_STRIP_MOLARITIES,
        water_name="H2O",
        water_molarity_pure=params.X3_STRIP_WATER_MOLARITY,
    )

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

    F19 = fs.new_process_stream("F19", phase="L", mol={})  # T101 liquor
    F20 = fs.new_process_stream("F20", phase="G", mol={})  # T101 vapour

    F21 = fs.new_process_stream("F21", phase="L", mol={})  # T102 liquor
    F22 = fs.new_process_stream("F22", phase="G", mol={})  # T102 vapour

    F23 = fs.new_process_stream("F23", phase="G", mol={})  # combined vapour


    # -------------------------------------------------------------------------
    # UNITS
    # -------------------------------------------------------------------------
    M102 = Mixer("M102_AcidMixer")
    M102.add_inlet("makeup", F2_M)
    M102.add_outlet("out", F2)
    _add_unit(M102)

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
    R101.add_inlet("fuel", condition_to_unit(F1, "R101_Dissolver"))
    R101.add_inlet("acid", condition_to_unit(F2, "R101_Dissolver"))
    R101.add_inlet("cleaning_water", condition_to_unit(F2_B, "R101_Dissolver"))
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
        nonextract_to_org=getattr(params, "X1_NONEXTRACT_TO_ORG", ["TBP", "Dodecane"]),
        N_max=params.X1_N_MAX,
        OA_max=params.X1_OA_MAX,
        N_balance_cap=getattr(params, "X1_N_BALANCE_CAP", 20),
    )
    X101.add_inlet("aq_in", condition_to_unit(F6, "X101_TBP_Extraction"))
    X101.add_inlet("org_in", condition_to_unit(F7, "X101_TBP_Extraction"))
    X101.add_outlet("org_out", F8)
    X101.add_outlet("aq_out", F17)
    _add_unit(X101)

    T101 = ContinuousEvaporatorAHA(
        "T101_Evaporator",
        T_set_K=params.EVAP_T_K,
        frac_vap=params.T101_FRAC_VAP,
        min_liq_mol_s=params.T101_MIN_LIQ_MOL_S,
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
    T101.add_inlet("in", condition_to_unit(F17, "T101_Evaporator"))
    T101.add_outlet("liq", F19)
    T101.add_outlet("vap", F20)
    _add_unit(T101)


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

    T102 = ContinuousEvaporatorAHA(
        "T102_Evaporator",
        T_set_K=params.EVAP_T_K,
        frac_vap=params.T102_FRAC_VAP,
        min_liq_mol_s=params.T102_MIN_LIQ_MOL_S,
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
    T102.add_inlet("in", condition_to_unit(F18, "T102_Evaporator"))
    T102.add_outlet("liq", F21)
    T102.add_outlet("vap", F22)
    _add_unit(T102)

    M201 = Mixer("M201_EvapOverheadMixer")
    M201.add_inlet("T101_vap", condition_to_unit(F20, "M201_EvapOverheadMixer"))
    M201.add_inlet("T102_vap", condition_to_unit(F22, "M201_EvapOverheadMixer"))
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
    X103.add_inlet("aq_in", F11)
    X103.add_outlet("org_out", F13)
    X103.add_outlet("aq_out", F12)
    _add_unit(X103)

    V133 = Splitter("V133_SolventPurge", frac_to_A=params.PURGE_FRACTION_SOLVENT_LOOP)
    V133.add_inlet("in", condition_to_unit(F13, "V133_SolventPurge"))
    V133.add_outlet("A", F14)  # purge
    V133.add_outlet("B", F15)  # recycle
    _add_unit(V133)

    # Freeze after wiring
    reg.freeze()

    # -------------------------------------------------------------------------
    # Seed recycle
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
    DEFAULT_TEARS = ["F15"]
    tears = list(tear_stream_names) if tear_stream_names is not None else DEFAULT_TEARS

    prev_snap = snapshot_streams(fs)
    converged = False

    last_OA1 = last_AO2 = last_AO3 = None
    last_Vdot2 = last_Vdot3 = None
    last_norg_req = None

    def run_until(target_unit_name: str) -> None:
        """
        Run units in unit_sequence in order until and including target_unit_name.
        Raises if target not found.
        """
        found = False
        for u in unit_sequence:
            run_unit(u)
            if u.name == target_unit_name:
                found = True
                break
        if not found:
            raise KeyError(f"Unit '{target_unit_name}' not found in unit_sequence")

    def run_range(start_unit_name: str, stop_before_unit_name: str) -> None:
        """
        Run units from start_unit_name (inclusive) up to but not including stop_before_unit_name.
        Useful for running conditioning between two known units.
        """
        started = False
        for u in unit_sequence:
            if u.name == start_unit_name:
                started = True
            if started:
                if u.name == stop_before_unit_name:
                    break
                run_unit(u)

    def get_conditioned_name(base_name: str, unit_name: str) -> str:
        """
        Given a base stream name and a unit, return the conditioned stream name if it exists.
        We create streams like '{base}_T' and '{base}_p' in condition_to_unit(), possibly both.
        The final 'out' is either base, base_T, base_p, or base_T_p (depending on need_T/need_p).
        We detect the actually-wired inlet stream from the unit itself whenever possible.
        """
        # Prefer looking up the actual inlet stream name from the unit ports:
        u = fs.units[unit_name]
        # Guess port: most units use 'in' or 'aq_in'/'org_in'
        for port in ("in", "aq_in", "org_in"):
            if port in u.inlets and u.inlets[port].name.startswith(base_name):
                return u.inlets[port].name
        return base_name  # fallback

    # Identify key conditioning units by name pattern so we can run them after resizing
    def run_conditioning_for_stream_to_unit(stream_name: str, unit_name: str) -> None:
        """
        Run HX/PC units that condition a specific stream into a specific unit.
        This is robust because those units were inserted into unit_sequence by condition_to_unit()
        and are named 'HX__{stream}__to__{unit}' / 'PC__{stream}__to__{unit}'.
        """
        hx_name = f"HX__{stream_name}__to__{unit_name}"
        pc_name = f"PC__{stream_name}__to__{unit_name}"

        # Run in the order they appear in unit_sequence (HX then PC if both exist)
        for u in unit_sequence:
            if u.name == hx_name or u.name == pc_name:
                run_unit(u)

    for it in range(1, max_iter + 1):

        # -----------------------------------------------------------------
        # (1) Run up to V101 so aqueous to X101 exists (F6 conditioned)
        # -----------------------------------------------------------------
        run_until("V101_Liquid_Buffer")

        # -----------------------------------------------------------------
        # (2) Size X101 organic makeup for current aqueous feed to X101
        #     (X101.design uses its inlet streams; make sure conditioning to X101 is up to date)
        # -----------------------------------------------------------------
        # Ensure conditioning on F6->X101 and F7->X101 will be updated once we run M114 later,
        # but X101.design_N_and_OA() does not depend on org_in composition, only A_tot.
        # (A_tot is from aq_in; it is already current after run_until(V101).)
        _, OA_design = X101.design_N_and_OA()

        # Enforce minimum O/A if provided (prevents unrealistically low organic flow)
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

        # -----------------------------------------------------------------
        # (3) Run through T101 so F8/F17/F19/F20 are current
        #     This executes M114, X101, any HX/PC created for their inlets, and T101.
        # -----------------------------------------------------------------
        run_until("T101_Evaporator")

        # -----------------------------------------------------------------
        # (4) Size X102 strip feed (F9) to hit AO2 for current organic feed (F8)
        # -----------------------------------------------------------------
        try:
            _, AO2_design = X102.design_N_and_AO()
        except ValueError as e:
            raise RuntimeError(
                f"X102 design failed: {e}\n"
                f"Try increasing params.X2_AO_MAX or relaxing X2 required recovery.\n"
                f"Current X2 D keys={list(params.X2_DISTRIBUTION_COEFFICIENTS.keys())}, "
                f"required={params.X2_REQUIRED_RECOVERY_TO_AQ}"
            ) from e

        AO2 = max(float(AO2_design), float(getattr(params, "X2_AO_MIN", 0.0)))

        O2_tot = X102.inlets["org_in"].total_molar_flow()
        A2_target = AO2 * O2_tot
        Vdot2 = A2_target / max(params.X2_STRIP_WATER_MOLARITY, EPS)

        # Optional: enforce a minimum volumetric strip flow if you want
        if hasattr(params, "X2_STRIP_VDOT_MIN_L_S"):
            Vdot2 = max(float(Vdot2), float(params.X2_STRIP_VDOT_MIN_L_S))

        last_AO2 = float(AO2)
        last_Vdot2 = float(Vdot2)

        set_aqueous_molarities(
            F9,
            Vdot_L_s=Vdot2,
            molarity=params.X2_STRIP_MOLARITIES,
            water_name="H2O",
            water_molarity_pure=params.X2_STRIP_WATER_MOLARITY,
        )

        # IMPORTANT: update the conditioned copy of F9 that feeds X102 (if any)
        run_conditioning_for_stream_to_unit("F9", "X102_AHA_Strip")

        # -----------------------------------------------------------------
        # (5) Run through M201 so X102, T102, M201 are current (F10/F18/F22/F23)
        # -----------------------------------------------------------------
        run_until("M201_EvapOverheadMixer")

        # -----------------------------------------------------------------
        # (6) Size X103 strip feed (F11) to hit AO3 for current organic feed to X103
        #     CRITICAL: ensure conditioning between X102 and X103 has been run so X103.org_in is current.
        # -----------------------------------------------------------------
        # Run any conditioning units inserted after M201 and before X103 (e.g., HX__F10__to__X103..., PC__...)
        run_range("M201_EvapOverheadMixer", "X103_Final_Strip")

        try:
            _, AO3_design = X103.design_N_and_AO()
        except ValueError as e:
            raise RuntimeError(
                f"X103 design failed: {e}\n"
                f"Try increasing params.X3_AO_MAX or relaxing X3 required recovery.\n"
                f"Current X3 D={getattr(params, 'X3_DISTRIBUTION_COEFFICIENTS', None)}, "
                f"required={getattr(params, 'X3_REQUIRED_RECOVERY_TO_AQ', None)}"
            ) from e

        # Enforce minimum A/O if provided (prevents unrealistically low strip flow)
        AO3 = max(float(AO3_design), float(getattr(params, "X3_AO_MIN", 0.0)))

        O3_tot = X103.inlets["org_in"].total_molar_flow()
        A3_target = AO3 * O3_tot
        Vdot3 = A3_target / max(params.X3_STRIP_WATER_MOLARITY, EPS)

        # Optional: enforce minimum strip flow
        if hasattr(params, "X3_STRIP_VDOT_MIN_L_S"):
            Vdot3 = max(float(Vdot3), float(params.X3_STRIP_VDOT_MIN_L_S))

        last_AO3 = float(AO3)
        last_Vdot3 = float(Vdot3)

        set_aqueous_molarities(
            F11,
            Vdot_L_s=Vdot3,
            molarity=params.X3_STRIP_MOLARITIES,
            water_name="H2O",
            water_molarity_pure=params.X3_STRIP_WATER_MOLARITY,
        )

        # If X103 aq_in is a conditioned copy of F11, update it now.
        # If you wired X103.aq_in directly to F11, this is harmless (units won't exist).
        run_conditioning_for_stream_to_unit("F11", "X103_Final_Strip")

        # -----------------------------------------------------------------
        # (7) Run from X103 through V133 so recycle is updated
        # -----------------------------------------------------------------
        run_until("V133_SolventPurge")

        # -----------------------------------------------------------------
        # (8) Convergence
        # -----------------------------------------------------------------
        curr_snap = snapshot_streams(fs)
        err = max_rel_change_all_streams(prev_snap, curr_snap, include_Tp=include_Tp)
        if err < tol:
            converged = True
            break

        # -----------------------------------------------------------------
        # (9) Relax tears
        # -----------------------------------------------------------------
        if tears:
            for sname in tears:
                if sname not in fs.streams:
                    raise KeyError(f"Tear stream '{sname}' not found in flowsheet.streams")
                relax_stream(sname, prev_snap, curr_snap, alpha=relax)
            prev_snap = snapshot_streams(fs)
        else:
            prev_snap = curr_snap

    # -------------------------------------------------------------------------
    # Post-convergence: recompute required N at ACTUAL OA/AO used
    # -------------------------------------------------------------------------
    def _xN_over_x0(E: float, N: int) -> float:
        if abs(E - 1.0) < 1e-12:
            return 1.0 / (N + 1.0)
        s = (E ** (N + 1) - 1.0) / (E - 1.0)
        return 1.0 / s

    def _min_N_for_target_recovery(E: float, recovery: float, *, N_max: int = 200) -> int:
        # Find smallest N such that (1 - xN/x0) >= recovery
        # Handle corner cases
        if recovery <= 0.0:
            return 1
        if recovery >= 1.0:
            # effectively impossible for finite N unless E huge; return cap
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
        """
        X101 extraction: E = D * OA_used (D = org/aq, OA = O/A).
        Recovery spec is to ORGANIC.
        """
        OA = float(getattr(X101, "OA_used", 0.0))
        D = getattr(X101, "D", {})  # should exist in your TBP class
        req = getattr(X101, "required_recovery", {})  # to organic

        per_species: dict[str, int] = {}
        N_req = 1

        for sp, r in req.items():
            r = _as_fraction(r)
            if sp not in D:
                continue  # or raise; up to you
            Di = float(D[sp])

            if Di <= 0.0 or OA <= 0.0:
                # If E ~ 0 then essentially no extraction -> requires "infinite" stages
                n_sp = 200
            else:
                E = Di * OA
                n_sp = _min_N_for_target_recovery(E, r, N_max=max(getattr(X101, "N_max", 50), 200))

            per_species[sp] = int(n_sp)
            N_req = max(N_req, int(n_sp))

        return int(N_req), per_species

    def recompute_N_for_strip_at_AO_used(unit) -> tuple[int, dict[str, int]]:
        """
        Stripping: E = (AO_used)/D  (D = org/aq, AO = A/O)
        Recovery spec is to AQUEOUS.
        """
        AO = float(getattr(unit, "AO_used", 0.0))
        D = getattr(unit, "D", {})
        req = getattr(unit, "required_recovery_to_aq", {})

        per_species: dict[str, int] = {}
        N_req = 1

        for sp, r in req.items():
            r = _as_fraction(r)
            if sp not in D:
                continue  # or raise
            Di = float(D[sp])

            if Di <= 0.0:
                # D ~ 0 => overwhelmingly aqueous => stripping essentially complete
                n_sp = 1
            elif AO <= 0.0:
                n_sp = 200
            else:
                E = AO / Di
                n_sp = _min_N_for_target_recovery(E, r, N_max=max(getattr(unit, "N_max", 50), 200))

            per_species[sp] = int(n_sp)
            N_req = max(N_req, int(n_sp))

        return int(N_req), per_species

    # Compute required stages at actual ratios (after final iteration)
    X101_N_req_at_OA_used, X101_N_req_by_species = recompute_N_for_X101_at_OA_used()
    X102_N_req_at_AO_used, X102_N_req_by_species = recompute_N_for_strip_at_AO_used(X102)
    X103_N_req_at_AO_used, X103_N_req_by_species = recompute_N_for_strip_at_AO_used(X103)

    # Optional diagnostic print
    print("\n--- Stage back-calculation at actual ratios ---")
    print(f"X101: OA_used={getattr(X101,'OA_used',None)} -> N_required={X101_N_req_at_OA_used} (spec species: {X101_N_req_by_species})")
    print(f"X102: AO_used={getattr(X102,'AO_used',None)} -> N_required={X102_N_req_at_AO_used} (spec species: {X102_N_req_by_species})")
    print(f"X103: AO_used={getattr(X103,'AO_used',None)} -> N_required={X103_N_req_at_AO_used} (spec species: {X103_N_req_by_species})")
    print("---------------------------------------------\n")

    # -------------------------------------------------------------------------
    # Report / store chosen parameters
    # -------------------------------------------------------------------------
    fs.design = {
        "converged": True,
        "iterations": it,
        "tol": float(tol),
        "relax": float(relax),
        "tear_streams": list(tears),

        "X101_OA_target_for_sizing": float(last_OA1) if last_OA1 is not None else None,
        "X101_n_org_req_total": float(last_norg_req) if last_norg_req is not None else None,
        "X102_AO_target_for_sizing": float(last_AO2) if last_AO2 is not None else None,
        "X102_Vdot_L_s": float(last_Vdot2) if last_Vdot2 is not None else None,
        "X103_AO_target_for_sizing": float(last_AO3) if last_AO3 is not None else None,
        "X103_Vdot_L_s": float(last_Vdot3) if last_Vdot3 is not None else None,

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
    print("----------------------------------------------\n")

    return fs
