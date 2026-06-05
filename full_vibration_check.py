"""
SCI P354 Section 7 — Full Composite Floor Vibration Check
==========================================================

Standalone script for verifying the vibration performance of a composite
steel-deck slab floor system using the simplified method in SCI P354 §7.

Usage
-----
1. Set the slab and bay geometry in the USER INPUTS section.
2. Run the script.  Results are printed to the console.

References
----------
    SCI P354 — Design of Floors for Vibration: A New Approach (2009 ed.)
    EN 1994-1-1: 2004 — Composite slabs
    BS 6472-1: 2008 — Evaluation of human exposure to vibration in buildings
"""

import math
import os
import sys

# ==============================================================================
# KNOWN DECK PROFILES
# ==============================================================================
# Typical trapezoidal and re-entrant deck properties for quick input.
# Each entry: (profile_type, h_p mm, e mm [centroid height above soffit],
#              conc_vol_trough m³/m²)
# Adjust or extend as needed.
KNOWN_DECKS = {
    # Tata Steel ComFlor range (trapezoidal unless noted)
    "CF46":    ("trapezoidal", 46,  24.4, 0.027),
    "CF51":    ("trapezoidal", 51,  22.5, 0.027),
    "CF60":    ("trapezoidal", 60,  28.8, 0.033),
    "CF70":    ("trapezoidal", 70,  30.3, 0.032),
    "CF80":    ("trapezoidal", 80,  39.0, 0.043),
    "CF100":   ("trapezoidal", 100, 49.8, 0.054),
    "CF225":   ("trapezoidal", 225, 95.0, 0.078),
    # Kingspan Multideck range (re-entrant)
    "MD50":    ("reentrant",   50,  20.0, 0.020),
    "MD60":    ("reentrant",   60,  25.0, 0.022),
}


# ==============================================================================
# STANDARD BEAM SECTIONS FOR ASSUMED BEAM STIFFNESSES
# ==============================================================================
# Assumed composite beam sections (SCI P354 §7 requires an assumption on beam
# section where actual sizes are not yet known).
#   Shallow deck (h_p ≤ 80 mm): 305×165×40 UB  — min b_f = 152 mm per ComFlor
#   Deep deck    (h_p > 80 mm): 305×305×97 UC  — min b_f = 285 mm per ComFlor
# (d_mm, I_x_cm4, A_cm2)
VIB_BEAM_SHALLOW = (303.8, 8503, 51.3)   # 305×165×40 UB
VIB_BEAM_DEEP    = (300.0, 23600, 125.0)  # 305×305×97 UC


# ==============================================================================
# CORE FUNCTIONS
# ==============================================================================

def vib_Is_slab(h_mm, h_p, profile_type, nwc=True):
    """Dynamic second moment of area I_s [mm⁴/m] of composite slab per unit width.

    Uses the Table 7.2 approximation from SCI P354:
        I_s = coeff × h_s^3.7  (NWC)   or   coeff × h_s^3.5  (LWC)
    where h_s is the overall slab depth in mm.

    Parameters
    ----------
    h_mm         : overall slab depth [mm]
    h_p          : deck profile height [mm]
    profile_type : "reentrant" or "trapezoidal"
    nwc          : True for Normal Weight Concrete (≥ 20 kN/m³); False for LWC
    """
    # Table 7.2  (h_p, coeff_NWC, exp_NWC, coeff_LWC, exp_LWC)
    T72_REENTRANT   = [(51,  0.37, 3.7, 0.65, 3.5)]
    T72_TRAPEZOIDAL = [(60,  0.23, 3.7, 0.40, 3.5),
                       (80,  0.19, 3.7, 0.37, 3.5),
                       (225, 0.05, 3.7, 0.12, 3.5)]

    table = (T72_REENTRANT if profile_type == "reentrant"
             else T72_TRAPEZOIDAL)

    if len(table) == 1:
        _, cn, en, cl, el = table[0]
    elif h_p <= table[0][0]:
        _, cn, en, cl, el = table[0]
    elif h_p >= table[-1][0]:
        _, cn, en, cl, el = table[-1]
    else:
        for i in range(len(table) - 1):
            if table[i][0] <= h_p <= table[i+1][0]:
                t = (h_p - table[i][0]) / (table[i+1][0] - table[i][0])
                cn = table[i][1] + t * (table[i+1][1] - table[i][1])
                en = table[i][2]
                cl = table[i][3] + t * (table[i+1][3] - table[i][3])
                el = table[i][4]
                break

    return (cn * h_mm ** en) if nwc else (cl * h_mm ** el)


def vib_composite_EI(h_mm, h_p, beam_d_mm, beam_I_cm4, beam_A_cm2,
                     slab_b_eff_m, E_dyn):
    """Composite beam EI [N·m²] using transformed-section method.

    Fully composite (no slip). Uses dynamic concrete modulus E_dyn.

    Parameters
    ----------
    h_mm         : overall slab depth [mm]
    h_p          : deck profile height [mm]
    beam_d_mm    : steel beam depth [mm]
    beam_I_cm4   : steel beam I_x [cm⁴]
    beam_A_cm2   : steel beam area [cm²]
    slab_b_eff_m : effective slab breadth contributing to composite action [m]
    E_dyn        : dynamic concrete modulus [N/m²]
    """
    E_s_Pa = 210.0e9
    E_s_mm = 210.0e3    # MPa
    n      = E_s_mm / (E_dyn / 1e6)   # dynamic modular ratio

    I_steel = beam_I_cm4 * 1e4   # mm⁴
    A_steel = beam_A_cm2 * 1e2   # mm²
    e_a     = beam_d_mm / 2.0    # centroid of symmetric section from soffit [mm]
    b_eff   = slab_b_eff_m * 1000.0   # mm

    h_c = max(h_mm - h_p, 0.0)        # concrete above deck [mm]
    A_c = b_eff * h_c
    z_c = beam_d_mm + h_p + h_c / 2.0  # centroid of concrete slab from beam soffit [mm]

    denom = A_steel + A_c / n
    z_NA  = (A_steel * e_a + (A_c / n) * z_c) / denom  # neutral axis from soffit [mm]

    I_comp = (I_steel + A_steel * (z_NA - e_a) ** 2
              + b_eff * h_c ** 3 / (12.0 * n)
              + (A_c / n) * (z_c - z_NA) ** 2)   # mm⁴

    return E_s_Pa * I_comp * 1e-12   # N·m²


def vib_EIb_beam(h_mm, h_p, b, L_y, E_dyn):
    """Composite secondary beam EI_b [N·m²].

    Beam section assumed: 305×165×40 UB (shallow, h_p ≤ 80 mm)
                      or  305×305×97 UC (deep,    h_p > 80 mm).
    Effective slab breadth = min(L_y/4, b).

    Parameters
    ----------
    h_mm  : overall slab depth [mm]
    h_p   : deck profile height [mm]
    b     : slab span = secondary beam spacing [m]
    L_y   : secondary beam span [m]
    E_dyn : dynamic concrete modulus [N/m²]
    """
    use_deep = (h_p > 80)
    d_mm, I_cm4, A_cm2 = (VIB_BEAM_DEEP if use_deep else VIB_BEAM_SHALLOW)
    b_eff = min(L_y / 4.0, b)
    return vib_composite_EI(h_mm, h_p, d_mm, I_cm4, A_cm2, b_eff, E_dyn)


def vib_EIp_beam(h_mm, h_p, L_y, L_x, E_dyn):
    """Composite primary beam EI_p [N·m²].

    Beam section assumed: 305×165×40 UB (primary always uses shallow section;
    only the secondary beam needs a wide bottom flange for deep deck profiles).
    Effective slab breadth = min(L_x/4, L_y).

    Parameters
    ----------
    h_mm  : overall slab depth [mm]
    h_p   : deck profile height [mm]
    L_y   : secondary beam span = primary beam spacing [m]
    L_x   : primary beam span [m]
    E_dyn : dynamic concrete modulus [N/m²]
    """
    d_mm, I_cm4, A_cm2 = VIB_BEAM_SHALLOW
    b_eff = min(L_x / 4.0, L_y)
    return vib_composite_EI(h_mm, h_p, d_mm, I_cm4, A_cm2, b_eff, E_dyn)


def check_vibration(
        # Slab properties
        h_mm, h_p, profile_type,
        # Bay geometry
        b, L_y, L_x, n_y=2, n_x=2,
        # Material
        nwc=True,
        # Floor mass
        g_slab_kPa=None, finishes_kPa=0.75, ceiling_services_kPa=0.75,
        imposed_kPa=2.5,
        # Beam stiffnesses (None → use assumed standard sections)
        EI_b=None, EI_p=None,
        # Occupancy
        zeta=0.03, R_limit=8.0,
        # Resonance build-up factor (None → compute via SCI P354 Eq 37)
        rho=None):
    """SCI P354 Section 7 simplified vibration check for composite floors.

    Returns a results dictionary with all intermediate values.

    Parameters
    ----------
    h_mm            : overall slab depth [mm]
    h_p             : deck profile height [mm]
    profile_type    : "reentrant" or "trapezoidal"
    b               : slab span = secondary beam spacing [m]
    L_y             : secondary beam span [m]
    L_x             : primary beam span [m]
    n_y             : number of secondary bays (capped at 4 per SCI P354)
    n_x             : number of primary bays (capped at 4 per SCI P354)
    nwc             : True = Normal Weight Concrete; False = Lightweight Concrete
    g_slab_kPa      : slab self-weight [kN/m²]  (None → estimated from h_mm + profile)
    finishes_kPa    : permanent load from floor finishes and screed [kN/m²]
    ceiling_services_kPa : permanent load from ceiling, services, fittings [kN/m²]
    imposed_kPa     : characteristic imposed load [kN/m²] — only 10% included per P354
    EI_b            : secondary beam composite stiffness [N·m²]  (None → auto)
    EI_p            : primary beam composite stiffness   [N·m²]  (None → auto)
    zeta            : critical damping ratio (SCI P354 Table 4.1)
                        0.03 — fitted office / retail
                        0.045 — residential (fitted)
                        0.05  — industrial
    R_limit         : limiting response factor (SCI P354 Table 5.3)
                        4  — residential
                        8  — offices / retail
                        16 — industrial workshops
    rho             : resonance build-up factor (None → computed via Eq 37)
    """

    g_acc = 9.81

    # ── Dynamic elastic modulus (SCI P354 §4.1.3) ────────────────────────────
    # NWC → 38 GPa; LWC → 22 GPa
    E_dyn = 38.0e9 if nwc else 22.0e9   # N/m²

    # ── Slab I_s per unit width (Table 7.2 power-law) ────────────────────────
    I_s_mm = vib_Is_slab(h_mm, h_p, profile_type, nwc)  # mm⁴/m
    I_s    = I_s_mm * 1.0e-12                            # m⁴/m (= m³)
    EI_s   = 210.0e9 * I_s                              # N·m²/m

    # ── Floor mass m [kg/m²] (SCI P354 §4.1.2) ───────────────────────────────
    # Permanent loads + 10% imposed; partitions excluded
    if g_slab_kPa is None:
        # Simple estimate: average concrete depth × density (2400 kg/m³)
        h_c_avg_mm = (h_mm - h_p) + h_p * 0.5   # approximate average concrete depth
        g_slab_kPa = h_c_avg_mm / 1000.0 * 24.0  # kN/m²
    m_floor = (g_slab_kPa + finishes_kPa + ceiling_services_kPa
               + 0.1 * imposed_kPa) * 1000.0 / g_acc   # kg/m²

    # ── Beam stiffnesses ──────────────────────────────────────────────────────
    if EI_b is None:
        EI_b = vib_EIb_beam(h_mm, h_p, b, L_y, E_dyn)
    if EI_p is None:
        EI_p = vib_EIp_beam(h_mm, h_p, L_y, L_x, E_dyn)

    # ── Deflections (Table 7.1, arrangement 3 — ≥ 3 bays, internal) ─────────
    mg_b = m_floor * g_acc * b   # N/m — UDL × tributary width
    delta_sec_mode = (mg_b / 384.0
                      * (5.0 * L_y ** 4 / EI_b
                         + b ** 3 / EI_s)) * 1000.0   # mm

    # Primary beam: UDL = m·g·L_y [N/m] over span L_x (SS)
    delta_prim_mode = (m_floor * g_acc * L_y * 5.0 * L_x ** 4
                       / (384.0 * EI_p)) * 1000.0   # mm

    delta_mm = max(delta_sec_mode, delta_prim_mode, 1e-6)
    governing_mode = "secondary" if delta_sec_mode >= delta_prim_mode else "primary"

    # ── Fundamental frequency (Eq 4) ─────────────────────────────────────────
    f0 = 18.0 / math.sqrt(delta_mm)   # Hz

    # ── Minimum frequency check (§7.2) ───────────────────────────────────────
    if f0 < 3.0:
        return {
            "passed":          False,
            "util":            3.0 / max(f0, 0.01),
            "f0_Hz":           f0,
            "delta_mm":        delta_mm,
            "delta_sec_mm":    delta_sec_mode,
            "delta_prim_mm":   delta_prim_mode,
            "governing_mode":  governing_mode,
            "m_floor_kg_m2":   m_floor,
            "note": f"FAIL: f0 = {f0:.2f} Hz < 3 Hz minimum (SCI P354 §7.2)",
        }

    # ── η factor (Table 7.3) ──────────────────────────────────────────────────
    if f0 < 5.0:
        eta = 0.5
    elif f0 <= 6.0:
        eta = 0.21 * f0 - 0.55
    else:
        eta = 0.71

    # ── Modal mass M (§7.3.1 Eq 45–46) ──────────────────────────────────────
    n_y_eff = min(n_y, 4)
    n_x_eff = min(n_x, 4)

    L_eff = (1.09 * (1.10 ** ((n_y_eff - 1) / 2.0))
             * (EI_b / (m_floor * b * f0 ** 2)) ** 0.25
             * L_y)
    L_eff = min(L_eff, n_y_eff * L_y)

    S = (eta * (1.15 ** ((n_x_eff - 1) / 2.0))
         * (EI_s / (m_floor * f0 ** 2)) ** 0.25)
    S = min(S, n_x_eff * L_x)

    M = m_floor * L_eff * S   # kg

    # ── Resonance build-up factor ρ (§6.5.2 Eq 37) ──────────────────────────
    if rho is None:
        f_p = 2.0    # Hz — typical walking pace
        v_w = 1.5    # m/s — typical walking speed
        L_p = L_y    # walking path = secondary beam span (conservative)
        rho = 1.0 - math.exp(-2.0 * math.pi * zeta * L_p * f_p / v_w)

    # ── Frequency weighting W (§7.6, BS 6472 Wg z-axis) ─────────────────────
    # Flat at 1.0 for 2–8 Hz; decreasing above 8 Hz
    W = 1.0 if f0 <= 8.0 else 8.0 / f0

    # ── Response acceleration (§7.5) ─────────────────────────────────────────
    Q    = 746.0   # N — average person weight (76 kg × 9.81)
    phi_e = 1.0    # mode shape at excitation (conservative)
    phi_r = 1.0    # mode shape at response (conservative)

    if f0 <= 10.0:
        # Low-frequency floor — resonant response (Eq 50)
        a_rms = (0.1 * Q * rho * phi_e * phi_r * W
                 / (math.sqrt(2.0) * zeta * M))
    else:
        # High-frequency floor — transient response (Eq 51)
        a_rms = (185.0 * Q * phi_e * phi_r * W
                 / (M * f0 ** 0.3 * 700.0 * math.sqrt(2.0)))

    # ── Response factor R (Eq 38) ─────────────────────────────────────────────
    R    = a_rms / 0.005
    util = R / R_limit

    return {
        "passed":          util <= 1.0,
        "util":            util,
        "f0_Hz":           f0,
        "R":               R,
        "R_limit":         R_limit,
        "M_kg":            M,
        "rho":             rho,
        "W":               W,
        "a_rms_mm_s2":     a_rms * 1000.0,
        "delta_mm":        delta_mm,
        "delta_sec_mm":    delta_sec_mode,
        "delta_prim_mm":   delta_prim_mode,
        "governing_mode":  governing_mode,
        "m_floor_kg_m2":   m_floor,
        "EI_b_Nm2":        EI_b,
        "EI_p_Nm2":        EI_p,
        "EI_s_Nm2_per_m":  EI_s,
        "eta":             eta,
        "L_eff_m":         L_eff,
        "S_m":             S,
        "note": (f"{'PASS' if util <= 1.0 else 'FAIL'}: "
                 f"R = {R:.1f} (limit = {R_limit:.0f});  "
                 f"f0 = {f0:.2f} Hz;  "
                 f"M = {M:.0f} kg;  "
                 f"δ = {delta_mm:.1f} mm ({governing_mode} mode)"),
    }


def print_results(res, label=""):
    """Print a formatted summary of check_vibration() results."""
    tag = f"  [{label}]" if label else ""
    status = "PASS" if res["passed"] else "FAIL"
    bar = "=" * 70
    print(bar)
    print(f"  SCI P354 §7 Vibration Check{tag}  →  {status}")
    print(bar)
    print(f"  Fundamental frequency     f0  = {res['f0_Hz']:.3f} Hz")
    print(f"  Governing deflection      δ   = {res['delta_mm']:.2f} mm  ({res['governing_mode']} beam mode)")
    print(f"    Secondary beam mode     δ_s = {res['delta_sec_mm']:.2f} mm")
    print(f"    Primary beam mode       δ_p = {res['delta_prim_mm']:.2f} mm")
    print(f"  Floor mass                m   = {res['m_floor_kg_m2']:.1f} kg/m²")
    if res["passed"] or res.get("R") is not None:
        print(f"  Modal mass                M   = {res.get('M_kg', float('nan')):.0f} kg")
        print(f"  Resonance factor          ρ   = {res.get('rho', float('nan')):.3f}")
        print(f"  Frequency weighting       W   = {res.get('W', float('nan')):.2f}")
        print(f"  RMS acceleration          a   = {res.get('a_rms_mm_s2', float('nan')):.3f} mm/s²")
        print(f"  Response factor           R   = {res.get('R', float('nan')):.2f}  (limit = {res.get('R_limit', float('nan')):.0f})")
        print(f"  Utilisation               U   = {res['util']:.3f}")
    print(f"  {res['note']}")
    print()


# ==============================================================================
# USER INPUTS — edit this section to run a check
# ==============================================================================
if __name__ == "__main__":

    # ── Slab ──────────────────────────────────────────────────────────────────
    DECK_NAME    = "CF60"          # key into KNOWN_DECKS, or set manually below
    H_MM         = 130             # overall slab depth [mm]

    # Override with custom deck properties if needed:
    #   PROFILE_TYPE = "trapezoidal"  # or "reentrant"
    #   H_P          = 60             # deck profile height [mm]
    profile_type, h_p, _, _ = KNOWN_DECKS[DECK_NAME]

    # ── Bay geometry ──────────────────────────────────────────────────────────
    B   = 3.0     # slab span = secondary beam spacing [m]
    L_Y = 9.0     # secondary beam span [m]
    L_X = 9.0     # primary beam span [m]
    N_Y = 4       # number of secondary beam bays (SCI P354 caps at 4)
    N_X = 4       # number of primary beam bays

    # ── Floor loading ─────────────────────────────────────────────────────────
    # Permanent loads (excluding slab self-weight, which is estimated automatically)
    FINISHES_KPA    = 0.75   # kN/m²  floor finishes + screed
    SERVICES_KPA    = 0.75   # kN/m²  ceiling + services + fittings
    IMPOSED_KPA     = 3.0    # kN/m²  characteristic imposed (only 10% enters mass)

    # ── Occupancy (damping and response factor limit) ─────────────────────────
    #
    #  Occupancy     zeta    R_limit   Note
    #  -----------   -----   -------   ------------------------------------------
    #  Residential   0.045     4       SCI P354 Table 4.1 / Table 5.3
    #  Office        0.030     8       fitted office
    #  Retail        0.030     8
    #  Industrial    0.050    16       workshops
    ZETA    = 0.03
    R_LIMIT = 8.0

    # ── Concrete type ─────────────────────────────────────────────────────────
    NWC = True   # True = Normal Weight Concrete; False = Lightweight

    # ── Beam stiffnesses (leave None to use assumed standard sections) ─────────
    EI_B = None   # secondary beam composite EI [N·m²]  (None → auto)
    EI_P = None   # primary beam composite EI   [N·m²]  (None → auto)

    # ==============================================================================
    # RUN CHECK
    # ==============================================================================
    print(f"\nDeck: {DECK_NAME}  h = {H_MM} mm  (h_p = {h_p} mm, {profile_type})")
    print(f"Bay: b = {B} m  L_y = {L_Y} m  L_x = {L_X} m  "
          f"n_y = {N_Y}  n_x = {N_X}")
    print(f"Loads: finishes = {FINISHES_KPA} kN/m²  services = {SERVICES_KPA} kN/m²  "
          f"imposed = {IMPOSED_KPA} kN/m²")
    print(f"Damping: zeta = {ZETA}  R_limit = {R_LIMIT}\n")

    res = check_vibration(
        h_mm=H_MM, h_p=h_p, profile_type=profile_type,
        b=B, L_y=L_Y, L_x=L_X, n_y=N_Y, n_x=N_X,
        nwc=NWC,
        finishes_kPa=FINISHES_KPA, ceiling_services_kPa=SERVICES_KPA,
        imposed_kPa=IMPOSED_KPA,
        EI_b=EI_B, EI_p=EI_P,
        zeta=ZETA, R_limit=R_LIMIT,
    )
    print_results(res, label=f"{DECK_NAME} h={H_MM}mm  b={B}m L_y={L_Y}m")

    # ── Sweep over secondary beam spacings ────────────────────────────────────
    print("\n--- Parametric sweep: secondary beam spacing b [m] ---\n")
    print(f"  {'b (m)':<8} {'f0 (Hz)':<10} {'δ (mm)':<10} {'R':<8} {'Util':<8} Status")
    print(f"  {'-'*56}")
    for b_val in [2.0, 2.5, 3.0, 3.5, 4.0, 4.5]:
        r = check_vibration(
            h_mm=H_MM, h_p=h_p, profile_type=profile_type,
            b=b_val, L_y=L_Y, L_x=L_X, n_y=N_Y, n_x=N_X,
            nwc=NWC,
            finishes_kPa=FINISHES_KPA, ceiling_services_kPa=SERVICES_KPA,
            imposed_kPa=IMPOSED_KPA,
            EI_b=EI_B, EI_p=EI_P,
            zeta=ZETA, R_limit=R_LIMIT,
        )
        status = "PASS" if r["passed"] else "FAIL"
        R_val  = r.get("R", float("nan"))
        print(f"  {b_val:<8.1f} {r['f0_Hz']:<10.3f} {r['delta_mm']:<10.2f} "
              f"{R_val:<8.2f} {r['util']:<8.3f} {status}")
