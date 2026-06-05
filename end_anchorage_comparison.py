"""end_anchorage_comparison.py

Compares the effect of adding headed-stud end anchorage (EN 1994-1-1 Cl 9.7.4)
to the embossment-only (τ_u) shear bond for all τ_u composite decks in the
database.  m-k-only decks are skipped (τ_u is required for the partial
connection method).

Three parts
-----------
Part 1 – Moment capacity table
    Compares M_Rd at each span for: embossments only / + end anchorage / full
    connection.  Fixed representative slab depth (h_c above deck = H_C).

Part 2 – Pass/fail against design M_Ed (office loading, 1-span)
    Loads the actual optimised section depths from the results CSV and checks
    whether each composite slab passes the longitudinal shear check with and
    without end-anchorage studs.  Shows GWP cost of the studs where needed.

Part 3 – Case study
    Detailed worked example for a single chosen deck × span, showing all
    intermediate values.

End anchorage formulae (EN 1994-1-1 Cl 9.7.4 / SCI P300):
    P_pb,Rd = k_φ · d_do · t · f_yp,d          (Eq 9.10, bearing on deck)
    k_φ     = 1 + a / d_do  ≤  6.0
    d_do    = 1.1 × d_shank
    a       = distance from stud centre to sheet end  (≥ 1.5 d_do)
    N_a     = (1000 / l_1) × P_pb,Rd            (per metre width, 1 stud/rib)
"""

import sqlite3
import math
import os

# ── Configuration ──────────────────────────────────────────────────────────────

DATABASE    = "database_270426.db"
RESULTS_CSV = "results_office_ENV.csv"       # 1-span office envelope results

SPANS    = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]   # m (Part 1)
H_C_REP  = 75.0    # mm — representative concrete above deck (Part 1 only)
B_MM     = 1000.0  # mm — per metre width

# ── Office loading (consistent with optimiser) ─────────────────────────────────
QK_KNM2         = 2.5    # kN/m²  imposed (excl. partitions)
PARTITION_KNM2  = 0.5    # kN/m²  (EN 1991-1-1 Cl 6.3.1.2)
FINISHES_KNM2   = 0.355  # kN/m²  (Spanplatte + Gipsfaserplatte + Steinwolle)
SERVICES_KNM2   = 0.75   # kN/m²  ceiling + services
GAMMA_G         = 1.35
GAMMA_Q         = 1.50
M_COEFF_SAG     = 1.0 / 8.0   # SS single-span sagging coefficient

# ── Concrete C25/30 ────────────────────────────────────────────────────────────
FCK_MPA  = 25.0
FCD_MPA  = FCK_MPA / 1.5
ECM_MPA  = 31000.0   # MPa

# ── Stud properties (EN 13918, 19 mm shank) ───────────────────────────────────
D_SHANK    = 19.0            # mm shank diameter (most common through-deck size)
D_DO       = 1.1 * D_SHANK   # weld-collar outer diameter = 20.9 mm
A_END_MM   = 40.0            # mm — stud centre to sheet end (practical minimum)
K_PHI      = min(1.0 + A_END_MM / D_DO, 6.0)   # ≈ 2.91
F_U_STUD   = 450.0           # N/mm²  (EN 13918)
H_SC       = 100.0           # mm  standard height
GAMMA_VS   = 1.25            # EN 1994-1-1 shear bond partial factor
GAMMA_V    = 1.25            # stud shear partial factor

# Legacy limits used in Part 1 only
A_MIN      = 1.5 * D_DO
K_PHI_MIN  = min(1.0 + A_MIN / D_DO, 6.0)   # = 2.5
K_PHI_MAX  = 6.0

# GWP of headed studs — structural steel proxy per Sansom (2014)
STUD_GWP   = 1.823   # kg CO₂e/kg
STEEL_DENSITY = 7850.0  # kg/m³


# ── Helpers ────────────────────────────────────────────────────────────────────

def _M_Rd(N_c, N_cf, h_mm, e_p, e, M_cRd):
    """Composite sagging moment resistance [kNm/m] at partial connection.
    EN 1994-1-1 Cl 9.7.4 — partial connection method (Eqs 9.8–9.9).
    """
    N_c  = min(N_c, N_cf)
    eta  = N_c / N_cf
    x_pl = N_c / (0.85 * FCD_MPA * B_MM)
    z    = h_mm - 0.5 * x_pl - e_p + (e_p - e) * eta
    M_pa = M_cRd
    M_pr = min(1.25 * M_pa * (1.0 - eta), M_pa)
    return N_c / 1000.0 * z / 1000.0 + M_pr   # kNm/m


def _P_pb_Rd(k_phi, t_mm, f_ypd):
    """End-anchorage bearing resistance per stud [N] (Eq 9.10)."""
    return k_phi * D_DO * t_mm * f_ypd


def _N_a(t_mm, f_ypd, l_1_mm):
    """End-anchorage force per metre width [N/m], 1 stud per rib."""
    P_bear   = _P_pb_Rd(K_PHI, t_mm, f_ypd)
    alpha_s  = 1.0 if H_SC / D_SHANK >= 4.0 else 0.2 * (H_SC / D_SHANK + 1.0)
    P_shear1 = 0.8 * F_U_STUD * math.pi * D_SHANK**2 / 4.0 / GAMMA_V
    P_shear2 = 0.29 * alpha_s * D_SHANK**2 * math.sqrt(FCK_MPA * ECM_MPA) / GAMMA_V
    P_stud   = min(P_shear1, P_shear2)
    P_pb     = min(P_bear, P_stud)
    return (1000.0 / l_1_mm) * P_pb


def _stud_gwp(l_1_mm):
    """GWP contribution of end-anchorage studs per m² of slab [kg CO₂e/m²]."""
    d_s   = D_SHANK / 1000.0
    h_sc  = H_SC / 1000.0
    V     = math.pi / 4.0 * d_s**2 * h_sc + math.pi / 4.0 * (2.5*d_s)**2 * (0.4*d_s)
    mass  = (1000.0 / l_1_mm) * V * STEEL_DENSITY   # kg/m² (1 stud/rib per end)
    return mass * STUD_GWP


# ── Database query ─────────────────────────────────────────────────────────────

conn = sqlite3.connect(DATABASE)
cur  = conn.cursor()
cur.execute("""
    SELECT Deck, t, h_p, A_pe, e, e_p, f_ypd, tau_uRk, M_cRd, l_1
    FROM   sheeting_prop
    WHERE  tau_uRk IS NOT NULL
    ORDER  BY h_p, Deck
""")
rows = cur.fetchall()
conn.close()
rows = [r for r in rows if r[7] == r[7]]   # drop NaN tau_uRk

if not rows:
    print("No decks with τ_u data found in the database.")
    raise SystemExit

# ══════════════════════════════════════════════════════════════════════════════
# PART 1 — Moment capacity table (fixed representative H_C)
# ══════════════════════════════════════════════════════════════════════════════

print()
print("=" * 110)
print("PART 1 — MOMENT CAPACITY COMPARISON  (representative h_c above deck"
      f" = {H_C_REP:.0f} mm, C25/30)")
print(f"  Stud: ø{D_SHANK:.0f} mm, d_do = {D_DO:.1f} mm, a = {A_END_MM:.0f} mm → k_φ = {K_PHI:.2f}")
print(f"  k_φ_min (a = 1.5·d_do = {A_MIN:.0f} mm): {K_PHI_MIN:.1f}   |   k_φ_max: {K_PHI_MAX:.1f}")
print(f"  Layout: 1 stud per rib per slab end")
print("=" * 110)

COL = (f"{'Deck':<26} {'L':>4}  {'η_emb':>6}  "
       f"{'M_emb':>7}  {'M_anch(min)':>11}  {'M_anch(max)':>11}  {'M_pl':>7}  "
       f"{'Δanch(min)%':>11}  {'Δanch(max)%':>11}  {'Δfull%':>7}")
SEP = "-" * len(COL)
print(COL); print(SEP)

pct_min_all, pct_max_all, pct_full_all = [], [], []
latex_p1_rows = []   # (deck_name, span, eta_emb, M_emb, M_anch_mn, M_anch_mx, M_pl, pct_min, pct_max, pct_full)

for (deck_name, t, h_p, A_pe, e, e_p, f_ypd, tau_uRk, M_cRd, l_1) in rows:
    h_mm  = h_p + H_C_REP
    N_cf  = A_pe * f_ypd
    n_pp  = B_MM / l_1
    p_min = _P_pb_Rd(K_PHI_MIN, t, f_ypd)
    p_max = _P_pb_Rd(K_PHI_MAX, t, f_ypd)
    M_pl  = _M_Rd(N_cf, N_cf, h_mm, e_p, e, M_cRd)

    first = True
    for span in SPANS:
        L_x  = span * 1000.0 / 2
        N_c0 = (tau_uRk / GAMMA_VS) * B_MM * L_x
        if N_c0 >= N_cf:
            continue

        eta_emb   = N_c0 / N_cf
        N_c_min   = min(N_c0 + n_pp * p_min, N_cf)
        N_c_max   = min(N_c0 + n_pp * p_max, N_cf)
        M_emb     = _M_Rd(N_c0,    N_cf, h_mm, e_p, e, M_cRd)
        M_anch_mn = _M_Rd(N_c_min, N_cf, h_mm, e_p, e, M_cRd)
        M_anch_mx = _M_Rd(N_c_max, N_cf, h_mm, e_p, e, M_cRd)
        pct_min   = (M_anch_mn - M_emb) / M_emb * 100.0
        pct_max   = (M_anch_mx - M_emb) / M_emb * 100.0
        pct_full  = (M_pl      - M_emb) / M_emb * 100.0
        pct_min_all.append(pct_min); pct_max_all.append(pct_max); pct_full_all.append(pct_full)
        latex_p1_rows.append((deck_name, span, eta_emb, M_emb, M_anch_mn, M_anch_mx, M_pl,
                               pct_min, pct_max, pct_full))

        label = deck_name if first else ""
        first = False
        print(f"{label:<26} {span:>4.1f}  {eta_emb:>6.3f}  "
              f"{M_emb:>7.2f}  {M_anch_mn:>11.2f}  {M_anch_mx:>11.2f}  {M_pl:>7.2f}  "
              f"{pct_min:>11.1f}  {pct_max:>11.1f}  {pct_full:>7.1f}")
    if not first:
        print()

print(SEP)
if pct_min_all:
    n = len(pct_min_all)
    print(f"\nSummary across {n} partially-connected deck × span combinations:\n")
    for label, data in [
        (f"End anchorage, k_φ = {K_PHI_MIN:.1f} (min)", pct_min_all),
        (f"End anchorage, k_φ = {K_PHI_MAX:.1f} (max)", pct_max_all),
        ("Full shear connection (upper bound)",           pct_full_all),
    ]:
        print(f"  {label:<44}  min={min(data):5.1f}%  max={max(data):5.1f}%  "
              f"mean={sum(data)/len(data):5.1f}%")
    print()

# ══════════════════════════════════════════════════════════════════════════════
# PART 2 — Pass/fail at design M_Ed using actual optimised section depths
# ══════════════════════════════════════════════════════════════════════════════

print()
print("=" * 120)
print("PART 2 — LONGITUDINAL SHEAR PASS/FAIL AT DESIGN M_Ed"
      "  (office, 1-span, ENV criterion)")
print(f"  Loading: qk = {QK_KNM2} + {PARTITION_KNM2} kN/m²  |  "
      f"finishes = {FINISHES_KNM2} kN/m²  |  services = {SERVICES_KNM2} kN/m²")
print(f"  ULS combination: {GAMMA_G}·Gk + {GAMMA_Q}·Qk  |  M_Ed = (1/8)·w_uls·L²")
print(f"  Stud: ø{D_SHANK:.0f} mm, k_φ = {K_PHI:.2f}, 1/rib, N_a as computed")
print("=" * 120)

HDR2 = (f"{'Deck':<26} {'L':>4}  {'h_mm':>5}  {'M_Ed':>6}  {'M_Rd_emb':>9}  "
        f"{'M_Rd_stud':>10}  {'N_a':>7}  "
        f"{'No studs':>9}  {'+ studs':>8}  {'GWP_stud':>9}")
print(HDR2); print("-" * len(HDR2))

# Build a lookup: deck family → (t, h_p, A_pe, e, e_p, f_ypd, tau_uRk, M_cRd, l_1)
deck_lu = {}
for (deck_name, t, h_p, A_pe, e, e_p, f_ypd, tau_uRk, M_cRd, l_1) in rows:
    # Group by deck family (strip gauge suffix)
    family = ' '.join(deck_name.split()[:2])
    if family not in deck_lu:
        deck_lu[family] = (t, h_p, A_pe, e, e_p, f_ypd, tau_uRk, M_cRd, l_1)

# Try to load actual optimised h values from results CSV
opt_h = {}   # (deck_family, span) -> h_mm
if os.path.exists(RESULTS_CSV):
    import csv
    with open(RESULTS_CSV, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get('slab_type') != 'comp_slab':
                continue
            try:
                h_val  = float(row['h_struct_m']) * 1000.0
                span   = float(row['span_m'])
                d_name = row.get('deck_name', '')
                fam    = ' '.join(d_name.split()[:2])
                if fam and not math.isnan(h_val):
                    key = (fam, span)
                    if key not in opt_h or h_val < opt_h[key]:
                        opt_h[key] = h_val
            except (ValueError, KeyError):
                pass

stud_helps_count = 0
stud_not_needed  = 0
stud_insufficient = 0
latex_p2_rows = []   # (family, span, h_mm, M_Ed, M_Rd_emb, M_Rd_stud, Na, status_emb, status_stud, gwp_st_or_nan)

for family, (t, h_p, A_pe, e, e_p, f_ypd, tau_uRk, M_cRd, l_1) in sorted(deck_lu.items()):
    N_cf   = A_pe * f_ypd
    Na     = _N_a(t, f_ypd, l_1)
    gwp_st = _stud_gwp(l_1)
    first  = True

    for span in [3.0, 4.0, 5.0, 6.0, 7.0, 8.0]:
        L_x    = (span / 2.0) * 1000.0        # mm — shear span (SS)
        N_c0   = (tau_uRk / GAMMA_VS) * B_MM * L_x

        # Use actual optimised h if available, otherwise estimate
        h_mm = opt_h.get((family, span))
        if h_mm is None:
            h_mm = h_p + H_C_REP   # fallback

        # ULS design moment (office, 1-span)
        gk_sdl  = (FINISHES_KNM2 + SERVICES_KNM2)   # kN/m²
        qk_tot  = QK_KNM2 + PARTITION_KNM2            # kN/m²
        # self-weight estimated from h_mm (conservative — solid depth)
        gk_sw   = 25.0 * h_mm / 1000.0               # kN/m²
        w_uls   = GAMMA_G * (gk_sw + gk_sdl) + GAMMA_Q * qk_tot
        M_Ed    = M_COEFF_SAG * w_uls * span**2        # kNm/m

        M_Rd_emb  = _M_Rd(N_c0,        N_cf, h_mm, e_p, e, M_cRd)
        M_Rd_stud = _M_Rd(min(N_c0 + Na, N_cf), N_cf, h_mm, e_p, e, M_cRd)

        pass_emb  = M_Ed <= M_Rd_emb
        pass_stud = M_Ed <= M_Rd_stud

        if pass_emb:
            stud_not_needed += 1
            status_emb  = "PASS"
            status_stud = "n/a"
            gwp_str     = "—"
        elif pass_stud:
            stud_helps_count += 1
            status_emb  = "FAIL"
            status_stud = "PASS"
            gwp_str     = f"{gwp_st:.3f}"
        else:
            stud_insufficient += 1
            status_emb  = "FAIL"
            status_stud = "FAIL"
            gwp_str     = "—"

        latex_p2_rows.append((family, span, h_mm, M_Ed, M_Rd_emb, M_Rd_stud,
                               Na / 1000.0, status_emb, status_stud,
                               gwp_st if not pass_emb and pass_stud else None))

        label = family if first else ""
        first = False
        print(f"{label:<26} {span:>4.1f}  {h_mm:>5.0f}  {M_Ed:>6.2f}  "
              f"{M_Rd_emb:>9.2f}  {M_Rd_stud:>10.2f}  {Na/1000:>7.1f}  "
              f"{status_emb:>9}  {status_stud:>8}  {gwp_str:>9}")
    if not first:
        print()

print("-" * len(HDR2))
total = stud_not_needed + stud_helps_count + stud_insufficient
print(f"\nSummary ({total} deck × span combinations):")
print(f"  Bond/friction passes without studs : {stud_not_needed:3d}"
      f" ({stud_not_needed/total*100:.0f}%)")
print(f"  Studs required and sufficient      : {stud_helps_count:3d}"
      f" ({stud_helps_count/total*100:.0f}%)")
print(f"  Studs required but still failing   : {stud_insufficient:3d}"
      f" ({stud_insufficient/total*100:.0f}%)")
gwp_vals = [_stud_gwp(r[9]) for r in rows]
print(f"  GWP cost of studs where needed     : "
      f"~{min(gwp_vals):.2f}-{max(gwp_vals):.2f} kg CO2e/m2")
print()

# ══════════════════════════════════════════════════════════════════════════════
# PART 3 — Case study: ComFlor 60 / C25/30 at 5 m span (typical office floor)
# ══════════════════════════════════════════════════════════════════════════════

CASE_FAMILY = "Cofrastra 60"
CASE_SPAN   = 4.0   # m — chosen because bond/friction alone fails but studs suffice
CASE_H_MM   = opt_h.get((CASE_FAMILY, CASE_SPAN))
if CASE_H_MM is None:
    CASE_H_MM = 135.0   # fallback if no results file

print()
print("=" * 80)
print(f"PART 3 — CASE STUDY: {CASE_FAMILY} (representative gauge), C25/30, L = {CASE_SPAN:.0f} m")
print("=" * 80)

case_deck = deck_lu.get(CASE_FAMILY)
if case_deck is None:
    print(f"  Deck family '{CASE_FAMILY}' not found in database (no τ_u data).")
else:
    t, h_p, A_pe, e, e_p, f_ypd, tau_uRk, M_cRd, l_1 = case_deck

    N_cf   = A_pe * f_ypd
    L_x    = (CASE_SPAN / 2.0) * 1000.0       # mm shear span
    N_c0   = (tau_uRk / GAMMA_VS) * B_MM * L_x
    Na     = _N_a(t, f_ypd, l_1)
    N_c_with = min(N_c0 + Na, N_cf)

    gk_sdl  = FINISHES_KNM2 + SERVICES_KNM2
    qk_tot  = QK_KNM2 + PARTITION_KNM2
    gk_sw   = 25.0 * CASE_H_MM / 1000.0
    w_uls   = GAMMA_G * (gk_sw + gk_sdl) + GAMMA_Q * qk_tot
    M_Ed    = M_COEFF_SAG * w_uls * CASE_SPAN**2

    M_Rd_emb  = _M_Rd(N_c0,    N_cf, CASE_H_MM, e_p, e, M_cRd)
    M_Rd_stud = _M_Rd(N_c_with, N_cf, CASE_H_MM, e_p, e, M_cRd)
    gwp_st    = _stud_gwp(l_1)

    eta_emb  = min(N_c0, N_cf) / N_cf
    eta_stud = N_c_with / N_cf

    alpha_s  = 1.0 if H_SC / D_SHANK >= 4.0 else 0.2 * (H_SC / D_SHANK + 1.0)
    P_shear1 = 0.8 * F_U_STUD * math.pi * D_SHANK**2 / 4.0 / GAMMA_V
    P_shear2 = 0.29 * alpha_s * D_SHANK**2 * math.sqrt(FCK_MPA * ECM_MPA) / GAMMA_V
    P_stud   = min(P_shear1, P_shear2)
    P_bear   = _P_pb_Rd(K_PHI, t, f_ypd)
    P_pb_Rd  = min(P_bear, P_stud)

    print(f"\n  Deck properties (representative gauge):")
    print(f"    h_p = {h_p:.0f} mm,  t = {t:.1f} mm,  f_ypd = {f_ypd:.0f} N/mm²")
    print(f"    A_pe = {A_pe:.0f} mm²/m,  e = {e:.1f} mm,  e_p = {e_p:.1f} mm")
    print(f"    τ_u,Rk = {tau_uRk:.0f} N/m²,  l_1 = {l_1:.0f} mm (rib spacing)")
    print(f"    M_cRd (bare deck) = {M_cRd:.2f} kNm/m")

    print(f"\n  Section depth (from optimiser): h = {CASE_H_MM:.0f} mm")
    print(f"  Full-connection force: N_cf = A_pe × f_ypd = {A_pe:.0f} × {f_ypd:.0f}"
          f" = {N_cf/1000:.0f} kN/m")

    print(f"\n  Shear span:  L_x = L/2 = {CASE_SPAN/2:.2f} m = {L_x:.0f} mm  (simply-supported)")
    print(f"  Shear bond force:")
    print(f"    N_c,bond = (τ_u,Rk/γ_vs) × b × L_x")
    print(f"             = ({tau_uRk:.0f}/{GAMMA_VS}) × {B_MM:.0f} × {L_x:.0f}")
    print(f"             = {N_c0/1000:.1f} kN/m  →  η = {eta_emb:.3f}")

    print(f"\n  End anchorage (Eq 9.10):")
    print(f"    d_do = 1.1 × {D_SHANK:.0f} = {D_DO:.1f} mm")
    print(f"    k_φ  = 1 + {A_END_MM:.0f}/{D_DO:.1f} = {K_PHI:.2f}  (≤ 6.0)")
    print(f"    P_bear = k_φ × d_do × t × f_ypd = {K_PHI:.2f} × {D_DO:.1f} × {t:.1f}"
          f" × {f_ypd:.0f} = {P_bear:.0f} N/stud")
    print(f"    P_shear (Cl 6.6.3.1):")
    print(f"      P_Rd,1 = 0.8 × {F_U_STUD:.0f} × π × {D_SHANK:.0f}²/4 / {GAMMA_V}"
          f" = {P_shear1:.0f} N/stud")
    print(f"      P_Rd,2 = 0.29 × α_s × {D_SHANK:.0f}² × √({FCK_MPA:.0f}×{ECM_MPA:.0f})"
          f" / {GAMMA_V} = {P_shear2:.0f} N/stud")
    print(f"    P_pb,Rd = min({P_bear:.0f}, {P_stud:.0f}) = {P_pb_Rd:.0f} N/stud")
    print(f"    N_a = (1000/{l_1:.0f}) × {P_pb_Rd:.0f} = {Na/1000:.1f} kN/m")

    print(f"\n  Design moment:  M_Ed = (1/8) × {w_uls:.2f} × {CASE_SPAN:.1f}² = {M_Ed:.2f} kNm/m")
    print(f"\n  Moment resistance comparison:")
    print(f"    Without end anchorage (embossments only):")
    print(f"      N_c = {N_c0/1000:.1f} kN/m  →  M_Rd = {M_Rd_emb:.2f} kNm/m"
          f"  →  {'PASS' if M_Ed <= M_Rd_emb else 'FAIL'}  "
          f"(util = {M_Ed/M_Rd_emb:.3f})")
    print(f"    With end anchorage (ø{D_SHANK:.0f} mm studs, 1/rib):")
    print(f"      N_c = {N_c0/1000:.1f} + {Na/1000:.1f} = {N_c_with/1000:.1f} kN/m"
          f"  →  η = {eta_stud:.3f}")
    print(f"      M_Rd = {M_Rd_stud:.2f} kNm/m"
          f"  →  {'PASS' if M_Ed <= M_Rd_stud else 'FAIL'}  "
          f"(util = {M_Ed/M_Rd_stud:.3f})")
    print(f"\n  Stud GWP cost: {gwp_st:.3f} kg CO₂e/m²"
          f"  (mass = {gwp_st/STUD_GWP:.3f} kg/m²)")
    if M_Ed > M_Rd_emb and M_Ed <= M_Rd_stud:
        print(f"\n  ► Studs enable a {(M_Rd_stud - M_Rd_emb)/M_Rd_emb*100:.1f}% increase in M_Rd")
        print(f"    at a GWP cost of {gwp_st:.3f} kg CO₂e/m² — a minor addition relative")
        print(f"    to the total structural GWP (~40–70 kg CO₂e/m² for composite slabs).")
    elif M_Ed <= M_Rd_emb:
        print(f"\n  ► Embossment bond alone is sufficient at this span — studs not needed.")
    else:
        print(f"\n  ► Even with studs the check fails at this span —")
        print(f"    a deeper section or different deck would be required.")

print()

# Capture case study values NOW before the LaTeX loops overwrite loop variables.
_cs = None
if case_deck is not None:
    _cs = dict(
        family=CASE_FAMILY, span=CASE_SPAN, h_mm=CASE_H_MM,
        t=t, h_p=h_p, A_pe=A_pe, e=e, e_p=e_p,
        f_ypd=f_ypd, tau_uRk=tau_uRk, M_cRd=M_cRd, l_1=l_1,
        N_cf=N_cf, L_x=L_x, N_c0=N_c0, Na=Na,
        eta_emb=eta_emb, eta_stud=eta_stud, N_c_with=N_c_with,
        P_bear=P_bear, P_shear1=P_shear1, P_shear2=P_shear2, P_pb_Rd=P_pb_Rd,
        w_uls=w_uls, M_Ed=M_Ed, M_Rd_emb=M_Rd_emb, M_Rd_stud=M_Rd_stud,
        gwp_st=gwp_st,
    )

# ══════════════════════════════════════════════════════════════════════════════
# LATEX OUTPUT
# ══════════════════════════════════════════════════════════════════════════════

LATEX_FILE = "end_anchorage_comparison.tex"

def _tex_pct(v):
    return f"{v:.1f}"

def _tex_status(s):
    if s == "PASS":
        return r"\textbf{PASS}"
    elif s == "FAIL":
        return r"\textcolor{red}{FAIL}"
    else:
        return r"\textit{n/a}"

# ── collect Part 2 into a compact per-family×span matrix ──────────────────────
# status codes: P = pass without studs, S = studs needed (pass with), F = fail both
p2_families = []
p2_matrix   = {}   # (family, span) -> (status_code, h_mm, M_Ed, M_Rd_emb, gwp_st)
p2_spans    = [3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
for (fam, span, h_mm, M_Ed, M_Rd_emb, M_Rd_stud, Na, se, ss, gwp) in latex_p2_rows:
    if fam not in p2_families:
        p2_families.append(fam)
    if se == "PASS":
        code = "P"
    elif ss == "PASS":
        code = "S"
    else:
        code = "F"
    p2_matrix[(fam, span)] = (code, h_mm, M_Ed, M_Rd_emb, gwp)

def _cell(fam, span):
    v = p2_matrix.get((fam, span))
    if v is None:
        return "---"
    code, h_mm, M_Ed, M_Rd_emb, gwp = v
    if code == "P":
        return r"\cellcolor{green!15}\textbf{P}"
    elif code == "S":
        return r"\cellcolor{yellow!40}\textbf{P*}"
    else:
        return r"\cellcolor{red!15}F"

with open(LATEX_FILE, "w", encoding="utf-8") as f:

    f.write("% Auto-generated by end_anchorage_comparison.py\n")
    f.write("% Requires: booktabs, colortbl, xcolor, longtable packages\n\n")

    # ── Preamble note ──────────────────────────────────────────────────────────
    f.write("\\subsection*{End-anchorage shear stud comparison}\n")
    f.write("\\label{app:end_anchorage}\n\n")
    f.write(
        "The following tables summarise the effect of headed-stud end anchorage "
        "(EN~1994-1-1 Cl~9.7.4) on the longitudinal shear capacity of all "
        r"$\tau_u$" "-method composite decks in the database. "
        f"Stud geometry: \\o\\,{D_SHANK:.0f}\\,mm shank, "
        f"$d_{{do}} = {D_DO:.1f}$\\,mm, $a = {A_END_MM:.0f}$\\,mm "
        f"$\\Rightarrow k_\\phi = {K_PHI:.2f}$, one stud per rib per slab end.\n\n"
    )

    # ── Part 1: Summary statistics table ──────────────────────────────────────
    f.write("\\subsubsection*{Part~1 --- Moment capacity increase from end anchorage}\n\n")
    f.write(
        f"Representative slab depth: $h_c = {H_C_REP:.0f}$\\,mm above the deck, "
        "C25/30 concrete. "
        f"Spans evaluated: {SPANS[0]:.1f}--{SPANS[-1]:.1f}\\,m. "
        f"Results cover {len(pct_min_all)} partially-connected deck"
        r"$\,\times\,$" "span combinations.\n\n"
    )
    f.write("\\begin{table}[h!]\n\\centering\\small\n")
    f.write("\\caption{Percentage increase in composite moment resistance $M_\\mathrm{Rd}$ "
            "relative to embossment-only shear bond.}\n")
    f.write("\\label{tab:ea_summary}\n")
    f.write("\\begin{tabular}{lccc}\n\\toprule\n")
    f.write("Approach & Min\\,(\\%) & Max\\,(\\%) & Mean\\,(\\%) \\\\\n\\midrule\n")
    for label, data in [
        (f"End anchorage, $k_\\phi = {K_PHI_MIN:.1f}$ (min, $a = 1.5\\,d_{{do}}$)", pct_min_all),
        (f"End anchorage, $k_\\phi = {K_PHI_MAX:.1f}$ (max)",                        pct_max_all),
        ("Full shear connection (upper bound)",                                        pct_full_all),
    ]:
        mn  = min(data);  mx  = max(data);  avg = sum(data) / len(data)
        f.write(f"{label} & {mn:.1f} & {mx:.1f} & {avg:.1f} \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")

    # ── Part 1: Full detail table (longtable, key decks only) ─────────────────
    f.write("\\subsubsection*{Part~1 --- Full moment capacity table}\n\n")
    f.write(
        "Table~\\ref{tab:ea_full} lists $M_\\mathrm{Rd}$ at each partially-connected "
        "span for each deck gauge, comparing embossments only, end anchorage "
        "(min and max $k_\\phi$), and full connection.\n\n"
    )
    f.write("\\begin{longtable}{llcccccrrrr}\n")
    f.write("\\caption{End-anchorage moment capacity comparison for all "
            r"$\tau_u$-method decks.} "
            "\\label{tab:ea_full}\\\\\n")
    f.write("\\toprule\n")
    f.write("Deck & $L$\\,(m) & $\\eta_\\mathrm{emb}$ & "
            "$M_\\mathrm{emb}$ & $M_\\mathrm{anch,min}$ & $M_\\mathrm{anch,max}$ & "
            "$M_\\mathrm{pl}$ & "
            "$\\Delta_\\mathrm{min}$\\,(\\%) & $\\Delta_\\mathrm{max}$\\,(\\%) & "
            "$\\Delta_\\mathrm{full}$\\,(\\%) \\\\\n")
    f.write(" & & & \\multicolumn{4}{c}{[kNm/m]} & & & \\\\\n")
    f.write("\\midrule\n\\endfirsthead\n")
    f.write("\\toprule\n")
    f.write("Deck & $L$\\,(m) & $\\eta_\\mathrm{emb}$ & "
            "$M_\\mathrm{emb}$ & $M_\\mathrm{anch,min}$ & $M_\\mathrm{anch,max}$ & "
            "$M_\\mathrm{pl}$ & "
            "$\\Delta_\\mathrm{min}$\\,(\\%) & $\\Delta_\\mathrm{max}$\\,(\\%) & "
            "$\\Delta_\\mathrm{full}$\\,(\\%) \\\\\n")
    f.write("\\midrule\n\\endhead\n")
    f.write("\\bottomrule\n\\endfoot\n")

    prev_deck = None
    for (dn, span, eta, M_emb, M_mn, M_mx, M_pl, pm, px, pf) in latex_p1_rows:
        if dn != prev_deck and prev_deck is not None:
            f.write("\\midrule\n")
        prev_deck = dn
        label = dn.replace("_", r"\_") if dn != prev_deck else ""
        f.write(f"{dn} & {span:.1f} & {eta:.3f} & "
                f"{M_emb:.2f} & {M_mn:.2f} & {M_mx:.2f} & {M_pl:.2f} & "
                f"{pm:.1f} & {px:.1f} & {pf:.1f} \\\\\n")

    f.write("\\end{longtable}\n\n")

    # ── Part 2: Pass/fail matrix table ────────────────────────────────────────
    f.write("\\subsubsection*{Part~2 --- Longitudinal shear pass/fail at design $M_\\mathrm{Ed}$}\n\n")
    f.write(
        "Table~\\ref{tab:ea_passfail} shows whether each deck family passes "
        "the longitudinal shear check at office loading "
        f"($q_k = {QK_KNM2}+{PARTITION_KNM2}$\\,kN/m\\textsuperscript{{2}}, "
        f"SDL $= {FINISHES_KNM2+SERVICES_KNM2:.3f}$\\,kN/m\\textsuperscript{{2}}) "
        "without studs~(P), only with studs~(P*), or not at all~(F). "
        "Section depths are taken from the optimiser results where available "
        f"(fallback: $h = h_p + {H_C_REP:.0f}$\\,mm). "
        "\\textbf{P*} indicates spans where end anchorage is both necessary and sufficient; "
        "the stud GWP cost at those spans is listed in "
        "Table~\\ref{tab:ea_gwp}.\n\n"
    )
    span_cols = p2_spans
    f.write("\\begin{table}[h!]\n\\centering\\small\n")
    f.write("\\caption{Longitudinal shear check result by deck family and span. "
            "\\textbf{P}~=~pass without studs; "
            "\\textbf{P*}~=~pass only with end-anchorage studs; "
            "F~=~fail even with studs.}\n")
    f.write("\\label{tab:ea_passfail}\n")
    col_spec = "l" + "c" * len(span_cols)
    f.write(f"\\begin{{tabular}}{{{col_spec}}}\n\\toprule\n")
    span_header = " & ".join(f"{s:.0f}\\,m" for s in span_cols)
    f.write(f"Deck family & {span_header} \\\\\n\\midrule\n")
    for fam in p2_families:
        cells = " & ".join(_cell(fam, s) for s in span_cols)
        f.write(f"{fam} & {cells} \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")

    # GWP cost where studs needed
    stud_needed = [(fam, span, gwp)
                   for (fam, span, h_mm, M_Ed, M_Rd_emb, M_Rd_stud, Na, se, ss, gwp)
                   in latex_p2_rows if gwp is not None]
    if stud_needed:
        f.write("\\begin{table}[h!]\n\\centering\\small\n")
        f.write("\\caption{GWP cost of end-anchorage studs at spans where they are "
                "required and sufficient (P* in Table~\\ref{tab:ea_passfail}).}\n")
        f.write("\\label{tab:ea_gwp}\n")
        f.write("\\begin{tabular}{llcc}\n\\toprule\n")
        f.write("Deck family & $L$\\,(m) & "
                "$\\mathrm{GWP}_\\mathrm{stud}$\\,(kg\\,CO\\textsubscript{2}e/m\\textsuperscript{2}) & "
                "As\\,\\%\\,of\\,total\\,GWP\\,range \\\\\n\\midrule\n")
        for (fam, span, gwp) in stud_needed:
            pct_lo = gwp / 70.0 * 100.0
            pct_hi = gwp / 40.0 * 100.0
            f.write(f"{fam} & {span:.0f} & {gwp:.2f} & "
                    f"{pct_lo:.0f}--{pct_hi:.0f}\\,\\% \\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")

    # ── Part 3: Case study ─────────────────────────────────────────────────────
    if _cs is not None:
        cs = _cs   # use captured values — loop variables above have been overwritten
        util_emb  = cs['M_Ed'] / cs['M_Rd_emb']
        util_stud = cs['M_Ed'] / cs['M_Rd_stud']
        pct_increase = (cs['M_Rd_stud'] - cs['M_Rd_emb']) / cs['M_Rd_emb'] * 100.0

        f.write("\\subsubsection*{Part~3 --- Case study: "
                f"{cs['family']}, C25/30, $L = {cs['span']:.0f}$\\,m}}\n\n")
        f.write(
            f"This example uses a {cs['family']} deck (representative gauge: "
            f"$t = {cs['t']:.1f}$\\,mm, $h_p = {cs['h_p']:.0f}$\\,mm, "
            f"$A_{{pe}} = {cs['A_pe']:.0f}$\\,mm\\textsuperscript{{2}}/m, "
            f"$f_{{ypd}} = {cs['f_ypd']:.0f}$\\,N/mm\\textsuperscript{{2}}, "
            f"$\\tau_{{u,Rk}} = {cs['tau_uRk']:.0f}$\\,N/m\\textsuperscript{{2}}, "
            f"rib spacing $l_1 = {cs['l_1']:.0f}$\\,mm) "
            f"at $L = {cs['span']:.0f}$\\,m simply-supported, "
            f"$h = {cs['h_mm']:.0f}$\\,mm, C25/30 ($f_{{cd}} = {FCD_MPA:.2f}$\\,MPa).\n\n"
        )
        f.write(
            "\\paragraph{Shear span and bond force.}\n"
            f"The shear span is $L_x = L/2 = {cs['span']/2:.2f}$\\,m. "
            "The embossment bond force per unit width is\n"
            "\\begin{equation}\n"
            "    N_{c,\\mathrm{bond}} = \\frac{\\tau_{u,Rk}}{\\gamma_{vs}}\\,b\\,L_x\n"
            f"    = \\frac{{{cs['tau_uRk']:.0f}}}{{1.25}}"
            f"\\times 1000\\times{cs['L_x']:.0f}\n"
            f"    = {cs['N_c0']/1000:.1f}\\text{{\\,kN/m}}\n"
            "\\end{equation}\n"
            f"giving a degree of connection $\\eta = {cs['eta_emb']:.3f}$ "
            f"(full connection: $N_{{cf}} = {cs['N_cf']/1000:.0f}$\\,kN/m).\n\n"
        )
        f.write(
            "\\paragraph{End-anchorage resistance (Eq~9.10).}\n"
            f"For a \\o\\,{D_SHANK:.0f}\\,mm stud: "
            f"$d_{{do}} = 1.1\\times{D_SHANK:.0f} = {D_DO:.1f}$\\,mm, "
            f"$k_\\phi = 1 + {A_END_MM:.0f}/{D_DO:.1f} = {K_PHI:.2f}$.\n"
            "\\begin{align}\n"
            f"    P_{{pb,Rd}} &= k_\\phi\\,d_{{do}}\\,t\\,f_{{ypd}} "
            f"= {K_PHI:.2f}\\times{D_DO:.1f}\\times{cs['t']:.1f}\\times{cs['f_ypd']:.0f} "
            f"= {cs['P_pb_Rd']:.0f}\\text{{\\,N/stud}} \\\\\n"
            f"    N_a &= \\frac{{1000}}{{{cs['l_1']:.0f}}}\\times{cs['P_pb_Rd']:.0f} "
            f"= {cs['Na']/1000:.1f}\\text{{\\,kN/m}}\n"
            "\\end{align}\n\n"
        )
        f.write(
            "\\paragraph{Design moment and resistance.}\n"
            f"The ULS design moment is $M_{{Ed}} = {cs['M_Ed']:.2f}$\\,kNm/m. "
            "The partial-connection resistance without and with end anchorage:\n"
            "\\begin{align}\n"
            f"    M_{{Rd,emb}}  &= {cs['M_Rd_emb']:.2f}\\text{{\\,kNm/m}}"
            f"\\quad\\Rightarrow\\quad \\eta_{{util}} = {util_emb:.3f}"
            f"\\quad \\textbf{{({'PASS' if cs['M_Ed']<=cs['M_Rd_emb'] else 'FAIL'})}} \\\\\n"
            f"    M_{{Rd,stud}} &= {cs['M_Rd_stud']:.2f}\\text{{\\,kNm/m}}"
            f"\\quad\\Rightarrow\\quad \\eta_{{util}} = {util_stud:.3f}"
            f"\\quad \\textbf{{({'PASS' if cs['M_Ed']<=cs['M_Rd_stud'] else 'FAIL'})}} \n"
            "\\end{align}\n"
            f"End anchorage increases $M_{{Rd}}$ by {pct_increase:.1f}\\,\\% "
            f"at a GWP cost of {cs['gwp_st']:.2f}\\,kg\\,CO\\textsubscript{{2}}e/m\\textsuperscript{{2}} "
            f"({cs['gwp_st']/STUD_GWP:.2f}\\,kg\\,steel/m\\textsuperscript{{2}}).\n\n"
        )

print(f"\nLaTeX output saved to '{LATEX_FILE}'")
