"""
Slab-alone natural frequency check for all optimised composite slab configurations.

Reads the three optimisation result CSVs, extracts every unique composite slab
cross-section, and checks the slab-alone natural frequency against the 4 Hz
office criterion.

Method:
  Treat the slab as a 1 m wide simply-supported strip spanning b between
  rigid supports (secondary beams assumed infinitely stiff).  This isolates
  the slab's own contribution and is what ComFlor's design tool reports.

  Slab stiffness: SCI P354 Table 7.2 power-law fit
      I_s = coeff * h_s^3.7  (NWC, E_c,dyn = 38 GPa)
  Deflection:
      delta = 5 * m * g * b^4 / (384 * EI_s)   [m]
  Natural frequency:
      f1 = 18 / sqrt(delta_mm)                  [Hz]   (SCI P354 Eq 4)

  Conservative: the full combined floor frequency (slab + beams) is always
  lower than the slab-alone frequency, so passing here is a necessary but not
  sufficient condition.  The slab-alone check identifies slabs that would
  struggle even with perfect support.

Threshold:
  f1 >= 4 Hz  — standard simple criterion for composite office floors
               (SCI P354 para 4.1.5; used by ComFlor design software)
  f1 >= 3 Hz  — absolute SCI P354 minimum (floors below this resonant
               with the 3rd harmonic of walking)

Loading (office case, consistent with optimiser):
  self-weight  = concrete (profiled) + deck sheeting, matching the main optimiser:
                   h_conc_avg = (h_mm - h_p) + conc_vol_trough*1000   [mm]
                   g_slab     = 24 kN/m3 * h_conc_avg/1000            [kN/m2]
                   g_deck     = deck.weight_kN                        [kN/m2]
                 (h_p, conc_vol_trough, weight_kN read from the database per
                  gauge-specific deck; 24 kN/m3 = NWC dry density)
  SDL          = 1.105 kN/m2  (finishes 0.355 from office floor layers + services 0.75)
  imposed      = 10% of Q_k = 0.25 kN/m2  (Q_k = 2.5; SCI P354 para 4.1.2)
"""

import math
import os
import sqlite3
import pandas as pd

base = r"C:\Users\ellie\OneDrive - Imperial College London\Year 4\Final Project\Implementation\Slimfloor Case"

# ==============================================================================
# LOADING CONSTANTS  (office case — match optimiser settings)
# ==============================================================================
SDL_kNm2        = 1.105  # finishes 0.355 (office floor layers) + services 0.75 [kN/m2]
Q_k_kNm2        = 2.5    # characteristic imposed (office, excl. partitions) [kN/m2]
rho_c_kNm3      = 24.0   # NWC dry density (matches optimiser dry_density_load) [kN/m3]
E_s_Pa          = 210.0e9  # steel Young's modulus        [Pa]
E_c_dyn_Pa      = 38.0e9   # dynamic concrete modulus NWC [Pa]
g_acc           = 9.81   # m/s2

# ==============================================================================
# DECK SELF-WEIGHT DATA  (from database, per gauge-specific deck — matches optimiser)
# ==============================================================================
# Pre-load h_p, conc_vol_trough and deck weight for every deck so the slab
# self-weight reproduces the main optimiser exactly:
#   h_conc_avg = (h_mm - h_p) + conc_vol_trough*1000   [mm]
#   g_slab     = rho_c_kNm3 * h_conc_avg/1000          [kN/m2]  (concrete, profiled)
#   g_deck     = weight_kN                              [kN/m2]
_DB_PATH = os.path.join(base, "database_270426.db")
DECK_DB = {}   # deck_name -> (h_p_mm, conc_vol_trough_m3m2, weight_kN)
with sqlite3.connect(_DB_PATH) as _con:
    for _name, _h_p, _cvt, _wkn in _con.execute(
            "SELECT Deck, h_p, conc_vol_trough, weight_kN FROM sheeting_prop"):
        DECK_DB[_name] = (_h_p, _cvt, _wkn)

# ==============================================================================
# DECK PROPERTIES  (h_p in mm, profile type for Table 7.2 coefficient)
# ==============================================================================
DECK_PROPS = {
    'ComFlor 46':    (46,  'reentrant'),
    'ComFlor 51':    (51,  'trapezoidal'),
    'ComFlor 60':    (60,  'trapezoidal'),
    'ComFlor 80':    (80,  'trapezoidal'),
    'ComFlor 210':   (210, 'trapezoidal'),
    'ComFlor 225':   (225, 'trapezoidal'),
    'Multideck 50':  (50,  'trapezoidal'),
    'Multideck 60':  (60,  'trapezoidal'),
    'Multideck 80':  (80,  'trapezoidal'),
    'Multideck 146': (146, 'trapezoidal'),
}

def get_deck_props(deck_name):
    """Return (h_p_mm, profile_type), stripping gauge suffix."""
    for key in DECK_PROPS:
        if deck_name.startswith(key):
            return DECK_PROPS[key]
    return None, None

# ==============================================================================
# SCI P354 TABLE 7.2 — SLAB SECOND MOMENT OF AREA
# ==============================================================================
def I_s_mm4_per_m(h_mm, h_p_mm, profile_type):
    """Dynamic I_s [mm4/m] per unit width.  SCI P354 Table 7.2 power-law fit:
        I_s = coeff * h_s^3.7   (NWC)
    Coefficients interpolated by h_p for trapezoidal decks.
    """
    if profile_type == 'reentrant':
        coeff = 0.37
    else:
        T = [(60, 0.23), (80, 0.19), (225, 0.05)]
        if h_p_mm <= T[0][0]:
            coeff = T[0][1]
        elif h_p_mm >= T[-1][0]:
            coeff = T[-1][1]
        else:
            for i in range(len(T) - 1):
                if T[i][0] <= h_p_mm <= T[i+1][0]:
                    t = (h_p_mm - T[i][0]) / (T[i+1][0] - T[i][0])
                    coeff = T[i][1] + t * (T[i+1][1] - T[i][1])
                    break
    return coeff * float(h_mm) ** 3.7

# ==============================================================================
# SLAB-ALONE FREQUENCY CHECK
# ==============================================================================
def slab_frequency(h_mm, b_m, deck_name):
    """Return (f1_Hz, delta_mm, I_s_mm4, EI_s_MNm2, m_floor_kg_m2) or None."""
    h_p_mm, profile_type = get_deck_props(deck_name)
    db = DECK_DB.get(deck_name)
    if h_p_mm is None or db is None:
        return None
    h_p_db, conc_vol_trough, deck_weight_kN = db

    I_s   = I_s_mm4_per_m(h_mm, h_p_mm, profile_type)      # mm4/m
    EI_s  = E_s_Pa * I_s * 1.0e-12                          # N.m2/m

    # Self-weight matches the main optimiser: profiled (average) concrete depth
    # plus deck sheeting (uses the same database h_p as the optimiser).
    h_c        = max(h_mm - h_p_db, 1.0)                    # mm concrete above deck
    h_conc_avg = h_c + conc_vol_trough * 1000.0            # mm equivalent solid depth
    g_slab     = rho_c_kNm3 * h_conc_avg / 1000.0          # kN/m2  concrete (profiled)
    g_self     = g_slab + deck_weight_kN                   # kN/m2  concrete + deck
    m_floor    = (g_self + SDL_kNm2 + 0.1 * Q_k_kNm2) * 1000.0 / g_acc  # kg/m2

    delta_m  = 5.0 * m_floor * g_acc * b_m**4 / (384.0 * EI_s)
    delta_mm = delta_m * 1000.0
    f1       = 18.0 / math.sqrt(max(delta_mm, 1e-9))

    return f1, delta_mm, I_s, EI_s / 1e6, m_floor

# ==============================================================================
# LOAD OPTIMISATION RESULTS
# ==============================================================================
records = []
for label, fn in [('1-span', 'results_office_ENV.csv'),
                  ('2-span', 'results_office_2span_ENV.csv'),
                  ('3-span', 'results_office_3span_ENV.csv')]:
    path = os.path.join(base, fn)
    df   = pd.read_csv(path)
    comp = df[(df['slab_type'] == 'comp_slab')].dropna(subset=['h_struct_m']).copy()
    comp['n_spans'] = label
    comp['h_mm']    = (comp['h_struct_m'] * 1000).round(0)
    records.append(comp)

all_comp = pd.concat(records, ignore_index=True)

# Drop duplicates: same deck / h_mm / span across the three files
uniq = (all_comp
        .drop_duplicates(subset=['deck_name', 'h_mm', 'span_m'])
        .sort_values(['deck_name', 'span_m', 'h_mm'])
        .reset_index(drop=True))

# ==============================================================================
# RUN CHECK
# ==============================================================================
rows = []
for _, row in uniq.iterrows():
    deck  = row['deck_name']
    h_mm  = float(row['h_mm'])
    b_m   = float(row['span_m'])

    result = slab_frequency(h_mm, b_m, deck)
    if result is None:
        rows.append({'deck': deck, 'h_mm': h_mm, 'b_m': b_m,
                     'I_s': None, 'EI_s': None, 'm': None,
                     'delta': None, 'f1': None,
                     'pass_4hz': None, 'pass_3hz': None, 'note': 'Unknown deck'})
        continue

    f1, delta_mm, I_s, EI_s_MN, m_floor = result
    rows.append({
        'deck':     deck,
        'h_mm':     h_mm,
        'b_m':      b_m,
        'I_s':      I_s,
        'EI_s':     EI_s_MN,
        'm':        m_floor,
        'delta':    delta_mm,
        'f1':       f1,
        'pass_4hz': f1 >= 4.0,
        'pass_3hz': f1 >= 3.0,
        'note':     '',
    })

res = pd.DataFrame(rows)

# ==============================================================================
# PRINT RESULTS
# ==============================================================================
HDR = (f"{'Deck':<22} {'h':>5} {'b':>5}  "
       f"{'Is (mm4/m)':>12} {'EIs (MN.m2/m)':>14} "
       f"{'m (kg/m2)':>10} {'d (mm)':>8} {'f1 (Hz)':>8}  4Hz  3Hz")
SEP = '-' * len(HDR)

print()
print('=' * len(HDR))
print('COMPOSITE SLAB NATURAL FREQUENCY — SLAB ALONE ON RIGID SUPPORTS')
print(f'SDL={SDL_kNm2} kN/m2  Q_k={Q_k_kNm2} kN/m2  '
      f'rho_c={rho_c_kNm3} kN/m3 (dry)  E_s=210 GPa  E_c,dyn=38 GPa')
print('Self-weight = profiled concrete (h_conc_avg) + deck, matching the optimiser')
print('I_s from SCI P354 Table 7.2  |  delta = 5*m*g*b^4 / (384*EI_s)  |  f1 = 18/sqrt(delta_mm)')
print('Thresholds: 4 Hz (office simple criterion)  |  3 Hz (absolute minimum)')
print('=' * len(HDR))

prev_deck = None
for _, row in res.iterrows():
    if row['deck'] != prev_deck:
        print(SEP)
        prev_deck = row['deck']

    if row['f1'] is None:
        print(f"{row['deck']:<22} {row['h_mm']:>5.0f} {row['b_m']:>5.1f}  "
              f"{'N/A':>12} {'N/A':>14} {'N/A':>10} {'N/A':>8} {'N/A':>8}  N/A  N/A")
        continue

    p4 = 'PASS' if row['pass_4hz'] else 'FAIL'
    p3 = 'PASS' if row['pass_3hz'] else 'FAIL'
    print(f"{row['deck']:<22} {row['h_mm']:>5.0f} {row['b_m']:>5.1f}  "
          f"{row['I_s']:>12.0f} {row['EI_s']:>14.2f} "
          f"{row['m']:>10.1f} {row['delta']:>8.3f} {row['f1']:>8.3f}  {p4}  {p3}")

print(SEP)

# Summary
valid = res[res['f1'].notna()]
n_total = len(valid)
n_pass4 = valid['pass_4hz'].sum()
n_fail4 = (~valid['pass_4hz']).sum()
n_pass3 = valid['pass_3hz'].sum()

print(f'\nConfigurations checked: {n_total}  (+ {res["f1"].isna().sum()} unknown decks)')
print(f'  4 Hz criterion:  {n_pass4} PASS  {n_fail4} FAIL')
print(f'  3 Hz criterion:  {n_pass3} PASS  {n_total - n_pass3} FAIL')

if n_fail4 > 0:
    print('\nFailing configurations (f1 < 4 Hz):')
    fails = valid[~valid['pass_4hz']].sort_values('f1')
    for _, r in fails.iterrows():
        print(f"  {r['deck']:<22}  h={r['h_mm']:>5.0f} mm  b={r['b_m']:.1f} m  "
              f"f1={r['f1']:.3f} Hz  delta={r['delta']:.2f} mm")
else:
    print('\n  All configurations pass the 4 Hz slab-alone criterion.')
print()

# Save CSV
out = os.path.join(base, 'vibration_slab_alone.csv')
res.to_csv(out, index=False)
print(f'Saved: {out}')
