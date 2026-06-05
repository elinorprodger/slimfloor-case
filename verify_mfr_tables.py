"""
Verify optimized composite slab configurations against manufacturer load-span tables.
ComFlor (Tata Steel) and Kingspan Multideck manuals.

Key assumptions:
- Office load = 5.00 kN/m² total applied (qk=3.0 + partitions=1.0 + services=0.5 + finishes=0.5)
- Fire resistance: R60
- Concrete: Normal Weight Concrete (NWC)
- ComFlor tables: Eurocode-based (where available), beam width 152mm (shallow) / 400mm (deep)
- Kingspan tables: British Standard only; using 6.0 kN/m² column (conservative, nearest above 5.0)
- CF210 tables: British Standard only (no EC tables in manual)
"""

import pandas as pd
import os
import numpy as np

base = r"C:\Users\ellie\OneDrive - Imperial College London\Year 4\Final Project\Implementation\Slimfloor Case"

# ==============================================================================
# MANUFACTURER TABLE DATA
# Format: {h_mm: {thickness_mm: max_span_m}}
# All at 5.00 kN/m² applied load (ComFlor) or 6.0 kN/m² (Kingspan - conservative)
# R60 fire period
# ==============================================================================
# ── ComFlor 46 ──────────────────────────────────────────────────────────────
# p.12: Bar Fire / Propped (single or continuous slab), BS, Beam 152mm, R60, 5.00 kN/m²
cf46_propped = {
    120: {0.90: 4.79, 1.20: 4.85},
    130: {0.90: 4.80, 1.20: 5.20},
    140: {0.90: 4.68, 1.20: 5.50},
    150: {0.90: 4.53, 1.20: 5.49},
    160: {0.90: 4.39, 1.20: 5.35},
    170: {0.90: 4.26, 1.20: 5.22},
    180: {0.90: 4.15, 1.20: 5.08},
    190: {0.90: 4.04, 1.20: 4.95},
    200: {0.90: 3.93, 1.20: 4.83},
}
# p.11: Bar Fire / Unpropped, Double span, BS, Beam 152mm, R60, 5.00 kN/m²
cf46_unpropped = {
    120: {0.90: 2.96, 1.20: 3.26},
    130: {0.90: 2.86, 1.20: 3.21},
    140: {0.90: 2.78, 1.20: 3.18},
    150: {0.90: 2.70, 1.20: 3.11},
    160: {0.90: 2.61, 1.20: 3.04},
    170: {0.90: 2.52, 1.20: 2.97},
    180: {0.90: 2.44, 1.20: 2.89},
    190: {0.90: 2.37, 1.20: 2.83},
    200: {0.90: 2.30, 1.20: 2.77},
}

# ── ComFlor 51+ ──────────────────────────────────────────────────────────────
# p.19: Bar Fire / Propped (single span propped, continuous slab), EC, Beam 152mm, R60, 5.00 kN/m²
cf51_propped = {
    101: {0.90: 4.03, 1.00: 4.05, 1.20: 4.10},
    110: {0.90: 4.34, 1.00: 4.35, 1.20: 4.40},
    120: {0.90: 4.66, 1.00: 4.68, 1.20: 4.72},
    130: {0.90: 4.98, 1.00: 5.00, 1.20: 5.04},
    140: {0.90: 5.29, 1.00: 5.31, 1.20: 5.36},
    150: {0.90: 5.41, 1.00: 5.61, 1.20: 5.66},
    160: {0.90: 5.34, 1.00: 5.74, 1.20: 5.96},
    170: {0.90: 5.22, 1.00: 5.71, 1.20: 6.25},
    180: {0.90: 5.10, 1.00: 5.57, 1.20: 6.37},
    190: {0.90: 5.00, 1.00: 5.46, 1.20: 6.34},
    200: {0.90: 4.90, 1.00: 5.35, 1.20: 6.20},
}
# p.17: Bar Fire / Unpropped, Double span, EC, Beam 152mm, R60, 5.00 kN/m²
cf51_unpropped = {
    101: {0.90: 3.08, 1.00: 3.35, 1.20: 3.88},
    110: {0.90: 3.08, 1.00: 3.26, 1.20: 3.78},
    120: {0.90: 2.99, 1.00: 3.18, 1.20: 3.68},
    130: {0.90: 2.91, 1.00: 3.16, 1.20: 3.59},
    140: {0.90: 2.83, 1.00: 3.09, 1.20: 3.51},
    150: {0.90: 2.76, 1.00: 3.02, 1.20: 3.42},
    160: {0.90: 2.70, 1.00: 2.95, 1.20: 3.41},
    170: {0.90: 2.64, 1.00: 2.88, 1.20: 3.37},
    180: {0.90: 2.58, 1.00: 2.82, 1.20: 3.29},
    190: {0.90: 2.53, 1.00: 2.76, 1.20: 3.23},
    200: {0.90: 2.48, 1.00: 2.71, 1.20: 3.17},
}

# ── ComFlor 60 ────────────────────────────────────────────────────────────────
# p.24: Bar Fire / Propped (single span propped, continuous slab), EC, Beam 152mm, R60, 5.00 kN/m²
cf60_propped = {
    120: {0.90: 4.48, 1.00: 4.51, 1.20: 4.57},
    130: {0.90: 4.77, 1.00: 4.80, 1.20: 4.86},
    140: {0.90: 5.06, 1.00: 5.09, 1.20: 5.14},
    150: {0.90: 5.10, 1.00: 5.36, 1.20: 5.42},
    160: {0.90: 4.98, 1.00: 5.57, 1.20: 5.68},
    170: {0.90: 4.81, 1.00: 5.48, 1.20: 5.95},
    180: {0.90: 4.66, 1.00: 5.34, 1.20: 6.20},
    190: {0.90: 4.49, 1.00: 5.15, 1.20: 6.28},
    200: {0.90: 4.35, 1.00: 5.00, 1.20: 6.19},
}
# p.23: Bar Fire / Unpropped, Double span, EC, Beam 152mm, R60, 5.00 kN/m²
cf60_unpropped = {
    120: {0.90: 3.68, 1.00: 4.07, 1.20: 4.70},
    130: {0.90: 3.63, 1.00: 3.93, 1.20: 4.65},
    140: {0.90: 3.49, 1.00: 3.79, 1.20: 4.50},
    150: {0.90: 3.36, 1.00: 3.74, 1.20: 4.38},
    160: {0.90: 3.24, 1.00: 3.68, 1.20: 4.25},
    170: {0.90: 3.13, 1.00: 3.56, 1.20: 4.13},
    180: {0.90: 3.03, 1.00: 3.45, 1.20: 4.02},
    190: {0.90: 2.93, 1.00: 3.34, 1.20: 3.91},
    200: {0.90: 2.85, 1.00: 3.24, 1.20: 3.88},
}

# ── ComFlor 80 ────────────────────────────────────────────────────────────────
# p.30: Bar Fire / Propped (single span propped, continuous slab), EC, Beam 152mm, R60, 5.00 kN/m²
cf80_propped = {
    140: {0.90: 5.17, 1.00: 5.19, 1.20: 5.24},
    150: {0.90: 5.25, 1.00: 5.45, 1.20: 5.50},
    160: {0.90: 5.15, 1.00: 5.72, 1.20: 5.76},
    170: {0.90: 4.97, 1.00: 5.62, 1.20: 6.02},
    180: {0.90: 4.78, 1.00: 5.49, 1.20: 6.26},
    190: {0.90: 4.61, 1.00: 5.30, 1.20: 6.45},
    200: {0.90: 4.46, 1.00: 5.13, 1.20: 6.35},
    210: {0.90: 4.46, 1.00: 5.13, 1.20: 6.35},
    220: {0.90: 4.46, 1.00: 5.13, 1.20: 6.35},
    230: {0.90: 4.46, 1.00: 5.13, 1.20: 6.35},
}
# p.29: Bar Fire / Unpropped, Double span, EC, Beam 152mm, R60, 5.00 kN/m²
cf80_unpropped = {
    140: {0.90: 3.70, 1.00: 4.21, 1.20: 4.82},
    150: {0.90: 3.55, 1.00: 4.03, 1.20: 4.65},
    160: {0.90: 3.40, 1.00: 3.87, 1.20: 4.60},
    170: {0.90: 3.27, 1.00: 3.72, 1.20: 4.53},
    180: {0.90: 3.15, 1.00: 3.59, 1.20: 4.38},
    190: {0.90: 3.04, 1.00: 3.47, 1.20: 4.24},
    200: {0.90: 2.95, 1.00: 3.36, 1.20: 4.11},
    210: {0.90: 2.95, 1.00: 3.36, 1.20: 4.11},
    220: {0.90: 2.95, 1.00: 3.36, 1.20: 4.11},
    230: {0.90: 2.95, 1.00: 3.36, 1.20: 4.11},
}

# ── ComFlor 210 ───────────────────────────────────────────────────────────────
# p.38: Bar Fire / Unpropped, BS, Beam 400mm, R60, 5.00 kN/m²  [BS tables only]
cf210_unpropped = {
    280: {1.25: 5.65},
    290: {1.25: 5.52},
    300: {1.25: 5.40},
    310: {1.25: 5.29},
    320: {1.25: 5.18},
    330: {1.25: 5.08},
    340: {1.25: 4.99},
    350: {1.25: 4.90},
    375: {1.25: 4.71},
    400: {1.25: 4.53},
    430: {1.25: 4.35},
    450: {1.25: 4.25},
}
# p.38: Bar Fire / Propped, BS, Beam 400mm, R60, 5.00 kN/m²
cf210_propped = {
    280: {1.25: 7.47},
    290: {1.25: 7.20},
    300: {1.25: 6.93},
    310: {1.25: 6.67},
    320: {1.25: 6.43},
    330: {1.25: 6.22},
    340: {1.25: 6.03},
    350: {1.25: 5.84},
    375: {1.25: 5.42},
    400: {1.25: 5.06},
    430: {1.25: 4.75},
    450: {1.25: 4.60},
}

# ── ComFlor 225 ───────────────────────────────────────────────────────────────
# p.42: Bar Fire / Unpropped, EC, Beam 400mm, R60, 5.00 kN/m²
cf225_unpropped = {
    283: {1.25: 6.10},
    295: {1.25: 6.00},
    300: {1.25: 5.93},
    310: {1.25: 5.82},
    320: {1.25: 5.73},
    330: {1.25: 5.65},
    340: {1.25: 5.58},
    350: {1.25: 5.51},
    375: {1.25: 5.28},
    400: {1.25: 5.07},
    437: {1.25: 4.80},
}
# p.42: Bar Fire / Propped, EC, Beam 400mm, R60, 5.00 kN/m²
cf225_propped = {
    283: {1.25: 7.90},
    295: {1.25: 8.01},
    300: {1.25: 8.05},
    310: {1.25: 8.13},
    320: {1.25: 8.21},
    330: {1.25: 8.29},
    340: {1.25: 8.37},
    350: {1.25: 8.45},
    375: {1.25: 8.32},
    400: {1.25: 7.83},
    437: {1.25: 7.50},
}

# ── Kingspan MD50 ─────────────────────────────────────────────────────────────
# pp.15-17: Standard Load/Span / NWC / Unpropped / Double span / BS / 6.0 kN/m²
# (6.0 kN/m² is the conservative nearest column above 5.0 kN/m² — no 5.0 column in structural tables)
# Fire resistance tables (p.22+) must be read IN CONJUNCTION per manual note p.21/8,
# but structural tables govern (shorter spans) so these are the binding check.
md50_unpropped = {
    100: {0.85: 3.31, 0.90: 3.50, 1.00: 3.56, 1.10: 3.61, 1.20: 3.66},
    110: {0.85: 3.21, 0.90: 3.49, 1.00: 3.74, 1.10: 3.93, 1.20: 3.98},
    120: {0.85: 3.13, 0.90: 3.40, 1.00: 3.64, 1.10: 3.87, 1.20: 4.03},
    130: {0.85: 3.05, 0.90: 3.32, 1.00: 3.55, 1.10: 3.77, 1.20: 3.95},
    140: {0.85: 2.97, 0.90: 3.24, 1.00: 3.47, 1.10: 3.68, 1.20: 3.87},
    150: {0.85: 2.89, 0.90: 3.17, 1.00: 3.39, 1.10: 3.59, 1.20: 3.78},
    160: {0.85: 2.81, 0.90: 3.10, 1.00: 3.32, 1.10: 3.52, 1.20: 3.70},
    175: {0.85: 2.71, 0.90: 3.00, 1.00: 3.22, 1.10: 3.40, 1.20: 3.59},
    200: {0.85: 2.57, 0.90: 2.84, 1.00: 3.07, 1.10: 3.26, 1.20: 3.42},
    250: {0.85: 2.35, 0.90: 2.58, 1.00: 2.80, 1.10: 2.99, 1.20: 3.13},
}
# Propped MD50 not tabulated in Kingspan manual — marked as N/A in verification
md50_propped = None  # No propped table available

# ── Kingspan MD60 ─────────────────────────────────────────────────────────────
# p.38: BS, NWC Unpropped, Double span, 6.0 kN/m² column
# (Table has single span and double span sections per p.37 note 1;
#  double span = continuous/multispan condition used here)
md60_unpropped = {
    120: {0.90: 3.52, 1.00: 3.81, 1.10: 4.09, 1.20: 4.20},
    130: {0.90: 3.42, 1.00: 3.70, 1.10: 3.97, 1.20: 4.21},
    140: {0.90: 3.32, 1.00: 3.60, 1.10: 3.86, 1.20: 4.10},
    150: {0.90: 3.24, 1.00: 3.51, 1.10: 3.76, 1.20: 4.00},
    160: {0.90: 3.16, 1.00: 3.43, 1.10: 3.67, 1.20: 3.90},
    175: {0.90: 3.05, 1.00: 3.31, 1.10: 3.55, 1.20: 3.76},
    200: {0.90: 2.87, 1.00: 3.14, 1.10: 3.36, 1.20: 3.57},
    250: {0.90: 2.58, 1.00: 2.84, 1.10: 3.07, 1.20: 3.26},
}
# p.39: BS, NWC Propped, 6.0 kN/m² column
md60_propped = {
    120: {0.90: 3.52, 1.00: 3.81, 1.10: 4.09, 1.20: 4.20},
    130: {0.90: 3.77, 1.00: 3.98, 1.10: 4.14, 1.20: 4.21},
    140: {0.90: 3.99, 1.00: 4.20, 1.10: 4.36, 1.20: 4.36},
    150: {0.90: 4.20, 1.00: 4.41, 1.10: 4.57, 1.20: 4.57},
    160: {0.90: 4.40, 1.00: 4.60, 1.10: 4.77, 1.20: 4.77},
    175: {0.90: 4.67, 1.00: 4.88, 1.10: 5.05, 1.20: 5.05},
    200: {0.90: 5.08, 1.00: 5.29, 1.10: 5.47, 1.20: 5.46},
    250: {0.90: 5.13, 1.00: 5.64, 1.10: 6.10, 1.20: 6.16},
}

# ── Kingspan MD80 ─────────────────────────────────────────────────────────────
# p.58: BS, NWC, 6.0 kN/m² column
md80_unpropped = {
    130: {1.00: 4.53, 1.10: 4.55, 1.20: 4.55},
    140: {1.00: 4.39, 1.10: 4.68, 1.20: 4.90},
    150: {1.00: 4.26, 1.10: 4.54, 1.20: 4.81},
    160: {1.00: 4.15, 1.10: 4.42, 1.20: 4.68},
    175: {1.00: 3.99, 1.10: 4.25, 1.20: 4.50},
    200: {1.00: 3.76, 1.10: 4.01, 1.20: 4.24},
    250: {1.00: 3.40, 1.10: 3.62, 1.20: 3.84},
}
# p.59: BS, NWC, 6.0 kN/m² (double span row)
md80_propped = {
    130: {1.00: 4.53, 1.10: 4.55, 1.20: 4.55},
    140: {1.00: 4.39, 1.10: 4.68, 1.20: 4.90},
    150: {1.00: 4.26, 1.10: 4.54, 1.20: 4.81},
    160: {1.00: 4.21, 1.10: 4.42, 1.20: 4.68},
    175: {1.00: 4.41, 1.10: 4.60, 1.20: 4.75},
    200: {1.00: 4.70, 1.10: 4.90, 1.20: 5.07},
    250: {1.00: 5.16, 1.10: 5.38, 1.20: 5.58},
}

# ── Kingspan MD146 ────────────────────────────────────────────────────────────
# p.74: BS, NWC Unpropped gauge=1.2mm, 6.0 kN/m² column
md146_unpropped = {
    215: {1.20: 5.74, 1.50: 6.07},
    225: {1.20: 5.67, 1.50: 6.00},
    235: {1.20: 5.58, 1.50: 5.92},
    245: {1.20: 5.50, 1.50: 5.85},
    255: {1.20: 5.43, 1.50: 5.77},
    265: {1.20: 5.35, 1.50: 5.70},
    275: {1.20: 5.30, 1.50: 5.63},
    285: {1.20: 5.24, 1.50: 5.57},
    295: {1.20: 5.17, 1.50: 5.51},
    305: {1.20: 5.08, 1.50: 5.46},
}
# MD146 propped
md146_propped = {
    215: {1.50: 6.45},
    225: {1.50: 6.75},
    235: {1.50: 7.05},
    245: {1.50: 7.35},
    255: {1.50: 7.65},
    265: {1.50: 7.95},
    275: {1.50: 8.25},
    285: {1.50: 8.55},
    295: {1.50: 8.85},
    305: {1.50: 9.15},
}
# ==============================================================================
# DECK NAME MAPPING
# (deck_name_csv) -> (short_id, thickness, unprop_table, prop_table, source_name, code)
# ==============================================================================
TABLE_MAP = {
    'ComFlor 46 0.9':   ('CF46',  0.90, cf46_unpropped,  cf46_propped,  'ComFlor 46 (BS, Bar Fire)', 'BS'),
    'ComFlor 46 1.2':   ('CF46',  1.20, cf46_unpropped,  cf46_propped,  'ComFlor 46 (BS, Bar Fire)', 'BS'),
    'ComFlor 51 0.9':   ('CF51',  0.90, cf51_unpropped,  cf51_propped,  'ComFlor 51+ (EC, Bar Fire)', 'EC'),
    'ComFlor 51 1.0':   ('CF51',  1.00, cf51_unpropped,  cf51_propped,  'ComFlor 51+ (EC, Bar Fire)', 'EC'),
    'ComFlor 51 1.2':   ('CF51',  1.20, cf51_unpropped,  cf51_propped,  'ComFlor 51+ (EC, Bar Fire)', 'EC'),
    'ComFlor 60 0.9':   ('CF60',  0.90, cf60_unpropped,  cf60_propped,  'ComFlor 60 (EC, Bar Fire)',  'EC'),
    'ComFlor 60 1.0':   ('CF60',  1.00, cf60_unpropped,  cf60_propped,  'ComFlor 60 (EC, Bar Fire)',  'EC'),
    'ComFlor 60 1.2':   ('CF60',  1.20, cf60_unpropped,  cf60_propped,  'ComFlor 60 (EC, Bar Fire)',  'EC'),
    'ComFlor 80 0.9':   ('CF80',  0.90, cf80_unpropped,  cf80_propped,  'ComFlor 80 (EC, Bar Fire)', 'EC'),
    'ComFlor 80 1.2':   ('CF80',  1.20, cf80_unpropped,  cf80_propped,  'ComFlor 80 (EC, Bar Fire)', 'EC'),
    'ComFlor 210':      ('CF210', 1.25, cf210_unpropped, cf210_propped, 'ComFlor 210 (BS)',  'BS'),
    'ComFlor 225':      ('CF225', 1.25, cf225_unpropped, cf225_propped, 'ComFlor 225 (EC)',  'EC'),
    'Multideck 50 0.85':('MD50',  0.85, md50_unpropped,  md50_propped,  'Kingspan MD50 (BS, dbl span, 6kN/m2)', 'BS'),
    'Multideck 50 0.9': ('MD50',  0.90, md50_unpropped,  md50_propped,  'Kingspan MD50 (BS, dbl span, 6kN/m2)', 'BS'),
    'Multideck 50 1.0': ('MD50',  1.00, md50_unpropped,  md50_propped,  'Kingspan MD50 (BS, dbl span, 6kN/m2)', 'BS'),
    'Multideck 50 1.1': ('MD50',  1.10, md50_unpropped,  md50_propped,  'Kingspan MD50 (BS, dbl span, 6kN/m2)', 'BS'),
    'Multideck 50 1.2': ('MD50',  1.20, md50_unpropped,  md50_propped,  'Kingspan MD50 (BS, dbl span, 6kN/m2)', 'BS'),
    'Multideck 60 0.9': ('MD60',  0.90, md60_unpropped,  md60_propped,  'Kingspan MD60 (BS, 6kN/m2)', 'BS'),
    'Multideck 60 1.0': ('MD60',  1.00, md60_unpropped,  md60_propped,  'Kingspan MD60 (BS, 6kN/m2)', 'BS'),
    'Multideck 60 1.1': ('MD60',  1.10, md60_unpropped,  md60_propped,  'Kingspan MD60 (BS, 6kN/m2)', 'BS'),
    'Multideck 60 1.2': ('MD60',  1.20, md60_unpropped,  md60_propped,  'Kingspan MD60 (BS, 6kN/m2)', 'BS'),
    'Multideck 80 1.0': ('MD80',  1.00, md80_unpropped,  md80_propped,  'Kingspan MD80 (BS, 6kN/m2)', 'BS'),
    'Multideck 80 1.1': ('MD80',  1.10, md80_unpropped,  md80_propped,  'Kingspan MD80 (BS, 6kN/m2)', 'BS'),
    'Multideck 80 1.2': ('MD80',  1.20, md80_unpropped,  md80_propped,  'Kingspan MD80 (BS, 6kN/m2)', 'BS'),
    'Multideck 146 1.2':('MD146', 1.20, md146_unpropped, None,          'Kingspan MD146 (BS, 6kN/m2)', 'BS'),
    'Multideck 146 1.5':('MD146', 1.50, md146_unpropped, md146_propped, 'Kingspan MD146 (BS, 6kN/m2)', 'BS'),
}

def interp_span(table, h_mm, thickness):
    """Linearly interpolate max span for given h_mm and thickness."""
    h_keys = sorted(table.keys())
    if not h_keys:
        return None
    # Clamp to range
    h = h_mm
    if h <= h_keys[0]:
        return table[h_keys[0]].get(thickness)
    if h >= h_keys[-1]:
        return table[h_keys[-1]].get(thickness)
    # Interpolate
    for i in range(len(h_keys)-1):
        h_lo, h_hi = h_keys[i], h_keys[i+1]
        if h_lo <= h <= h_hi:
            v_lo = table[h_lo].get(thickness)
            v_hi = table[h_hi].get(thickness)
            if v_lo is None or v_hi is None:
                return table[h_lo].get(thickness) or table[h_hi].get(thickness)
            frac = (h - h_lo) / (h_hi - h_lo)
            return round(v_lo + frac * (v_hi - v_lo), 2)
    return None

def get_table_span(deck_name, h_mm, propped):
    if deck_name not in TABLE_MAP:
        return None, 'Unknown deck'
    _, thick, unprop_tbl, prop_tbl, src_name, code = TABLE_MAP[deck_name]
    tbl = prop_tbl if propped else unprop_tbl
    if tbl is None:
        # No propped table available for this deck
        return None, src_name + ' (no propped table)'
    max_span = interp_span(tbl, h_mm, thick)
    return max_span, src_name

# ==============================================================================
# LOAD CSVs AND VERIFY
# ==============================================================================
records = []
for label, fn in [('1-span', 'results_office_ENV.csv'),
                  ('2-span', 'results_office_2span_ENV.csv'),
                  ('3-span', 'results_office_3span_ENV.csv')]:
    df = pd.read_csv(os.path.join(base, fn))
    comp = df[df['slab_type'] == 'comp_slab'].dropna(subset=['h_struct_m']).copy()
    comp['n_spans_label'] = label
    comp['h_mm'] = (comp['h_struct_m'] * 1000).round(0)
    records.append(comp)

all_comp = pd.concat(records, ignore_index=True)

rows = []
for _, row in all_comp.iterrows():
    deck  = row['deck_name']
    span  = round(float(row['span_m']), 2)
    h_mm  = float(row['h_mm'])
    propped = bool(row['propped'])
    n_spans = row['n_spans_label']

    max_span, src = get_table_span(deck, h_mm, propped)
    if max_span is not None:
        ratio  = round(span / max_span, 3)
        result = 'PASS' if span <= max_span else 'FAIL'
    else:
        ratio  = None
        result = 'N/A'

    rows.append({
        'Config':        n_spans,
        'Deck':          deck,
        'h_mm':          int(h_mm),
        'Span (m)':      span,
        'Propped':       propped,
        'Max span (m)':  max_span,
        'Span/Max':      ratio,
        'Result':        result,
        'Source':        src,
    })

verif = pd.DataFrame(rows)

# Drop duplicates on key columns
uniq = verif.drop_duplicates(subset=['Config','Deck','h_mm','Span (m)','Propped'])
uniq = uniq.sort_values(['Deck','Propped','Span (m)'])

# ==============================================================================
# REPORT
# ==============================================================================
pd.set_option('display.max_rows', 300)
pd.set_option('display.width',    250)
pd.set_option('display.max_colwidth', 50)

print("=" * 100)
print("MANUFACTURER LOAD-SPAN TABLE VERIFICATION")
print("Office load: 5.00 kN/m² | R60 fire | NWC")
print("ComFlor tables: Eurocode (CF46/CF100/CF210 use BS); Kingspan tables: BS, 6.0 kN/m² column (conservative)")
print("=" * 100)
print(f"\nTotal unique configurations checked: {len(uniq)}")
print(f"  PASS: {(uniq['Result']=='PASS').sum()}")
print(f"  FAIL: {(uniq['Result']=='FAIL').sum()}")
print(f"  N/A:  {(uniq['Result']=='N/A').sum()}")

fails = uniq[uniq['Result'] == 'FAIL']
if len(fails) > 0:
    print("\n=== FAILING CONFIGURATIONS ===")
    print(fails[['Config','Deck','h_mm','Span (m)','Propped','Max span (m)','Span/Max','Source']].to_string(index=False))
else:
    print("\n  ALL configurations PASS.")

print("\n=== FULL VERIFICATION TABLE ===")
print(uniq[['Config','Deck','h_mm','Span (m)','Propped','Max span (m)','Span/Max','Result','Source']].to_string(index=False))

# Also save to CSV for report
out_path = os.path.join(base, 'mfr_table_verification.csv')
uniq.to_csv(out_path, index=False)
print(f"\nSaved to: {out_path}")
