# =============================================================================
#  TWO-SPAN CONTINUOUS SLAB OPTIMISATION
# -----------------------------------------------------------------------------
#  Optimises composite steel-concrete, RC and timber slab cross-sections for
#  minimum GWP across a range of spans.
#
#  RC solid and composite slabs are analysed as two-span continuous beams.
#  Timber slabs and RC ribbed stay simply-supported (continuous not implemented).
#
# =============================================================================

import struct_analysis  # structural analysis
import plot_datasets    # standardized plotting of results
import matplotlib.pyplot as plt
import os

# #############################################################################
# #                                                                         # #
# #                            USER INPUTS                                  # #
# #            Everything you may want to change is in this block.          # #
# #                                                                         # #
# #############################################################################

# ── ANALYSIS RANGE ─────────────────────────────────────────────────────────
n_spans  = 2                           # spans for RC / composite (timber stays at 1)
lengths  = [3, 4, 5, 6, 7, 8, 9, 10]   # spans to evaluate (x-axis of plot) [m]
idx_vrc  = 1                           # index in `lengths` of the verified span
max_iter = 100                         # iterations per optimisation (higher = better)

# ── DESIGN CASE ────────────────────────────────────────────────────────────
# use_case : occupancy type — sets imposed/permanent loads and floor build-up
#   "residential"  → EN 1991-1-1 Cat. A   qk = 2.0 kN/m²
#   "office"       → EN 1991-1-1 Cat. B   qk = 3.0 kN/m²  + 1.0 kN/m² partitions
#   "retail"       → EN 1991-1-1 Cat. D1  qk = 4.0 kN/m²
#   "industrial"   → EN 1991-1-1 Cat. E1  qk = 7.5 kN/m²
use_case = "office"

criteria = ["ENV"]   # design envelope to satisfy (e.g. "ULS", "SLS", "FIRE")
optima   = ["GWP"]    # quantity to minimise

# ── FIRE REQUIREMENT ───────────────────────────────────────────────────────
# fire      : required fire resistance period — 'R30'/'R60'/'R90'/'R120'.
#             Use 'R0' to skip the fire check entirely.
# cover_top : nominal cover, top slab surface to mesh bar face (mm).
#             EN 1992-1-1: c_nom = c_min,dur + 10 mm  (XC1 interior → 20 mm).
#             Fire axis distance = cover_top + d_s/2; EN 1994-1-2 minimums:
#             R30 → 10 mm, R60 → 20 mm, R90 → 30 mm, R120 → 40 mm (axis dist).
fire      = 'R60'
cover_top = 20.0

# ── COMPOSITE SLAB OPTIONS ─────────────────────────────────────────────────
database_name = "database_270426.db"

# manufacturer : restrict composite deck products by manufacturer.
#   None                          → every deck in the database (all manufacturers)
#   ["Tata Steel"]                → only Tata Steel ComFlor decks
#   ["Kingspan"]                  → only Kingspan Multideck decks
#   ["ArcelorMittal"]             → only ArcelorMittal Cofraplus / Cofrastra decks
#   ["Tata Steel", "Kingspan"]    → Tata Steel and Kingspan (exclude ArcelorMittal)
manufacturer = ["Tata Steel", "Kingspan"]

# concrete_grades : grades to iterate through for the composite slab.
#   ["C25/30"]                     → C25/30 only (fastest; default)
#   ["C25/30", "C30/37"]           → both grades
#   ["C20/25", "C25/30", "C30/37"] → every grade in database
#   None                           → falls back to ["C25/30"]
concrete_grades = ["C25/30", "C30/37"]

# ── RUN CONTROL ────────────────────────────────────────────────────────────
# REPLOT_ONLY : skip the (slow) optimisation and regenerate plots from the
#               saved CSV. Set True to re-format plots without re-running.
REPLOT_ONLY = False

# #############################################################################
# #         Below this point: definitions and logic — no need to edit.       #
# #############################################################################

plt.close('all')  # clear any leftover figures from previous runs

# ── PER-OCCUPANCY DEFINITIONS ─────────────────────────────────────────────
# floor_layers   : finish layers on top of the slab (used for all slab types)
# g2k_services   : additional permanent load beyond floor structure (services/fittings)
#                  used for all slab types (ceiling_services for composite, g2k for RC/timber)
# partition_load : movable partition allowance (kN/m²), treated as variable action for all types

_use_cases = {
    "residential": {
        # EN 1991-1-1 Cat. A – dwellings / apartments
        "qk":             2.0,    # kN/m²  EN 1991-1-1 Table 6.2 Cat. A upper bound
        "psi_fi":         0.5,    # EN 1990 Table A1.3 — Cat. A
        "partition_load": 0.0,
        "g2k_services":   0.75,   # kN/m²  Einbauten (services/fittings beyond floor layers)
        "floor_layers": [
            ["'Parkett 2-Schicht werkversiegelt, 11 mm'", False, False],
            ["'Unterlagsboden Zement, 85 mm'",            False, False],
            ["'Glaswolle'",                               0.03,  False],
        ],
    },
    "office": {
        # EN 1991-1-1 Cat. B – office areas
        # Raised platform floor Plus replaces screed + parquet; services routed underneath
        "qk":             2.5,    # kN/m²  UK NA Cat. B recommended value
        "psi_fi":         0.5,    # EN 1990 Table A1.3 — Cat. B
        "partition_load": 0.5,    # kN/m²  lightweight demountable partitions (EN 1991-1-1 Cl. 6.3.1.2 ≤1.0 kN/m wall)
        "g2k_services":   0.75,   # kN/m²  services
        "floor_layers": [
            ["'Spanplatte'", False, False],
            ["'Gipsfaserplatte'", 0.019, False],
            ["'Steinwolle'", 0.03,  False],
        ],
    },
    "retail": {
        # EN 1991-1-1 Cat. D1 – general retail
        "qk":             4.0,
        "psi_fi":         0.7,    # EN 1990 Table A1.3 — Cat. D
        "partition_load": 0.0,
        "g2k_services":   0.5,    # kN/m²  services/fittings allowance
        "floor_layers": [
            ["'Parkett 2-Schicht werkversiegelt, 11 mm'", False, False],
            ["'Unterlagsboden Zement, 85 mm'", False, False],
            ["'Glaswolle'", 0.03, False],
        ],
    },
    "industrial": {
        # EN 1991-1-1 Cat. E1 – storage / industrial
        "qk":             7.5,
        "psi_fi":         0.9,    # EN 1990 Table A1.3 — Cat. E
        "partition_load": 0.0,
        "g2k_services":   0.5,    # kN/m²  services/equipment allowance
        "floor_layers": [
            ["'Unterlagsboden Zement, 85 mm'", False, False],
            ["'Steinwolle'", 0.03,  False]],
    },
}

if use_case not in _use_cases:
    raise ValueError(f"Unknown use_case '{use_case}'. "
                     f"Choose from: {list(_use_cases.keys())}")

_cfg = _use_cases[use_case]

# ── FLOOR BUILD-UPS ───────────────────────────────────────────────────────

# composite slab and RC solid: occupancy finish layers
bodenaufbau_comp = struct_analysis.FloorStruc(_cfg["floor_layers"], database_name)
bodenaufbau_rc   = struct_analysis.FloorStruc(_cfg["floor_layers"], database_name)

# RC ribbed (slim profile): finish + thin gravel layer
bodenaufbau_rc_rib = struct_analysis.FloorStruc(
    _cfg["floor_layers"] + [["'Kies gebrochen'", 0.06, False]], database_name)

# timber solid: finish + gravel for impact-sound insulation
bodenaufbau_wd_solid = struct_analysis.FloorStruc(
    _cfg["floor_layers"] + [["'Kies gebrochen'", 0.12, False]], database_name)

# timber ribbed: finish + fire-protection assembly (REI60, Lignum 4.1 Table 433-2 Col G):
# Gipsfaserplatte (30mm) + Glaswolle (30mm) + Kies (120mm) + Steinwolle (180mm within rib)
h_ins = 0.18
bodenaufbau_wd_rib = struct_analysis.FloorStruc(
    _cfg["floor_layers"] + [
        ["'Gipsfaserplatte'", 0.03, False],
        ["'Glaswolle'",       0.03, False],
        ["'Kies gebrochen'",  0.12, False],
        ["'Steinwolle'",      h_ins, False],
    ], database_name)
bodenaufbau_wd_rib.h -= h_ins   # steinwolle sits within rib depth, not additive

# ── LOAD VARIABLES ────────────────────────────────────────────────────────

# finishes_load for composite slab derived from floor layer weights in database
# FloorStruc.gk_area is already in N/m² (weight column in N/m³ × layer height in m)
finishes_load = bodenaufbau_comp.gk_area             # N/m²

# g2k / qk in N/m² (used by all slab types)
# Partition load (EN 1991-1-1 Cl. 6.3.1.2) is treated as a variable action for
# all slab types so the comparison is on equal terms with the composite slab.
qk  = (_cfg["qk"] + _cfg["partition_load"]) * 1e3   # N/m²  imposed + partition
g2k = _cfg["g2k_services"] * 1e3                     # N/m²  additional permanent (services/fittings)

# ── REQUIREMENTS ──────────────────────────────────────────────────────────
req = struct_analysis.Requirements(fire=fire, cover_top=cover_top)

# composite slab loads — all in N/m² (SI), matching g2k / qk used by RC and timber
comp_slab_loads = {
    "imposed_load":      _cfg["qk"] * 1e3,             # N/m²  imposed (excl. partitions)
    "finishes_load":     finishes_load,                 # N/m²  floor layers from database
    "ceiling_services":  g2k,                           # N/m²  services/fittings (= g2k)
    "construction_load": 750.0,                         # N/m²  same for all occupancies
    "partition_load":    _cfg["partition_load"] * 1e3,  # N/m²
    # fire design — period and cover from requirements; psi_fi from occupancy
    "R_fi":      req.t_fire,
    "cover_top": req.cover_top,
    "psi_fi":    _cfg["psi_fi"],
}

print(f"Use case      : {use_case}  [{n_spans}-span continuous for RC/composite]")
print(f"  qk (imposed)  = {_cfg['qk'] * 1e3:.0f} N/m²")
print(f"  partitions    = {_cfg['partition_load'] * 1e3:.0f} N/m²  (all slab types, EN 1991-1-1 Cl. 6.3.1.2)")
print(f"  qk (total)    = {qk:.0f} N/m²  (all slab types)")
print(f"  finishes      = {finishes_load:.0f} N/m²  (from database layer weights, composite)")
print(f"  g2k_services  = {g2k:.0f} N/m²  (all slab types: ceiling/services permanent load)")
print(f"  Floor layers  : {[l[0] for l in _cfg['floor_layers']]}")
print(f"  Fire check    : R{req.t_fire}  cover_top={req.cover_top:.0f}mm  psi_fi={_cfg['psi_fi']}")
print()

# ── OPTIMISATION ───────────────────────────────────────────────────────────
results_csv = f"results_{use_case}_2span_{'_'.join(criteria)}.csv"

best_vrc = None   # best composite slab at verification span (live object)

if not REPLOT_ONLY:
    # ── delete old CSV so this run starts with a clean file ───────────────
    if os.path.exists(results_csv):
        os.remove(results_csv)

    # ── TIMBER SOLID (rectangular) — simply-supported ─────────────────────
    mat_names = ["'Glue_laminated_timber'", "'3- and 5-ply wood'", "'Solid_structural_timber'"]
    plot_datasets.run_dataset(
        lengths, database_name, criteria, optima, bodenaufbau_wd_solid, req,
        "wd_rec", mat_names, g2k, qk, max_iter, idx_vrc,
        n_spans=1,
        results_csv=results_csv, use_case=use_case)

    # ── TIMBER RIBBED — simply-supported ──────────────────────────────────
    mat_names = ["'Glue_laminated_timber'", "'Solid_structural_timber'"]
    plot_datasets.run_dataset(
        lengths, database_name, criteria, optima, bodenaufbau_wd_rib, req,
        "wd_rib", mat_names, g2k, qk, max_iter, idx_vrc,
        n_spans=1,
        results_csv=results_csv, use_case=use_case)

    # ── RC SOLID (rectangular) — two-span continuous ───────────────────────
    mat_names = ["'ready_mixed_concrete'"]
    plot_datasets.run_dataset(
        lengths, database_name, criteria, optima, bodenaufbau_rc, req,
        "rc_rec", mat_names, g2k, qk, max_iter, idx_vrc,
        n_spans=n_spans,
        results_csv=results_csv, use_case=use_case)

    # ── RC RIBBED — simply-supported (continuous not yet implemented) ──────
    mat_names = ["'ready_mixed_concrete'"]
    plot_datasets.run_dataset(
        lengths, database_name, criteria, optima, bodenaufbau_rc_rib, req,
        "rc_rib", mat_names, g2k, qk, max_iter, idx_vrc,
        n_spans=1,
        results_csv=results_csv, use_case=use_case)

    # ── COMPOSITE SLAB — two-span continuous ──────────────────────────────
    mat_names = ["'ready_mixed_concrete'"]
    _, best_vrc = plot_datasets.run_dataset(
        lengths, database_name, criteria, optima, bodenaufbau_comp, req,
        "comp_slab", mat_names, g2k, qk, max_iter, idx_vrc,
        deck=None, comp_slab_loads=comp_slab_loads,
        manufacturer=manufacturer, concrete_grades=concrete_grades,
        n_spans=n_spans,
        results_csv=results_csv, use_case=use_case)

    print(f"\nResults saved to '{results_csv}'")

# ── PLOT FROM CSV (always runs) ────────────────────────────────────────────
data_max = plot_datasets.plot_from_csv(results_csv, lengths)

# cross-section plot of best composite slab (only when we have live objects)
# if best_vrc is not None:
#     plot_datasets.plot_section(best_vrc.section)
#     plt.title(
#         f"Best composite slab at L = {lengths[idx_vrc]:.1f} m  [{n_spans}-span]\n"
#         f"Deck: {best_vrc.section.deck.mech_prop}  |  "
#         f"h = {best_vrc.section.h_mm:.0f} mm  |  "
#         f"GWP = {best_vrc.section.co2:.1f} kg-CO₂-eq/m²",
#         fontsize=9)

# ── PLOT FORMATTING ───────────────────────────────────────────────────────
_title = (f"Two-span continuous — {use_case.capitalize()} ({'/'.join(criteria)})  "
          f"|  qk={qk:.0f} N/m², finishes={finishes_load:.0f} N/m²\n"
          f"RC solid and composite: two-span continuous  |  RC ribbed, timber solid, timber ribbed: simply-supported")

# idx→(fig, subplot): 0=h_struct,1=h_tot, 2=gwp_struct,3=gwp_tot, 4=mass_struct,5=mass_tot
_fig_info = {
    1: (["h$_{struct}$", "[m]"],              ["h$_{tot}$", "[m]"]),
    2: (["GWP$_{struct}$", "[kg-CO$_2$-eq/m²]"], ["GWP$_{tot}$", "[kg-CO$_2$-eq/m²]"]),
    3: (["Mass$_{struct}$", "[kg/m²]"],       ["Mass$_{tot}$", "[kg/m²]"]),
}
_plot_dir = os.path.join("plots", use_case, '_'.join(criteria))
os.makedirs(_plot_dir, exist_ok=True)
_fig_filenames = {
    1: os.path.join(_plot_dir, f"plot_{use_case}_2span_{'_'.join(criteria)}_height.png"),
    2: os.path.join(_plot_dir, f"plot_{use_case}_2span_{'_'.join(criteria)}_gwp.png"),
    3: os.path.join(_plot_dir, f"plot_{use_case}_2span_{'_'.join(criteria)}_mass.png"),
}

for fig_n, (info_left, info_right) in _fig_info.items():
    fig = plt.figure(fig_n)
    fig.set_size_inches(10, 4)
    for sub_pos, info in enumerate([info_left, info_right], start=1):
        idx = (fig_n - 1) * 2 + (sub_pos - 1)
        plt.subplot(1, 2, sub_pos)
        plt.xlabel('l [m]', fontsize=12)
        plt.ylabel(info[0] + " " + info[1], fontsize=12)
        pair_max = max(data_max[(fig_n - 1) * 2], data_max[(fig_n - 1) * 2 + 1])
        plt.axis((min(lengths), max(lengths), 0, pair_max * 1.1))
        plt.grid()
    plt.suptitle(_title, fontsize=9)
    plt.tight_layout(rect=[0, 0.10, 1, 0.93])

plot_datasets.add_dataset_legend(has_propped=True)

for fig_n, fname in _fig_filenames.items():
    plt.figure(fig_n)
    plt.savefig(fname, dpi=150, bbox_inches='tight')
    print(f"Plot saved to '{fname}'")
plt.show()
