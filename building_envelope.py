"""
Building envelope comparison — storey count analysis.

For each use case, computes the number of storeys achievable in a fixed building
height using the optimum (minimum GWP) slab configuration at each span.

Two comparisons are produced for each use case:
  1. Within slab type : how storey count changes with span for each type
  2. Across slab types: for the same span, which type yields the most storeys

Formula:
    floor-to-floor (ftf) = h_tot + ceiling_void + clear_height   [m]
    n_storeys = building_height / ftf                              [-]

Where:
    h_tot        = h_struct + floor build-up above slab (finishes, screed,
                   acoustic insulation) — taken directly from the CSV.
                   This captures the full floor zone consumed above the
                   structural element, which differs significantly between
                   types (composite ~67 mm, timber solid ~187 mm).
    ceiling_void = space below the slab for MEP routing and hung ceiling.
    clear_height = usable floor-to-ceiling height.

n_storeys is reported as a continuous (non-integer) value so fractional
differences between slab types remain visible in the comparison. In practice,
the achievable storey count is the floor of n_storeys.

Multi-span composite (2-span, 3-span) results are included as separate
sub-rows under comp_slab where those CSV files are available.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

base = r"C:\Users\ellie\OneDrive - Imperial College London\Year 4\Final Project\Implementation\Slimfloor Case"

# ── matplotlib style — matches plot_datasets.py ───────────────────────────────
plt.rcParams.update({
    'font.family':    'Times New Roman',
    'font.size':      11,
    'axes.grid':      True,
    'grid.color':     '#cccccc',
    'grid.linewidth': 0.6,
    'axes.axisbelow': True,
})

# ==============================================================================
# BUILDING PARAMETERS PER USE CASE
# ==============================================================================
#   bldg_h       : total building height [m]
#   clear_h      : floor-to-ceiling clear height [m]
#   ceiling_void : space below slab for MEP routing and hung ceiling [m]
#                  (floor finishes above slab are already in h_tot from CSV)

USE_CASES = {
    "office": {
        "bldg_h":       40.0,
        "clear_h":       2.7,
        "ceiling_void":  0.40,
        "csvs": [
            ("SS",     "results_office_ENV.csv"),
            ("2-span", "results_office_2span_ENV.csv"),
            ("3-span", "results_office_3span_ENV.csv"),
        ],
    },
    "residential": {
        "bldg_h":       30.0,
        "clear_h":       2.5,
        "ceiling_void":  0.20,   # increased from 0.15 m — minimum for service pipes
        "csvs": [
            ("SS",     "results_residential_ENV.csv"),
            ("2-span", "results_residential_2span_ENV.csv"),
            ("3-span", "results_residential_3span_ENV.csv"),
        ],
    },
    "retail": {
        "bldg_h":       20.0,
        "clear_h":       4.0,
        "ceiling_void":  0.35,
        "csvs": [
            ("SS",     "results_retail_ENV.csv"),
            ("2-span", "results_retail_2span_ENV.csv"),
            ("3-span", "results_retail_3span_ENV.csv"),
        ],
    },
    "industrial": {
        "bldg_h":       12.0,
        "clear_h":       6.0,
        "ceiling_void":  0.25,
        "csvs": [
            ("SS", "results_industrial_ENV.csv"),
        ],
    },
}

# ==============================================================================
# SLAB TYPE LABELS AND PLOT STYLE
# ==============================================================================
SLAB_LABELS = {
    "comp_slab": "Composite",
    "rc_rec":    "RC Solid",
    "rc_rib":    "RC Ribbed",
    "wd_rec":    "Timber Solid",
    "wd_rib":    "Timber Ribbed",
}

SLAB_ORDER = ["comp_slab", "rc_rec", "rc_rib", "wd_rec", "wd_rib"]

USE_CASE_TITLES = {
    "office":      "Office",
    "residential": "Residential",
    "retail":      "Retail",
    "industrial":  "Industrial",
}

# Plot style per slab type — colours match plot_datasets.py
SLAB_STYLE = {
    "comp_slab": {"color": "steelblue",   "ls": "-",   "marker": "o"},
    "rc_rec":    {"color": "green",       "ls": "-",   "marker": "s"},
    "rc_rib":    {"color": "limegreen",   "ls": "-",   "marker": "^"},
    "wd_rec":    {"color": "saddlebrown", "ls": "-",   "marker": "D"},
    "wd_rib":    {"color": "sandybrown",  "ls": "-",   "marker": "v"},
}

# Composite sub-span line styles
COMP_SPAN_STYLE = {
    "SS":     {"ls": "-",   "marker": "o"},
    "2-span": {"ls": "--",  "marker": "s"},
    "3-span": {"ls": "-.",  "marker": "^"},
}


# ==============================================================================
# HELPERS
# ==============================================================================
def _is_propped(val):
    return str(val).strip().lower() in ('true', '1')


def load_best(csv_path, span_label):
    """Return DataFrame of best (min GWP) config per (slab_type, span_m).

    For composite slabs: unpropped configurations are preferred; propped is
    used only when no unpropped feasible section exists at that span.
    This is consistent with the treatment in the plotting scripts.
    """
    if not os.path.exists(csv_path):
        return None
    df = pd.read_csv(csv_path)
    df = df[(df["optimum"] == "GWP") & (df["criterion"] == "ENV")].dropna(
        subset=["h_struct_m", "gwp_tot", "h_tot_m"]
    )

    rows = []
    for (slab_type, span_m), group in df.groupby(["slab_type", "span_m"]):
        if slab_type == 'comp_slab':
            unp = group[~group['propped'].apply(_is_propped)]
            pool = unp if not unp.empty else group
        else:
            pool = group
        rows.append(pool.loc[pool['gwp_tot'].idxmin()])
    best = pd.DataFrame(rows).sort_values(["slab_type", "span_m"]).copy()
    best["span_label"] = span_label
    mat = (best["material"]
           .fillna("")
           .str.replace(r"^error:.*$", "—", regex=True)
           .str.replace("'", "", regex=False))
    best["section"] = best["deck_name"].fillna(mat)
    return best


def n_storeys_calc(h_tot, clear_h, ceiling_void, bldg_h):
    ftf = h_tot + clear_h + ceiling_void
    return bldg_h / ftf, ftf


def _latex_escape(s):
    """Escape characters that are special in LaTeX."""
    return (str(s)
            .replace('&', r'\&')
            .replace('%', r'\%')
            .replace('_', r'\_')
            .replace('#', r'\#')
            .replace('~', r'\textasciitilde{}')
            .replace('^', r'\textasciicircum{}'))


def _build_latex(use_case, cfg, data, span_labels_comp):
    """Return a LaTeX string with two tables for one use case."""
    bldg_h       = cfg["bldg_h"]
    clear_h      = cfg["clear_h"]
    ceiling_void = cfg["ceiling_void"]
    title        = USE_CASE_TITLES.get(use_case, use_case.capitalize())
    all_spans    = sorted(data["span_m"].unique())
    comp_data    = data[data["slab_type"] == "comp_slab"]
    other_data   = data[data["slab_type"] != "comp_slab"]

    lines = []
    lines.append(f"\n% {'='*60}")
    lines.append(f"% {title.upper()}")
    lines.append(f"% {'='*60}\n")

    # ── Table 1: storey count matrix ──────────────────────────────────────────
    n_cols  = len(all_spans)
    col_fmt = "l" + "r" * n_cols
    span_hdr = " & ".join(f"{int(s)}\\,m" for s in all_spans)

    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\small")
    caption = (f"Achievable storey count vs span --- {title} "
               f"(building height = {bldg_h:.0f}\\,m, "
               f"clear height = {clear_h:.1f}\\,m, "
               f"ceiling void = {ceiling_void:.2f}\\,m). "
               r"Values are continuous; $\lfloor n \rfloor$ gives the achievable integer count.")
    lines.append(f"\\caption{{{caption}}}")
    lines.append(f"\\label{{tab:envelope_{use_case}_storey}}")
    lines.append(f"\\begin{{tabular}}{{{col_fmt}}}")
    lines.append(r"\toprule")
    lines.append(f"Slab type & {span_hdr} \\\\")
    lines.append(r"\midrule")

    # Composite rows
    for slbl in span_labels_comp:
        sub = comp_data[comp_data["span_label"] == slbl]
        if sub.empty:
            continue
        row_vals = []
        for sp in all_spans:
            match = sub[sub["span_m"] == sp]
            row_vals.append(f"{match.iloc[0]['n']:.2f}" if len(match) else "---")
        lines.append(f"Composite ({_latex_escape(slbl)}) & " + " & ".join(row_vals) + r" \\")

    lines.append(r"\midrule")

    # Other slab types
    for stype in SLAB_ORDER[1:]:
        sub = other_data[other_data["slab_type"] == stype]
        if sub.empty:
            continue
        row_vals = []
        for sp in all_spans:
            match = sub[sub["span_m"] == sp]
            row_vals.append(f"{match.iloc[0]['n']:.2f}" if len(match) else "---")
        lbl = _latex_escape(SLAB_LABELS.get(stype, stype))
        lines.append(f"{lbl} & " + " & ".join(row_vals) + r" \\")

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")
    lines.append("")

    # ── Table 2: section detail at each span ──────────────────────────────────
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\small")
    caption2 = (f"Optimal section details and floor-to-floor height --- {title}. "
                r"Ranked by storey count (descending) at each span. "
                r"$h$ = structural depth, ftf = floor-to-floor height.")
    lines.append(f"\\caption{{{caption2}}}")
    lines.append(f"\\label{{tab:envelope_{use_case}_detail}}")
    lines.append(r"\begin{tabular}{llrrr}")
    lines.append(r"\toprule")
    lines.append(r"Slab type & Section & $h$ (mm) & ftf (m) & Storeys \\")
    lines.append(r"\midrule")

    for sp in all_spans:
        sp_data   = data[data["span_m"] == sp]
        comp_sp   = sp_data[sp_data["slab_type"] == "comp_slab"]
        others_sp = sp_data[sp_data["slab_type"] != "comp_slab"]

        rows = []
        for slbl in span_labels_comp:
            sub = comp_sp[comp_sp["span_label"] == slbl]
            if not sub.empty:
                r = sub.iloc[0]
                rows.append((f"Composite ({slbl})",
                             str(r["deck_name"]),
                             r["h_struct_m"] * 1000, r["ftf"], r["n"]))
        for stype in SLAB_ORDER[1:]:
            sub = others_sp[others_sp["slab_type"] == stype]
            if not sub.empty:
                r = sub.iloc[0]
                rows.append((SLAB_LABELS[stype], str(r["section"]),
                             r["h_struct_m"] * 1000, r["ftf"], r["n"]))

        if not rows:
            continue
        rows.sort(key=lambda x: x[4], reverse=True)

        lines.append(f"\\multicolumn{{5}}{{l}}{{\\textit{{Span = {int(sp)}\\,m}}}} \\\\")
        for lbl, sect, hmm, ftf, n in rows:
            sect_clean = _latex_escape(sect) if sect and sect != 'nan' else '---'
            lines.append(f"{_latex_escape(lbl)} & {sect_clean} & "
                         f"{hmm:.0f} & {ftf:.3f} & {n:.2f} \\\\")
        lines.append(r"\midrule")

    # Remove last \midrule before \bottomrule
    if lines and lines[-1] == r"\midrule":
        lines[-1] = r"\bottomrule"
    else:
        lines.append(r"\bottomrule")

    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")
    lines.append("")

    return "\n".join(lines)


# ==============================================================================
# MAIN LOOP — PRINT TABLES + COLLECT DATA FOR PLOTTING
# ==============================================================================
plot_data   = {}   # {use_case: DataFrame}
latex_parts = []   # accumulate LaTeX for all use cases

for use_case, cfg in USE_CASES.items():

    bldg_h       = cfg["bldg_h"]
    clear_h      = cfg["clear_h"]
    ceiling_void = cfg["ceiling_void"]
    fixed        = clear_h + ceiling_void

    frames = []
    for span_label, fn in cfg["csvs"]:
        path = os.path.join(base, fn)
        df   = load_best(path, span_label)
        if df is not None:
            frames.append(df)

    if not frames:
        print(f"\n{'='*80}")
        print(f"USE CASE: {use_case.upper()}  (no CSV results found — skipping)")
        print(f"{'='*80}")
        continue

    data = pd.concat(frames, ignore_index=True)
    data["ftf"] = data["h_tot_m"] + fixed
    data["n"]   = bldg_h / data["ftf"]
    data["slab_label"] = data["slab_type"].map(SLAB_LABELS)

    plot_data[use_case] = data
    all_spans = sorted(data["span_m"].unique())
    comp_data  = data[data["slab_type"] == "comp_slab"]
    other_data = data[data["slab_type"] != "comp_slab"]
    span_labels_comp = sorted(comp_data["span_label"].unique(),
                               key=lambda x: {"SS": 0, "2-span": 1, "3-span": 2}.get(x, 9))
    latex_parts.append(_build_latex(use_case, cfg, data, span_labels_comp))

    print(f"\n{'='*80}")
    print(f"USE CASE: {use_case.upper()}")
    print(f"  Building height = {bldg_h} m  |  Clear height = {clear_h} m  "
          f"|  Ceiling void = {ceiling_void} m")
    print(f"  Fixed floor-to-floor component = {fixed:.2f} m  "
          f"(excl. floor slab + build-up above)")
    print(f"  n_storeys is continuous (non-integer) — floor() gives achievable count")
    print(f"{'='*80}")

    # ── Comparison 1: within slab type ────────────────────────────────────────
    print("\n--- 1. Within slab type: storey count vs span ---\n")
    hdr_spans = "".join(f"{s:>7.0f}m" for s in all_spans)
    print(f"  {'Slab type / condition':<28}{hdr_spans}")
    print(f"  {'-'*28}{'-'*8*len(all_spans)}")

    for slbl in span_labels_comp:
        sub = comp_data[comp_data["span_label"] == slbl]
        row_str = f"  {'Composite (' + slbl + ')':<28}"
        for sp in all_spans:
            match = sub[sub["span_m"] == sp]
            row_str += f"{match.iloc[0]['n']:>7.2f} " if len(match) else f"{'--':>7} "
        print(row_str)

    for stype in SLAB_ORDER[1:]:
        sub = other_data[other_data["slab_type"] == stype]
        if sub.empty:
            continue
        row_str = f"  {SLAB_LABELS.get(stype, stype):<28}"
        for sp in all_spans:
            match = sub[sub["span_m"] == sp]
            row_str += f"{match.iloc[0]['n']:>7.2f} " if len(match) else f"{'--':>7} "
        print(row_str)

    # ── Comparison 2: across slab types at each span ───────────────────────────
    print("\n--- 2. Across slab types: storey count at each span ---\n")

    for sp in all_spans:
        sp_data   = data[data["span_m"] == sp]
        comp_sp   = sp_data[sp_data["slab_type"] == "comp_slab"]
        others_sp = sp_data[sp_data["slab_type"] != "comp_slab"]

        rows = []
        for slbl in span_labels_comp:
            sub = comp_sp[comp_sp["span_label"] == slbl]
            if not sub.empty:
                r = sub.iloc[0]
                rows.append(("Composite (" + slbl + ")",
                             r["deck_name"], r["h_struct_m"] * 1000, r["ftf"], r["n"]))
        for stype in SLAB_ORDER[1:]:
            sub = others_sp[others_sp["slab_type"] == stype]
            if not sub.empty:
                r = sub.iloc[0]
                rows.append((SLAB_LABELS[stype], r["section"],
                             r["h_struct_m"] * 1000, r["ftf"], r["n"]))

        if not rows:
            continue
        rows.sort(key=lambda x: x[4], reverse=True)

        print(f"  Span = {sp:.0f} m")
        print(f"    {'Type':<26} {'Section':<28} {'h (mm)':>7} "
              f"{'ftf (m)':>8} {'Storeys':>8}")
        print(f"    {'-'*78}")
        for lbl, sect, hmm, ftf, n in rows:
            print(f"    {lbl:<26} {str(sect):<28} {hmm:>7.0f} {ftf:>8.3f} {n:>8.2f}")
        print()

    # ── Summary ────────────────────────────────────────────────────────────────
    print(f"--- Summary ---")
    print(f"  Max storeys: {data['n'].max():.2f}  "
          f"({data.loc[data['n'].idxmax(), 'slab_label']}  "
          f"h={data.loc[data['n'].idxmax(), 'h_struct_m']*1000:.0f} mm  "
          f"span={data.loc[data['n'].idxmax(), 'span_m']:.0f} m)")
    print(f"  Min storeys: {data['n'].min():.2f}  "
          f"({data.loc[data['n'].idxmin(), 'slab_label']}  "
          f"h={data.loc[data['n'].idxmin(), 'h_struct_m']*1000:.0f} mm  "
          f"span={data.loc[data['n'].idxmin(), 'span_m']:.0f} m)")
    print(f"  Range: {data['n'].max() - data['n'].min():.2f} storeys across "
          f"all types and spans\n")


# ==============================================================================
# WRITE LaTeX FILE
# ==============================================================================
latex_header = r"""% building_envelope_results.tex
% Auto-generated by building_envelope.py — do not edit by hand.
% Requires packages: booktabs, caption
% Usage: \input{building_envelope_results.tex}
%
% n_storeys is continuous (non-integer).
% The achievable integer storey count is floor(n_storeys).

"""

latex_path = os.path.join(base, 'building_envelope_results.tex')
with open(latex_path, 'w', encoding='utf-8') as f:
    f.write(latex_header)
    f.write("\n".join(latex_parts))
print(f'Saved {latex_path}')


# ==============================================================================
# PLOTS
# ==============================================================================

def _plot_envelope(ax, data, cfg, use_case):
    """Plot n_storeys vs span for all slab types on ax."""
    bldg_h       = cfg["bldg_h"]
    clear_h      = cfg["clear_h"]
    ceiling_void = cfg["ceiling_void"]

    comp_data  = data[data["slab_type"] == "comp_slab"]
    other_data = data[data["slab_type"] != "comp_slab"]

    span_labels_comp = sorted(comp_data["span_label"].unique(),
                               key=lambda x: {"SS": 0, "2-span": 1, "3-span": 2}.get(x, 9))

    # Composite — one line per span condition
    base_style = SLAB_STYLE["comp_slab"]
    for slbl in span_labels_comp:
        sub = comp_data[comp_data["span_label"] == slbl].sort_values("span_m")
        if sub.empty:
            continue
        s = COMP_SPAN_STYLE[slbl]
        ax.plot(sub["span_m"], sub["n"],
                color=base_style["color"], linestyle=s["ls"], marker=s["marker"],
                markersize=5, linewidth=1.5,
                label=f"Composite ({slbl})")

    # Other slab types
    for stype in SLAB_ORDER[1:]:
        sub = other_data[other_data["slab_type"] == stype].sort_values("span_m")
        if sub.empty:
            continue
        st = SLAB_STYLE[stype]
        ax.plot(sub["span_m"], sub["n"],
                color=st["color"], linestyle=st["ls"], marker=st["marker"],
                markersize=5, linewidth=1.5,
                label=SLAB_LABELS[stype])

    # Integer storey reference lines
    n_min = data["n"].min()
    n_max = data["n"].max()
    for ni in range(int(np.floor(n_min)), int(np.ceil(n_max)) + 1):
        ax.axhline(ni, color='grey', linewidth=0.5, linestyle=':', zorder=0)

    ax.set_title(USE_CASE_TITLES.get(use_case, use_case))
    ax.set_xlabel(r'Span $l$ [m]')
    ax.set_ylabel(r'Achievable storeys (continuous)')
    all_spans = sorted(data["span_m"].unique())
    ax.set_xticks(all_spans)


# ── Figure A: one subplot per use case (2 × 2 grid) ─────────────────────────
cases = [uc for uc in USE_CASES if uc in plot_data]
n_cases = len(cases)
ncols = 2
nrows = (n_cases + 1) // 2

figA, axes = plt.subplots(nrows, ncols, figsize=(14, 5 * nrows),
                           sharex=False, sharey=False)
axes = np.array(axes).flatten()

for idx, use_case in enumerate(cases):
    ax  = axes[idx]
    cfg = USE_CASES[use_case]
    _plot_envelope(ax, plot_data[use_case], cfg, use_case)

# Hide unused subplots
for idx in range(n_cases, len(axes)):
    axes[idx].set_visible(False)

# Shared legend below the figure
handles, labels = axes[0].get_legend_handles_labels()
figA.legend(handles, labels, loc='lower center',
            ncol=min(len(labels), 4), fontsize=9, frameon=True,
            bbox_to_anchor=(0.5, 0.01))
figA.suptitle('Storey count vs span — minimum GWP slab per type', y=1.01)
figA.tight_layout(rect=[0, 0.07, 1, 1])
outA = os.path.join(base, 'fig_envelope_storey_count.png')
figA.savefig(outA, dpi=200, bbox_inches='tight')
print(f'Saved {outA}')


# ── Figure B: composite vs RC solid only, all use cases overlaid ─────────────
# Shows clearly how composite compares to RC solid across occupancies
figB, axB = plt.subplots(figsize=(10, 6))
axB.set_title('Composite (SS) vs RC Solid — storey count across occupancies')

OCC_COLORS = {
    "office":      "steelblue",
    "residential": "mediumpurple",
    "retail":      "green",
    "industrial":  "saddlebrown",
}

for use_case in cases:
    data = plot_data[use_case]
    color = OCC_COLORS.get(use_case, 'grey')
    label = USE_CASE_TITLES.get(use_case, use_case)

    # Composite SS
    comp_ss = (data[(data["slab_type"] == "comp_slab") & (data["span_label"] == "SS")]
               .sort_values("span_m"))
    if not comp_ss.empty:
        axB.plot(comp_ss["span_m"], comp_ss["n"],
                 color=color, linestyle='-', marker='o', markersize=5,
                 linewidth=1.5, label=f'{label} — composite')

    # RC Solid
    rc = data[data["slab_type"] == "rc_rec"].sort_values("span_m")
    if not rc.empty:
        axB.plot(rc["span_m"], rc["n"],
                 color=color, linestyle='--', marker='s', markersize=5,
                 linewidth=1.5, label=f'{label} — RC solid')

all_spans_B = sorted(set(
    sp for d in plot_data.values() for sp in d["span_m"].unique()
))
axB.set_xlabel(r'Span $l$ [m]')
axB.set_ylabel(r'Achievable storeys (continuous)')
axB.set_xticks(all_spans_B)
axB.legend(ncol=2, fontsize=9, frameon=True, loc='upper right')
figB.tight_layout()
outB = os.path.join(base, 'fig_envelope_comp_vs_rc.png')
figB.savefig(outB, dpi=200, bbox_inches='tight')
print(f'Saved {outB}')


# ── Figure C: office only — all composite span configs vs other types ─────────
if "office" in plot_data:
    figC, axC = plt.subplots(figsize=(10, 6))
    axC.set_title('Office — storey count: composite span configurations vs all slab types')
    _plot_envelope(axC, plot_data["office"], USE_CASES["office"], "office")
    axC.legend(ncol=2, fontsize=9, frameon=True, loc='upper right')
    figC.tight_layout()
    outC = os.path.join(base, 'fig_envelope_office_detail.png')
    figC.savefig(outC, dpi=200, bbox_inches='tight')
    print(f'Saved {outC}')

plt.show()
