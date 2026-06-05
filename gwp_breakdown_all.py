"""
gwp_breakdown_all.py  —  GWP breakdown for ALL slab types.

Reads a span-optimisation CSV (e.g. results_office_2span_ENV.csv) and produces
two figures, each with one subplot per slab type (2 × 3 grid):

  Fig 1  —  Structural GWP breakdown  [kg-CO₂-eq/m²]
             Composite bars: split into Concrete / Deck / Mesh / Trough bars
             RC & timber bars: structural GWP as a single block

  Fig 2  —  Total GWP  =  structural  +  floor build-up  [kg-CO₂-eq/m²]
             Same structural stacking as Fig 1, with floor layers added on top.
             Floor layers stacked and coloured by material name so the same
             material (e.g. "Raised Access Floor Plus") shares a colour across
             all slab types.

Both figures share a common y-axis range so slab types are directly comparable.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import numpy as np
import os

# ── CONFIGURATION ──────────────────────────────────────────────────────────────

# CSV file to analyse — edit path as needed
CSV_FILE = "results_industrial_2span_ENV.csv"

# Spans to include in the plot (None = all found in CSV)
SPANS_TO_SHOW = None

# Save PNG figures alongside the CSV?
SAVE_FIGS = True

# Figure DPI
DPI = 150

# ── SLAB TYPE DISPLAY ──────────────────────────────────────────────────────────
SLAB_CONFIG = [
    ("comp_slab", "Composite slab",  "steelblue"),
    ("rc_rec",    "RC solid",        "green"),
    ("rc_rib",    "RC ribbed",       "limegreen"),
    ("wd_rec",    "Timber solid",    "saddlebrown"),
    ("wd_rib",    "Timber ribbed",   "sandybrown"),
]

# ── COMPOSITE COMPONENT COLOURS ────────────────────────────────────────────────
COMP_COMPONENTS = ["Concrete", "Steel deck", "Mesh", "Trough bars"]
COMP_COLS       = ["gwp_concrete", "gwp_deck", "gwp_mesh", "gwp_trough_bars"]
COMP_COLOURS    = ["#4C72B0", "#DD8452", "#55A868", "#C44E52"]

# ── FLOOR LAYER COLOURS ────────────────────────────────────────────────────────
FLOOR_PALETTE = [
    "#aec7e8", "#ffbb78", "#c7e9c0", "#ff9896",
    "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d",
    "#9edae5", "#e7ba52",
]

# ── HELPERS ────────────────────────────────────────────────────────────────────

def _best_rows(df: pd.DataFrame) -> dict:
    """
    For each slab type and span, return the row with the lowest gwp_struct.
    Returns {slab_type: {span: pd.Series or None}}.
    """
    result = {}
    for slab_type, _, _ in SLAB_CONFIG:
        sub = df[df["slab_type"] == slab_type].copy()
        result[slab_type] = {}
        # coerce gwp_struct to numeric — drops error rows silently
        sub["_gwp"] = pd.to_numeric(sub["gwp_struct"], errors="coerce")
        sub = sub.dropna(subset=["_gwp"])
        for span in sorted(sub["span_m"].unique()):
            grp = sub[sub["span_m"] == span]
            if grp.empty:
                result[slab_type][span] = None
            else:
                result[slab_type][span] = grp.loc[grp["_gwp"].idxmin()]
    return result


def _floor_layers(row: pd.Series) -> list:
    """
    Return a list of (layer_name, gwp_value) for all non-zero floor layers.
    Layer names are cleaned (single quotes, underscores removed).
    """
    layers = []
    for i in range(8):
        name_col = f"floor_L{i}_name"
        gwp_col  = f"gwp_floor_L{i}"
        if name_col not in row.index or gwp_col not in row.index:
            break
        gwp = pd.to_numeric(row[gwp_col], errors="coerce")
        if pd.isna(gwp) or gwp <= 0:
            continue
        name = str(row[name_col]).strip().strip("'")
        layers.append((name, float(gwp)))
    return layers


def _build_floor_color_map(best: dict) -> dict:
    """Assign a consistent colour to every unique floor layer name."""
    all_names = set()
    for spans_dict in best.values():
        for row in spans_dict.values():
            if row is None:
                continue
            for name, _ in _floor_layers(row):
                all_names.add(name)
    return {name: FLOOR_PALETTE[i % len(FLOOR_PALETTE)]
            for i, name in enumerate(sorted(all_names))}


def _draw_subplot(ax, slab_type, slab_label, slab_colour,
                  best_rows_for_type: dict, spans: list,
                  mode: str, floor_cmap: dict) -> float:
    """
    Draw one subplot for a given slab type.
    mode : 'struct' | 'total'
    Returns max y value (for shared axis scaling).
    """
    x     = np.arange(len(spans))
    width = 0.55
    max_y = 0.0

    for xi, span in enumerate(spans):
        row = best_rows_for_type.get(span)
        if row is None:
            continue

        bottom = 0.0

        # ── structural part ───────────────────────────────────────────────
        if slab_type == "comp_slab":
            for col, colour, comp in zip(COMP_COLS, COMP_COLOURS, COMP_COMPONENTS):
                val = float(row[col]) if col in row.index and pd.notna(row.get(col)) else 0.0
                if val <= 0:
                    continue
                ax.bar(xi, val, width, bottom=bottom, color=colour,
                       edgecolor="white", linewidth=0.3)
                # label components large enough to read
                if val > 2.5:
                    ax.text(xi, bottom + val / 2, f"{val:.0f}",
                            ha="center", va="center", fontsize=6, color="white", fontweight="bold")
                bottom += val
        else:
            gwp_s = float(row["gwp_struct"])
            ax.bar(xi, gwp_s, width, bottom=bottom, color=slab_colour,
                   edgecolor="white", linewidth=0.3)
            if gwp_s > 5:
                ax.text(xi, bottom + gwp_s / 2, f"{gwp_s:.0f}",
                        ha="center", va="center", fontsize=6.5, color="white", fontweight="bold")
            bottom += gwp_s

        # ── floor build-up (total mode only) ─────────────────────────────
        if mode == "total":
            layers = _floor_layers(row)
            for lname, lgwp in layers:
                colour = floor_cmap.get(lname, "#cccccc")
                ax.bar(xi, lgwp, width, bottom=bottom, color=colour,
                       edgecolor="white", linewidth=0.3, hatch="//",
                       alpha=0.85)
                bottom += lgwp

        max_y = max(max_y, bottom)

    ax.set_xticks(x)
    ax.set_xticklabels([f"{s:.0f}" for s in spans], fontsize=10)
    ax.set_title(slab_label, fontsize=10, pad=4)
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.grid(axis="y", which="major", linewidth=0.4)
    ax.set_axisbelow(True)
    return max_y


def _make_legend_handles(mode: str, floor_cmap: dict,
                          slab_types_present: list) -> tuple:
    """Return (handles, labels) for the figure legend."""
    handles, labels = [], []

    # composite components
    for comp, colour in zip(COMP_COMPONENTS, COMP_COLOURS):
        handles.append(mpatches.Patch(facecolor=colour, edgecolor="white", linewidth=0.3))
        labels.append(comp)

    # non-composite slab types
    for slab_type, label, colour in SLAB_CONFIG:
        if slab_type == "comp_slab":
            continue
        if slab_type not in slab_types_present:
            continue
        handles.append(mpatches.Patch(facecolor=colour, edgecolor="white", linewidth=0.3))
        labels.append(label)

    # floor layers (total mode only)
    if mode == "total" and floor_cmap:
        for lname, colour in floor_cmap.items():
            handles.append(mpatches.Patch(facecolor=colour, edgecolor="grey",
                                          linewidth=0.3, hatch="//", alpha=0.85))
            labels.append(lname)

    return handles, labels


# ── STANDALONE COMPOSITE FIGURE ───────────────────────────────────────────────

def plot_composite_standalone(df: pd.DataFrame, best: dict, floor_cmap: dict,
                               spans: list, csv_file: str,
                               plot_dir: str = ".") -> None:
    """
    Two standalone figures for the composite slab only (no subplots):
      Fig A  —  Structural GWP breakdown  (absolute + percentage side by side)
      Fig B  —  Total GWP  (structural components + floor layers)
    """
    rows = best.get("comp_slab", {})
    valid_spans = [s for s in spans if rows.get(s) is not None]
    if not valid_spans:
        print("  [skip] no composite slab data for standalone plot")
        return

    x     = np.arange(len(valid_spans))
    width = 0.55
    base  = os.path.splitext(os.path.basename(csv_file))[0]
    use_case = df["use_case"].iloc[0].capitalize() if "use_case" in df.columns else ""
    n_sp_val = df["n_spans"].dropna().unique()
    n_sp_str = f"{int(n_sp_val[0])}-span" if len(n_sp_val) == 1 else "multi-span"

    # ── Fig A: absolute + percentage, side by side ────────────────────────────
    fig_a, (ax_abs, ax_pct) = plt.subplots(1, 2, figsize=(13, 5))

    values = {comp: [] for comp in COMP_COMPONENTS}
    deck_names = []
    for span in valid_spans:
        row = rows[span]
        for comp, col in zip(COMP_COMPONENTS, COMP_COLS):
            val = float(row[col]) if col in row.index and pd.notna(row.get(col)) else 0.0
            values[comp].append(max(val, 0.0))
        deck_names.append(str(row.get("deck_name", "")).strip())

    for ax, is_pct in [(ax_abs, False), (ax_pct, True)]:
        bottom = np.zeros(len(valid_spans))
        totals = np.array([sum(values[c][i] for c in COMP_COMPONENTS)
                           for i in range(len(valid_spans))])

        for comp, colour in zip(COMP_COMPONENTS, COMP_COLOURS):
            vals = np.array(values[comp])
            data = vals / totals * 100 if is_pct else vals
            bars = ax.bar(x, data, width, bottom=bottom,
                          label=comp, color=colour,
                          edgecolor="white", linewidth=0.4)
            threshold = 3.0 if is_pct else 2.0
            for rect, v in zip(bars, data):
                if v > threshold:
                    ax.text(rect.get_x() + rect.get_width() / 2,
                            rect.get_y() + rect.get_height() / 2,
                            f"{v:.0f}{'%' if is_pct else ''}",
                            ha="center", va="center", fontsize=7,
                            color="white", fontweight="bold")
            bottom += data

        # total GWP + deck name above each bar (absolute panel only)
        if not is_pct:
            for i, (total, dn) in enumerate(zip(totals, deck_names)):
                lbl = f"{total:.1f}"
                if dn:
                    lbl += f"\n{dn}"
                ax.text(i, bottom[i] + 0.5, lbl,
                        ha="center", va="bottom", fontsize=6.5)

        ax.set_xticks(x)
        ax.set_xticklabels([f"{s:.0f} m" for s in valid_spans], fontsize=10)
        ax.set_xlabel("Span length", fontsize=12)
        if is_pct:
            ax.set_ylabel("Share of GWP$_{struct}$  [%]", fontsize=12)
            ax.set_ylim(0, 108)
            ax.legend(loc="upper right", fontsize=9)
        else:
            ax.set_ylabel("GWP$_{struct}$  [kg-CO$_2$-eq/m²]", fontsize=12)
            ax.legend(loc="upper left", fontsize=9)
        ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
        ax.grid(axis="y", which="major", linewidth=0.5)
        ax.set_axisbelow(True)

    fig_a.suptitle(
        f"Composite slab — structural GWP breakdown  ({use_case}, {n_sp_str})",
        fontsize=10)
    fig_a.tight_layout(rect=[0, 0, 1, 0.95])
    if SAVE_FIGS:
        fname = os.path.join(plot_dir, f"{base}_gwp_comp_struct.png")
        fig_a.savefig(fname, dpi=DPI, bbox_inches="tight")
        print(f"Saved: {fname}")

    # ── Fig B: total GWP (structural components + floor layers) ─────────────
    fig_b, ax_tot = plt.subplots(figsize=(8, 5))
    bottom = np.zeros(len(valid_spans))

    # structural components
    for comp, col, colour in zip(COMP_COMPONENTS, COMP_COLS, COMP_COLOURS):
        vals = np.array(values[comp])
        bars = ax_tot.bar(x, vals, width, bottom=bottom,
                          label=comp, color=colour,
                          edgecolor="white", linewidth=0.4)
        for rect, v in zip(bars, vals):
            if v > 2.0:
                ax_tot.text(rect.get_x() + rect.get_width() / 2,
                            rect.get_y() + rect.get_height() / 2,
                            f"{v:.0f}", ha="center", va="center",
                            fontsize=7, color="white", fontweight="bold")
        bottom += vals

    # floor layers stacked on top
    # collect all layer names in order of first appearance
    seen_layers = []
    for span in valid_spans:
        for lname, _ in _floor_layers(rows[span]):
            if lname not in seen_layers:
                seen_layers.append(lname)

    for lname in seen_layers:
        layer_vals = []
        for span in valid_spans:
            layer_gwp = next((g for n, g in _floor_layers(rows[span]) if n == lname), 0.0)
            layer_vals.append(layer_gwp)
        layer_vals = np.array(layer_vals)
        colour = floor_cmap.get(lname, "#cccccc")
        bars = ax_tot.bar(x, layer_vals, width, bottom=bottom,
                          label=lname, color=colour, hatch="//",
                          edgecolor="grey", linewidth=0.3, alpha=0.85)
        bottom += layer_vals

    # total GWP label above each bar
    for i, (tot, dn) in enumerate(zip(bottom, deck_names)):
        lbl = f"{tot:.1f}"
        if dn:
            lbl += f"\n{dn}"
        ax_tot.text(i, tot + 0.5, lbl,
                    ha="center", va="bottom", fontsize=6.5)

    ax_tot.set_xticks(x)
    ax_tot.set_xticklabels([f"{s:.0f} m" for s in valid_spans], fontsize=10)
    ax_tot.set_xlabel("Span length", fontsize=12)
    ax_tot.set_ylabel("GWP$_{tot}$  [kg-CO$_2$-eq/m²]", fontsize=12)
    ax_tot.set_title(
        f"Composite slab — total GWP breakdown  ({use_case}, {n_sp_str})",
        fontsize=10)
    ax_tot.legend(loc="upper left", fontsize=9, ncol=2)
    ax_tot.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax_tot.grid(axis="y", which="major", linewidth=0.5)
    ax_tot.set_axisbelow(True)
    fig_b.tight_layout()

    if SAVE_FIGS:
        fname = os.path.join(plot_dir, f"{base}_gwp_comp_total.png")
        fig_b.savefig(fname, dpi=DPI, bbox_inches="tight")
        print(f"Saved: {fname}")


# ── MAIN ───────────────────────────────────────────────────────────────────────

def run(csv_file: str = CSV_FILE) -> None:
    plt.rcParams.update({"font.family": "Times New Roman"})
    if not os.path.exists(csv_file):
        print(f"CSV not found: {csv_file}")
        return

    df = pd.read_csv(csv_file)

    # determine spans
    all_spans = sorted(df["span_m"].dropna().unique().tolist())
    spans = SPANS_TO_SHOW if SPANS_TO_SHOW is not None else all_spans
    spans = [s for s in spans if s in all_spans]

    best = _best_rows(df)
    floor_cmap = _build_floor_color_map(best)

    # output folder: plots/<use_case>/GWP breakdown/
    _use_case_raw = df["use_case"].iloc[0] if "use_case" in df.columns else "unknown"
    _plot_dir = os.path.join("plots", _use_case_raw, "GWP breakdown")
    os.makedirs(_plot_dir, exist_ok=True)

    # which slab types actually have data?
    present = [st for st, _, _ in SLAB_CONFIG
               if any(v is not None for v in best.get(st, {}).values())]

    # subplot grid: 3 rows × 2 cols (last cell used for legend)
    n_cols = 2
    n_rows = 3
    fig_w, fig_h = 10, 13

    for mode, mode_label, fig_suffix in [
        ("struct", "Structural GWP",              "struct"),
        ("total",  "Total GWP  (structural + floor build-up)", "total"),
    ]:
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_w, fig_h),
                                 sharey=False)
        axes_flat = axes.flatten()

        # hide the 6th (unused) subplot
        for ai in range(len(SLAB_CONFIG), len(axes_flat)):
            axes_flat[ai].set_visible(False)

        all_max_y = 0.0

        for ai, (slab_type, slab_label, slab_colour) in enumerate(SLAB_CONFIG):
            ax = axes_flat[ai]
            rows_for_type = best.get(slab_type, {})
            max_y = _draw_subplot(
                ax, slab_type, slab_label, slab_colour,
                rows_for_type, spans, mode, floor_cmap)
            all_max_y = max(all_max_y, max_y)

            # y-axis label only on left column
            if ai % n_cols == 0:
                ax.set_ylabel("GWP  [kg-CO$_2$-eq/m²]", fontsize=12)
            # x-axis label only on bottom row
            if ai >= (n_rows - 1) * n_cols:
                ax.set_xlabel("Span  [m]", fontsize=12)

        # apply shared y scale with 10 % headroom
        y_top = all_max_y * 1.10
        for ai in range(len(SLAB_CONFIG)):
            axes_flat[ai].set_ylim(0, y_top)

        # global legend in the 6th (empty) subplot position
        legend_ax = axes_flat[len(SLAB_CONFIG)]
        legend_ax.set_visible(True)
        legend_ax.axis("off")
        handles, labels = _make_legend_handles(mode, floor_cmap, present)
        legend_ax.legend(handles, labels, loc="center", fontsize=8.5,
                         frameon=True, framealpha=0.9,
                         title="GWP components", title_fontsize=9)

        # CSV name in suptitle
        use_case = df["use_case"].iloc[0] if "use_case" in df.columns else ""
        n_sp_val = df["n_spans"].dropna().unique()
        n_sp_str = f"{int(n_sp_val[0])}-span" if len(n_sp_val) == 1 else "multi-span"
        fig.suptitle(
            f"{mode_label}  —  {use_case.capitalize()} ({n_sp_str})",
            fontsize=10)
        fig.tight_layout(rect=[0, 0, 1, 0.95])

        if SAVE_FIGS:
            base = os.path.splitext(os.path.basename(csv_file))[0]
            fname = os.path.join(_plot_dir, f"{base}_gwp_{fig_suffix}.png")
            fig.savefig(fname, dpi=DPI, bbox_inches="tight")
            print(f"Saved: {fname}")

    # standalone composite figure
    plot_composite_standalone(df, best, floor_cmap, spans, csv_file, _plot_dir)

    plt.show()


if __name__ == "__main__":
    run(CSV_FILE)
