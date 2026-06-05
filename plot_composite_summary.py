"""
plot_composite_summary.py

Produces three publication figures from pre-computed CSV results:
  Fig 1 — Office composite slab GWP vs span (1-, 2-, 3-span configurations)
  Fig 2 — Composite vs RC flat slab GWP vs span, all occupancies (single span)
  Fig 3 — Manufacturer load-span table verification dot plot (L / L_max)

Run directly:  python plot_composite_summary.py
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry import Polygon

BASE = os.path.dirname(os.path.abspath(__file__))

# ── matplotlib style — matches plot_datasets.py ───────────────────────────────
plt.rcParams.update({
    'font.family':    'Times New Roman',
    'font.size':      11,
    'axes.grid':      True,
    'grid.color':     '#cccccc',
    'grid.linewidth': 0.6,
    'axes.axisbelow': True,
})

# Colours from plot_datasets.py
C_COMP   = 'steelblue'     # composite slab (unpropped)
C_PROPPED= 'mediumpurple'  # composite slab (propped / 2-span)
C_RC_SOL = 'green'         # RC solid
C_RC_RIB = 'limegreen'     # RC ribbed
C_WD_SOL = 'saddlebrown'   # Timber solid
C_WD_RIB = 'sandybrown'    # Timber ribbed

# ── helpers ───────────────────────────────────────────────────────────────────

def _is_propped(val):
    return str(val).strip().lower() in ('true', '1')


def _agg_per_span(df, spans, slab_type, col='gwp_tot'):
    """Min, max, mean of *col* per span for a given slab type.
    For comp_slab: unpropped rows are preferred over propped."""
    mins, maxs, means = [], [], []
    for s in spans:
        rows = df[(df['slab_type'] == slab_type) & (df['span_m'] == s)].dropna(subset=[col])
        if slab_type == 'comp_slab':
            unp = rows[~rows['propped'].apply(_is_propped)]
            rows = unp if not unp.empty else rows
        if rows.empty:
            mins.append(np.nan); maxs.append(np.nan); means.append(np.nan)
        else:
            mins.append(rows[col].min())
            maxs.append(rows[col].max())
            means.append(rows[col].mean())
    return np.array(mins), np.array(maxs), np.array(means)


def _fill_band(ax, spans, lo, hi, color, alpha=0.20):
    """Polygon shaded band — alpha and style match plot_datasets.py."""
    valid = ~(np.isnan(lo) | np.isnan(hi))
    if valid.sum() < 2:
        return
    xs = np.array(spans)[valid]
    coords = list(zip(xs, hi[valid])) + list(zip(xs[::-1], lo[valid][::-1]))
    x, y = Polygon(coords).exterior.xy
    ax.fill(x, y, alpha=alpha, facecolor=color, edgecolor=color, linewidth=1.5)


def _line(ax, spans, means, color, ls, marker, label, lw=1.5, ms=6):
    """Mean line — linewidth matches plot_datasets.py."""
    valid = ~np.isnan(means)
    xs = np.array(spans)[valid]
    ys = means[valid]
    if len(xs):
        ax.plot(xs, ys, color=color, linestyle=ls, marker=marker,
                markersize=ms, linewidth=lw, label=label)


# ==============================================================================
# FIGURE 1 — Office composite slab: GWP by span configuration
# ==============================================================================

SPANS1 = [3, 4, 5, 6, 7, 8, 9]

df_1sp = pd.read_csv(os.path.join(BASE, 'results_office_ENV.csv'))
df_2sp = pd.read_csv(os.path.join(BASE, 'results_office_2span_ENV.csv'))
df_3sp = pd.read_csv(os.path.join(BASE, 'results_office_3span_ENV.csv'))

mn1, mx1, av1 = _agg_per_span(df_1sp, SPANS1, 'comp_slab')
mn2, mx2, av2 = _agg_per_span(df_2sp, SPANS1, 'comp_slab')
mn3, mx3, av3 = _agg_per_span(df_3sp, SPANS1, 'comp_slab')

# Blue / orange / green — maximally distinct, colourblind-friendly
C1 = 'steelblue'   # single span
C2 = 'darkorange'  # two-span
C3 = 'seagreen'    # three-span

fig1, ax1 = plt.subplots(figsize=(10, 6))
ax1.set_title('Office — composite slab GWP by span configuration')

_fill_band(ax1, SPANS1, mn1, mx1, C1)
_fill_band(ax1, SPANS1, mn2, mx2, C2)
_fill_band(ax1, SPANS1, mn3, mx3, C3)

_line(ax1, SPANS1, av1, C1, '-',   'o', 'single span')
_line(ax1, SPANS1, av2, C2, '--',  's', 'two-span continuous')
_line(ax1, SPANS1, av3, C3, '-.',  '^', 'three-span continuous')

ax1.set_xlabel(r'Span $l$ [m]')
ax1.set_ylabel(r'GWP$_{tot}$ [kg-CO$_2$-eq/m$^2$]')
ax1.set_xticks(SPANS1)
ax1.legend(title='Composite slab', frameon=True, loc='upper left')
fig1.tight_layout()
fig1.savefig(os.path.join(BASE, 'fig_composite_span_configs.png'), dpi=200, bbox_inches='tight')
print('Saved fig_composite_span_configs.png')


# ==============================================================================
# FIGURE 2 — Composite vs RC flat slab: all occupancies, single span
# ==============================================================================

SPANS2 = [3, 4, 5, 6, 7, 8, 9, 10]

# Occupancy colours drawn from the same optimizer palette
OCCUPANCIES = [
    ('office',       'results_office_ENV.csv',       C_COMP,    'o', 'Office'),
    ('residential',  'results_residential_ENV.csv',  C_PROPPED, 's', 'Residential'),
    ('retail',       'results_retail_ENV.csv',        C_RC_SOL, '^', 'Retail'),
    ('industrial',   'results_industrial_ENV.csv',    C_WD_SOL, 'D', 'Industrial'),
]

fig2, (ax2l, ax2r) = plt.subplots(1, 2, figsize=(14, 6), sharey=False)
fig2.suptitle('Composite vs RC flat slab across occupancies (single span)')

for _, fname, color, marker, label in OCCUPANCIES:
    df = pd.read_csv(os.path.join(BASE, fname))

    # structural GWP (left subplot)
    mn_cs, mx_cs, av_cs = _agg_per_span(df, SPANS2, 'comp_slab', col='gwp_struct')
    mn_rs, mx_rs, av_rs = _agg_per_span(df, SPANS2, 'rc_rec',    col='gwp_struct')
    _fill_band(ax2l, SPANS2, mn_cs, mx_cs, color, alpha=0.20)
    _fill_band(ax2l, SPANS2, mn_rs, mx_rs, color, alpha=0.20)
    _line(ax2l, SPANS2, av_cs, color, '-',  marker, f'{label} — composite')
    _line(ax2l, SPANS2, av_rs, color, '--', marker, f'{label} — RC flat')

    # total GWP (right subplot)
    mn_ct, mx_ct, av_ct = _agg_per_span(df, SPANS2, 'comp_slab', col='gwp_tot')
    mn_rt, mx_rt, av_rt = _agg_per_span(df, SPANS2, 'rc_rec',    col='gwp_tot')
    _fill_band(ax2r, SPANS2, mn_ct, mx_ct, color, alpha=0.20)
    _fill_band(ax2r, SPANS2, mn_rt, mx_rt, color, alpha=0.20)
    _line(ax2r, SPANS2, av_ct, color, '-',  marker, '_nolegend_')
    _line(ax2r, SPANS2, av_rt, color, '--', marker, '_nolegend_')

for ax, ylabel in [(ax2l, r'GWP$_{struct}$ [kg-CO$_2$-eq/m$^2$]'),
                   (ax2r, r'GWP$_{tot}$ [kg-CO$_2$-eq/m$^2$]')]:
    ax.set_xlabel(r'Span $l$ [m]')
    ax.set_ylabel(ylabel)
    ax.set_xticks(SPANS2)

ax2l.legend(ncol=2, fontsize=9, frameon=True, loc='upper left')
fig2.tight_layout()
fig2.savefig(os.path.join(BASE, 'fig_composite_vs_rc_occupancies.png'), dpi=200, bbox_inches='tight')
print('Saved fig_composite_vs_rc_occupancies.png')


# ==============================================================================
# FIGURE 3 — Manufacturer load-span table verification
# ==============================================================================

verif = pd.read_csv(os.path.join(BASE, 'mfr_table_verification.csv'))
verif = verif[verif['Result'].isin(['PASS', 'FAIL'])].copy()

# Deck display order: Multideck (80→60→146), ComFlor shallow (80→60→51→46), deep (225→210)
DECK_ORDER = [
    'Multideck 80 1.2', 'Multideck 80 1.1', 'Multideck 80 1.0',
    'Multideck 60 1.2', 'Multideck 60 1.1', 'Multideck 60 1.0', 'Multideck 60 0.9',
    'Multideck 146 1.5', 'Multideck 146 1.2',
    'ComFlor 80 1.2', 'ComFlor 80 0.9',
    'ComFlor 60 1.2', 'ComFlor 60 1.0', 'ComFlor 60 0.9',
    'ComFlor 51 1.2', 'ComFlor 51 1.0', 'ComFlor 51 0.9',
    'ComFlor 46 1.2', 'ComFlor 46 0.9',
    'ComFlor 225', 'ComFlor 210',
]

# Keep only decks present in the data; append any unexpected ones at the bottom
present = verif['Deck'].unique().tolist()
ordered = [d for d in DECK_ORDER if d in present]
ordered += sorted(set(present) - set(DECK_ORDER))

y_pos = {deck: i for i, deck in enumerate(ordered)}
n_decks = len(ordered)

fig3, ax3 = plt.subplots(figsize=(10, max(6, n_decks * 0.42)))

PASS_COLOR = '#4878d0'   # blue
FAIL_COLOR = '#d65f5f'   # red
ALPHA      = 0.55
SIZE       = 28

for _, row in verif.iterrows():
    deck = row['Deck']
    if deck not in y_pos:
        continue
    y   = y_pos[deck]
    x   = row['Span/Max']
    col = PASS_COLOR if row['Result'] == 'PASS' else FAIL_COLOR
    ax3.scatter(x, y, color=col, s=SIZE, alpha=ALPHA, zorder=3, linewidths=0)

# Manufacturer limit line
ax3.axvline(1.0, color='black', linestyle='--', linewidth=1.2, zorder=4,
            label=r'Manufacturer limit ($L/L_{max}=1$)')

ax3.set_yticks(range(n_decks))
ax3.set_yticklabels(ordered, fontsize=9)
ax3.set_xlabel(r'Design span / manufacturer maximum span, $L\,/\,L_{max}$')
ax3.set_xlim(left=max(0, verif['Span/Max'].min() - 0.05))
ax3.invert_yaxis()

# Custom legend
from matplotlib.lines import Line2D
legend_handles = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=PASS_COLOR,
           markersize=8, alpha=ALPHA, label='Pass (within manufacturer envelope)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=FAIL_COLOR,
           markersize=8, alpha=ALPHA, label='Fail (design span exceeds maximum)'),
    Line2D([0], [0], color='black', linestyle='--', linewidth=1.2,
           label=r'Manufacturer limit ($L/L_{max}=1$)'),
]
ax3.legend(handles=legend_handles, loc='upper right', ncol=1, fontsize=9, frameon=True)

fig3.tight_layout()
fig3.savefig(os.path.join(BASE, 'fig_mfr_verification.png'), dpi=200, bbox_inches='tight')
print('Saved fig_mfr_verification.png')

# ==============================================================================
# FIGURES 4–6: Criteria comparison for composite slabs (office)
# One figure per span configuration (single / two-span / three-span).
# Each figure: two subplots — gwp_struct (left) and gwp_tot (right).
# Four overlaid lines: Envelope (ENV), ULS-only, SLS-only, Fire-only.
# Shaded band = material spread (min–max) for each criterion.
# ==============================================================================

SPAN_CONFIGS = [
    ('SS',     'single span',         'results_office_ENV.csv',
               'results_office_ULS.csv',
               'results_office_SLS1.csv',
               'results_office_FIRE.csv'),
    ('2-span', 'two-span continuous', 'results_office_2span_ENV.csv',
               'results_office_2span_ULS.csv',
               'results_office_2span_SLS1.csv',
               'results_office_2span_FIRE.csv'),
    ('3-span', 'three-span continuous', 'results_office_3span_ENV.csv',
               'results_office_3span_ULS.csv',
               'results_office_3span_SLS1.csv',
               'results_office_3span_FIRE.csv'),
]

# label, display name, colour, linestyle, marker
CRITERIA_STYLE = [
    ('ENV',  'Envelope (ENV)', C_COMP,      '-',   'o'),
    ('ULS',  'ULS only',       'darkorange', '--',  's'),
    ('SLS1', 'SLS only',       'seagreen',   '-.',  '^'),
    ('FIRE', 'Fire only',      'firebrick',  ':',   'D'),
]

SPANS_CRIT = [3, 4, 5, 6, 7, 8, 9, 10]

for span_key, span_title, f_env, f_uls, f_sls, f_fire in SPAN_CONFIGS:
    crit_files = [
        ('ENV',  f_env),
        ('ULS',  f_uls),
        ('SLS1', f_sls),
        ('FIRE', f_fire),
    ]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(14, 5), sharey=False)
    fig.suptitle(
        f'Office — composite slab: design criteria comparison ({span_title})')

    for crit_key, fname in crit_files:
        df = pd.read_csv(os.path.join(BASE, fname))
        _, display, color, ls, marker = next(
            s for s in CRITERIA_STYLE if s[0] == crit_key)

        mn_s, mx_s, av_s = _agg_per_span(df, SPANS_CRIT, 'comp_slab',
                                          col='gwp_struct')
        mn_t, mx_t, av_t = _agg_per_span(df, SPANS_CRIT, 'comp_slab',
                                          col='gwp_tot')

        _fill_band(axL, SPANS_CRIT, mn_s, mx_s, color, alpha=0.15)
        _fill_band(axR, SPANS_CRIT, mn_t, mx_t, color, alpha=0.15)

        _line(axL, SPANS_CRIT, av_s, color, ls, marker, display)
        _line(axR, SPANS_CRIT, av_t, color, ls, marker, '_nolegend_')

    for ax, ylabel in [
        (axL, r'GWP$_{struct}$ [kg-CO$_2$-eq/m$^2$]'),
        (axR, r'GWP$_{tot}$  [kg-CO$_2$-eq/m$^2$]'),
    ]:
        ax.set_xlabel(r'Span $l$ [m]')
        ax.set_ylabel(ylabel)
        ax.set_xticks(SPANS_CRIT)

    axL.legend(frameon=True, loc='upper left', fontsize=9)
    fig.tight_layout()

    fname_out = os.path.join(
        BASE, f'fig_criteria_comparison_{span_key.replace("-", "")}.png')
    fig.savefig(fname_out, dpi=200, bbox_inches='tight')
    print(f'Saved {fname_out}')

plt.show()
