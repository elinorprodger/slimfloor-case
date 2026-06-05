import struct_analysis  # file with code for structural analysis
import struct_optimization  # file with code for structural optimization
import sqlite3  # import modul for SQLite
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from shapely.geometry import Polygon
from scipy.interpolate import interp1d
from scipy.spatial import ConvexHull
import pandas as pd
import os


# Mapping from check-name prefix → canonical CSV column name.
# Order matters: more specific prefixes must appear before shorter ones.
_CHECK_COL_MAP = [
    ("Construction ULS - Bending",         "util_constr_uls_bending"),
    ("Construction ULS - Shear",           "util_constr_uls_shear"),
    ("Construction SLS - Deflection",      "util_constr_sls_deflection"),
    ("Constr SLS - Internal Support",      "util_constr_sls_support"),
    ("Composite ULS - Bending",            "util_comp_uls_bending"),
    ("Composite ULS - Hogging",            "util_comp_uls_hogging"),
    ("Composite ULS - Longitudinal Shear", "util_comp_uls_long_shear"),
    ("Composite ULS - Vertical Shear",     "util_comp_uls_vert_shear"),
    ("Composite SLS - Deflection",         "util_comp_sls_deflection"),
    ("SLS - Vibration",                    "util_sls_vibration"),
    ("D1 -",                               "util_fire_d1_insulation"),
    ("D2 -",                               "util_fire_d2_sagging"),
    ("D3 -",                               "util_fire_d3_hogging"),
    ("D4 -",                               "util_fire_d4_heff"),
    ("D5 -",                               "util_fire_d5_field"),
]

# All possible check-util column names in a fixed order (ensures consistent CSV columns
# across runs even when some checks are absent for a given section/span).
_ALL_CHECK_COLS = [col for _, col in _CHECK_COL_MAP]


def _check_col_name(check_name):
    """Return the canonical column name for a check, or a sanitised fallback."""
    for prefix, col in _CHECK_COL_MAP:
        if check_name.startswith(prefix):
            return col
    import re
    return "util_" + re.sub(r"[^a-z0-9]+", "_", check_name.lower()).strip("_")


def _extract_check_utils(checks):
    """Return a dict {col_name: util} for every check in *checks*.

    Initialised with NaN for every column in _ALL_CHECK_COLS so the CSV
    always has a consistent set of columns regardless of which checks ran.
    """
    _nan = float('nan')
    result = {col: _nan for col in _ALL_CHECK_COLS}
    for c in checks:
        col = _check_col_name(c.get("name", ""))
        util = c.get("util", _nan)
        # keep the maximum if the same canonical column is hit more than once
        # (e.g. two D3 variants map to the same column)
        if col not in result or not (result[col] == result[col]) or util > result[col]:
            result[col] = util
    return result


class _CompSlabMember:
    """Lightweight wrapper giving composite slab results a Member1D-like interface."""
    def __init__(self, section, floorstruc):
        self.section = section
        self.floorstruc = floorstruc

# PLOT DATASETS OF MEMBERS WITH DEFINED CROSS_SECTIONS AND VARIED MATERIALS
# ----------------------------------------------------------------------------------------------------------------------
def _build_member_list(lengths, database_name, criteria, optima, floorstruc, requirements,
                       crsec_type, mat_names, g2k, qk, max_iter, idx_vrfctn,
                       deck, comp_slab_loads, manufacturer, concrete_grades, n_spans=1):
    """Run DB queries, create sections and optimise for all span lengths.
    Returns (member_list, legend, idx_vrfctn) — no plotting side-effects."""

    if idx_vrfctn == -1:
        idx_vrfctn = random.randint(0, len(lengths)-1)

    # GENERATE INITIAL CROSS-SECTIONS
    # Search database (table products, attribute material) for products
    # get prod_id of relevant materials from database and create initial cross-section for each product
    to_plot = []
    connection = sqlite3.connect(database_name)
    cursor = connection.cursor()
    for mat_name in mat_names:
        # Wählt alle EPDs vom Material "mat-name" (z.B. ready mixed concrete), welche sich gem. Spalte Statistik zwischen dem 10% und 90% Quantil befindet. Wo Source = Betonsortenrechenr, Ecoinvent oder KBOB ist, wird die Zeile nicht gewählt.
        # For composite slabs only: restrict concrete products to the grades
        # listed in `concrete_grades` (default ["C25/30"]) to cut optimisation
        # runtime. Every EPD product for each listed grade is kept so the
        # optimiser still picks the lowest-GWP one.
        # RC rectangular/ribbed slabs continue to iterate through every grade.
        if crsec_type == "comp_slab" and mat_name.strip("'").lower() == "ready_mixed_concrete":
            grades = concrete_grades if concrete_grades is not None else ["C25/30"]
            # build a safe SQL IN-list: escape single quotes in grade strings
            quoted = ", ".join("'" + g.replace("'", "''") + "'" for g in grades)
            grade_filter = " AND MECH_PROP IN (" + quoted + ")"
        else:
            grade_filter = ""
        inquiry = ("""
                SELECT PRO_ID FROM products
                WHERE DENSITY IS NOT NULL
                AND MECH_PROP IS NOT NULL
                AND Statistik = 1
                AND "SOURCE" NOT LIKE '%Betonsortenrechner%'
                AND "SOURCE" NOT LIKE '%Ecoinvent%'
                AND "SOURCE" NOT LIKE '%KBOB%'
                AND "MATERIAL" LIKE """ + mat_name + grade_filter
        )
        # inquiry = ("SELECT PRO_ID FROM products WHERE"
        #            " material=" + mat_name)
        cursor.execute(inquiry)
        result = cursor.fetchall()
        for i, prod_id in enumerate(result):
            prod_id_str = "'" + str(prod_id[0]) + "'"
            inquiry = ("""
                    SELECT MECH_PROP FROM products
                    WHERE  PRO_ID LIKE """ + prod_id_str
            )
            # inquiry = ("SELECT mech_prop FROM products WHERE"
            #            " PRO_ID=" + prod_id_str)
            cursor.execute(inquiry)
            result = cursor.fetchall()
            mech_prop = "'" + result[0][0] + "'"

            if crsec_type == "wd_rec":
                # create a Wood material object
                timber = struct_analysis.Wood(mech_prop, database_name, prod_id_str)
                timber.get_design_values()
                # create initial wooden rectangular cross-section
                section_0 = struct_analysis.RectangularWood(timber, 1.0, 0.1, xi=0.02)
                # add section to content-definition of plot-line
                line_i = [section_0, floorstruc]
                to_plot.append(line_i)

            elif crsec_type == "rc_rec":
                # create a Concrete material object
                concrete = struct_analysis.ReadyMixedConcrete(mech_prop, database_name, prod_id=prod_id_str)
                concrete.get_design_values()
                # search database for rebar material of type B500B with lowest and highes emissions
                # exclude not epd sources from the data.
                # only take values, which are inside an 80% confidence interval
                inquiry = ("""
                            SELECT PRO_ID FROM products
                            WHERE Total_GWP = (SELECT MIN(Total_GWP) FROM products
                                                WHERE "MATERIAL" LIKE '%Steel_reinforcing_bar%'
                                                AND DENSITY IS NOT NULL
                                                AND MECH_PROP IS NOT NULL
                                                AND Statistik = 1
                                                AND "SOURCE" NOT LIKE '%Betonsortenrechner%'
                                                AND "SOURCE" NOT LIKE '%Ecoinvent%'
                                                AND "SOURCE" NOT LIKE '%KBOB%')
                            OR Total_GWP = (SELECT MAX(Total_GWP) FROM products
                                                WHERE "MATERIAL" LIKE '%Steel_reinforcing_bar%'
                                                AND DENSITY IS NOT NULL
                                                AND MECH_PROP IS NOT NULL
                                                AND Statistik = 1
                                                AND "SOURCE" NOT LIKE '%Betonsortenrechner%'
                                                AND "SOURCE" NOT LIKE '%Ecoinvent%'
                                                AND "SOURCE" NOT LIKE '%KBOB%')
                            """
                           )
                cursor.execute(inquiry)
                result = cursor.fetchall()
                prod_id_low = result[0]
                prod_id_low_str = "'" + str(prod_id_low[0]) + "'"
                prod_id_high = result[1]
                prod_id_high_str = "'" + str(prod_id_high[0]) + "'"
                # create a rebar material objects with mech prop B500B and low rsp high emission values
                rebar_low_em = struct_analysis.SteelReinforcingBar("'B500B'", database_name, prod_id=prod_id_low_str)
                rebar_high_em = struct_analysis.SteelReinforcingBar("'B500B'", database_name, prod_id=prod_id_high_str)
                # create initial cross-sections
                section_00 = struct_analysis.RectangularConcrete(concrete, rebar_low_em, 1.0, 0.20,
                                                                0.014, 0.15, 0.01, 0.15, 0.01, 0.15, 0.01, 0.15,
                                                                0, 0.15, 2)
                section_01 = struct_analysis.RectangularConcrete(concrete, rebar_high_em, 1.0, 0.20,
                                                                 0.014, 0.15, 0.01, 0.15, 0.01, 0.15, 0.01, 0.15,
                                                                 0, 0.15, 2)

                # add sections to content-definition of plot-line
                line_i0 = [section_00, floorstruc]
                line_i1 = [section_01, floorstruc]
                to_plot.extend([line_i0, line_i1])

            elif crsec_type == "rc_rib":
                # create a Concrete material object
                concrete = struct_analysis.ReadyMixedConcrete(mech_prop, database_name, prod_id=prod_id_str)
                concrete.get_design_values()
                # search database for rebar material of type B500B with lowest and highes emissions
                inquiry = ("""
                                            SELECT PRO_ID FROM products
                                            WHERE Total_GWP = (SELECT MIN(Total_GWP) FROM products
                                                                WHERE "MATERIAL" LIKE '%Steel_reinforcing_bar%'
                                                                AND DENSITY IS NOT NULL
                                                                AND MECH_PROP IS NOT NULL
                                                                AND Statistik = 1
                                                                AND "SOURCE" NOT LIKE '%Betonsortenrechner%'
                                                                AND "SOURCE" NOT LIKE '%Ecoinvent%'
                                                                AND "SOURCE" NOT LIKE '%KBOB%')
                                            OR Total_GWP = (SELECT MAX(Total_GWP) FROM products
                                                                WHERE "MATERIAL" LIKE '%Steel_reinforcing_bar%'
                                                                AND DENSITY IS NOT NULL
                                                                AND MECH_PROP IS NOT NULL
                                                                AND Statistik = 1
                                                                AND "SOURCE" NOT LIKE '%Betonsortenrechner%'
                                                                AND "SOURCE" NOT LIKE '%Ecoinvent%'
                                                                AND "SOURCE" NOT LIKE '%KBOB%')
                                            """
                           )
                cursor.execute(inquiry)
                result = cursor.fetchall()
                prod_id_low = result[0]
                prod_id_low_str = "'" + str(prod_id_low[0]) + "'"
                prod_id_high = result[1]
                prod_id_high_str = "'" + str(prod_id_high[0]) + "'"

                # create a rebar material objects with mech prop B500B and low rsp high emission values
                rebar_low_em = struct_analysis.SteelReinforcingBar("'B500B'", database_name, prod_id=prod_id_low_str)
                rebar_high_em = struct_analysis.SteelReinforcingBar("'B500B'", database_name, prod_id=prod_id_high_str)

                # create initial cross-sections
                section_00 = struct_analysis.RibbedConcrete(concrete, rebar_low_em, 4, 1.0, 0.15, 0.3, 0.18, 0.01, 0.15, 0.01, 0.15, 0.02, 2, 0.01, 0.15, 2)
                section_01 = struct_analysis.RibbedConcrete(concrete, rebar_high_em, 4, 1.0, 0.15, 0.3, 0.18, 0.01, 0.15,
                                                            0.01, 0.15, 0.02, 2, 0.01, 0.15, 2)
                # add sections to content-definition of plot-line
                line_i0 = [section_00, floorstruc]
                line_i1 = [section_01, floorstruc]
                to_plot.extend([line_i0, line_i1])

            elif crsec_type == "wd_rib":
                # create a Wood material object
                timber1 = struct_analysis.Wood(mech_prop, database_name, prod_id=prod_id_str)  # create a Wood material object
                timber1.get_design_values()

                # search database for timber board material (CLT) with lowest and highes emissions
                inquiry = ("""
                                                            SELECT PRO_ID FROM products
                                                            WHERE Total_GWP = (SELECT MIN(Total_GWP) FROM products
                                                                                WHERE "MATERIAL" LIKE '%3- and 5-ply wood%'
                                                                                AND DENSITY IS NOT NULL
                                                                                AND MECH_PROP IS NOT NULL
                                                                                AND Statistik = 1
                                                                                AND "SOURCE" NOT LIKE '%Betonsortenrechner%'
                                                                                AND "SOURCE" NOT LIKE '%Ecoinvent%'
                                                                                AND "SOURCE" NOT LIKE '%KBOB%')
                                                            OR Total_GWP = (SELECT MAX(Total_GWP) FROM products
                                                                                WHERE "MATERIAL" LIKE '%3- and 5-ply wood%'
                                                                                AND DENSITY IS NOT NULL
                                                                                AND MECH_PROP IS NOT NULL
                                                                                AND Statistik = 1
                                                                                AND "SOURCE" NOT LIKE '%Betonsortenrechner%'
                                                                                AND "SOURCE" NOT LIKE '%Ecoinvent%'
                                                                                AND "SOURCE" NOT LIKE '%KBOB%')
                                                            """
                           )
                cursor.execute(inquiry)
                result = cursor.fetchall()
                prod_id_low = result[0]
                prod_id_low_str = "'" + str(prod_id_low[0]) + "'"
                prod_id_high = result[1]
                prod_id_high_str = "'" + str(prod_id_high[0]) + "'"

                # create a timber material objects in timber board (CLT, C24) with and low rsp high emission values
                clt_low_em = struct_analysis.Wood("'C24'", database_name, prod_id=prod_id_low_str)
                clt_low_em.get_design_values()
                clt_high_em = struct_analysis.Wood("'C24'", database_name, prod_id=prod_id_high_str)
                clt_high_em.get_design_values()

                # create initial cross-sections
                section_00 = struct_analysis.RibWood(timber1, clt_low_em, clt_low_em, 4, 0.12, 0.22, 0.625,
                                                     0.027, 0.027)
                section_01 = struct_analysis.RibWood(timber1, clt_high_em, clt_high_em, 4, 0.12, 0.22, 0.625,
                                                    0.027, 0.027)

                # add sections to content-definition of plot-line
                line_i0 = [section_00, floorstruc]
                line_i1 = [section_01, floorstruc]
                to_plot.extend([line_i0, line_i1])

            elif crsec_type == "comp_slab":
                # create a Concrete material object
                concrete = struct_analysis.ReadyMixedConcrete(mech_prop, database_name, prod_id=prod_id_str)
                concrete.get_design_values()
                loads = comp_slab_loads or {}

                # determine which deck products to use
                if deck is not None:
                    decks_to_use = [deck]
                else:
                    # query deck products from the sheeting_prop table,
                    # optionally restricted to one or more manufacturers
                    if manufacturer is None:
                        cursor.execute("SELECT DISTINCT Deck FROM sheeting_prop")
                    elif isinstance(manufacturer, (list, tuple)):
                        placeholders = ", ".join("?" * len(manufacturer))
                        cursor.execute(
                            f"SELECT DISTINCT Deck FROM sheeting_prop WHERE Manufacturer IN ({placeholders})",
                            list(manufacturer))
                    else:
                        cursor.execute("SELECT DISTINCT Deck FROM sheeting_prop WHERE Manufacturer = ?",
                                       (manufacturer,))
                    deck_names = cursor.fetchall()
                    decks_to_use = [struct_analysis.Decking(mech_prop=name[0], database=database_name)
                                    for name in deck_names]

                for dk in decks_to_use:
                    section_0 = struct_analysis.CompositeSlab(dk, concrete, 130,
                                                              database=database_name,
                                                              n_spans=n_spans,
                                                              **loads)
                    line_i = [section_0, floorstruc]
                    to_plot.append(line_i)

            else:
                print("cross-section type is not defined inside function plot_dataset()")



    # ANALYSIS AND OPTIMIZATION OF CROSS-SECTIONS
    member_list = []
    legend = []
    # create plot data
    for i in to_plot:
        for criterion in criteria:
            for optimum in optima:
                members = []
                for length in lengths:
                    section0 = i[0]
                    floorstruc = i[1]

                    if section0.section_type == "comp_slab":
                        # Unpropped series (with propped fallback where unpropped is infeasible)
                        opt_section = struct_optimization.opt_comp_slab(
                            section0, span=length, to_opt=optimum,
                            criterion=criterion, max_iter=max_iter, propped=False)
                        opt_section.run_all_checks(length, criterion)
                        if not opt_section.all_passed:
                            opt_section = struct_optimization.opt_comp_slab(
                                section0, span=length, to_opt=optimum,
                                criterion=criterion, max_iter=max_iter, propped=True)
                            opt_section.run_all_checks(length, criterion)
                        opt_member = _CompSlabMember(opt_section, floorstruc) if opt_section.all_passed else None
                        members.append(opt_member)
                    else:
                        if n_spans == 1:
                            sys = struct_analysis.BeamSimpleSup(length)
                        elif n_spans == 2:
                            sys = struct_analysis.BeamTwoSpan(length)
                        else:
                            sys = struct_analysis.BeamContinuousSupEl(length)
                        member0 = struct_analysis.Member1D(section0, sys, floorstruc, requirements, g2k, qk)
                        opt_section = struct_optimization.get_optimized_section(member0, criterion, optimum, max_iter)
                        opt_member = struct_analysis.Member1D(opt_section, sys, floorstruc, requirements, g2k, qk)
                        # search for an alternative solution for rectangular concrete section with lower minimal h and fill in floorstructure
                        if section0.section_type == "rc_rec":
                            # create floor structure for slim reinforced concrete cross-section
                            bodenaufbau_rcdecke_slim = [["'Parkett 2-Schicht werkversiegelt, 11 mm'", False, False],
                                                        ["'Unterlagsboden Zement, 85 mm'", False, False],
                                                        ["'Glaswolle'", 0.03, False], ["'Kies gebrochen'", 0.06, False]]
                            floorstruc_alt = struct_analysis.FloorStruc(bodenaufbau_rcdecke_slim, database_name)
                            member0_alt = struct_analysis.Member1D(section0, sys, floorstruc_alt, requirements, g2k, qk)
                            opt_section_alt = struct_optimization.get_optimized_section(member0_alt, criterion, optimum, max_iter, h_min=0.12)
                            opt_member_alt = struct_analysis.Member1D(opt_section_alt, sys, floorstruc_alt, requirements, g2k, qk)
                            # update opt_member, if alternative solution has lower GWP
                            if opt_member_alt.co2 < opt_member.co2:
                                opt_member = opt_member_alt
                        members.append(opt_member)
                if i[0].section_type == "comp_slab":
                    material_lg = i[0].concrete.mech_prop + " + " + i[0].deck.mech_prop
                elif i[0].section_type[0:2] == "rc":
                    material_lg = i[0].concrete_type.mech_prop + " + " + i[0].rebar_type.mech_prop
                elif i[0].section_type == "wd":
                    material_lg = i[0].wood_type.mech_prop
                elif i[0].section_type == "wd_rib":
                    material_lg = i[0].wood_type_1.mech_prop
                else:
                    material_lg = "error: section material is not defined"
                member_list.append(members)
                legend.append([i[0].section_type, material_lg, criterion, optimum])

    return member_list, legend, idx_vrfctn


# Maximum number of floor build-up layer slots written to the CSV.
# Set large enough to cover the deepest build-up in any script
# (timber ribbed + extra acoustic/fire layers = up to 7 layers).
_MAX_FLOOR_LAYERS = 8


def _rows_from_member_list(member_list, legend, lengths, use_case=None):
    """Convert member_list + legend into a list of dicts suitable for a pandas DataFrame.

    For composite slabs the row also includes section details, individual check
    utilisations, and the GWP component breakdown.  Non-composite rows carry NaN
    for composite-only columns so the schema is uniform across all section types.

    Floor build-up GWP is written to gwp_floor_L0 … gwp_floor_L{N-1} (fixed-width,
    NaN for unused slots) with matching floor_L0_name … floor_L{N-1}_name columns
    so layers can be identified when comparing across different occupancies.
    """
    _nan = float('nan')
    rows = []
    for members, (sec_typ, mat, cri, opt) in zip(member_list, legend):
        for span, mem in zip(lengths, members):
            sec = mem.section if mem is not None else None
            is_comp = sec is not None and hasattr(sec, 'deck')

            # composite slab: extract governing check, deflection info, section details
            gov_util, gov_name = _nan, ''
            d_tot, d_var, lim_tot, lim_var, d_util = _nan, _nan, _nan, _nan, _nan
            bar_dia     = _nan
            check_utils = {col: _nan for col in _ALL_CHECK_COLS}
            concrete_grade = ''
            mesh_type      = ''
            mesh_layers    = _nan
            n_spans_val    = _nan
            h_c_mm         = _nan

            if is_comp:
                bar_dia = getattr(sec, 'bar_dia_trough', _nan)
                checks  = getattr(sec, 'checks', [])

                for c in checks:
                    if c.get('name', '').startswith('D5'):
                        continue  # D5 is warning-only, excluded from governing
                    u = c.get('util', 0) or 0
                    if u > (gov_util if gov_util == gov_util else 0):
                        gov_util = u
                        gov_name = c.get('name', '')

                for c in checks:
                    if c.get('name') == 'Composite SLS - Deflection':
                        d_tot   = c.get('delta_total_mm', _nan)
                        d_var   = c.get('delta_var_mm',   _nan)
                        lim_tot = c.get('limit_total_mm', _nan)
                        lim_var = c.get('limit_var_mm',   _nan)
                        d_util  = c.get('util',           _nan)
                        break

                check_utils = _extract_check_utils(checks)

                # section design details
                concrete_grade = getattr(sec.concrete, 'mech_prop', '').strip("'")
                mesh_type   = sec.mesh.mech_prop if sec.mesh is not None else ''
                mesh_layers = getattr(sec, 'mesh_layers', _nan)
                n_spans_val = getattr(sec, 'n_spans', _nan)
                h_c_mm      = getattr(sec, 'h_c', _nan)

            # structural mass per m²: composite uses .mass property; others use g0k/9.81
            sec_mass = (sec.mass if is_comp else sec.g0k / 9.81) if sec is not None else _nan

            # floor build-up layer GWP (all section types; NaN for unused slots)
            floor_layers = mem.floorstruc.layers if mem is not None else []
            floor_gwp   = {f'gwp_floor_L{i}':      (floor_layers[i].co2
                                                      if i < len(floor_layers) else _nan)
                           for i in range(_MAX_FLOOR_LAYERS)}
            floor_names = {f'floor_L{i}_name':      (floor_layers[i].name.strip("'")
                                                      if i < len(floor_layers) else '')
                           for i in range(_MAX_FLOOR_LAYERS)}

            row = {
                # ── run metadata ─────────────────────────────────────────────
                'use_case':        use_case or '',
                'slab_type':       sec_typ,
                'material':        mat,
                'criterion':       cri,
                'optimum':         opt,
                'n_spans':         n_spans_val,
                # ── geometry ─────────────────────────────────────────────────
                'span_m':          span,
                'h_struct_m':      sec.h if sec is not None else _nan,
                'h_tot_m':         (sec.h + mem.floorstruc.h) if mem is not None else _nan,
                'h_c_mm':          h_c_mm,
                # ── environmental / mass ──────────────────────────────────────
                'gwp_struct':      sec.co2 if sec is not None else _nan,
                'gwp_tot':         (sec.co2 + mem.floorstruc.co2) if mem is not None else _nan,
                'mass_struct':     sec_mass,
                'mass_tot':        (sec_mass + mem.floorstruc.gk_area / 9.81) if mem is not None else _nan,
                # ── composite section details ─────────────────────────────────
                'deck_name':       sec.deck.mech_prop if is_comp else '',
                'concrete_grade':  concrete_grade,
                'mesh_type':       mesh_type,
                'mesh_layers':     mesh_layers,
                'propped':         bool(getattr(sec, 'propped', False)) if is_comp else False,
                'bar_dia':         bar_dia,
                # ── governing check ───────────────────────────────────────────
                'gov_util':        gov_util,
                'gov_name':        gov_name,
                # ── deflection ────────────────────────────────────────────────
                'd_tot_mm':        d_tot,
                'd_var_mm':        d_var,
                'lim_tot_mm':      lim_tot,
                'lim_var_mm':      lim_var,
                'd_util':          d_util,
                # ── GWP component breakdown (composite only; NaN for others) ──
                'gwp_concrete':    sec.co2_concrete    if is_comp else _nan,
                'gwp_deck':        sec.co2_deck        if is_comp else _nan,
                'gwp_mesh':        sec.co2_mesh        if is_comp else _nan,
                'gwp_trough_bars': sec.co2_trough_bars if is_comp else _nan,
            }
            # ── individual check utilisations (composite only) ────────────────
            row.update(check_utils)
            # ── floor build-up layer GWP (all section types) ──────────────────
            row.update(floor_gwp)
            row.update(floor_names)
            rows.append(row)
    return rows


def _save_to_csv(rows, path):
    """Append rows to a CSV file, writing the header only when the file does not yet exist."""
    df = pd.DataFrame(rows)
    write_header = not os.path.exists(path)
    df.to_csv(path, mode='a', header=write_header, index=False)


def plot_dataset(lengths, database_name, criteria, optima, floorstruc, requirements, crsec_type, mat_names,
                 g2k=0.75, qk=2.0, max_iter=100, idx_vrfctn=-1, deck=None, comp_slab_loads=None,
                 manufacturer=None, concrete_grades=None, results_csv=None, n_spans=1, use_case=None):
    member_list, legend, idx_vrfctn = _build_member_list(
        lengths, database_name, criteria, optima, floorstruc, requirements,
        crsec_type, mat_names, g2k, qk, max_iter, idx_vrfctn,
        deck, comp_slab_loads, manufacturer, concrete_grades, n_spans=n_spans)

    # optionally save results to CSV (append; caller is responsible for clearing the
    # file at the start of a fresh run)
    if results_csv is not None:
        _save_to_csv(_rows_from_member_list(member_list, legend, lengths, use_case=use_case), results_csv)

    # CREATE DATA OF ENVELOPE AREA OF DATASET
    # Helper: extract values from member lists, returning NaN for infeasible (None) members
    _nan = float('nan')
    _valid = lambda vals: [v for v in vals if v == v]

    def _sec_mass(mem):
        sec = mem.section
        return sec.mass if hasattr(sec, 'mass') else sec.g0k / 9.81

    def _agg(raw):
        mn  = [min((v for v in vals if v == v), default=_nan) for vals in zip(*raw)]
        mx  = [max((v for v in vals if v == v), default=_nan) for vals in zip(*raw)]
        avg = [sum(_valid(vals)) / len(_valid(vals)) if _valid(vals) else _nan for vals in zip(*raw)]
        return mn, mx, avg

    raw_h        = [[mem.section.h if mem is not None else _nan                        for mem in sl] for sl in member_list]
    raw_h_tot    = [[mem.section.h + mem.floorstruc.h if mem is not None else _nan     for mem in sl] for sl in member_list]
    raw_co2      = [[mem.section.co2 if mem is not None else _nan                      for mem in sl] for sl in member_list]
    raw_co2_tot  = [[mem.section.co2 + mem.floorstruc.co2 if mem is not None else _nan for mem in sl] for sl in member_list]
    raw_mass     = [[_sec_mass(mem) if mem is not None else _nan                        for mem in sl] for sl in member_list]
    raw_mass_tot = [[_sec_mass(mem) + mem.floorstruc.gk_area / 9.81 if mem is not None else _nan for mem in sl] for sl in member_list]

    h_min,    h_max,    h_mean    = _agg(raw_h)
    ht_min,   ht_max,   ht_mean   = _agg(raw_h_tot)
    co2_min,  co2_max,  co2_mean  = _agg(raw_co2)
    ct_min,   ct_max,   ct_mean   = _agg(raw_co2_tot)
    m_min,    m_max,    m_mean    = _agg(raw_mass)
    mt_min,   mt_max,   mt_mean   = _agg(raw_mass_tot)

    # idx→(fig, subplot): 0=h_struct,1=h_tot, 2=gwp_struct,3=gwp_tot, 4=mass_struct,5=mass_tot
    values_min  = [h_min,  ht_min,  co2_min, ct_min,  m_min,  mt_min]
    values_max  = [h_max,  ht_max,  co2_max, ct_max,  m_max,  mt_max]
    values_mean = [h_mean, ht_mean, co2_mean,ct_mean, m_mean, mt_mean]

    def _activate(idx):
        """Switch to correct figure and subplot for index idx."""
        plt.figure(idx // 2 + 1)
        plt.subplot(1, 2, idx % 2 + 1)

    # PLOT DATASET TO FIGURE
    plt.rcParams.update({'font.family': 'Times New Roman'})
    for fig_n in (1, 2, 3):
        plt.figure(fig_n)
    data_max = [0] * 6
    vrfctn_members = [[], []]
    for i, members in enumerate(member_list):
        plotdata = [[], [], [], [], [], []]
        for j, mem in enumerate(members):
            if mem is not None:
                plotdata[0].append(mem.section.h)
                plotdata[1].append(mem.section.h + mem.floorstruc.h)
                plotdata[2].append(mem.section.co2)
                plotdata[3].append(mem.section.co2 + mem.floorstruc.co2)
                plotdata[4].append(_sec_mass(mem))
                plotdata[5].append(_sec_mass(mem) + mem.floorstruc.gk_area / 9.81)
            else:
                for k in range(6):
                    plotdata[k].append(_nan)
            if j == idx_vrfctn and mem is not None:
                vrfctn_members[0].append(mem)
                vrfctn_members[1].append(i)
        sec_typ, mat, cri, opt = legend[i]
        # set line color
        if sec_typ == "rc_rec":       color = 'green'
        elif sec_typ == "wd_rec":     color = 'saddlebrown'
        elif sec_typ == "rc_rib":     color = 'limegreen'
        elif sec_typ == "wd_rib":     color = 'sandybrown'
        elif sec_typ == "comp_slab":  color = 'steelblue'
        else:                         color = "k"
        # set linestyle
        if cri == "ULS":    linestyle = "--"
        elif cri == "SLS1": linestyle = (0, (3, 1, 1, 1))
        elif cri == "SLS2": linestyle = ":"
        elif cri == "ENV":  linestyle = "-"
        else:               linestyle = (0, (1, 10))
        # set linewidth
        if opt == "h":     linewidth = 0.5
        elif opt == "GWP": linewidth = 1.0
        else:              linewidth = 0.1
        for idx, data in enumerate(plotdata):
            valid_data = [v for v in data if v == v]
            if valid_data:
                data_max[idx] = max(data_max[idx], max(valid_data))

    # Draw envelope fill and mean line once per plot_dataset call.
    for idx in range(6):
        _activate(idx)
        valid_max = [(l, v) for l, v in zip(lengths, values_max[idx]) if v == v]
        valid_min = [(l, v) for l, v in zip(lengths, values_min[idx]) if v == v]
        if len(valid_max) >= 2 and len(valid_min) >= 2:
            coords = valid_max + valid_min[::-1]
            polygon = Polygon(coords)
            x, y = polygon.exterior.xy
            plt.fill(x, y, alpha=0.20, facecolor=color, edgecolor=color, linewidth=1.5)
        # mean line – skip NaN spans
        valid_l_mean = [l for l, v in zip(lengths, values_mean[idx]) if v == v]
        valid_mean   = [v for v in values_mean[idx] if v == v]
        if valid_l_mean:
            plt.plot(valid_l_mean, valid_mean, color=color, linestyle=linestyle, linewidth=1.5)

    # For comp_slab: overlay spans where the minimum-GWP section is propped.
    comp_series = [ml for ml, leg in zip(member_list, legend) if leg[0] == 'comp_slab']
    if comp_series:
        _val_fns = [
            lambda m: m.section.h,
            lambda m: m.section.h + m.floorstruc.h,
            lambda m: m.section.co2,
            lambda m: m.section.co2 + m.floorstruc.co2,
            lambda m: _sec_mass(m),
            lambda m: _sec_mass(m) + m.floorstruc.gk_area / 9.81,
        ]
        for idx, vfn in enumerate(_val_fns):
            _activate(idx)
            p_lengths, p_means = [], []
            for j, length in enumerate(lengths):
                cands = [ml[j] for ml in comp_series if j < len(ml) and ml[j] is not None]
                if not cands:
                    continue
                # Prefer unpropped; fall back to propped only if no unpropped candidate.
                _unp = [m for m in cands if not getattr(m.section, 'propped', False)]
                best = min(_unp if _unp else cands, key=lambda m: m.section.co2)
                if getattr(best.section, 'propped', False):
                    val = vfn(best)
                    if val == val:
                        p_lengths.append(length)
                        p_means.append(val)
            if p_lengths:
                plt.plot(p_lengths, p_means, color='mediumpurple', linestyle=linestyle,
                         linewidth=2.0, marker='o', markersize=5, label='_nolegend_')

    return data_max, vrfctn_members


def run_dataset(lengths, database_name, criteria, optima, floorstruc, requirements, crsec_type, mat_names,
                g2k=0.75, qk=2.0, max_iter=100, idx_vrfctn=-1, deck=None, comp_slab_loads=None,
                manufacturer=None, concrete_grades=None, results_csv=None, n_spans=1,
                use_case=None):
    """Run optimisation for one slab type without producing any plots.

    Saves results to *results_csv* (appended) if a path is given.
    Returns (rows, vrfctn_best) where:
      rows         – list of dicts (same schema as the CSV)
      vrfctn_best  – for comp_slab: the lowest-GWP _CompSlabMember at idx_vrfctn;
                     for all other types: None
    """
    if idx_vrfctn == -1:
        import random
        idx_vrfctn = random.randint(0, len(lengths) - 1)

    member_list, legend, idx_vrfctn = _build_member_list(
        lengths, database_name, criteria, optima, floorstruc, requirements,
        crsec_type, mat_names, g2k, qk, max_iter, idx_vrfctn,
        deck, comp_slab_loads, manufacturer, concrete_grades, n_spans=n_spans)

    rows = _rows_from_member_list(member_list, legend, lengths, use_case=use_case)

    if results_csv is not None:
        _save_to_csv(rows, results_csv)

    # print per-span best for composite slabs (same table as the inline version)
    if crsec_type == "comp_slab":
        n_combos = sum(1 for leg in legend if leg[0] == 'comp_slab')

        # ── helper: governing check (highest util, D5 excluded as warning-only) ──
        _NAME_MAP = [
            ("Construction ULS - Bending",                          "Constr ULS bending"),
            ("Construction ULS - Shear",                            "Constr ULS shear"),
            ("Construction SLS - Deflection",                       "Constr SLS deflection"),
            ("Composite ULS - Longitudinal Shear (m-k)",            "Comp ULS long shear"),
            ("Composite ULS - Longitudinal Shear (partial)",        "Comp ULS long shear"),
            ("Composite ULS - Longitudinal Shear",                  "Comp ULS long shear"),
            ("Composite ULS - Vertical Shear",                      "Comp ULS vert shear"),
            ("Composite ULS - Hogging",                             "Comp ULS hogging"),
            ("Composite SLS - Deflection",                          "Comp SLS deflection"),
            ("D1 - Thermal insulation",                             "Fire D1 insulation"),
            ("D2 - Sagging moment",                                 "Fire D2 moment"),
            ("D3 - Hogging moment",                                 "Fire D3 hogging"),
            ("D4 - Minimum effective thickness",                    "Fire D4 h_eff"),
            ("SLS - Vibration",                                     "SLS vibration"),
        ]

        def _governing(section):
            checks = getattr(section, "checks", [])
            relevant = [c for c in checks if not c["name"].startswith("D5")]
            if not relevant:
                return 0.0, "—"
            gov = max(relevant, key=lambda c: c["util"])
            short = gov["name"]
            for full, abbr in _NAME_MAP:
                if gov["name"].startswith(full):
                    short = abbr
                    break
            return gov["util"], short

        def _deflection_info(section):
            """Return (delta_total_mm, delta_var_mm, limit_total_mm, limit_var_mm, util)
            from the Composite SLS - Deflection check, or None if not found."""
            checks = getattr(section, "checks", [])
            for c in checks:
                if c["name"] == "Composite SLS - Deflection":
                    return (c.get("delta_total_mm"), c.get("delta_var_mm"),
                            c.get("limit_total_mm"), c.get("limit_var_mm"),
                            c["util"])
            return None, None, None, None, None

        def _best_candidate(candidates):
            """Return the lowest-GWP unpropped member; fall back to lowest-GWP
            propped member only if no unpropped candidate exists."""
            unpropped = [m for m in candidates if not getattr(m.section, "propped", False)]
            pool = unpropped if unpropped else candidates
            return min(pool, key=lambda m: m.section.co2)

        # ── GWP breakdown diagnostic: top-5 candidates per span ─────────────────
        print(f"\nComposite slab — GWP breakdown, top 5 candidates per span:")
        hdr = (f"  {'L':>4}  {'Deck':<28} {'h':>5} {'GWP':>7}"
               f"  {'Conc':>6}  {'Deck':>6}  {'Mesh':>6}  {'Bars':>6}  {'Governing check'}")
        print(hdr)
        print("  " + "-" * (len(hdr) - 2))
        primary_series_dbg = [ml for ml, leg in zip(member_list, legend) if leg[0] == 'comp_slab']
        for j, length in enumerate(lengths):
            cands = [ml[j] for ml in primary_series_dbg if ml[j] is not None]
            cands_sorted = sorted(cands, key=lambda m: m.section.co2)[:5]
            for k, m in enumerate(cands_sorted):
                s = m.section
                _, gov = _governing(s)
                conc = getattr(s, 'co2_concrete',    float('nan'))
                deck = getattr(s, 'co2_deck',        float('nan'))
                mesh = getattr(s, 'co2_mesh',        float('nan'))
                bars = getattr(s, 'co2_trough_bars', float('nan'))
                l_lbl = f"{length:.1f}" if k == 0 else ""
                print(f"  {l_lbl:>4}  {s.deck.mech_prop:<28} {s.h_mm:>5.0f} {s.co2:>7.1f}"
                      f"  {conc:>6.1f}  {deck:>6.1f}  {mesh:>6.1f}  {bars:>6.1f}  {gov}")
            if cands_sorted:
                print()
        # ─────────────────────────────────────────────────────────────────────

        print(f"\nComposite slab — best section per span  ({n_combos} deck×concrete combinations):")
        print(f"  {'L [m]':<8} {'Deck profile':<32} {'h [mm]':>7} {'Bar ø [mm]':>10} "
              f"{'GWP [kg-CO2/m²]':>16}  {'Util':>5}  {'Governing check':<22}"
              f"  {'d_tot [mm]':>10}  {'d_lim [mm]':>10}  {'d_var [mm]':>10}  {'d_varlim[mm]':>12}  {'d util':>6}")
        print(f"  {'-'*8} {'-'*32} {'-'*7} {'-'*10} {'-'*16}  {'-'*5}  {'-'*22}"
              f"  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*12}  {'-'*6}")
        primary_series = [ml for ml, leg in zip(member_list, legend) if leg[0] == 'comp_slab']
        for j, length in enumerate(lengths):
            candidates = [ml[j] for ml in primary_series if ml[j] is not None]
            if candidates:
                best = _best_candidate(candidates)
                bar_dia = best.section.bar_dia_trough
                bar_str = f"{bar_dia:.0f}" if bar_dia is not None else "—"
                util, gov_name = _governing(best.section)
                d_tot, d_var, lim_tot, lim_var, d_util = _deflection_info(best.section)
                d_tot_s   = f"{d_tot:.2f}"   if d_tot   is not None else "—"
                d_var_s   = f"{d_var:.2f}"   if d_var   is not None else "—"
                lim_tot_s = f"{lim_tot:.1f}" if lim_tot is not None else "—"
                lim_var_s = f"{lim_var:.1f}" if lim_var is not None else "—"
                d_util_s  = f"{d_util:.3f}"  if d_util  is not None else "—"
                propped_flag = " [P]" if getattr(best.section, "propped", False) else ""
                print(f"  {length:<8.1f} {best.section.deck.mech_prop:<32s} "
                      f"{best.section.h_mm:>7.0f} {bar_str:>10} "
                      f"{best.section.co2:>16.1f}  {util:>5.3f}  {gov_name:<22}"
                      f"  {d_tot_s:>10}  {lim_tot_s:>10}  {d_var_s:>10}  {lim_var_s:>10}  {d_util_s:>6}"
                      f"{propped_flag}")
            else:
                print(f"  {length:<8.1f} {'— no feasible section —':<32}")
        print()

    # find the lowest-GWP comp_slab member at the verification span
    # Prefer unpropped; only use propped if no unpropped candidate exists.
    vrfctn_best = None
    if crsec_type == "comp_slab":
        primary_series = [ml for ml, leg in zip(member_list, legend) if leg[0] == 'comp_slab']
        candidates = [ml[idx_vrfctn] for ml in primary_series
                      if idx_vrfctn < len(ml) and ml[idx_vrfctn] is not None]
        if candidates:
            unpropped = [m for m in candidates if not getattr(m.section, "propped", False)]
            pool = unpropped if unpropped else candidates
            vrfctn_best = min(pool, key=lambda m: m.section.co2)

    return rows, vrfctn_best


def add_dataset_legend(has_propped=False, fig_nums=(1, 2, 3)):
    """Add a figure-level legend to each of the three comparison figures.

    Call this after all plot_dataset / plot_from_csv calls but before
    tight_layout / show.  The legend is placed below the subplots so the
    caller should use tight_layout(rect=[0, 0.08, 1, 1]) (or similar) to
    leave room for it.

    Parameters
    ----------
    has_propped : bool
        If True, include a 'Composite slab — propped' entry in mediumpurple.
    fig_nums : iterable of int
        Figure numbers to attach the legend to (default: 1, 2, 3).
    """
    import matplotlib.patches as mpatches
    import matplotlib.lines as mlines

    def _entry(color, label):
        """Filled band patch + solid mean line combined into one handle."""
        patch = mpatches.Patch(facecolor=color, alpha=0.20,
                               edgecolor=color, linewidth=1.5)
        line  = mlines.Line2D([], [], color=color, linewidth=1.5)
        return (patch, line), label

    # Build entries by material group so legend columns read:
    #   Col 1: Composite  |  Col 2: RC  |  Col 3: Timber
    # matplotlib fills legends row-by-row, so interleave the groups.

    h_comp_unp, _ = _entry('steelblue',   'Composite slab — unpropped')
    h_rc_sol,   _ = _entry('green',       'RC solid')
    h_rc_rib,   _ = _entry('limegreen',   'RC ribbed')
    h_wd_sol,   _ = _entry('saddlebrown', 'Timber solid')
    h_wd_rib,   _ = _entry('sandybrown',  'Timber ribbed')

    if has_propped:
        patch = mpatches.Patch(facecolor='mediumpurple', alpha=0.20,
                               edgecolor='mediumpurple', linewidth=1.5)
        line  = mlines.Line2D([], [], color='mediumpurple', linewidth=2.0,
                              linestyle='-', marker='o', markersize=5)
        h_comp_pro = (patch, line)
        # 6 entries, ncol=3 → row-major layout gives:
        # Col 1: Composite unpropped, Composite propped
        # Col 2: RC solid,            RC ribbed
        # Col 3: Timber solid,        Timber ribbed
        handles = [h_comp_unp, h_rc_sol,  h_wd_sol,
                   h_comp_pro, h_rc_rib,  h_wd_rib]
        labels  = ['Composite slab — unpropped', 'RC solid',   'Timber solid',
                   'Composite slab — propped',   'RC ribbed',  'Timber ribbed']
    else:
        # 5 entries — place composite first, then RC pair, then timber pair
        handles = [h_comp_unp, h_rc_sol, h_wd_sol, h_rc_rib, h_wd_rib]
        labels  = ['Composite slab — unpropped', 'RC solid', 'Timber solid',
                   'RC ribbed', 'Timber ribbed']

    for fig_n in fig_nums:
        fig = plt.figure(fig_n)
        fig.legend(
            handles, labels,
            handler_map={tuple: _TupleHandler()},
            loc='lower center',
            ncol=3,
            fontsize=9,
            frameon=True,
            bbox_to_anchor=(0.5, 0.01),
        )


class _TupleHandler(plt.matplotlib.legend_handler.HandlerBase):
    """Render a (patch, line) tuple as an overlaid patch + line."""
    def create_artists(self, legend, orig_handle, xdescent, ydescent,
                       width, height, fontsize, trans):
        patch, line = orig_handle
        # patch spanning full width
        from matplotlib.patches import FancyBboxPatch
        p = plt.matplotlib.patches.Rectangle(
            [xdescent, ydescent], width, height,
            facecolor=patch.get_facecolor(),
            edgecolor=patch.get_edgecolor(),
            linewidth=patch.get_linewidth(),
            alpha=patch.get_alpha(),
            transform=trans)
        # line across the centre
        import matplotlib.lines as mlines
        l = mlines.Line2D([xdescent, xdescent + width],
                          [ydescent + height / 2, ydescent + height / 2],
                          color=line.get_color(),
                          linewidth=line.get_linewidth(),
                          transform=trans)
        return [p, l]


def print_comp_slab_table_from_csv(csv_path, lengths):
    """Print the composite slab best-section table from a saved CSV.

    Reproduces the same output as the table printed during run_dataset so it
    is available when REPLOT_ONLY = True.  Requires the CSV to have been
    generated after the extended _rows_from_member_list (which saves gov_util,
    gov_name, bar_dia, deflection columns).
    """
    df = pd.read_csv(csv_path)
    comp = df[df['slab_type'] == 'comp_slab']
    if comp.empty:
        return

    def _is_propped(val):
        return str(val).strip().lower() in ('true', '1')

    n_combos = comp.groupby('span_m').size().max() if not comp.empty else 0
    print(f"\nComposite slab — best section per span  ({n_combos} deck×concrete combinations):")
    print(f"  {'L [m]':<8} {'Deck profile':<32} {'h [mm]':>7} {'Bar ø [mm]':>10} "
          f"{'GWP [kg-CO2/m²]':>16}  {'Util':>5}  {'Governing check':<22}"
          f"  {'d_tot [mm]':>10}  {'d_lim [mm]':>10}  {'d_var [mm]':>10}"
          f"  {'d_varlim[mm]':>12}  {'d util':>6}")
    print(f"  {'-'*8} {'-'*32} {'-'*7} {'-'*10} {'-'*16}  {'-'*5}  {'-'*22}"
          f"  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*12}  {'-'*6}")

    for span in lengths:
        rows = comp[comp['span_m'] == span].dropna(subset=['gwp_struct'])
        if rows.empty:
            print(f"  {span:<8.1f} {'— no feasible section —':<32}")
            continue
        # prefer lowest-GWP unpropped; fall back to lowest-GWP propped
        unpropped = rows[~rows['propped'].apply(_is_propped)]
        pool = unpropped if not unpropped.empty else rows
        best = pool.loc[pool['gwp_struct'].idxmin()]

        h_mm     = best['h_struct_m'] * 1000
        bar_dia  = best.get('bar_dia', float('nan'))
        bar_s    = f"{bar_dia:.0f}" if bar_dia == bar_dia else '—'
        util     = best.get('gov_util', float('nan'))
        util_s   = f"{util:.3f}" if util == util else '—'
        gov      = best.get('gov_name', '') or '—'
        d_tot    = best.get('d_tot_mm',  float('nan'))
        d_var    = best.get('d_var_mm',  float('nan'))
        lim_tot  = best.get('lim_tot_mm',float('nan'))
        lim_var  = best.get('lim_var_mm',float('nan'))
        d_util   = best.get('d_util',    float('nan'))
        d_tot_s  = f"{d_tot:.2f}"  if d_tot  == d_tot  else '—'
        d_var_s  = f"{d_var:.2f}"  if d_var  == d_var  else '—'
        l_tot_s  = f"{lim_tot:.1f}"if lim_tot== lim_tot else '—'
        l_var_s  = f"{lim_var:.1f}"if lim_var== lim_var else '—'
        d_util_s = f"{d_util:.3f}" if d_util == d_util  else '—'
        prop_s   = ' [P]' if _is_propped(best['propped']) else ''

        print(f"  {span:<8.1f} {best['deck_name']:<32s} "
              f"{h_mm:>7.0f} {bar_s:>10} "
              f"{best['gwp_struct']:>16.1f}  {util_s:>5}  {gov:<22}"
              f"  {d_tot_s:>10}  {l_tot_s:>10}  {d_var_s:>10}"
              f"  {l_var_s:>12}  {d_util_s:>6}{prop_s}")
    print()


def plot_from_csv(csv_path, lengths):
    """Load optimisation results from *csv_path* and produce the standard 2×2
    comparison plot on figure 1.  Returns data_max (4-element list) for axis
    scaling in the caller."""
    _nan = float('nan')
    df = pd.read_csv(csv_path)

    plt.rcParams.update({'font.family': 'Times New Roman'})
    for fig_n in (1, 2, 3):
        plt.figure(fig_n)
    data_max = [0] * 6

    def _color(slab_type):
        return {'rc_rec': 'green', 'wd_rec': 'saddlebrown',
                'rc_rib': 'limegreen', 'wd_rib': 'sandybrown',
                'comp_slab': 'steelblue'}.get(slab_type, 'k')

    def _linestyle(criterion):
        return '-'

    def _activate_csv(idx):
        plt.figure(idx // 2 + 1)
        plt.subplot(1, 2, idx % 2 + 1)

    cols = ['h_struct_m', 'h_tot_m', 'gwp_struct', 'gwp_tot', 'mass_struct', 'mass_tot']

    def _is_propped(val):
        """Robustly parse propped flag regardless of CSV dtype (bool/int/string)."""
        return str(val).strip().lower() in ('true', '1')

    # Pre-compute which spans have NO feasible unpropped comp_slab section.
    # Purple is only shown when every candidate at that span required propping.
    comp_df = df[df['slab_type'] == 'comp_slab']
    propped_spans = set()
    if not comp_df.empty and 'propped' in comp_df.columns:
        for span in lengths:
            span_rows = comp_df[comp_df['span_m'] == span].dropna(subset=['gwp_struct'])
            if span_rows.empty:
                continue
            has_unpropped = span_rows['propped'].apply(_is_propped).eq(False).any()
            if not has_unpropped:
                propped_spans.add(span)

    for (slab_type, criterion, optimum), group in df.groupby(
            ['slab_type', 'criterion', 'optimum']):
        color     = _color(slab_type)
        linestyle = _linestyle(criterion)

        mins, maxs, means = {c: [] for c in cols}, {c: [] for c in cols}, {c: [] for c in cols}
        for span in lengths:
            span_rows = group[group['span_m'] == span]
            # For comp_slab: exclude propped spans from the blue (unpropped) band.
            if slab_type == 'comp_slab' and span in propped_spans:
                for col in cols:
                    mins[col].append(_nan)
                    maxs[col].append(_nan)
                    means[col].append(_nan)
                continue
            # At mixed spans (some decks unpropped-feasible, some only propped),
            # unpropped is always preferred, so the blue band must use unpropped
            # rows only — a lower-GWP propped section must not pull the band down.
            if slab_type == 'comp_slab':
                span_rows = span_rows[~span_rows['propped'].apply(_is_propped)]
            for col in cols:
                feasible = span_rows[col].dropna().tolist()
                if feasible:
                    mins[col].append(min(feasible))
                    maxs[col].append(max(feasible))
                    means[col].append(sum(feasible) / len(feasible))
                else:
                    mins[col].append(_nan)
                    maxs[col].append(_nan)
                    means[col].append(_nan)

        for idx, col in enumerate(cols):
            _activate_csv(idx)
            valid_max = [(l, v) for l, v in zip(lengths, maxs[col]) if v == v]
            valid_min = [(l, v) for l, v in zip(lengths, mins[col]) if v == v]
            if len(valid_max) >= 2 and len(valid_min) >= 2:
                coords = valid_max + valid_min[::-1]
                polygon = Polygon(coords)
                x, y = polygon.exterior.xy
                plt.fill(x, y, alpha=0.20, facecolor=color, edgecolor=color, linewidth=1.5)
            valid_l = [l for l, v in zip(lengths, means[col]) if v == v]
            valid_m = [v for v in means[col] if v == v]
            if valid_l:
                plt.plot(valid_l, valid_m, color=color, linestyle=linestyle, linewidth=1.5)
            valid_maxvals = [v for v in maxs[col] if v == v]
            if valid_maxvals:
                data_max[idx] = max(data_max[idx], max(valid_maxvals))

    # Purple overlay: shaded band + mean line for propped spans.
    # Prepend the last unpropped span as a bridge point so the purple band
    # connects seamlessly to the blue band with no gap at the transition.
    if propped_spans and not comp_df.empty:
        # Use the same linestyle as the rest of the composite slab data
        comp_crit = comp_df['criterion'].iloc[0] if 'criterion' in comp_df.columns else 'ENV'
        propped_linestyle = _linestyle(comp_crit)

        col_map = [('h_struct_m', 0), ('h_tot_m', 1), ('gwp_struct', 2), ('gwp_tot', 3),
                   ('mass_struct', 4), ('mass_tot', 5)]

        sorted_propped = sorted(propped_spans)
        # find the span in lengths immediately before the first propped span
        first_propped_idx = next(
            (i for i, l in enumerate(lengths) if l == sorted_propped[0]), None)
        bridge_span = lengths[first_propped_idx - 1] \
            if first_propped_idx is not None and first_propped_idx > 0 else None

        for col, idx in col_map:
            _activate_csv(idx)
            p_lengths, p_mins, p_maxs, p_means = [], [], [], []

            # bridge point: use unpropped comp_slab data at last unpropped span
            if bridge_span is not None:
                bridge_rows = comp_df[comp_df['span_m'] == bridge_span].dropna(subset=[col])
                unpropped_bridge = bridge_rows[
                    ~bridge_rows['propped'].apply(_is_propped)]
                if not unpropped_bridge.empty:
                    feasible = unpropped_bridge[col].tolist()
                    p_lengths.append(bridge_span)
                    p_mins.append(min(feasible))
                    p_maxs.append(max(feasible))
                    p_means.append(sum(feasible) / len(feasible))

            for span in sorted_propped:
                span_rows = comp_df[comp_df['span_m'] == span].dropna(subset=[col])
                feasible = span_rows[col].tolist()
                if not feasible:
                    continue
                p_lengths.append(span)
                p_mins.append(min(feasible))
                p_maxs.append(max(feasible))
                p_means.append(sum(feasible) / len(feasible))

            # update data_max so the purple band is not clipped by the y-axis
            if p_maxs:
                data_max[idx] = max(data_max[idx], max(p_maxs))

            if len(p_lengths) >= 2:
                coords = list(zip(p_lengths, p_maxs)) + list(zip(reversed(p_lengths), reversed(p_mins)))
                polygon = Polygon(coords)
                x, y = polygon.exterior.xy
                plt.fill(x, y, alpha=0.20, facecolor='mediumpurple',
                         edgecolor='mediumpurple', linewidth=1.5)
            if p_lengths:
                plt.plot(p_lengths, p_means, color='mediumpurple', linestyle=propped_linestyle,
                         linewidth=1.5, label='_nolegend_')

    return data_max



# PLOT GEOMETRY OF SECTIONS
# ----------------------------------------------------------------------------------------------------------------------
def plot_section(section):
    # Create a figure and axis
    if section.section_type == "rc_rec":  # Rectangular Reinforced Concrete Cross-Section
        fig, ax, offset = plot_rectangle_with_dimensions(section.b, section.h, 'green', 'x')
        plot_rebars_long(ax, section, offset)
        # add stirrups to plot (if stirrups are defined)
        if section.bw_bg[0] > 0 and section.bw_bg[2] > 0:
            plot_stirrups(ax, section, offset)
        legend = (f'{section.concrete_type.mech_prop}, prod_ID:{section.concrete_type.prod_id} \n'
                  f'{section.rebar_type.mech_prop}, prod_ID:{section.rebar_type.prod_id} \n'
                  f'di_xo / s_xo = {section.bw[1][0]:.3f} / {section.bw[1][1]} \n'
                  f'di_xu / s_xu = {section.bw[0][0]:.3f} / {section.bw[0][1]} \n'
                  f'di_stir / s_stir / n = {section.bw_bg[0]} / {section.bw_bg[1]} / {section.bw_bg[2] }\n'
                  f'c_nom = {100*section.c_nom:.1f} cm \n'
                  f'x/d = {section.x_p/section.d:.2f} \n'
                  f'GWP = {section.co2:.0f} kg/m^2')
    elif section.section_type == "wd_rec":  # Rectangular Wooden Cross-Section
        fig, ax, offset = plot_rectangle_with_dimensions(section.b, section.h, 'brown', '/')
        legend = (f'{section.wood_type.mech_prop}, prod_ID:{section.wood_type.prod_id} \n'
                  f'GWP = {section.co2:.0f} kg/m^2')

    elif section.section_type == "rc_rib": #Betonrippenquerschnitte
        fig, ax, offset = plot_rib_with_dimensions(section.b, section.b_w, section.h, section.h_f, 'green', 'x')
        legend = (f'{section.concrete_type.mech_prop}, prod_ID:{section.concrete_type.prod_id} \n'
              f'length = {section.l0} \n'
              f'{section.rebar_type.mech_prop}, prod_ID:{section.rebar_type.prod_id} \n'
              f'di_r = {section.bw_r[0]:.3f} \n'
              f'di_xu / s_xu = {section.bw[0][0]:.3f} / {section.bw[0][1]} \n'
              #f'di_stir / s_stir / n = {section.bw_bg[0]} / {section.bw_bg[1]} / {section.bw_bg[2]}\n'
              f'c_nom = {100 * section.c_nom:.1f} cm \n'
              f'x/d = {section.x_p / section.d:.2f} \n'
              f'h, hf, hw, b, bw = {section.h:.2f}, {section.h_f:.2f}, {section.h_w:.2f}, {section.b:.2f}, {section.b_w:.2f} \n'
              f'GWP = {section.co2:.0f} kg/m^2')

    elif section.section_type == "wd_rib": #Betonrippenquerschnitte
        fig, ax, offset = plot_wd_rib_with_dimensions(section.b, section.h, section.a, section.t2, section.t3, 'brown', 'x')
        legend = (f'{section.wood_type_1.mech_prop}, prod_ID:{section.wood_type_1.prod_id} \n'
              f'length = {section.l0} \n'
              f'h, b, a, t2, t3 = {section.h:.2f}, {section.b:.2f}, {section.a:.2f}, {section.t2:.2f}, {section.t3:.2f} \n'
              f'GWP = {section.co2:.0f} kg/m^2')

    elif section.section_type == "comp_slab":
        fig, ax, offset = plot_comp_slab_with_dimensions(section)
        mesh_name = section.mesh.mech_prop if section.mesh else "none"
        bar_dia = getattr(section, 'bar_dia_trough', None)
        if bar_dia is not None:
            bar_str = (f'Trough bar: ø{bar_dia:.0f} mm  '
                       f'(A_s={section.A_s_trough:.0f} mm²/m)\n'
                       f'  trough bars: {section.co2_trough_bars:.1f}')
        else:
            bar_str = 'Trough bar: none'
        legend = (f'Deck: {section.deck.mech_prop}\n'
                  f'Concrete: {section.concrete.mech_prop}\n'
                  f'h = {section.h_mm:.0f} mm,  h_c = {section.h_c:.0f} mm\n'
                  f'Mesh: {mesh_name}\n'
                  f'{bar_str}\n'
                  f'GWP = {section.co2:.1f} kg-CO2-eq/m²\n'
                  f'  concrete: {section.co2_concrete:.1f}\n'
                  f'  deck: {section.co2_deck:.1f}\n'
                  f'  mesh: {section.co2_mesh:.1f}')

    else:
        print("no plot for specified section_type defined jet")
        fig, ax = plt.subplots()
        legend = f'no plot for section_type "{section.section_type}" defined jet'
    fig.text(0.01, 0.99, legend, ha='left', va='top', fontsize=9, color='black',
             bbox=dict(facecolor='lightgrey', edgecolor='black', boxstyle='round,pad=0.2'))


# Function to plot a rectangular cross-section with given dimensions
def plot_rectangle_with_dimensions(width, height, color='black', hatch='*', offset=0.1):
    # Create a figure and axis
    fig, ax = plt.subplots()

    # Define the rectangle with hatching (lower-left corner at (x, y), width, and height)
    rect = patches.Rectangle((offset, offset), width, height, linewidth=1, edgecolor=color, facecolor='none',
                             hatch=hatch, fill=False)

    # Add the rectangle to the plot
    ax.add_patch(rect)

    # Add dimension annotations
    ax.annotate(f'b = {width:.2f} m', xy=(offset + width / 2, 0.05), xytext=(offset + width / 2, 0.06), ha='center')
    ax.annotate(f'h = {height:.2f} m', xy=(0.02, offset + height / 2), xytext=(0.01, offset + height / 2),
                va='center', rotation='vertical')

    # Draw arrows for dimensions
    # ax.annotate('', xy=(0.1, 0.05), xytext=(0.1 + width, 0.05), arrowprops=dict(arrowstyle='|-|', color='black'))
    # ax.annotate('', xy=(0.05, 0.1), xytext=(0.05, 0.1 + height), arrowprops=dict(arrowstyle='|-|', color='black'))

    # Hide the x and y axes
    ax.axis('off')

    # Set the aspect of the plot to be equal
    ax.set_aspect('equal')

    # Set the limits of the plot
    ax.set_xlim(0, width+5*offset)
    ax.set_ylim(0, height+4*offset)

    return fig, ax, offset

def plot_rib_with_dimensions(b, bw, h, hf, color='black', hatch='*', offset=0.1):
    # Create a figure and axis
    fig, ax = plt.subplots()

    # Define the rectangle with hatching (lower-left corner at (x, y), width, and height)
    rect_flange = patches.Rectangle((offset, offset+h-hf), 2*b, hf, linewidth=1, edgecolor=color, facecolor='none',
                             hatch=hatch, fill=False)
    rect_rib1 = patches.Rectangle((offset+b/2-bw/2, offset), bw, h-hf, linewidth=1, edgecolor=color, facecolor='none',
                             hatch=hatch, fill=False)
    rect_rib2 = patches.Rectangle((offset + 3*b / 2 - bw / 2, offset), bw, h - hf, linewidth=1, edgecolor=color,
                                 facecolor='none',
                                 hatch=hatch, fill=False)

    # Add the rectangle to the plot
    ax.add_patch(rect_flange)
    ax.add_patch(rect_rib1)
    ax.add_patch(rect_rib2)

    # # Add dimension annotations
    # ax.annotate(f'b = {width:.2f} m', xy=(offset + width / 2, 0.05), xytext=(offset + width / 2, 0.06), ha='center')
    # ax.annotate(f'h = {height:.2f} m', xy=(0.02, offset + height / 2), xytext=(0.01, offset + height / 2),
    #             va='center', rotation='vertical')

    # Draw arrows for dimensions
    # ax.annotate('', xy=(0.1, 0.05), xytext=(0.1 + width, 0.05), arrowprops=dict(arrowstyle='|-|', color='black'))
    # ax.annotate('', xy=(0.05, 0.1), xytext=(0.05, 0.1 + height), arrowprops=dict(arrowstyle='|-|', color='black'))

    # Hide the x and y axes
    ax.axis('off')

    # Set the aspect of the plot to be equal
    ax.set_aspect('equal')

    # Set the limits of the plot
    ax.set_xlim(0, b + 10 * offset)
    ax.set_ylim(0, h + 4 * offset)

    return fig, ax, offset

def plot_wd_rib_with_dimensions(b, h, a, t2, t3, color='black', hatch='--', offset=0.1):
    # Create a figure and axis
    fig, ax = plt.subplots()

    # Define the rectangle with hatching (lower-left corner at (x, y), width, and height)
    rect_flange2 = patches.Rectangle((offset, offset), 2*a, t2, linewidth=1, edgecolor=color, facecolor='none',
                             hatch='--', fill=False)
    rect_flange3 = patches.Rectangle((offset, offset + t2+h), 2*a, t3, linewidth=1, edgecolor=color, facecolor='none',
                                     hatch='--', fill=False)
    rect_rib1 = patches.Rectangle((offset+a/2, offset+t2), b, h, linewidth=1, edgecolor=color, facecolor='none',
                             hatch='-', fill=False)
    rect_rib2 = patches.Rectangle((offset+3*a/2, offset+t2), b, h, linewidth=1, edgecolor=color, facecolor='none',
                             hatch='-', fill=False)


    # Add the rectangle to the plot
    ax.add_patch(rect_flange2)
    ax.add_patch(rect_flange3)
    ax.add_patch(rect_rib1)
    ax.add_patch(rect_rib2)

    # Add dimension annotations
    ax.annotate(f'b_Rippe = {b:.2f} m', xy=(offset + b, 0.05), xytext=(offset + b / 2, 0.06), ha='center')
    ax.annotate(f'b = {b:.2f} m', xy=(offset + b/2, 0.05), xytext=(offset + 3*b / 2, -0.16), ha='center')
    ax.annotate(f'h = {h:.2f} m', xy=(0.02, offset + h / 2), xytext=(0.01, offset + h / 2),
                va='center', rotation='vertical')

    # #Draw arrows for dimensions
    # ax.annotate('', xy=(0.1, 0.05), xytext=(0.1 + width, 0.05), arrowprops=dict(arrowstyle='|-|', color='black'))
    # ax.annotate('', xy=(0.05, 0.1), xytext=(0.05, 0.1 + height), arrowprops=dict(arrowstyle='|-|', color='black'))

    # Hide the x and y axes
    ax.axis('off')

    # Set the aspect of the plot to be equal
    ax.set_aspect('equal')

    # Set the limits of the plot
    ax.set_xlim(0, 2*a + 2 * offset)
    ax.set_ylim(0, h+t2+t3 + 4 * offset)

    return fig, ax, offset

def plot_comp_slab_with_dimensions(section, offset=0.005):
    """Plot a composite slab cross-section showing trapezoidal deck profile and concrete topping.
    All dimensions are converted from mm to m for plotting (consistent with other section plots)."""
    fig, ax = plt.subplots()

    # convert key dimensions to metres for plotting
    h    = section.h_mm / 1000.0   # total depth
    h_p  = section.deck.h_p / 1000.0   # deck profile height
    h_c  = section.h_c / 1000.0   # concrete above deck
    t    = section.deck.t / 1000.0   # sheet thickness (visual only)
    b    = 0.3  # width of section to draw (show ~300 mm strip)

    # ── concrete topping (rectangle above deck) ─────────────────────────────
    rect_conc = patches.Rectangle((offset, offset + h_p), b, h_c,
                                  linewidth=1, edgecolor='grey', facecolor='lightgrey',
                                  hatch='..', label='Concrete')
    ax.add_patch(rect_conc)

    # ── trapezoidal deck profile ─────────────────────────────────────────────
    # draw a simplified repeating trapezoidal profile across the width
    n_ribs = max(int(b / 0.15), 1)  # approximate number of ribs in drawn width
    rib_pitch = b / n_ribs
    top_frac = 0.55    # fraction of pitch that is the top flange
    bot_frac = 0.35    # fraction of pitch that is the bottom trough

    for i in range(n_ribs):
        x0 = offset + i * rib_pitch
        # trapezoidal shape: top-left, top-right, bottom-right, bottom-left
        margin = rib_pitch * (1 - top_frac) / 2
        top_l = x0 + margin
        top_r = x0 + rib_pitch - margin
        bot_margin = rib_pitch * (1 - bot_frac) / 2
        bot_l = x0 + bot_margin
        bot_r = x0 + rib_pitch - bot_margin

        trap = patches.Polygon(
            [(top_l, offset + h_p),
             (top_r, offset + h_p),
             (bot_r, offset),
             (bot_l, offset)],
            closed=True, linewidth=1, edgecolor='steelblue', facecolor='lightsteelblue',
            hatch='//')
        ax.add_patch(trap)

    # ── mesh reinforcement dots ──────────────────────────────────────────────
    if section.mesh is not None:
        mesh_y = offset + h_p + h_c / 2   # mid-height of concrete topping
        mesh_d = section.mesh.d_s / 1000.0 if section.mesh.d_s else 0.005
        # space dots across width
        n_dots = max(int(b / 0.05), 3)
        for j in range(n_dots):
            cx = offset + (j + 0.5) * b / n_dots
            dot = plt.Circle((cx, mesh_y), mesh_d / 2, color='blue')
            ax.add_patch(dot)

    # ── trough reinforcement bars (fire design) ──────────────────────────────
    if getattr(section, 'bar_dia_trough', None) is not None:
        bar_y = offset + section.d_bar_trough / 1000.0   # mm → m
        bar_r = (section.bar_dia_trough / 1000.0) / 2    # radius in m
        for i in range(n_ribs):
            cx = offset + i * rib_pitch + rib_pitch / 2  # centre of each trough gap
            dot = plt.Circle((cx, bar_y), bar_r, color='orangered', zorder=5)
            ax.add_patch(dot)

    # ── dimension annotations ────────────────────────────────────────────────
    dim_x = offset + b + offset * 2
    ax.annotate('', xy=(dim_x, offset), xytext=(dim_x, offset + h),
                arrowprops=dict(arrowstyle='<->', color='black', lw=1.2))
    ax.text(dim_x + offset, offset + h / 2, f'h={section.h_mm:.0f}',
            va='center', fontsize=9)

    dim_x2 = dim_x + offset * 8
    ax.annotate('', xy=(dim_x2, offset + h_p), xytext=(dim_x2, offset + h),
                arrowprops=dict(arrowstyle='<->', color='black', lw=1.2))
    ax.text(dim_x2 + offset, offset + h_p + h_c / 2, f'h_c={section.h_c:.0f}',
            va='center', fontsize=9)

    ax.annotate('', xy=(dim_x2, offset), xytext=(dim_x2, offset + h_p),
                arrowprops=dict(arrowstyle='<->', color='black', lw=1.2))
    ax.text(dim_x2 + offset, offset + h_p / 2, f'h_p={section.deck.h_p:.0f}',
            va='center', fontsize=9)

    # axes
    ax.axis('off')
    ax.set_aspect('equal')
    ax.set_xlim(0, b + 15 * offset)
    ax.set_ylim(-offset, h + 4 * offset)

    return fig, ax, offset


def plot_rebars_long(ax, section, offset, color='blue'):
    # get rebar positions
    rebar_positions = get_rebar_positions(section)
    # plot rebars
    for (x, y, di) in rebar_positions:
        rebar = plt.Circle((x+offset, y+offset), di/2, color=color)
        ax.add_patch(rebar)

def plot_stirrups(ax, section, offset, color='blue'):
    # get stirrup positions
    stirrup_positions = get_stirrup_positions(section)
    # plot stirrups
    for (x1, y1, x2, y2) in stirrup_positions:
        stirrup = plt.Line2D([x1+offset, x2+offset], [y1+offset, y2+offset], color=color, linewidth=2)
        ax.add_line(stirrup)

def get_rebar_positions(section):
    # create x and y coordinates of lower longitudinal reinforcement
    y_u = section.h - section.d
    s_xu = section.bw[0][1]
    x_u = [section.b/2]
    while max(x_u) + s_xu < section.b:
        x_u.append(max(x_u) + s_xu)
        x_u.append(min(x_u) - s_xu)

    # create x and y coordinates of upper longitudinal reinforcement
    y_o = section.ds
    s_xo = section.bw[1][1]
    x_o = [section.b/2]
    while max(x_o) + s_xo < section.b:
        x_o.append(max(x_o) + s_xo)
        x_o.append(min(x_o) - s_xo)

    # assemble rebar positions and dimensions relative to cross-section left lower edge
    di_xu = section.bw[0][0]
    di_xo = section.bw[1][0]
    rebar_positions = []
    for xi in x_u:
        rebar_i = (xi, y_u, di_xu)
        rebar_positions.append(rebar_i)
    for xi in x_o:
        rebar_i = (xi, y_o, di_xo)
        rebar_positions.append(rebar_i)
    return rebar_positions


def get_stirrup_positions(section):
    # create coordinates of stirrup edge points
    n_stirrup = section.bw_bg[2]
    di_stirrup = section.bw_bg[0]
    edge_dist = section.c_nom + di_stirrup/2
    # create x-coordinates of vertical stirrups
    x_linspace = np.linspace(edge_dist, section.b-edge_dist, n_stirrup)
    # create y-coordinates of stirrups-endings
    y_u = edge_dist
    y_o = section.h - edge_dist
    # create lines between edge points
    stirrup_positions = []
    for idx, xi in enumerate(x_linspace):
        # create vertical lines
        line_vert = (xi, y_u, xi, y_o)
        stirrup_positions.append(line_vert)
        if idx < x_linspace.size-1:
            # create horizontal lines
            line_u = (xi, y_u, x_linspace[idx+1], y_u)
            stirrup_positions.append(line_u)
            line_o = (xi, y_o, x_linspace[idx+1], y_o)
            stirrup_positions.append(line_o)
    return stirrup_positions


# PLOT DATASETS OF CROSS_SECTION WITH VARIED MATERIALS IN M-CHI PLOT AND PLOT THE OPTIMIZED SECTIONS FOR VALIDATION
# ----------------------------------------------------------------------------------------------------------------------
def plot_section_dataset(database_name, crsec_type, mat_names, ax, gwp_budget=50):
    # GENERATE INITIAL CROSS-SECTIONS
    # Search database (table products, attribute material) for products,
    # get prod_id of relevant materials from database and create initial cross-section for each product
    to_plot = []
    connection = sqlite3.connect(database_name)
    cursor = connection.cursor()
    mat_nr = -1
    x_values = []
    y_values = []
    for mat_name in mat_names:
        inquiry = ("SELECT PRO_ID FROM products WHERE"
                   " MATERIAL=" + mat_name)
        cursor.execute(inquiry)
        result = cursor.fetchall()
        for i, prod_id in enumerate(result):
            mat_nr += 1  # number for annotation in plot
            # create materials for wooden cross-sections, derive corresponding design values
            prod_id_str = "'" + str(prod_id[0]) + "'"
            inquiry = ("SELECT MECH_PROP FROM products WHERE"
                       " PRO_ID=" + prod_id_str)
            cursor.execute(inquiry)
            result = cursor.fetchall()
            mech_prop = "'" + result[0][0] + "'"
            if crsec_type == "wd_rec":
                # create a Wood material object
                timber = struct_analysis.Wood(mech_prop, database_name, prod_id_str)
                timber.get_design_values()
                # create initial wooden rectangular cross-section
                section_0 = struct_analysis.RectangularWood(timber, 1.0, 0.12)
                color = "tab:brown"
            elif crsec_type == "rc_rec":
                # create a Concrete material object
                concrete = struct_analysis.ReadyMixedConcrete(mech_prop, database_name, prod_id=prod_id_str)
                concrete.get_design_values()
                # create a Rebar material object
                rebar = struct_analysis.SteelReinforcingBar("'B500B'", database_name)
                # create initial wooden rectangular cross-section
                section_0 = struct_analysis.RectangularConcrete(concrete, rebar, 1.0, 0.12,
                                                                0.014, 0.15, 0.01, 0.15,
                                                                0, 0.15, 2)
                color = "tab:green"
            elif crsec_type == "rc_rib":
                # create a Concrete material object
                concrete = struct_analysis.ReadyMixedConcrete(mech_prop, database_name, prod_id=prod_id_str)
                concrete.get_design_values()
                # create a Rebar material object
                rebar = struct_analysis.SteelReinforcingBar("'B500B'", database_name)
                # create initial wooden rectangular cross-section
                section_0 = struct_analysis.RibbedConcrete(concrete, rebar, 4, 1.0, 0.14, 0.3, 0.18, 0.01, 0.15, 0.01, 0.15, 0.02, 2, 0.01, 0.15, 2)
                color = 'mediumseagreen'
## XXXXXXXXXXX neuen Querschnittstyp initialisieren
            else:
                print("cross-section type is not defined inside function plot_dataset()")
                section_0 = []
                color = "tab:grey"


            #
            # maximizing Mu by varying the geometry within the max allowed gwp-budget.
            opt_section = struct_optimization.get_opt_sec(section_0, gwp_budget)

            # add M-Chi relationship to plot
            x_values, y_values = plot_m_chi(opt_section, ax, mat_nr, x_values, y_values)

            # plot cross-section
            plot_section(opt_section)
            plt.title(f'#{mat_nr}')

    # ## Add envelope of dataset to plot
    # # Combine x and y into a single array for ConvexHull
    # x = [item for sublist in x_values for item in sublist]
    # y = [item for sublist in y_values for item in sublist]
    # points = np.column_stack((x, y))
    # # points = [item for sublist in points_nested for item in sublist]
    # print(points.shape)
    # hull = ConvexHull(points)
    #
    # # Plot the envelope area (convex hull)
    # for simplex in hull.simplices:
    #     ax.plot(points[simplex, 0], points[simplex, 1], 'r-')
#
#     # Fill the hull area
#     hull_vertices = points[hull.vertices]
#     ax.fill(hull_vertices[:, 0], hull_vertices[:, 1], color, alpha=0.3, label="Envelope")

    ## Add envelope of dataset to plot
    # Define a common x-axis for interpolation
    common_x = sorted(set(np.concatenate(x_values)))
    # Interpolate y-values onto the common x-axis
    interpolated_y_values_min = []
    interpolated_y_values_max = []
    # y_test = []
    for x, y in zip(x_values, y_values):
        # define interpolation functions with different rules, when x is out of range
        f1 = interp1d(x, y, kind='linear', bounds_error=False, fill_value="extrapolate")
        f2 = interp1d(x, y, kind='linear', bounds_error=False, fill_value=(0, 0))
        a = []
        b = []
        for xi in common_x:
            if xi >= 0:  # apply interpolation functions withe rules for positive x_values
                a.append(float(f1(xi)))
                b.append(float(f2(xi)))
            else:  # apply interpolation functions withe rules for negative x_values
                b.append(float(f1(xi)))
                a.append(float(f2(xi)))
        interpolated_y_values_min.append(a)
        interpolated_y_values_max.append(b)

    y_lower = np.min(interpolated_y_values_min, axis=0) if interpolated_y_values_min else np.array([])
    y_upper = np.max(interpolated_y_values_max, axis=0) if interpolated_y_values_max else np.array([])

    # Fill the area between the envelope
    ax.fill_between(common_x, y_lower, y_upper, color=color, alpha=0.3, label='Envelope')



# plot m-chi relationship for a defined cross-section
def plot_m_chi(section, ax, i, x_values, y_values):
    # add M-Chi relationship to plot
    if section.section_type[0:2] == "wd":
        color = "tab:brown"  # color for wood
        chi_u_p, chi_u_n = section.mu_max/section.ei1, section.mu_min/section.ei1
        x = [chi_u_n, 0, chi_u_p]
        y = [section.mu_min, 0, section.mu_max]
    elif section.section_type[0:2] == "rc":
        color = "tab:green"  # color for reinforced concrete
        chi_r1_p, chi_r1_n = section.mr_p / section.ei1, section.mr_n / section.ei1
        chi_r2_p, chi_r2_n = section.mr_p / section.ei2, section.mr_n / section.ei2
        chi_y_p, chi_y_n = section.mu_max / section.ei2, section.mu_min / section.ei2
        chi_u_p, chi_u_n = section.concrete_type.ec2d/section.x_p, -section.concrete_type.ec2d/section.x_n
        if section.mr_p <= section.mu_max and section.mr_n >= section.mu_min:  # ductile behavior in both directions
            x = [chi_u_n, chi_y_n, chi_r2_n, chi_r1_n, 0, chi_r1_p, chi_r2_p, chi_y_p, chi_u_p]
            y = [section.mu_min, section.mu_min, section.mr_n, section.mr_n, 0, section.mr_p, section.mr_p,
                 section.mu_max, section.mu_max]
        elif section.mr_p <= section.mu_max and section.mr_n < section.mu_min:  # ductile behavior only for pos. moments
            x = [chi_r1_n, chi_r1_n, 0, chi_r1_p, chi_r2_p, chi_y_p, chi_u_p]
            y = [0, section.mr_n, 0, section.mr_p, section.mr_p,
                 section.mu_max, section.mu_max]
        elif section.mr_p > section.mu_max and section.mr_n >= section.mu_min:  # ductile behavior only for neg. moments
            x = [chi_u_n, chi_y_n, chi_r2_n, chi_r1_n, 0, chi_r1_p, chi_r1_p]
            y = [section.mu_min, section.mu_min, section.mr_n, section.mr_n, 0, section.mr_p, 0]
        else:  # no ductile behavior
            x = [chi_r1_n, chi_r1_n, 0, chi_r1_p, chi_r1_p]
            y = [0, section.mr_n, 0, section.mr_p, 0]
    else:
        print("M-Chi plot is not defined for sections of type " + section.section_type)
        x, y = 0, 0
        color = "tab:black"
    y_mod = [yi/1e3 for yi in y]  # modify unit of moment from Nm/m to kNm/m
    ax.plot(x, y_mod, color=color)  # plot m-Chi relationship with unit of moment: kNm/m
    ax.annotate(f'#{i}', xy=(x[-1], y_mod[-1]), xytext=(x[-1]*1.1, y_mod[-1]),
                 arrowprops=dict(facecolor='black', shrink=0.2, width=0.2, headwidth=2, headlength=4),
                 fontsize=9, color='black', va='center')
    x_values.append(x)
    y_values.append(y_mod)
    return x_values, y_values
    ## Define and plot envelope


# plt.annotate(f'#{i}', xy=(ver_x, ver_y),
#                          xytext=(ver_x + 0.05*lengths[-1], ver_y),
#                          arrowprops=dict(facecolor='black', shrink=0.2, width=0.2, headwidth=2, headlength=4),
#                          fontsize=9, color='black', va='center')