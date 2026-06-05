# File enthält Code für die Strukturanalyse (Bauteil- und Querschnittsanalyse)
# units: [m], [kg], [s], [N], [CHF]

# Abgebildete Materialien:
# - Beton
# - Betonstahl
# - Holz
#
# Abgebildete Querschnitte 1D:
# - Betonrechteck-QS
# - Holzrechteck-QS
# - Beton-Rippen-QS
# - Holz-Hohlkasten-QS
#
# Abgebildete Statische Systeme 1D:
# - Einfacher Balken
# - Durchlaufträger (in Bearbeitung)
#
#Statische Systeme 2D:
# - Platte 4-seitig gelagert für vordefinierte Spannweiten
#
# Weitere Klassen:
# - Bauteil 1D
# - Bodenaufbauschicht
# - Bodenaufbau
# - Rechteckquerschnitte
# - Anforderungen

import sqlite3  # import modul for SQLite
import math
import numpy as np
from scipy.optimize import minimize

#DEFINITONS OF MATERIAL PROPERTIES--------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
class Decking:
    def __init__(self, mech_prop, database, prod_id="undef"):
        self.mech_prop = mech_prop
        connection = sqlite3.connect(database)
        cursor = connection.cursor()
        # get mechanical properties from database
        inquiry = ("SELECT DECK_ID, t, h_p, weight_kg, weight_kN, A_p, A_pe, e, e_p, theta, "
                   "s_w, num_webs, I_p, M_cRd, M_tRd, f_ypd, m_p, k_p, tau_uRk, W_pl, W_elf, "
                   "GWP_UNIT, GWP, conc_vol_trough, l_1, l_2, l_3 "
                   "FROM sheeting_prop WHERE Deck=?")
        cursor.execute(inquiry, (mech_prop,))
        result = cursor.fetchall()
        self.E = 2.1E+11 #N/m2
        (self.DECK_ID, self.t, self.h_p, self.weight_kg, self.weight_kN, self.A_p, self.A_pe,
         self.e, self.e_p, self.theta, self.s_w, self.num_webs, self.I_p, self.M_cRd, self.M_tRd,
         self.f_ypd, self.m_p, self.k_p, self.tau_uRk, self.W_pl, self.W_elf,
         self.GWP_UNIT, self.GWP, self.conc_vol_trough,
         self.l_1, self.l_2, self.l_3) = result[0]
        # Rib profile type for EN 1994-1-2 Annex D (re-entrant: bottom wider than top)
        self.profile_type = "reentrant" if (self.l_2 > self.l_3) else "trapezoidal"
        connection.close()

class SteelReinforcingMesh:
    def __init__(self, mech_prop, database, prod_id="undef"):
        self.mech_prop = mech_prop
        connection = sqlite3.connect(database)
        cursor = connection.cursor()

        inquiry = "SELECT strength_tens, E_modulus, area_per_metre, bar_size FROM material_prop WHERE name=?"
        cursor.execute(inquiry, (mech_prop,))
        result = cursor.fetchall()
        self.fsk, self.Es, self.A_s, self.d_s = result[0]

        result = cursor.fetchall()
        self.fsd     = self.get_design_values()
        connection.close()

    def get_design_values(self, gamma_s=1.15):
        return self.fsk / gamma_s

class Wood:
    # defines properties of wooden material
    def __init__(self, mech_prop, database, prod_id="undef"):  # retrieve basic mechanical data from database
        self.mech_prop = mech_prop
        connection = sqlite3.connect(database)
        cursor = connection.cursor()
        # get mechanical properties from database
        inquiry = ("SELECT strength_bend, strength_shea, E_modulus, density_load, burn_rate FROM material_prop WHERE"
                   " name=" + mech_prop)
        cursor.execute(inquiry)
        result = cursor.fetchall()
        self.fmk, self.fvd, self.Emmean, self.weight, self.burn_rate = result[0]
        # get GWP properties from database
        if prod_id == "undef":  # no specific product is defined, chose first product entry with required mechanical
            # properties in database
            inquiry = "SELECT PRO_ID, DENSITY, Total_GWP, Cost FROM products WHERE MECH_PROP=" + mech_prop
        else:
            inquiry = "SELECT PRO_ID, DENSITY, Total_GWP, Cost FROM products WHERE PRO_ID=" + prod_id
        cursor.execute(inquiry)
        result = cursor.fetchall()
        # self.prod_id, self.density, self.GWP, self.cost, self.cost2 = result[0]
        self.prod_id, density, GWP, self.cost = result[0]
        self.GWP = GWP/1e3  # transform unit from [kg-Co2-eq/t] to [kg-Co2-eq/kg]
        self.density = float(density)
        self.cost = 0
        self.cost2 = 0
        self.fmd = self.get_design_values()

    def get_design_values(self, gamma_m=1.7, eta_m=1, eta_t=1, eta_w=1):  # calculate design values
        if self.mech_prop[1:3] == "GL":
            gamma_m = 1.5  # SIA 265, 2.2.5: reduzierter Sicherheitsbeiwert für BSH
        fmd = self.fmk * eta_m * eta_t * eta_w / gamma_m  # SIA 265, 2.2.2, Formel (3)
        return fmd


class ReadyMixedConcrete:
    # defines properties of concrete material
    def __init__(self, mech_prop, database, dmax=32,
                 prod_id="undef"):  # retrieve basic mechanical data from database (self, table,
        self.mech_prop = mech_prop
        connection = sqlite3.connect(database)
        cursor = connection.cursor()
        # get mechanical properties from database
        inquiry = ("SELECT strength_comp, strength_tens, E_modulus, density_load, dry_density_load, wet_density_load FROM material_prop WHERE name="
                   + mech_prop)
        cursor.execute(inquiry)
        result = cursor.fetchall()
        self.fck, self.fctm, self.Ecm, self.weight, self.dry_weight, self.wet_weight = result[0]
        # get GWP properties from database
        if prod_id == "undef":  # no specific product is defined, chose first product entry with required mechanical
            # properties in database
            # inquiry = ("""
            #         SELECT PRO_ID, density, Total_GWP, cost, cost2 FROM products WHERE "material [string]" LIKE """ + mech_prop
            #            )
            inquiry = ("""SELECT PRO_ID, DENSITY, Total_GWP, Cost FROM products WHERE MECH_PROP LIKE """
                       + mech_prop)
        else:
            # inquiry = ("""SELECT PRO_ID, density, Total_GWP, cost, cost2 FROM products WHERE PRO_ID LIKE """ + prod_id
            #            )
            inquiry = ("""SELECT PRO_ID, DENSITY, Total_GWP, Cost FROM products WHERE PRO_ID LIKE """
                       + prod_id)
        cursor.execute(inquiry)
        result = cursor.fetchall()
        # self.prod_id, self.density, self.GWP, self.cost, self.cost2 = result[0]
        self.prod_id, density, GWP, self.cost = result[0]
        self.GWP = GWP/1e3  # transform unit from [kg-Co2-eq/t] to [kg-Co2-eq/kg]
        self.density = float(density)
        self.cost = 0
        self.cost2 = 0
        self.dmax = dmax
        self.fcd, self.tcd, self.ec2d = self.get_design_values()

    def get_design_values(self, gamma_c=1.5, eta_t=1):  # calculate design values
        eta_fc = min((30e6 / self.fck) ** (1 / 3), 1)  # SIA 262, 4.2.1.2, Formel (26)
        fcd = self.fck * eta_fc * eta_t / gamma_c  # SIA 262, 2.3.2.3, Formel (2)
        tcd = 0.3 * eta_t * 1e6 * (self.fck * 1e-6) ** 0.5 / gamma_c  # SIA 262, 2.3.2.4, Formel (3)
        ec2d = 0.003  # SIA 262, 4.2.4, Tabelle 8
        return fcd, tcd, ec2d

class SteelReinforcingBar:
    # defines properties of reinforcement  material
    def __init__(self, mech_prop, database, prod_id="undef"):
        # retrieve basic mechanical data from database (self, table, database name)
        self.mech_prop = mech_prop
        connection = sqlite3.connect(database)
        cursor = connection.cursor()
        # get mechanical properties from database
        inquiry = "SELECT strength_tens, E_modulus FROM material_prop WHERE name=" + mech_prop
        cursor.execute(inquiry)
        result = cursor.fetchall()
        self.fsk, self.Es = result[0]
        # get GWP properties from database
        if prod_id == "undef":  # no specific product is defined, chose first product entry with required mechanical
            # properties in database
            inquiry = "SELECT PRO_ID, DENSITY, Total_GWP, Cost FROM products WHERE MECH_PROP=" + mech_prop
        else:
            inquiry = "SELECT PRO_ID, DENSITY, Total_GWP, Cost FROM products WHERE PRO_ID=" + prod_id
        cursor.execute(inquiry)
        result = cursor.fetchall()
        #self.prod_id, density, self.GWP, self.cost = result[0]
        self.prod_id, density, GWP, self.cost = result[0]
        self.GWP = GWP/1e3  # transform unit from [kg-Co2-eq/t] to [kg-Co2-eq/kg]
        self.density = float(density)
        self.cost = 0
        self.fsd = self.get_design_values()

    def get_design_values(self, gamma_s=1.15):  # calculate design values
        fsd = self.fsk / gamma_s  # SIA 262, 2.3.2.5, Formel (4)
        return fsd


#-----------------------------------------------------------------------------------------------------------------------
# COMPOSITE SLAB CROSS-SECTION (EN 1994-1-1)
#-----------------------------------------------------------------------------------------------------------------------
class CompositeSlab:
    """Composite steel-concrete slab with profiled sheeting.
    Design checks per EN 1994-1-1 (simply supported, unpropped).
    """
    MESH_TYPES    = ["A142", "A193", "A252", "A393"]  # lightest to heaviest
    STEEL_DENSITY = 7850.0  # kg/m³

    # ── EN 1994-1-2 Annex D fire tables ──────────────────────────────────────

    # EN 1994-1-2 Table 3.2  k_y,θ = f_ay,θ / f_ay for structural steel (deck)
    _KY_THETA_DECK = [
        (  20, 1.00), ( 100, 1.00), ( 200, 1.00), ( 300, 1.00), ( 400, 1.00),
        ( 500, 0.78), ( 600, 0.47), ( 700, 0.23), ( 800, 0.11), ( 900, 0.06),
        (1000, 0.04), (1100, 0.02), (1200, 0.00),
    ]

    # EN 1994-1-2 Table 3.2  k_y,θ for hot-rolled reinforcing bars (B500B)
    # EN 1994-1-2 Cl 3.2.3: hot-rolled bars use Table 3.2 k_y,θ (same as structural steel)
    _KY_THETA_BAR = [
        (  20, 1.00), ( 100, 1.00), ( 200, 1.00), ( 300, 1.00), ( 400, 1.00),
        ( 500, 0.78), ( 600, 0.47), ( 700, 0.23), ( 800, 0.11), ( 900, 0.06),
        (1000, 0.04), (1100, 0.02), (1200, 0.00),
    ]

    # Table D.1  Coefficients for thermal insulation formula D.1
    # t_i = a0 + a1*h1 + a2*Phi + a3*(A/Lr) + a4*(1/l3) + a5*(A/Lr)*(1/l3)
    # h1 = h_c = concrete thickness above the ribs [mm]
    _FIRE_D1_COEFF = {
        "NWC": (-28.8,  1.55, -12.6,  0.33, -735, 48.0),
        "LWC": (-79.2,  2.18,  -2.44, 0.56, -542, 52.3),
    }

    # Table D.2  Coefficients for deck temperature formula D.4
    # theta_a = b0 + b1*(1/l3) + b2*(A/Lr) + b3*Phi + b4*Phi^2
    _FIRE_D2_COEFF = {
        "NWC": {
            60:  {"lower": ( 951, -1197, -2.32,   86.4, -150.7),
                  "web":   ( 661,  -833, -2.96,  537.7, -351.9),
                  "upper": ( 340, -3269, -2.62, 1148.4, -679.8)},
            90:  {"lower": (1018,  -839, -1.55,   65.1, -108.1),
                  "web":   ( 816,  -959, -2.21,  464.9, -340.2),
                  "upper": ( 618, -2786, -1.79,  767.9, -472.0)},
            120: {"lower": (1063,  -679, -1.13,   46.7,  -82.8),
                  "web":   ( 925,  -949, -1.82,  344.2, -267.4),
                  "upper": ( 770, -2460, -1.67,  592.6, -379.0)},
        },
        "LWC": {
            30:  {"lower": ( 800, -1326, -2.65,  114.5, -181.2),
                  "web":   ( 483,  -286, -2.26,  439.6, -244.0),
                  "upper": ( 331, -2284, -1.54,  488.8, -131.7)},
            60:  {"lower": ( 955,  -622, -1.32,   47.7,  -81.1),
                  "web":   ( 761,  -558, -1.67,  426.5, -303.0),
                  "upper": ( 607, -2261, -1.02,  664.5, -410.0)},
            90:  {"lower": (1019,  -478, -0.91,   32.7,  -60.8),
                  "web":   ( 906,  -654, -1.36,  287.8, -230.3),
                  "upper": ( 789, -1847, -0.99,  469.5, -313.0)},
            120: {"lower": (1062,  -399, -0.65,   19.8,  -43.7),
                  "web":   ( 989,  -629, -1.07,  186.1, -152.6),
                  "upper": ( 903, -1561, -0.92,  305.2, -197.2)},
        },
    }

    # Table D.4  Coefficients for determination of the limiting temperature (formula D.7)
    # theta_lim = d0 + d1*N_s + d2*(A/Lr) + d3*Phi + d4*(1/l3)
    # N_s = normal force in hogging reinforcement [N]
    _FIRE_D4_COEFF = {
        "NWC": {
             60: ( 867, -1.9e-4, -8.75, -123, -1378),
             90: (1055, -2.2e-4, -9.91, -154, -1990),
            120: (1144, -2.2e-4, -9.71, -166, -2155),
        },
        "LWC": {
             30: ( 524, -1.6e-4, -3.43,  -80,  -392),
             60: (1030, -2.6e-4,-10.95, -181, -1834),
             90: (1159, -2.5e-4,-10.88, -208, -2233),
            120: (1213, -2.5e-4,-10.09, -214, -2320),
        },
    }

    # Table D.6  Minimum effective thickness h_eff [mm] as a function of fire resistance
    # Required: h_eff >= h_eff_min - h3, where h3 = screed thickness on top (mm)
    _FIRE_D6_MIN_HEFF = {
        30: 60, 60: 80, 90: 100, 120: 120, 180: 150, 240: 175,
    }

    # Table D.7  Field of application geometry limits (mm)
    # h_p = h1 (rib height), h_c = h2 (concrete above ribs)
    # Keys match the 'Type' field in the database ('reentrant' or 'trapezoidal')
    _FIRE_D7_LIMITS = {
        "reentrant":   {"l1": ( 77.0, 135.0), "l2": (110.0, 150.0),
                        "l3": ( 38.5,  97.5), "h_p": ( 50.0, 130.0), "h_c": ( 30.0,  60.0)},
        "trapezoidal": {"l1": ( 80.0, 155.0), "l2": ( 32.0, 132.0),
                        "l3": ( 40.0, 115.0), "h_p": ( 50.0, 125.0), "h_c": ( 50.0, 100.0)},
    }

    # Table D.3  Coefficients for rebar-in-rib temperature formula D.5
    # theta_s = c0 + c1*(u3/h1) + c2*sqrt(z) + c3*(A/Lr) + c4*alpha + c5*(1/l3)
    # u3  = distance from bar centre to lower flange [mm]
    # h1  = rib height h_p [mm]
    # z   = position factor from formula D.6: 1/z = 1/sqrt(u1) + 1/sqrt(u2) + 1/sqrt(u3)
    #       u1, u2 = shortest distances from bar centre to each web [mm]
    # alpha = web inclination angle [degrees]
    # A/Lr [mm], l3 [mm] → theta_s [°C]
    _FIRE_D3_COEFF = {
        "NWC": {
            60:  (1191, -250, -240, -5.01, 1.04, -925),
            90:  (1342, -256, -235, -5.30, 1.39, -1267),
            120: (1387, -238, -227, -4.79, 1.68, -1326),
        },
        "LWC": {
            30:  ( 809, -135, -243, -0.70, 0.48,  -315),
            60:  (1336, -242, -292, -6.11, 1.63,  -900),
            90:  (1381, -240, -269, -5.46, 2.24,  -918),
            120: (1397, -230, -253, -4.44, 2.47,  -906),
        },
    }

    # ── Annex D static helpers ────────────────────────────────────────────────

    @staticmethod
    def _fire_ky_theta(theta, table):
        """Linear interpolation of a strength-reduction table [(T, k), ...]."""
        if theta <= table[0][0]:  return table[0][1]
        if theta >= table[-1][0]: return table[-1][1]
        for i in range(len(table) - 1):
            t0, k0 = table[i];  t1, k1 = table[i + 1]
            if t0 <= theta <= t1:
                return k0 + (k1 - k0) * (theta - t0) / (t1 - t0)
        return 0.0

    @staticmethod
    def _fire_h_eff(h_c, h_p, l1, l2, l3):
        """Effective slab thickness h_eff [mm] — EN 1994-1-2 formulae D.15a/b."""
        if l3 > 2.0 * l1:
            return h_c
        ratio = (l1 + l2) / (l1 + l3)
        if h_p / h_c <= 1.5:
            return h_c + 0.5 * h_p * ratio   # D.15a
        else:
            return h_c * (1.0 + 0.75 * ratio)  # D.15b

    @staticmethod
    def _fire_ALr(h_p, l2, l3):
        """Rib geometry factor A/Lr [mm] — formula D.2."""
        A  = h_p * (l2 + l3) / 2.0
        Lr = l2 + 2.0 * math.sqrt(h_p**2 + ((l2 - l3) / 2.0)**2)
        return A / Lr

    @staticmethod
    def _fire_view_factor(h_p, l1, l2, l3):
        """Upper-flange view factor Phi [-] — formula D.3."""
        d1  = (l1 - l2) / 2.0
        phi = (math.sqrt(h_p**2 + (l3 + d1)**2) - math.sqrt(h_p**2 + d1**2)) / l3
        return max(0.0, min(1.0, phi))

    @staticmethod
    def _fire_t_i(h1, phi, ALr, l3, concrete_type="NWC"):
        """Insulation fire resistance time t_i [min] — formula D.1.
        h1 = h_c = concrete thickness above the ribs [mm]."""
        a = CompositeSlab._FIRE_D1_COEFF.get(concrete_type, CompositeSlab._FIRE_D1_COEFF["NWC"])
        a0, a1, a2, a3, a4, a5 = a
        inv_l3 = 1.0 / l3
        return a0 + a1*h1 + a2*phi + a3*ALr + a4*inv_l3 + a5*ALr*inv_l3

    @staticmethod
    def _fire_deck_temperature(l3, ALr, phi, R_fi, part, concrete_type="NWC"):
        """Deck steel temperature [°C] — formula D.4.  Returns None if no data."""
        ctype = CompositeSlab._FIRE_D2_COEFF.get(concrete_type,
                                                  CompositeSlab._FIRE_D2_COEFF["NWC"])
        row = ctype.get(R_fi)
        if row is None:
            return None
        b0, b1, b2, b3, b4 = row[part]
        return b0 + b1*(1.0/l3) + b2*ALr + b3*phi + b4*phi**2

    @staticmethod
    def _fire_z_factor(u1, u2, u3):
        u1 = max(u1, 5.0)
        u2 = max(u2, 5.0)
        u3 = max(u3, 5.0)
        return 1.0 / (1.0 / math.sqrt(u1) + 1.0 / math.sqrt(u2) + 1.0 / math.sqrt(u3))

    @staticmethod
    def _fire_rebar_temperature(u3, h_p, z, ALr, alpha, l3, R_fi, concrete_type="NWC"):
        """Temperature of additional reinforcement bar in rib [°C] — formula D.5.

        theta_s = c0 + c1*(u3/h_p) + c2*sqrt(z) + c3*(A/Lr) + c4*alpha + c5*(1/l3)

        u3    = distance from bar centre to lower flange [mm]
        h_p   = rib height [mm]
        z     = position factor from _fire_z_factor() [mm^0.5]
        ALr   = rib geometry factor A/Lr [mm]
        alpha = web inclination angle [degrees]
        l3    = upper flange width [mm]
        Returns None if no coefficients are available for the given R_fi.
        """
        ctype = CompositeSlab._FIRE_D3_COEFF.get(concrete_type,
                                                  CompositeSlab._FIRE_D3_COEFF["NWC"])
        row = ctype.get(R_fi)
        if row is None:
            return None
        c0, c1, c2, c3, c4, c5 = row
        return c0 + c1*(u3/h_p) + c2*math.sqrt(z) + c3*ALr + c4*alpha + c5*(1.0/l3)

    def __init__(self, deck, concrete, h_mm, database=None,
                 imposed_load=2000.0, finishes_load=1000.0, ceiling_services=750.0,
                 construction_load=750.0, partition_load=0.0,
                 gamma_c=1.5, gamma_g=1.35, gamma_q=1.5, gamma_m0=1.0,
                 RH=50.0, t0=28.0, cement_class="N",
                 R_fi=60, cover_top=20.0, psi_fi=0.5,
                 n_spans=1, propped=False):
        self.section_type = "comp_slab"
        self.deck = deck
        self.concrete = concrete
        self.database = database
        self.h = h_mm / 1000.0       # m (consistent with other section types)
        self.h_mm = h_mm
        self.h_c = max(h_mm - deck.h_p, 1.0)   # mm – concrete above deck profile (floor at 1 mm)
        self.d_p = max(h_mm - deck.e, 1.0)      # mm – effective depth to deck centroid (floor at 1 mm)
        self.b = 1.0                 # unit width, m
        self.b_mm = 1000.0           # unit width, mm

        # equivalent average concrete depth (topping + trough fill)
        # trough volume from manufacturer data (m³/m²), converted to mm
        self.h_conc_avg = self.h_c + deck.conc_vol_trough * 1000.0   # mm

        # loads (N/m²) — SI throughout; internal calcs convert to kN/m² where needed
        self.imposed_load = imposed_load
        self.finishes_load = finishes_load
        self.ceiling_services = ceiling_services
        self.construction_load = construction_load
        self.partition_load = partition_load

        # partial factors
        self.gamma_c = gamma_c
        self.gamma_g = gamma_g
        self.gamma_q = gamma_q
        self.gamma_m0 = gamma_m0

        # creep parameters (EN 1992-1-1 Annex B)
        self.RH = RH              # relative humidity (%), typically 50% indoors
        self.t0 = t0              # age at loading (days); 28 days ≈ when SDL (finishes) are applied
        self.cement_class = cement_class  # "S", "N", or "R"

        # propped construction flag (EN 1994-1-1 Cl 9.3.1)
        self.propped = propped

        # fire parameters (EN 1994-1-2 Annex D)
        self.R_fi      = R_fi       # required fire resistance period (min): 30/60/90/120
        self.cover_top = cover_top  # nominal concrete cover from top slab surface to mesh bar face (mm)
                                    # axis distance = cover_top + d_s/2
        self.psi_fi    = psi_fi     # fire combination factor for variable loads (EN 1990)

        # static coefficients — depend on number of continuous spans
        self.n_spans = n_spans

        # Construction stage: deck sheet treated as simply-supported for every span,
        # regardless of multi-span layout.  Deck sheets are lapped or butt-jointed
        # over intermediate supports in practice and are not moment-continuous;
        # manufacturer construction-stage span tables also assume SS.  The composite
        # stage still uses elastic continuous analysis (below).
        # SS coefficients: M_sag = wL²/8, V = wL/2, δ = 5wL⁴/384EI — same for all n_spans.
        self.M_COEFF_CONSTR_SAG = 1.0 / 8.0
        self.M_COEFF_CONSTR_HOG = 0.0          # no hogging under SS assumption
        self.V_COEFF_CONSTR     = 0.5
        self.K_DEFL_CONSTR      = 5.0 / 384.0

        # Composite stage coefficients — elastic continuous beam analysis.
        # EN 1994-1-1 Cl 9.4.2(1) permits linear elastic analysis for ULS.
        # The previous simply-supported simplification (Cl 9.4.2(5)) has been
        # replaced by elastic continuous coefficients so that continuity is fully
        # exploited for sagging bending and longitudinal shear.
        #
        # Sagging moment coefficients (elastic, UDL on all spans):
        #   1-span: 1/8  (simply-supported — no continuity)
        #   2-span: 9/128 ≈ 0.0703  (midspan of each span)
        #   3-span: 0.080  (midspan of interior span, governs)
        #
        # m-k equivalent isostatic span factor (EN 1994-1-1 Cl 9.7.3(6)):
        #   For continuous slabs the m-k check uses a reduced equivalent span:
        #     internal spans → 0.8L,  external spans → 0.9L.
        #   External spans (0.9L) always govern over internal (0.8L) because
        #   V_Ed × L_s ∝ L_equiv² and 0.9² > 0.8².  A single factor of 0.9 is
        #   therefore used for all multi-span cases (both 2-span and 3-span).
        # L_X_FACTOR: distance from the critical sagging cross-section to the
        # nearest support, as a fraction of span L.  Used in the partial
        # connection method (EN 1994-1-1 Cl 9.7.4) where L_x is that distance.
        # Derived from elastic analysis — location where shear = 0 (max sagging):
        #   1-span: midspan → L_x = L/2 = 0.500L
        #   2-span: shear = 0 at 3L/8 from pin end → L_x = 3/8 = 0.375L
        #   3-span: end span governs (0.08wL²); shear = 0 at 0.4L from end → L_x = 0.4L
        #           (interior span sagging = 0.025wL² — much less, not critical)
        if n_spans == 1:
            self.M_COEFF_SAG  = 1.0 / 8.0
            self.M_COEFF_HOG  = 0.0
            self.V_COEFF      = 0.5
            self.K_DEFL       = 5.0 / 384.0
            self.L_EQUIV_MK   = 1.0          # no reduction for SS span
            self.L_X_FACTOR   = 0.5          # midspan
        elif n_spans == 2:
            self.M_COEFF_SAG  = 9.0 / 128.0  # elastic 2-span midspan moment
            self.M_COEFF_HOG  = 1.0 / 8.0    # elastic 2-span support moment
            self.V_COEFF      = 5.0 / 8.0
            self.K_DEFL       = 2.0 / 369.0
            self.L_EQUIV_MK   = 0.9          # Cl 9.7.3(6): external spans
            self.L_X_FACTOR   = 3.0 / 8.0   # = 0.375L (shear = 0 at 3L/8 from pin)
        else:
            self.M_COEFF_SAG  = 0.080        # elastic 3-span end-span moment (governs)
            self.M_COEFF_HOG  = 0.100        # elastic 3-span support moment
            self.V_COEFF      = 0.600
            self.K_DEFL       = 1.0 / 145.0
            self.L_EQUIV_MK   = 0.9          # Cl 9.7.3(6): external spans govern
            self.L_X_FACTOR   = 0.4          # end span: shear = 0 at 0.4L from end
        # alias kept for backward compatibility with checks that reference M_COEFF
        self.M_COEFF = self.M_COEFF_SAG

        # dead weight of composite slab [N/m] – includes concrete in troughs
        self.g0k = (concrete.dry_weight / 1000 * self.h_conc_avg / 1000.0 + deck.weight_kN) * 1000.0 * self.b

        # select minimum mesh reinforcement (EN 1994-1-1 Cl 9.8.1)
        # 0.2% Ac for unpropped; 0.4% Ac for propped construction
        self.mesh = None
        if database is not None:
            As_min = (0.004 if propped else 0.002) * self.b_mm * self.h_c  # mm²/m
            for mesh_name in self.MESH_TYPES:
                mesh = SteelReinforcingMesh(mech_prop=mesh_name, database=database)
                if mesh.A_s >= As_min:
                    self.mesh = mesh
                    break
            if self.mesh is None:
                self.mesh = SteelReinforcingMesh(mech_prop=self.MESH_TYPES[-1], database=database)

        # 1 = single layer; 2 = double layer (set by _chk_comp_uls_hogging when needed)
        self.mesh_layers = 1

        # GWP of composite slab [kg-CO2-eq/m] (per metre of beam length, consistent with other sections)
        # concrete contribution – includes trough fill (concrete.GWP is kg-CO2-eq/kg)
        co2_concrete = self.b * (self.h_conc_avg / 1000.0) * concrete.density * concrete.GWP  # [kg-CO2-eq/m]

        # deck contribution – adapt to GWP_UNIT stored in database
        gwp_unit = str(deck.GWP_UNIT).strip().lower() if deck.GWP_UNIT else ""
        if "m2" in gwp_unit or "m²" in gwp_unit:
            co2_deck = deck.GWP * self.b                                     # [kg-CO2-eq/m]
        elif "t" in gwp_unit or "ton" in gwp_unit:
            co2_deck = deck.weight_kg * deck.GWP / 1000.0 * self.b           # [kg-CO2-eq/m]
        elif "kg" in gwp_unit:
            co2_deck = deck.weight_kg * deck.GWP * self.b                    # [kg-CO2-eq/m]
        else:
            co2_deck = deck.GWP * self.b                                     # fallback: assume per m²

        # reinforcement GWP (B500B / Steel_reinforcing_bar) — shared by mesh and trough bars
        _rebar_gwp = 0.0
        _bar_fsk   = 500e6   # Pa — nominal fsk for B500B (fallback if DB unavailable)
        if database is not None:
            _rebar     = SteelReinforcingBar("'B500B'", database)
            _rebar_gwp = _rebar.GWP   # kg-CO2-eq/kg
            _bar_fsk   = _rebar.fsk   # Pa

        # mesh contribution
        co2_mesh = 0.0
        if self.mesh is not None:
            mesh_mass = self.mesh.A_s * 1e-6 * self.STEEL_DENSITY * self.b   # kg/m
            co2_mesh  = mesh_mass * _rebar_gwp                                # [kg-CO2-eq/m]

        # trough bar reinforcement — diameter selected by fire check (_chk_fire_resistance);
        # initialised to zero here; fire check sets these attributes and updates self.co2.
        self.bar_dia_trough  = None   # mm, set by fire check
        self.A_s_trough      = 0.0   # mm²/m
        self.d_bar_trough    = 0.0   # mm from bottom face of slab to bar centroid
        self.co2_trough_bars = 0.0
        self._bar_fsk        = _bar_fsk    # Pa — for fire strength-reduction
        self._rebar_gwp      = _rebar_gwp  # kg-CO2-eq/kg — for GWP update after fire check

        self.co2_concrete = co2_concrete
        self.co2_deck     = co2_deck
        self.co2_mesh     = co2_mesh
        self.co2 = co2_concrete + co2_deck + co2_mesh + self.co2_trough_bars  # [kg-CO2-eq/m]

    @property
    def mass(self):
        """Structural mass per unit floor area [kg/m²]: concrete + deck + mesh + trough bars."""
        conc_mass = (self.concrete.dry_weight / 9.81) * (self.h_conc_avg / 1000.0)
        deck_mass = self.deck.weight_kg
        mesh_mass = (self.mesh.A_s * self.mesh_layers * 1e-6 * self.STEEL_DENSITY) if self.mesh is not None else 0.0
        bar_mass  = self.A_s_trough * 1e-6 * self.STEEL_DENSITY
        return conc_mass + deck_mass + mesh_mass + bar_mass

    # ── design checks (span-dependent) ───────────────────────────────────────
    def run_all_checks(self, span, criterion="ENV"):
        """Run EN 1994-1-1 and EN 1994-1-2 checks for given *span* (m).
        Returns max utilisation across all relevant checks.

        criterion : "ULS"  → ULS checks only
                    "SLS1"/"SLS2" → SLS checks only
                    "FIRE" → fire check only
                    "ENV"  → all checks including fire (if R_fi > 0)
        """
        # _r_redist: composite stage moment redistribution factor (0.0–0.30).
        # Set by _chk_comp_uls_hogging when elastic hogging exceeds mesh capacity;
        # consumed by _chk_comp_uls_bending and _chk_comp_uls_long_shear.
        # Reset to 0.0 at the start of every run so checks are independent.
        self._r_redist = 0.0
        # Tracks whether plastic redistribution was used at construction stage.
        # If True, the deck has been plastified over the internal support and
        # cannot contribute to composite hogging resistance (Cl 9.7.3(2)P).
        # For propped construction props are at midspan so the deck sub-spans
        # are simply-supported — no hogging develops at the permanent internal
        # supports, so this is always False.
        self._constr_redist_applied = False

        checks = []
        if self.propped:
            # Props placed at midspan divide each span into two L/2 simply-supported
            # sub-spans during casting.  The deck still carries the wet concrete load
            # but over L/2, not L.  Run construction checks with SS coefficients at
            # the propped sub-span; temporarily override the stored coefficients so
            # the check methods use the correct values, then restore them.
            # No hogging develops at permanent supports during construction (SS
            # sub-spans), so _chk_constr_sls_support is omitted.
            _saved_constr = (self.M_COEFF_CONSTR_SAG, self.M_COEFF_CONSTR_HOG,
                             self.V_COEFF_CONSTR, self.K_DEFL_CONSTR, self.n_spans)
            try:
                self.M_COEFF_CONSTR_SAG = 1.0 / 8.0    # SS sub-span sagging
                self.M_COEFF_CONSTR_HOG = 0.0           # SS — no hogging
                self.V_COEFF_CONSTR     = 0.5           # SS shear coefficient
                self.K_DEFL_CONSTR      = 5.0 / 384.0  # SS deflection coefficient
                self.n_spans            = 1             # suppress hogging branch
                checks.append(self._chk_constr_uls_bending(span / 2.0))
                checks.append(self._chk_constr_uls_shear(span / 2.0))
                checks.append(self._chk_constr_sls_deflection(span / 2.0))
            finally:
                (self.M_COEFF_CONSTR_SAG, self.M_COEFF_CONSTR_HOG,
                 self.V_COEFF_CONSTR, self.K_DEFL_CONSTR, self.n_spans) = _saved_constr
        else:
            # Unpropped: deck spans the full distance between permanent supports.
            # Construction stage uses SS assumption (deck lapped over supports) so
            # there is no hogging check and no internal-support interaction check.
            checks.append(self._chk_constr_uls_bending(span))
            checks.append(self._chk_constr_uls_shear(span))
            checks.append(self._chk_constr_sls_deflection(span))

        # Hogging check runs BEFORE sagging/shear: if redistribution is needed
        # it sets self._r_redist so the sagging checks use the increased demand.
        if self.n_spans > 1:
            checks.append(self._chk_comp_uls_hogging(span))

        checks.append(self._chk_comp_uls_bending(span))
        checks.append(self._chk_comp_uls_long_shear(span))
        checks.append(self._chk_comp_uls_vert_shear(span))
        delta_lock = self._locked_in_deflection(span)
        checks.append(self._chk_comp_sls_deflection(span, delta_lock))
        # fire checks — only when a fire period has been specified (R_fi > 0)
        if self.R_fi and self.R_fi > 0:
            for fire_chk in self._chk_fire_resistance(span):
                checks.append(fire_chk)

        if criterion == "ULS":
            relevant = [c for c in checks if "ULS" in c["name"]]
        elif criterion in ("SLS1", "SLS2"):
            relevant = [c for c in checks if "SLS" in c["name"]]
        elif criterion == "FIRE":
            relevant = [c for c in checks if "FIRE" in c["name"]]
        else:   # ENV — all checks
            relevant = checks

        self.checks = checks
        # D5 (field of application) is a warning-only check: it flags that the
        # Annex D method is technically outside its validation range, but the
        # calculation still proceeds conservatively.  Exclude it from all_passed
        # so that geometrically valid structural results are not discarded.
        structural = [c for c in relevant if not c["name"].startswith("D5")]
        self.all_passed = all(c["passed"] for c in structural)
        self.max_util = max((c["util"] for c in structural), default=0.0)
        return self.max_util

    # ── construction stage ───────────────────────────────────────────────────
    def _constr_loads(self):
        G_k1a = self.deck.weight_kN                                     # kN/m²
        Q_k1c = self.concrete.wet_weight / 1000 * self.h_conc_avg / 1000.0  # kN/m²
        Q_k1b = self.construction_load / 1e3                            # N/m² → kN/m²
        Q_k1a = max(0.75, 0.10 * Q_k1c)                                # kN/m²
        return G_k1a, Q_k1a, Q_k1b, Q_k1c

    @staticmethod
    def _patch_M_mid(q, a, L):
        """Max midspan moment from UDL *q* over patch length *a* centred at midspan
        on a simply-supported beam of span *L*."""
        a = min(a, L)
        return q * a / 2.0 * (L / 2.0 - a / 4.0)

    @staticmethod
    def _patch_V_max(q, a, L):
        """Max support reaction from UDL *q* over patch length *a* placed at the
        support on a simply-supported beam of span *L*."""
        a = min(a, L)
        return q * a * (2.0 * L - a) / (2.0 * L)

    @staticmethod
    def _patch_delta_mid(q, a, L, EI):
        """Midspan deflection from UDL *q* over patch length *a* centred at midspan
        on a simply-supported beam of span *L*, stiffness *EI*."""
        a = min(a, L)
        # closed-form: delta = q/(384*EI) * (a * (8*L^3 - 4*L*a^2 + a^3) * 2 ... )
        # simpler: use superposition of full span minus two end portions
        # For a symmetric patch of length a centred on span L:
        c = (L - a) / 2.0  # unloaded length each side
        if c <= 0:
            return 5.0 * q * L ** 4 / (384.0 * EI)
        # midspan deflection for UDL q from x=c to x=L-c (symmetric about midspan)
        # Using integration: delta_mid = q/(24*EI) * ( c*(L-c)/2*(L^2+c*(L-c)) - c^4/4 ... )
        # Easier to compute from known formula for partial UDL.
        # delta_mid = q*a*(5*a^4 - 24*a^2*L^2 + 16*L^4) / (384*EI*... )
        # Actually use the exact closed-form for symmetric load:
        # delta_mid = q / (384 EI) * (5*a^4/16 - 6*a^2*L^2/4 + ... )
        # Safest approach: numerical integration via beam equation
        # For a simply-supported beam, partial UDL from c to L-c:
        # R = q*a/2, M(x) for x in [0,c]: R*x
        # Integrate using Macaulay / virtual work:
        # delta_mid = q*a/(384*EI) * (8*L**3 - 4*L*a**2 + a**3)  (symmetric load)
        return q * a / (384.0 * EI) * (8.0 * L ** 3 - 4.0 * L * a ** 2 + a ** 3)

    def _chk_constr_uls_bending(self, span):
        G_k1a, Q_k1a, Q_k1b, Q_k1c = self._constr_loads()
        L = span  # m
        a_patch = 3.0  # m – 3 m × 3 m working area (SCI P300 Fig. 3.1)

        # Ponding load (SCI P300 Fig. 3.1 / EN 1994-1-1 Cl 9.3.2):
        # If the initial deflection under deck SW + wet concrete exceeds h/10,
        # add an extra variable action = wet density × 0.7 × δ_initial (kN/m²).
        L_mm = L * 1000.0
        EI_defl = self.deck.E / 1e6 * self.deck.I_p * 1e4  # N/mm² × mm⁴
        delta_init = self.K_DEFL_CONSTR * (G_k1a + Q_k1c) * L_mm ** 4 / EI_defl
        Q_pond = 0.0
        if delta_init > self.h_mm / 10.0:
            Q_pond = (self.concrete.wet_weight / 1000.0) * (0.7 * delta_init / 1000.0)

        # Loading model (SCI P300 Cl. 3.2.1):
        #   Q_k1b = 0.75 kN/m²  general construction load acting over full span
        #   Q_k1a = max(0.75, 0.10 × Q_k1c)  ADDITIONAL patch load over 3 m working area
        #   Q_pond = ponding load (variable action), included where δ_init > h/10
        #
        # Case 1: conservative upper bound — treat patch intensity Q_k1a as full-span
        #         UDL (maximises moment/shear for spans ≤ 3 m where patch = full span).
        w_uls_1 = self.gamma_g * G_k1a + self.gamma_q * (Q_k1b + Q_k1a + Q_k1c + Q_pond)
        M_Ed_sag_1 = self.M_COEFF_CONSTR_SAG * w_uls_1 * L ** 2

        # Case 2: general UDL over full span + patch at mid-span for max sagging moment.
        #         For spans > 3 m this is less critical than Case 1 and acts as a
        #         cross-check; for spans ≤ 3 m both cases are identical (a_patch → L).
        w_base = self.gamma_g * G_k1a + self.gamma_q * (Q_k1b + Q_k1c + Q_pond)
        q_patch = self.gamma_q * Q_k1a
        M_Ed_sag_2 = self.M_COEFF_CONSTR_SAG * w_base * L ** 2 + self._patch_M_mid(q_patch, a_patch, L)

        M_Ed_sag = max(M_Ed_sag_1, M_Ed_sag_2)
        # Deck is treated as simply-supported at construction stage (SS assumption):
        # no hogging check, no redistribution.  _constr_redist_applied remains False,
        # so the deck is always eligible to contribute to composite hogging (Cl 9.7.3(2)P).
        util_sag = M_Ed_sag / self.deck.M_cRd
        return {"name": "Construction ULS - Bending", "util": util_sag, "passed": util_sag <= 1.0}

    def _chk_constr_uls_shear(self, span):
        G_k1a, Q_k1a, Q_k1b, Q_k1c = self._constr_loads()
        L = span
        a_patch = 3.0

        # Ponding load (same calculation as in _chk_constr_uls_bending)
        L_mm = L * 1000.0
        EI_defl = self.deck.E / 1e6 * self.deck.I_p * 1e4
        delta_init = self.K_DEFL_CONSTR * (G_k1a + Q_k1c) * L_mm ** 4 / EI_defl
        Q_pond = 0.0
        if delta_init > self.h_mm / 10.0:
            Q_pond = (self.concrete.wet_weight / 1000.0) * (0.7 * delta_init / 1000.0)

        # Case 1: uniform
        w_uls_1 = self.gamma_g * G_k1a + self.gamma_q * (Q_k1a + Q_k1b + Q_k1c + Q_pond)
        V_Ed_1 = self.V_COEFF_CONSTR * w_uls_1 * L

        # Case 2: general UDL + patch placed at support for max end shear
        w_base = self.gamma_g * G_k1a + self.gamma_q * (Q_k1b + Q_k1c + Q_pond)
        q_patch = self.gamma_q * Q_k1a
        V_Ed_2 = self.V_COEFF_CONSTR * w_base * L + self._patch_V_max(q_patch, a_patch, L)

        V_Ed = max(V_Ed_1, V_Ed_2)

        h_w = self.deck.h_p - self.deck.t
        lambda_w = 0.346 * self.deck.s_w / self.deck.t * math.sqrt(self.deck.f_ypd / (self.deck.E / 1e6))
        if lambda_w <= 0.83:
            f_bv = 0.58 * self.deck.f_ypd
        elif lambda_w < 1.40:
            f_bv = 0.48 * self.deck.f_ypd / lambda_w
        else:
            f_bv = 0.67 * self.deck.f_ypd / lambda_w ** 2

        V_bRd_web = (h_w * self.deck.t * f_bv
                     / math.sin(math.radians(self.deck.theta))
                     / self.gamma_m0 / 1000)                    # kN per web
        if self.deck.num_webs is None:
            return {"name": "Construction ULS - Shear",
                    "util": float('inf'), "passed": False,
                    "note": "num_webs missing in database"}
        V_bRd = V_bRd_web * self.deck.num_webs                  # kN/m width
        util = V_Ed / V_bRd
        return {"name": "Construction ULS - Shear", "util": util, "passed": util <= 1.0}

    def _chk_constr_sls_deflection(self, span):
        """EN 1994-1-1 Cl 9.3.2 + UK NA Cl 4/NA2.15 (from BS 5950-4).

        Per SCI P359 Section 3.5.2, the two options are paired — the limit and
        the load model must match:

        Option 1 — ponding load OMITTED (initial δ ≤ h/10):
            Load  = deck SW + wet concrete only (no ponding extra weight)
            Check = δ_initial ≤ min(L/180, 20 mm)

        Option 2 — ponding load INCLUDED (initial δ > h/10):
            Extra load = wet density × 0.7 × δ_initial (variable action, SCI P359 Cl 3.2.1)
            Recalculate δ_ponded under the increased load
            Check = δ_ponded ≤ min(L/130, 30 mm)

        The lenient L/130 / 30 mm limit applies only when the ponding load has
        been included in the deflection calculation — not to δ_initial.
        Using δ_initial with the L/130 limit would combine the shortcut load
        model with the lenient limit, which is non-conservative (SCI P359).
        """
        G_k1a, Q_k1a, Q_k1b, Q_k1c = self._constr_loads()
        L = span * 1000.0  # mm
        EI = self.deck.E / 1e6 * self.deck.I_p * 1e4  # N/mm² × mm⁴

        # Characteristic load: deck SW + wet concrete only (EN 1994-1-1 Cl 9.6(2))
        w_sls = G_k1a + Q_k1c
        delta_initial = self.K_DEFL_CONSTR * w_sls * L ** 4 / EI

        if delta_initial <= self.h_mm / 10.0:
            # Option 1: ponding not significant — check initial δ against tighter limit
            delta_check = delta_initial
            limit = min(L / 180.0, 20.0)
        else:
            # Option 2: ponding significant — add 0.7×δ_initial as extra depth,
            # recalculate deflection under ponded load, check against lenient limit
            w_pond = (self.concrete.wet_weight / 1000.0) * (0.7 * delta_initial / 1000.0)
            delta_check = self.K_DEFL_CONSTR * (w_sls + w_pond) * L ** 4 / EI
            limit = min(L / 130.0, 30.0)

        util = delta_check / limit
        return {"name": "Construction SLS - Deflection", "util": util, "passed": util <= 1.0}

    def _chk_constr_sls_support(self, span):
        """EN 1993-1-3 Cl 7.2 + Cl 6.1.11 — internal support interaction check.

        Required whenever plastic global analysis (moment redistribution) is
        assumed at construction stage ULS (EN 1993-1-3 Cl 6.1.4.1(8)).

        At the internal support, using SLS (characteristic) loads and γ_M,ser = 1.0:

            M_Ed / M_c,Rd + F_Ed / R_w,Rd  ≤  0.9   (Eq. 6.28c, limit per Cl 7.2)

        where:
            M_Ed      = elastic hogging moment at SLS (kNm/m)
            M_c,Rd    = hogging moment resistance = M_tRd with γ_M,ser = 1.0 (kNm/m)
            F_Ed      = total reaction at internal support at SLS (kN/m)
            R_w,Rd    = web crippling resistance per unit width (kN/m),
                        from EN 1993-1-3 Cl 6.1.7.3, Eq. (6.18):
                        Category 2, α = 0.15, l_a = 10 mm (for β_v ≈ 0.5)
        """
        if self.deck.M_tRd is None:
            return {"name": "Constr SLS - Internal Support (Cl 7.2)",
                    "util": float('inf'), "passed": False,
                    "note": "M_tRd not in database"}
        if self.deck.num_webs is None:
            return {"name": "Constr SLS - Internal Support (Cl 7.2)",
                    "util": float('inf'), "passed": False,
                    "note": "num_webs missing in database"}

        G_k1a, Q_k1a, Q_k1b, Q_k1c = self._constr_loads()
        L = span  # m

        # Ponding load (characteristic) if significant
        L_mm = L * 1000.0
        EI_defl = self.deck.E / 1e6 * self.deck.I_p * 1e4
        delta_init = self.K_DEFL_CONSTR * (G_k1a + Q_k1c) * L_mm ** 4 / EI_defl
        Q_pond = 0.0
        if delta_init > self.h_mm / 10.0:
            Q_pond = (self.concrete.wet_weight / 1000.0) * (0.7 * delta_init / 1000.0)

        # SLS (characteristic) construction load — no γ factors
        w_sls = G_k1a + Q_k1b + Q_k1c + Q_pond  # kN/m²

        # Elastic hogging moment at SLS
        M_Ed_sls = self.M_COEFF_CONSTR_HOG * w_sls * L ** 2  # kNm/m

        # Internal support reaction at SLS:
        #   2-span UDL: R = 10/8 × w × L = 1.25 wL
        #   3-span UDL: R = 11/10 × w × L = 1.10 wL  (interior support)
        R_coeff = 1.25 if self.n_spans == 2 else 1.10
        F_Ed_sls = R_coeff * w_sls * L  # kN/m

        # Hogging moment resistance at SLS (γ_M,ser = 1.0)
        M_c_Rd = self.deck.M_tRd  # kNm/m

        # Web crippling resistance — EN 1993-1-3 Cl 6.1.7.3, Eq. (6.18)
        # Category 2 (internal support, sheeting): α = 0.15 (Eq. 6.20c)
        # β_v ≈ 0.5 for equal spans → interpolation gives l_a = 10 mm (Eq. 6.19b–c)
        # r = 4 mm assumed: manufacturer data for all profiled decks in the
        # database shows inner corner radii of 3–4 mm; 4 mm is adopted as a
        # representative conservative-side assumption.
        t      = self.deck.t            # mm
        f_yb   = self.deck.f_ypd        # N/mm²  (f_yb ≈ f_ypd for γ_M0 = 1.0)
        E_mm   = self.deck.E / 1e6      # N/mm²
        phi    = self.deck.theta        # web angle to horizontal (°)
        l_a    = 10.0                   # mm (Category 2, β_v ≈ 0.5)
        alpha  = 0.15                   # Category 2 sheeting (Eq. 6.20c)
        r      = 4.0                    # mm — inner corner radius (assumed, see note above)

        R_w_per_web = (alpha * t ** 2 * math.sqrt(f_yb * E_mm)
                       * (1.0 - 0.1 * math.sqrt(r / t))            # EN 1993-1-3 Eq. (6.18)
                       * (0.5 + math.sqrt(0.02 * l_a / t))
                       * (2.4 + (phi / 90.0) ** 2))                # N per web (γ_M,ser = 1.0)

        R_w_Rd = R_w_per_web * self.deck.num_webs / 1000.0  # kN/m

        # Interaction check (Cl 7.2): combined ratio must not exceed 0.9
        interaction = M_Ed_sls / M_c_Rd + F_Ed_sls / R_w_Rd
        util = interaction / 0.9
        return {"name": "Constr SLS - Internal Support (Cl 7.2)",
                "util": util, "passed": util <= 1.0}

    def _locked_in_deflection(self, span):
        """Deflection locked into the slab when the concrete hardens (unpropped only).

        When the initial deflection exceeds h/10, ponding is significant: extra
        wet concrete fills the depression and hardens in place, increasing the
        permanent locked-in deflection.

        Per SCI P359 Cl 3.2.1 / EN 1994-1-1 Cl 9.3.2, the extra concrete depth
        is 0.7 × δ_initial (where δ_initial is the deflection under deck SW +
        wet concrete at nominal thickness).  This is a single non-iterative step:
        the 0.7 factor is explicitly chosen in the code to avoid circular iteration.
        """
        if self.propped:
            # Props support the deck during casting — no deflection is locked in
            # before composite action.  Self-weight is instead carried by the
            # composite section after prop removal (see _chk_comp_sls_deflection).
            return 0.0
        L = span * 1000.0
        w_sw = self.concrete.wet_weight / 1000 * self.h_conc_avg / 1000.0 + self.deck.weight_kN
        EI = self.deck.E / 1e6 * self.deck.I_p * 1e4  # N/mm² × mm⁴
        delta_initial = self.K_DEFL_CONSTR * w_sw * L ** 4 / EI
        if delta_initial > self.h_mm / 10.0:
            w_pond = (self.concrete.wet_weight / 1000.0) * (0.7 * delta_initial / 1000.0)
            return self.K_DEFL_CONSTR * (w_sw + w_pond) * L ** 4 / EI
        return delta_initial

    # ── composite stage ──────────────────────────────────────────────────────
    def _comp_uls_loads(self):
        """ULS design loads on the composite section (kN/m²).

        Propped construction: props carry the wet concrete + deck SW during
        casting; after prop removal the composite section resists *all*
        permanent loads (slab SW + deck SW + SDL) plus imposed.

        Unpropped construction: the deck alone carries slab SW and deck SW
        before composite action is established.  The composite section only
        sees post-construction loads (SDL + imposed).  Including slab SW here
        would double-count a load that was already carried by the deck at
        construction stage — equivalent to the split already made in the SLS
        deflection check (_chk_comp_sls_deflection).
        """
        if self.propped:
            # All permanent loads go on the composite section
            gk = (self.concrete.dry_weight / 1000.0 * self.h_conc_avg / 1000.0
                  + self.deck.weight_kN
                  + (self.finishes_load + self.ceiling_services) / 1e3)
        else:
            # Unpropped: slab SW + deck SW already carried by deck at constr. stage
            gk = (self.finishes_load + self.ceiling_services) / 1e3
        qk = (self.imposed_load + self.partition_load) / 1e3          # N/m² → kN/m²
        return self.gamma_g * gk + self.gamma_q * qk, gk + qk

    def _chk_comp_uls_bending(self, span):
        """EN 1994-1-1 Cl 9.7.2 — sagging bending resistance (full shear connection).

        Plastic stress block method.  Applied regardless of which longitudinal
        shear design method (m-k or partial connection) is used.
        """
        w_uls, _ = self._comp_uls_loads()
        # If hogging redistribution was applied, the surplus is shared equally
        # to adjacent sagging spans: M_Ed_sag increases by r × M_hog / 2.
        # EN 1994-1-1 Cl 9.4.2(3).
        M_Ed = (self.M_COEFF_SAG + self._r_redist * self.M_COEFF_HOG / 2.0) * w_uls * span ** 2

        fcd_MPa = self.concrete.fcd / 1e6                       # Pa → N/mm²
        N_c_steel = self.deck.A_pe * self.deck.f_ypd            # N/m — steel limit
        N_c_conc  = 0.85 * fcd_MPa * self.b_mm * self.h_c      # N/m — concrete limit
        N_cf = min(N_c_steel, N_c_conc)
        x_pl = N_cf / (0.85 * fcd_MPa * self.b_mm)             # mm — compression block depth
        z    = self.h_mm - self.deck.e - 0.5 * x_pl             # mm — lever arm (gross centroid, EN 1994-1-1 Eq 9.9 at full connection)
        M_Rd = N_cf * z / 1e6                                   # kNm/m

        util = M_Ed / M_Rd if M_Rd > 0 else float('inf')
        return {"name": "Composite ULS - Bending (Cl 9.7.2)",
                "util": util, "passed": util <= 1.0}

    def _chk_comp_uls_long_shear(self, span):
        w_uls, _ = self._comp_uls_loads()
        # M_Ed includes any redistribution from hogging (see _chk_comp_uls_hogging).
        # V_Ed is based on elastic reactions and is unaffected by moment redistribution.
        M_Ed = (self.M_COEFF_SAG + self._r_redist * self.M_COEFF_HOG / 2.0) * w_uls * span ** 2
        V_Ed = self.V_COEFF * w_uls * span

        has_mk = (self.deck.m_p is not None and self.deck.k_p is not None
                  and not math.isnan(self.deck.m_p) and not math.isnan(self.deck.k_p))
        has_tau = (self.deck.tau_uRk is not None and not math.isnan(self.deck.tau_uRk))

        GAMMA_VS = 1.25

        if has_mk:
            # EN 1994-1-1 Cl 9.7.3(4): V_Ed is the maximum design vertical shear
            # force from structural analysis — i.e. the actual reaction, NOT a
            # reduced value.  For a 2-span elastic beam V_Ed = 5/8 wL (= V_COEFF × w × span).
            # EN 1994-1-1 Cl 9.7.3(6): for continuous slabs an equivalent isostatic
            # span L_equiv = L_EQUIV_MK × span is used "for the determination of the
            # RESISTANCE" (V_l,Rd) only — external spans 0.9L, internal spans 0.8L.
            # L_s (shear span) is derived from L_equiv per Table 9.1 (L_s = L_equiv/4).
            # V_Ed is NOT replaced by 0.5 × w × L_equiv; only L_s and V_l,Rd use L_equiv.
            L_equiv = self.L_EQUIV_MK * span          # m — equivalent isostatic span (for resistance only)
            L_s = L_equiv * 1000.0 / 4                # mm — shear span per Table 9.1 / Cl 9.7.3(6)
            V_Ed_mk = self.V_COEFF * w_uls * span     # kN/m — actual max shear from structural analysis
            V_l_Rd = ((self.b_mm * self.d_p / GAMMA_VS)
                      * (self.deck.m_p * self.deck.A_p / (self.b_mm * L_s) + self.deck.k_p)
                      / 1000.0)
            util = V_Ed_mk / V_l_Rd if V_l_Rd > 0 else float('inf')
            return {"name": "Composite ULS - Longitudinal Shear (m-k)",
                    "util": util, "passed": util <= 1.0}

        elif has_tau:
            # L_x = distance from critical sagging section to nearest support
            # (EN 1994-1-1 Cl 9.7.4(8)).  L_X_FACTOR gives the correct fraction
            # of span for each span type based on elastic analysis (0.5, 0.375, 0.4
            # for 1-, 2-, 3-span respectively — see __init__ for derivation).
            L_x = self.L_X_FACTOR * span * 1000.0   # mm
            N_c_bond = (self.deck.tau_uRk / GAMMA_VS) * self.b_mm * L_x  # N/m — bond/friction
            N_cf = self.deck.A_pe * self.deck.f_ypd                        # N/m — full connection

            def _partial_moment_resistance(N_c_total):
                """Return M_Rd (kNm/m) for a given total longitudinal force N_c (N/m)."""
                N_c_capped = min(N_cf, N_c_total)
                fcd_MPa = self.concrete.fcd / 1e6  # Pa → N/mm²
                # x_pl from partial force N_c per Eq (9.8); z per Eq (9.9) — denominator is
                # A_pe×f_yp,d (not N_cf) as stated explicitly in EN 1994-1-1 Cl 9.7.3(8).
                x_pl = N_c_capped / (0.85 * fcd_MPa * self.b_mm)  # mm — Eq (9.8)
                z = (self.h_mm - 0.5 * x_pl - self.deck.e_p
                     + (self.deck.e_p - self.deck.e)
                       * N_c_capped / (self.deck.A_pe * self.deck.f_ypd))  # Eq (9.9)
                M_pa = self.deck.M_cRd
                M_pr = min(1.25 * M_pa * (1 - N_c_capped / (self.deck.A_pe * self.deck.f_ypd)),
                           M_pa)
                return (N_c_capped / 1000.0 * z / 1000.0) + M_pr  # kNm/m

            # ── Step 1: check without end anchorage ──────────────────────────
            M_Rd = _partial_moment_resistance(N_c_bond)
            util = M_Ed / M_Rd if M_Rd > 0 else float('inf')
            if util <= 1.0:
                return {"name": "Composite ULS - Longitudinal Shear (partial)",
                        "util": util, "passed": True}

            # Bond/friction alone insufficient — check fails (no end anchorage in optimisation)
            return {"name": "Composite ULS - Longitudinal Shear (partial)",
                    "util": util, "passed": False,
                    "note": "bond/friction insufficient"}
        else:
            return {"name": "Composite ULS - Longitudinal Shear",
                    "util": float('inf'), "passed": False}

    def _chk_comp_uls_vert_shear(self, span):
        w_uls, _ = self._comp_uls_loads()
        V_Ed = self.V_COEFF * w_uls * span
        C_Rd_c = 0.18 / self.gamma_c
        k = min(1.0 + math.sqrt(200.0 / self.d_p), 2.0)
        rho_l = min(self.deck.A_pe / (self.b_mm * self.d_p), 0.02)
        V_Rd_c = (C_Rd_c * k * (100.0 * rho_l * self.concrete.fck * 1e-6) ** (1/3)
                  * self.b_mm * self.d_p / 1000.0)
        v_min = 0.035 * k ** 1.5 * (self.concrete.fck * 1e-6) ** 0.5
        V_Rd = max(V_Rd_c, v_min * self.b_mm * self.d_p / 1000.0)
        util = V_Ed / V_Rd
        return {"name": "Composite ULS - Vertical Shear", "util": util, "passed": util <= 1.0}

    def _chk_comp_uls_hogging(self, span):
        """EN 1994-1-1 Cl 9.7.3 — hogging ULS at internal support (n_spans > 1).

        Top mesh acts in tension; concrete above the deck provides the compression
        zone (deck ignored in hogging per EN 1994-1-1 Cl 9.7.3(5)).  If the minimum
        mesh is insufficient the method steps up through MESH_TYPES and upgrades
        self.mesh, updating GWP — mirroring how the fire check upgrades trough bars.
        """
        if self.mesh is None or self.database is None:
            return {"name": "Composite ULS - Hogging", "util": float('inf'), "passed": False}

        # ULS hogging demand (kN⋅m per m width)
        w_uls, _ = self._comp_uls_loads()
        M_Ed_hog = self.M_COEFF_HOG * w_uls * span ** 2

        # Design strengths
        fsd_mesh = 500.0 / 1.15          # N/mm² — B500A mesh (standard)
        fcd      = self.concrete.fcd / 1e6  # Pa → N/mm²

        # Find current mesh index, then try current and heavier meshes
        current_idx = 0
        for idx, name in enumerate(self.MESH_TYPES):
            if name == self.mesh.mech_prop:
                current_idx = idx
                break

        M_Rd_hog      = 0.0
        chosen_mesh   = self.mesh
        chosen_layers = 1

        # EN 1994-1-1 Cl 9.7.3(2)P: deck may contribute to hogging compression
        # only if (a) the sheet is continuous over the support and (b) no plastic
        # moment redistribution was used at the construction stage.
        # Condition (a) is satisfied whenever n_spans > 1 (deck is lapped
        # continuously).  Condition (b) is tracked by _constr_redist_applied.
        include_deck = (not self._constr_redist_applied
                        and self.deck.M_tRd is not None)

        def _hogging_resistance(mesh, n_layers):
            """Return M_Rd_hog (kNm/m) for a given mesh and number of layers."""
            d_s = mesh.d_s if (mesh.d_s or 0) > 0 else 6.0
            A_s = mesh.A_s * n_layers
            # Centroid from slab top: single = cover + d_s/2; double = cover + d_s
            axis = self.cover_top + d_s / 2.0 * n_layers
            d_eff = self.h_mm - axis  # mm from slab bottom
            N_s = A_s * fsd_mesh
            if include_deck:
                N_deck = min(self.deck.A_pe * self.deck.f_ypd, N_s)
                N_conc = max(N_s - N_deck, 0.0)
                x_pl   = N_conc / (0.85 * self.b_mm * fcd) if N_conc > 0 else 0.0
                # Lever arm uses gross centroid e (not effective e_p).
                # e_p is the centroid of the effective section for the construction
                # stage (sagging, top in compression).  In hogging, the deck carries
                # compression at the BOTTOM — the compression resultant acts near the
                # soffit, closer to e (gross centroid, 134 mm) than to e_p (185 mm).
                # Using e_p here severely under-estimates the lever arm for deep decks.
                return (N_deck * (d_eff - self.deck.e)
                        + N_conc * (d_eff - 0.5 * x_pl)) / 1e6
            else:
                x_pl = N_s / (0.85 * self.b_mm * fcd)
                return N_s * max(d_eff - 0.5 * x_pl, 0.0) / 1e6

        # ── pass 1: single-layer meshes (lightest sufficient) ────────────────
        for mesh_name in self.MESH_TYPES[current_idx:]:
            mesh = (self.mesh if mesh_name == self.mesh.mech_prop
                    else SteelReinforcingMesh(mech_prop=mesh_name, database=self.database))
            M_Rd_hog    = _hogging_resistance(mesh, 1)
            chosen_mesh  = mesh
            chosen_layers = 1
            if M_Rd_hog >= M_Ed_hog:
                break

        # ── pass 2: double-layer meshes (only if single-layer insufficient) ──
        if M_Rd_hog < M_Ed_hog:
            for mesh_name in self.MESH_TYPES:
                mesh = SteelReinforcingMesh(mech_prop=mesh_name, database=self.database)
                M_Rd_hog    = _hogging_resistance(mesh, 2)
                chosen_mesh  = mesh
                chosen_layers = 2
                if M_Rd_hog >= M_Ed_hog:
                    break

        # Apply mesh upgrade (type or layer count) if needed
        if chosen_mesh.mech_prop != self.mesh.mech_prop or chosen_layers != self.mesh_layers:
            mesh_mass        = chosen_mesh.A_s * chosen_layers * 1e-6 * self.STEEL_DENSITY * self.b
            self.co2_mesh    = mesh_mass * self._rebar_gwp
            self.co2         = (self.co2_concrete + self.co2_deck
                                + self.co2_mesh + self.co2_trough_bars)
            self.mesh        = chosen_mesh
            self.mesh_layers = chosen_layers

        layers_note   = " (double layer)" if chosen_layers == 2 else ""
        deck_note = ("deck included in compression (Cl 9.7.3(2)P)" if include_deck
                     else "deck excluded (construction redistribution applied)")

        if M_Rd_hog >= M_Ed_hog:
            # No redistribution needed — elastic hogging satisfied by mesh.
            self._r_redist = 0.0
            util   = M_Ed_hog / M_Rd_hog if M_Rd_hog > 0 else float('inf')
            return {"name": "Composite ULS - Hogging", "util": util, "passed": True,
                    "note": f"{deck_note}{layers_note}"}

        # Heaviest mesh (single or double) cannot satisfy elastic hogging demand.
        # Try moment redistribution per EN 1994-1-1 Cl 9.4.2(3): support moments
        # may be reduced by up to 30%, with corresponding increases to span moments.
        # The increased span demand is picked up by _chk_comp_uls_bending and
        # _chk_comp_uls_long_shear via self._r_redist.
        r_needed = 1.0 - (M_Rd_hog / M_Ed_hog)   # fraction that must be shed
        if r_needed <= 0.30:
            # Redistribution within Eurocode limit — hogging now just passes.
            self._r_redist = r_needed
            return {"name": "Composite ULS - Hogging",
                    "util": 1.0, "passed": True,
                    "note": f"redistribution r={r_needed:.2%} applied (Cl 9.4.2(3)){layers_note}"}
        else:
            # Even 30% redistribution is insufficient — hogging genuinely fails.
            self._r_redist = 0.30   # apply maximum; sagging checks will reflect this
            M_Ed_hog_redist = M_Ed_hog * (1.0 - 0.30)
            util = M_Ed_hog_redist / M_Rd_hog if M_Rd_hog > 0 else float('inf')
            return {"name": "Composite ULS - Hogging",
                    "util": util, "passed": False,
                    "note": f"max 30% redistribution applied but hogging still fails{layers_note}"}

    def _creep_coefficient(self):
        """Creep coefficient phi(inf, t0) per EN 1992-1-1 Annex B.
        Returns phi_0 for t = infinity (beta_c = 1.0)."""
        fck_MPa = self.concrete.fck / 1e6
        fcm = fck_MPa + 8.0  # EN 1992-1-1 Table 3.1

        # notional size h0 = 2*Ac/u (mm)  [EN 1992-1-1 Eq (B.6)]
        # Assumption: steel deck seals the bottom face, so only the top face
        # is exposed to drying → u = b (not 2b). This gives h0 = 2*h_c.
        h0 = 2.0 * (self.b_mm * self.h_c) / self.b_mm  # = 2*h_c

        # adjusted loading age for cement class (B.9)
        alpha_cement = {"S": -1, "N": 0, "R": 1}.get(self.cement_class, 0)
        t0_adj = max(self.t0 * (9.0 / (2.0 + self.t0 ** 1.2) + 1) ** alpha_cement, 0.5)

        # phi_RH (B.3a / B.3b)
        RH = self.RH
        if fcm <= 35.0:
            phi_RH = 1.0 + (1.0 - RH / 100.0) / (0.1 * h0 ** (1.0 / 3.0))
        else:
            alpha_1 = (35.0 / fcm) ** 0.7
            alpha_2 = (35.0 / fcm) ** 0.2
            phi_RH = (1.0 + (1.0 - RH / 100.0) / (0.1 * h0 ** (1.0 / 3.0)) * alpha_1) * alpha_2

        # beta(fcm) (B.4)
        beta_fcm = 16.8 / math.sqrt(fcm)

        # beta(t0) (B.5)
        beta_t0 = 1.0 / (0.1 + t0_adj ** 0.20)

        phi_0 = phi_RH * beta_fcm * beta_t0
        return phi_0

    def _composite_I(self, modular_ratio=None):
        """Composite transformed second moment of area (mm⁴/m).
        If *modular_ratio* is given it is used directly; otherwise
        the long-term modular ratio for permanent loads is computed
        from the creep coefficient per EN 1992-1-1 Annex B."""
        n_0 = self.deck.E / self.concrete.Ecm
        if modular_ratio is not None:
            n = modular_ratio
        else:
            phi_t = self._creep_coefficient()
            n_L = n_0 * (1.0 + 1.1 * phi_t)        # EN 1994-1-1 Cl 5.4.2.2
            n = n_0 / 3.0 + 2.0 * n_L / 3.0         # SCI P359: weighted average
        b_c = self.b_mm / n
        A_c = b_c * self.h_c
        y_c = self.deck.h_p + self.h_c / 2.0
        A_pe = self.deck.A_pe                  # effective area (mm²/m)
        y_p = self.deck.e                      # centroid of effective area from bottom (mm)
        y_na = (A_c * y_c + A_pe * y_p) / (A_c + A_pe)
        I_c = b_c * self.h_c ** 3 / 12.0 + A_c * (y_c - y_na) ** 2
        I_deck = self.deck.I_p * 1e4 + A_pe * (y_p - y_na) ** 2  # I_p: cm⁴ → mm⁴
        return I_c + I_deck

    def _chk_comp_sls_deflection(self, span, delta_locked_in):
        L = span * 1000.0
        n_0 = self.deck.E / self.concrete.Ecm
        phi_t = self._creep_coefficient()
        n_L = n_0 * (1.0 + 1.1 * phi_t)
        E_mm = self.deck.E / 1e6  # N/m² → N/mm²

        # variable load: use short-term modular ratio (no creep)
        I_short = self._composite_I(modular_ratio=n_0)
        qk = (self.imposed_load + self.partition_load) / 1e3   # N/m² → kN/m²
        delta_var = self.K_DEFL * qk * L ** 4 / (E_mm * I_short)

        # permanent load on composite section: SDL (finishes + ceiling/services)
        # long-term modular ratio for creep
        I_long = self._composite_I(modular_ratio=n_L)
        gk_sdl = (self.finishes_load + self.ceiling_services) / 1e3   # N/m² → kN/m²
        delta_sdl = self.K_DEFL * gk_sdl * L ** 4 / (E_mm * I_long)

        if self.propped:
            # Propped: props carry wet concrete load; on prop removal the composite
            # section resists the full self-weight.  The L/300 total limit applies
            # to (self-weight on composite + SDL + variable).
            gk_sw = (self.concrete.dry_weight / 1000.0 * self.h_conc_avg / 1000.0
                     + self.deck.weight_kN)   # kN/m²
            delta_sw = self.K_DEFL * gk_sw * L ** 4 / (E_mm * I_long)
            delta_total = delta_sw + delta_sdl + delta_var
        else:
            # Unpropped: the locked-in construction deflection is permanent sag in
            # the finished floor and must be included in the total deflection check.
            # The construction-stage L/180 check (EN 1994-1-1 Cl 9.3.2) limits
            # the deflection during casting; the composite-stage L/300 check is a
            # separate serviceability criterion for the finished floor.  Including
            # delta_locked_in here is not double-counting — the two limits serve
            # different purposes.
            delta_total = delta_locked_in + delta_sdl + delta_var

        limit_total = L / 300.0
        limit_var = min(L / 350.0, 20.0)  # EN 1994-1-1 Cl 9.8.2(4) / SCI P359
        util = max(delta_total / limit_total, delta_var / limit_var)
        return {"name": "Composite SLS - Deflection", "util": util,
                "passed": delta_total <= limit_total and delta_var <= limit_var,
                "delta_total_mm": delta_total,
                "delta_var_mm":   delta_var,
                "limit_total_mm": limit_total,
                "limit_var_mm":   limit_var}

    def _chk_fire_resistance(self, span):
        """EN 1994-1-2 Annex D — fire resistance checks for composite slabs.

        Standard fire ISO 834.  Returns a list of
        {"name": str, "passed": bool, "util": float} dicts, one per clause:

          D1 – Thermal insulation         (formula D.1, t_i ≥ R_fi)
          D2 – Sagging moment             (deck temps formula D.4, mesh temp from
                                           mesh-temperature table, trough bar temp
                                           formula D.5; M_Rd,fi ≥ M_Ed,fi)
          D3 – Hogging moment at internal support (n_spans > 1 only).
               Reduced cross-section per D.8–D.14: isotherm penetration b into
               the solid zone derived in one forward pass — N_s from the full
               h_c section → θ_lim (D.7/Table D.4) → z inverted from D.5 with
               u₃/h_p = 0.75, θ_s = θ_lim (Annex D.3 note 6) → b from D.11.
               Effective compression depth = h_c − b.  Final check: M_Rd,fi
               (reduced section) ≥ M_Ed,fi.  Falls back to direct moment check
               for R30/NWC (no D.4 data).
          D4 – Minimum effective thickness (h_eff ≥ Table D.6 minimum)
          D5 – Field of application        (Table D.7 geometry limits; checked first)

        Note: the ambient composite hogging ULS check (_chk_comp_uls_hogging) uses
        elastic continuous coefficients (Cl 9.4.2(1)) and may upgrade the mesh
        beyond the Cl 9.8.1 minimum.  D3 additionally checks hogging capacity in
        fire (n_spans > 1) and may further upgrade the mesh when fire loads govern.
        """
        R_fi   = self.R_fi
        psi_fi = self.psi_fi
        checks = []

        # ── rib geometry from deck ──────────────────────────────────────────
        l1  = float(self.deck.l_1)   # rib pitch [mm]
        l2  = float(self.deck.l_2)   # bottom rib width [mm]
        l3  = float(self.deck.l_3)   # top rib / upper-flange width [mm]
        h_p = float(self.deck.h_p)   # rib height [mm]
        h_c = float(self.h_c)        # concrete above ribs [mm]
        profile = self.deck.profile_type   # 're-entrant' or 'trapezoidal'

        # concrete type for coefficient tables (NWC / LWC)
        dry_w = self.concrete.dry_weight   # N/m³  (e.g. 24 000 N/m³ for NWC C25/30)
        concrete_type = "LWC" if dry_w < 20000 else "NWC"

        # ── FIRE D.5: Field of application (Table D.7) ──────────────────────
        limits = self._FIRE_D7_LIMITS.get(profile)
        if limits is None:
            checks.append({
                "name":   "D5 - Field of application (Table D.7)",
                "util":   float('inf'),
                "passed": False,
                "note":   f"Unknown profile type '{profile}' — not 'reentrant' or 'trapezoidal'",
            })
            return checks
        params = {
            "l1":  l1,
            "l2":  l2,
            "l3":  l3,
            "h_p": h_p,
            "h_c": h_c,
        }
        in_range = True
        worst_ratio = 0.0
        out_params  = []
        for key, val in params.items():
            lo, hi = limits[key]
            if not (lo <= val <= hi):
                in_range = False
                out_params.append(f"{key}={val:.0f}mm (limits {lo:.0f}–{hi:.0f}mm)")
            # Utilisation: distance from nearest limit / half-range (0 = mid, 1 = at limit)
            half = (hi - lo) / 2.0
            mid  = (hi + lo) / 2.0
            ratio = abs(val - mid) / half if half > 0 else 0.0
            worst_ratio = max(worst_ratio, ratio)

        passed_foa = in_range
        checks.append({
            "name":   "D5 - Field of application (Table D.7)",
            "util":   worst_ratio,
            "passed": passed_foa,
            "note":   ("Outside field of application: " + "; ".join(out_params))
                       if not passed_foa else "Within field of application",
        })
        # Note: a failed D.5 is a warning — the method is technically inapplicable,
        # but we continue with the other checks rather than stopping.

        # ── FIRE D.1: Thermal insulation (formula D.1) ──────────────────────
        ALr   = self._fire_ALr(h_p, l2, l3)                    # mm
        phi   = self._fire_view_factor(h_p, l1, l2, l3)        # -
        t_i   = self._fire_t_i(h_c, phi, ALr, l3, concrete_type)  # min  (h1 = h_c)

        util_ti  = R_fi / t_i if t_i > 0 else float('inf')
        passed_ti = t_i >= R_fi
        checks.append({
            "name":   "D1 - Thermal insulation (formula D.1)",
            "util":   util_ti,
            "passed": passed_ti,
            "note":   (f"t_i={t_i:.1f} min (req {R_fi} min); "
                       f"h1(=h_c)={h_c:.1f} mm, Phi={phi:.3f}, A/Lr={ALr:.1f} mm"),
        })

        # ── FIRE D4: Minimum effective thickness (Table D.6) ────────────────
        # h3 = screed thickness above slab (conservatively taken as 0 here)
        h_eff     = self._fire_h_eff(h_c, h_p, l1, l2, l3)
        h3        = 0.0   # screed layer above slab [mm]; 0 = conservative
        h_eff_min = self._FIRE_D6_MIN_HEFF.get(R_fi)
        if h_eff_min is None:
            checks.append({
                "name":   "D4 - Minimum effective thickness (Table D.6)",
                "util":   float('inf'),
                "passed": False,
                "note":   f"No Table D.6 entry for R_fi={R_fi} min",
            })
        else:
            h_eff_req  = h_eff_min - h3
            util_heff  = h_eff_req / h_eff if h_eff > 0 else float('inf')
            passed_eff = h_eff >= h_eff_req
            checks.append({
                "name":   "D4 - Minimum effective thickness (Table D.6)",
                "util":   util_heff,
                "passed": passed_eff,
                "note":   (f"h_eff={h_eff:.1f} mm (req {h_eff_req:.0f} mm = {h_eff_min:.0f} - h3={h3:.0f} mm); "
                           f"h_c={h_c:.1f} mm, h_p={h_p:.1f} mm"),
            })

        # ── FIRE D.2: Sagging moment capacity M_Rd,fi ≥ M_Ed,fi ────────────
        # Mesh contribution is ignored (crack control only); resistance comes
        # from the deck (if cool enough) and trough bars.

        # Fire moment demand — EN 1990 accidental combination.
        # Propped: composite section carries slab SW + deck SW + SDL + ψ₁Q.
        # Unpropped: slab SW and deck SW were carried by the deck during
        # construction; composite section carries SDL + ψ₁Q only.
        if self.propped:
            G_k = (self.concrete.dry_weight / 1000.0 * self.h_conc_avg / 1000.0
                   + self.deck.weight_kN
                   + (self.finishes_load + self.ceiling_services) / 1e3)
        else:
            G_k = (self.finishes_load + self.ceiling_services) / 1e3
        Q_k = (self.imposed_load + self.partition_load) / 1e3          # N/m² → kN/m²
        M_Ed_fi = self.M_COEFF_SAG * (G_k + psi_fi * Q_k) * span**2   # kNm/m

        # Deck temperatures — lower, web, upper flange (formula D.4, Table D.2 coefficients)
        theta_lower = self._fire_deck_temperature(l3, ALr, phi, R_fi, "lower", concrete_type)
        theta_upper = self._fire_deck_temperature(l3, ALr, phi, R_fi, "upper", concrete_type)

        # Deck contribution — include if lower flange temperature ≤ 350 °C
        N_deck_fi = 0.0
        ky_deck   = self._fire_ky_theta(theta_lower, self._KY_THETA_DECK)
        f_yp_fi   = self.deck.f_ypd * ky_deck   # N/mm²  (f_ypd already design value)
        N_deck_fi = self.deck.A_pe * f_yp_fi     # N/m

        # ── Select minimum trough bar diameter satisfying fire moment ─────────
        fck = self.concrete.fck
        if fck > 1000:    # stored as Pa → convert to N/mm²
            fck = fck / 1e6

        _cover_bottom = 70.0 if h_p > 200.0 else 25.0   # mm to bar centroid start
        TROUGH_BAR_DIAMETERS = [None, 8, 10, 12, 16, 20, 25, 32]   # mm; None = no bar

        chosen_bar_dia  = None
        chosen_N_bar    = 0.0
        chosen_d_bar    = 0.0
        chosen_A_s_bar  = 0.0
        chosen_bar_note = ""
        final_z_fi      = 0.0
        final_M_Rd_fi   = 0.0

        for bar_dia in TROUGH_BAR_DIAMETERS:
            N_bar_fi = 0.0
            bar_note = ""
            d_bar    = 0.0
            A_s_bar  = 0.0

            if bar_dia is not None:
                d_bar    = _cover_bottom + bar_dia / 2.0
                l3_pitch = float(self.deck.l_3)
                A_s_bar  = (math.pi / 4.0 * bar_dia**2) * (1000.0 / l3_pitch)
                u3_bar   = d_bar
                trough_w = l2 + (l1 - l2) * min(u3_bar / h_p, 1.0)
                u1_bar   = max(abs(trough_w) / 2.0, 5.0)
                u3_safe  = max(u3_bar, 5.0)
                z_pos    = self._fire_z_factor(u1_bar, u1_bar, u3_safe)
                theta_bar = self._fire_rebar_temperature(
                    u3_safe, h_p, z_pos, ALr,
                    float(self.deck.theta), l3, R_fi, concrete_type)
                if theta_bar is not None:
                    ky_bar   = self._fire_ky_theta(theta_bar, self._KY_THETA_BAR)
                    fsk_bar  = self._bar_fsk / 1e6
                    N_bar_fi = A_s_bar * fsk_bar * ky_bar
                    bar_note = (f"; bar ø{bar_dia:.0f}mm: "
                                f"T={theta_bar:.0f} degC, ky={ky_bar:.3f}, "
                                f"N={N_bar_fi/1e3:.1f} kN/m")
                else:
                    bar_note = f"; bar ø{bar_dia:.0f}mm: no D.5 data for this R_fi/concrete"

            N_total  = N_deck_fi + N_bar_fi
            x_fi     = N_total / (0.8 * fck * 1000.0) if N_total > 0 else 0.0
            z_deck   = (self.deck.e - 0.4 * x_fi) if N_deck_fi > 0 else 0.0
            z_bar    = ((float(self.h_mm) - d_bar) - 0.4 * x_fi
                        if (N_bar_fi > 0 and d_bar > 0) else 0.0)

            if N_total <= 0:
                continue   # no resistance for this combination, try next diameter

            M_Rd = 0.0
            if z_deck > 0:
                M_Rd += N_deck_fi * z_deck / 1e6
            if z_bar > 0:
                M_Rd += N_bar_fi * z_bar / 1e6

            # store as final values (overwritten each iteration; used if no bar passes)
            z_fi = (N_deck_fi * z_deck + N_bar_fi * z_bar) / N_total if N_total > 0 else 0.0
            final_z_fi    = z_fi
            final_M_Rd_fi = M_Rd
            chosen_bar_note = bar_note

            if M_Rd >= M_Ed_fi:
                chosen_bar_dia  = bar_dia
                chosen_N_bar    = N_bar_fi
                chosen_d_bar    = d_bar
                chosen_A_s_bar  = A_s_bar
                break   # smallest passing diameter found

        # ── Update section with chosen bar; refreshes GWP ──────
        if chosen_bar_dia is not None:
            bar_mass = chosen_A_s_bar * 1e-6 * self.STEEL_DENSITY * self.b
            self.bar_dia_trough  = chosen_bar_dia
            self.A_s_trough      = chosen_A_s_bar
            self.d_bar_trough    = chosen_d_bar
            self.co2_trough_bars = bar_mass * self._rebar_gwp
        else:
            self.bar_dia_trough  = None
            self.A_s_trough      = 0.0
            self.d_bar_trough    = 0.0
            self.co2_trough_bars = 0.0
        self.co2 = (self.co2_concrete + self.co2_deck
                    + self.co2_mesh + self.co2_trough_bars)

        # ── Final check result ───────────────────────────────────────────────────
        if final_M_Rd_fi == 0.0:
            checks.append({
                "name":   "D2 - Sagging moment M_Rd,fi",
                "util":   float('inf'),
                "passed": False,
                "note":   "No fire moment resistance (deck excluded and no bar data available)",
            })
            return checks

        util_m   = M_Ed_fi / final_M_Rd_fi if final_M_Rd_fi > 0 else float('inf')
        passed_m = util_m <= 1.0
        checks.append({
            "name":   "D2 - Sagging moment M_Rd,fi",
            "util":   util_m,
            "passed": passed_m,
            "note":   (f"M_Ed,fi={M_Ed_fi:.2f} kNm/m, M_Rd,fi={final_M_Rd_fi:.2f} kNm/m"),
        })

        # ── FIRE D.3: Hogging moment at internal support (continuous slabs) ──
        # EN 1994-1-2 Annex D.3: reduced cross-section + moment capacity check.
        #
        # Forward pass (no iteration, per Annex D.3 note 6):
        #   1. N_s from full h_c section (b = 0)
        #   2. θ_lim from formula D.7 (Table D.4) using N_s
        #   3. Invert formula D.5 with u₃/h_p = 0.75, θ_s = θ_lim → z
        #   4. b (isotherm penetration into solid zone) from D.11 using z
        #   5. Reduced lever arm with h_c − b; check M_Rd,fi ≥ M_Ed,fi
        #
        # Top mesh is on the unexposed face → η_s = 1.0, γ_M,fi = 1.0 (EN 1992-1-2 §2.4.2).
        if self.n_spans > 1 and self.mesh is not None:
            # Fire moment demand — same G_k / Q_k as D.2
            M_Ed_fi_hog = self.M_COEFF_HOG * (G_k + psi_fi * Q_k) * span ** 2  # kNm/m

            fsk_fi = 500.0   # N/mm²  η_s = 1.0 (cool top face)
            fck_fi = fck     # N/mm²  η_c = 1.0 (top compression zone is cool)

            A_s_m = self.mesh.A_s * self.mesh_layers                        # mm²/m
            d_s_m = self.mesh.d_s if (self.mesh.d_s or 0) > 0 else 6.0  # mm bar dia

            # Rib web inclination α = arctan(2·h_p / (l₁ − l₂))  [degrees]
            alpha_deg = (math.degrees(math.atan2(2.0 * h_p, l1 - l2))
                         if l1 > l2 else 90.0)
            sin_alpha = math.sin(math.radians(alpha_deg))

            # Coefficient rows
            row_d7 = CompositeSlab._FIRE_D4_COEFF.get(concrete_type, {}).get(R_fi)
            row_d5 = CompositeSlab._FIRE_D3_COEFF.get(concrete_type, {}).get(R_fi)

            # ── Step 1: N_s from full (unreduced) section ───────────────────────
            # Use h_c with no isotherm reduction to get an initial N_s, then
            # apply formula D.7 once to obtain θ_lim.  Per Annex D.3 note (6),
            # z is solved directly from D.5 using this θ_lim — no iteration needed.
            # Double mesh centroid is cover + d_s (two layers stacked); single = cover + d_s/2
            axis_dist = self.cover_top + d_s_m / 2.0 * self.mesh_layers
            x_h0  = A_s_m * fsk_fi / (0.8 * self.b_mm * fck_fi)
            d_eff0 = h_c - axis_dist
            z_h0   = max(d_eff0 - 0.4 * x_h0, 0.0)
            N_s0   = M_Ed_fi_hog * 1e6 / z_h0 if z_h0 > 0 else float('inf')

            # ── Step 2: θ_lim from formula D.7 (Table D.4) ─────────────────────
            b_iso     = 0.0
            theta_lim = None
            if row_d7 is not None and math.isfinite(N_s0):
                d0, d1, d2, d3, d4 = row_d7
                theta_lim = d0 + d1*N_s0 + d2*ALr + d3*phi + d4*(1.0/float(l3))

            # ── Step 3: z from D.5, then b from D.11 ────────────────────────────
            # Invert D.5 with u₃/h_p = 0.75 and θ_s = θ_lim to get z directly.
            if theta_lim is not None and row_d5 is not None:
                c0, c1, c2, c3, c4, c5 = row_d5
                sqrt_z = (theta_lim - c0 - 0.75*c1
                          - c3*ALr - c4*alpha_deg - c5*(1.0/l3)) / c2
                if sqrt_z > 0.0:
                    z_code = sqrt_z ** 2          # position factor [mm^0.5]
                    # D.11:  a = (1/z − 1/√h_p)² · l₁ · sinα
                    inv_diff = 1.0 / z_code - 1.0 / math.sqrt(h_p)
                    a_iso    = inv_diff ** 2 * l1 * sin_alpha
                    if a_iso > 0.0:
                        # c = −8(1+√(1+a)) for a ≥ 8;  +8(1+√(1+a)) for a < 8
                        c_iso = ((-8.0 if a_iso >= 8.0 else +8.0)
                                 * (1.0 + math.sqrt(1.0 + a_iso)))
                        disc  = a_iso**2 - 4.0*a_iso + c_iso
                        if disc >= 0.0:
                            b_iso = max(
                                0.5 * l1 * sin_alpha * (1.0 - math.sqrt(disc) / a_iso),
                                0.0)

            # ── Step 4: reduced section lever arm and final check ───────────────
            d_eff_h = (h_c - b_iso) - axis_dist

            if d_eff_h <= 0.0:
                checks.append({
                    "name":   "D3 - Hogging moment (reduced section, formula D.7)",
                    "util":   float('inf'),
                    "passed": False,
                    "note":   (f"Isotherm consumes solid zone: "
                               f"h_c−b = {h_c:.0f}−{b_iso:.1f} = {h_c - b_iso:.1f} mm"
                               f" ≤ axis dist {axis_dist:.0f} mm"),
                })
            elif row_d7 is None:
                # No D.7 data (e.g. R30/NWC) → direct moment check, reduced section
                x_h = A_s_m * fsk_fi / (0.8 * self.b_mm * fck_fi)
                z_h_fin = max(d_eff_h - 0.4 * x_h, 0.0)
                M_Rd_fi_hog = A_s_m * fsk_fi * z_h_fin / 1e6
                util_h   = M_Ed_fi_hog / M_Rd_fi_hog if M_Rd_fi_hog > 0 else float('inf')
                passed_h = M_Rd_fi_hog >= M_Ed_fi_hog
                checks.append({
                    "name":   "D3 - Hogging moment (reduced section)",
                    "util":   util_h,
                    "passed": passed_h,
                    "note":   (f"No D.7 data for R{R_fi}/{concrete_type}; "
                               f"M_Ed,fi={M_Ed_fi_hog:.2f} kNm/m, "
                               f"M_Rd,fi={M_Rd_fi_hog:.2f} kNm/m "
                               f"[{self.mesh.mech_prop}, b={b_iso:.1f} mm, "
                               f"h_c−b={h_c - b_iso:.0f} mm, d_eff={d_eff_h:.0f} mm]"),
                })
            else:
                x_h = A_s_m * fsk_fi / (0.8 * self.b_mm * fck_fi)
                z_h_fin = max(d_eff_h - 0.4 * x_h, 0.0)
                # D.3(2/3): remaining section (above θ_lim isotherm) at room temperature.
                # Check is M_Rd,fi ≥ M_Ed,fi — same as fallback path.
                M_Rd_fi_hog = A_s_m * fsk_fi * z_h_fin / 1e6  # kNm/m
                util_h   = M_Ed_fi_hog / M_Rd_fi_hog if M_Rd_fi_hog > 0 else float('inf')
                passed_h = M_Rd_fi_hog >= M_Ed_fi_hog
                checks.append({
                    "name":   "D3 - Hogging moment (reduced section, formula D.7)",
                    "util":   util_h,
                    "passed": passed_h,
                    "note":   (f"M_Ed,fi={M_Ed_fi_hog:.2f}, M_Rd,fi={M_Rd_fi_hog:.2f} kNm/m; "
                               f"b_iso={b_iso:.1f} mm, h_c−b={h_c - b_iso:.0f} mm, "
                               f"d_eff={d_eff_h:.0f} mm, z={z_h_fin:.0f} mm "
                               f"[{self.mesh.mech_prop}]"),
                })

        return checks



#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
class Section:
    # contains fundamental section properties like section type weight, resistance and stiffness
    def __init__(self, section_type):
        self.section_type = section_type
        # The following properties are defined in the specific cross-section classes. However, it could make sense to
        # provide them in this more general parent class.
        #
        # properties:
        # self.mu_max = float
        # self.mu_min = float
        # self.vu = float
        # self.qs_class_n = int
        # self.qs_class_p = int
        # self.g0k = float
        # self.ei1 = float
        # self.co2 = float
        # self.cost = float


class SupStrucRectangular(Section):
    # defines cross-section dimensions and has methods to calculate static properties of rectangular,
    # non-cracked sections
    def __init__(self, section_type, b, h, phi=0):  # create a rectangular object
        super().__init__(section_type)
        self.b = b  # width [m]
        self.h = h  # height [m]
        self.a_brutt = self.calc_area()
        self.iy = self.calc_moment_of_inertia()
        self.phi = phi

    def calc_area(self):
        #  in: width b [m], height h [m]
        #  out: area [m^2]
        a_brutt = self.b * self.h
        return a_brutt

    def calc_moment_of_inertia(self):
        #  in: width b [m], height h [m]
        #  out: second moment of inertia Iy [m^4]
        iy = self.b * self.h ** 3 / 12
        return iy

    def calc_strength_elast(self, fy, ty):
        #  in: yielding strength fy [Pa], shear strength ty [Pa]
        #  out: elastic bending resistance mu_el [Nm], elastic shear resistance vu_el [N]
        mu_el = self.iy * fy * 2 / self.h
        vu_el = self.b * self.h * ty / 1.5
        return mu_el, vu_el

    def calc_strength_plast(self, fy, ty):
        #  in: yielding strength fy [Pa], shear strength ty [Pa]
        #  out: plastic bending resistance mu_pl [Nm], plastic shear resistance vu_pl [N]
        mu_pl = self.b * self.h ** 2 * fy / 4
        vu_pl = self.b * self.h * ty
        return mu_pl, vu_pl

    def calc_weight(self, spec_weight):
        #  in: specific weight [N/m^3]
        #  out: weight of cross section per m length [N/m]
        w = spec_weight * self.a_brutt # w: spec_weight
        return w

#........................................................................
class RectangularWood(SupStrucRectangular, Section):
    # defines properties of rectangular, wooden cross-section
    def __init__(self, wood_type, b, h, phi=0.6, xi=0.01, ei_b=0.0):  # create a rectangular timber object
        section_type = "wd_rec"
        super().__init__(section_type, b, h, phi)
        self.wood_type = wood_type
        mu_el, vu_el = self.calc_strength_elast(wood_type.fmd, wood_type.fvd)
        self.mu_max, self.mu_min = [mu_el, -mu_el]  #Readme: Why is this needed for wood? -> is not needed for wood.
        # However, as the same resistance values should be provided for all cross-sections, I defined them for both
        # directions for wood too
        self.vu_p, self.vu_n = vu_el, vu_el
        self.qs_class_n, self.qs_class_p = [3, 3]  # Required cross-section class: 1:PP, 2:EP, 3:EE
        self.g0k = self.calc_weight(wood_type.weight) # dead weight of cross section [N/m]
        self.ei1 = self.wood_type.Emmean * self.iy  # elastic stiffness [Nm^2]
        self.co2 = self.a_brutt * self.wood_type.GWP * self.wood_type.density  # [kg_CO2_eq/m]
        self.cost = self.a_brutt * self.wood_type.cost # [CHF] #TODO: Not implemented yet!
        self.ei_b = ei_b  # stiffness perpendicular to direction of span [Nm^2]
        self.xi = xi  # damping factor, preset value see: HBT, Page 47 (higher value for some buildups possible)

    @staticmethod
    def fire_resistance(member):
        bnds = [(0, 240)]   #Randbedingungen für Definition Brand - mind. 0 min max. 240 min
        t0 = 60     #Brandeinwirkungsdauer
        max_t = minimize(RectangularWood.fire_minimizer, t0, args=[member], bounds=bnds)    #Brandwiderstanddauer → maximale Brandeinwirkungszeit
        t_max = max_t.x[0]
        return t_max

    @staticmethod
    def fire_minimizer(t, args):
        member = args[0]
        rem_sec = RectangularWood.remaining_section(member.section, member.fire, t)
        mu_fire = 1.8 * rem_sec.mu_max
        vu_fire = 1.8 * rem_sec.vu_p  # SIA 265 (51)
        qd_fire = member.psi[2] * member.qk + member.gk
        qd_fire_zul = min(mu_fire / (max(member.system.alpha_m) * member.system.l_tot ** 2),
                          vu_fire / (max(member.system.alpha_v) * member.system.l_tot))
        to_opt = abs(qd_fire - qd_fire_zul)
        return to_opt

    @staticmethod
    def remaining_section(section, fire, t=60, dred=0.007):
        betan = section.wood_type.burn_rate
        dcharn = betan * t
        d_ef = dcharn + dred
        h_fire = max(section.h - d_ef * (fire[0] + fire[2]))
        b_fire = max(section.b - d_ef * (fire[1] + fire[3]), 0)
        rem_sec = RectangularWood(section.wood_type, b_fire, h_fire)
        return rem_sec

# ........................................................................
class RectangularConcrete(SupStrucRectangular):
    # defines properties of rectangular, reinforced concrete cross-section
    def __init__(self, concrete_type, rebar_type, b, h, di_xu, s_xu, di_xo, s_xo, di_yu=0.01, s_yu=0.15, di_yo=0.01, s_yo=0.15, di_bw=0.0, s_bw=0.15, n_bw=0,
                 phi=2.0, c_nom=0.02, xi=0.02, jnt_srch=0.15):
        # create a rectangular concrete object
        section_type = "rc_rec"
        super().__init__(section_type, b, h, phi)
        self.concrete_type = concrete_type
        self.rebar_type = rebar_type
        self.c_nom = c_nom #Bewehrungsüberdeckung
        self.bw = [[di_xu, s_xu], [di_xo, s_xo], [di_yu, s_yu],[di_yo, s_yo]] #Definition Biegebewehrung 4-Lagig. x-Richtung ist dabei die Haupttragrichtung, di = Durchmesser, s = Abstand, u = untere Lagen (positives Biegemoment), o = obere Lagen (negatives Biegemoment)
        self.bw_bg = [di_bw, s_bw, n_bw] #Definition Querkraftbewehrung
        mr = self.b * self.h ** 2 / 6 * 1.3 * self.concrete_type.fctm  #cracking moment
        self.mr_p, self.mr_n = mr, -mr #mr_p: positives Rissmoment, mr_n: negatives Rissmoment
        [self.d, self.ds] = self.calc_d() #Statische Höhe. d für positive Biegung (untere Lagen), ds für negative Biegung (obere Lagen)
        #TODO: x und y Richtung Berücksichtigen
        [self.mu_max, self.x_p, self.as_p, self.qs_class_p] = self.calc_mu('pos')
        [self.mu_min, self.x_n, self.as_n, self.qs_class_n] = self.calc_mu('neg')
        self.roh, self.rohs = self.as_p / self.d, self.as_n / self.ds
        [self.vu_p, self.vu_n, self.as_bw] = self.calc_shear_resistance()
        self.g0k = self.calc_weight(concrete_type.weight)
        a_s_stat = self.as_p + self.as_n + self.as_bw  # rebar area without reinforcement joint surcharge
        self.joint_surcharge = jnt_srch  # joint surcharge
        #TODO: a_s_tot: Hat erst 1. & 4. Lage drin
        a_s_tot = a_s_stat * (1 + self.joint_surcharge)  # rebar area without reinforcement joint surcharge
        #TODO: Im GWP vom Gesamtquerschnitt müssen alle 4 Bewehrungslagen berücksichtigt werden! Prüfen ob alle 4 Lagen berücksichtigt werden
        co2_rebar = a_s_tot * self.rebar_type.GWP * self.rebar_type.density  # [kg_CO2_eq/m]
        co2_concrete = (self.a_brutt - a_s_tot) * self.concrete_type.GWP * self.concrete_type.density  # [kg_CO2_eq/m]
        self.ei1 = self.concrete_type.Ecm * self.iy  # elastic stiffness concrete (uncracked behaviour) [Nm^2]
        self.co2 = (co2_rebar + co2_concrete)
        self.mass = ((self.a_brutt - a_s_tot) * self.concrete_type.density
                     + a_s_tot * self.rebar_type.density)  # [kg/m²] concrete + rebar
        self.cost = (a_s_tot * self.rebar_type.cost + (self.a_brutt - a_s_tot) * self.concrete_type.cost
                     + self.concrete_type.cost2)
        self.ei_b = self.ei1
        self.xi = xi  # XXXXXXX preset value is an assumption. Has to be verified with literature. XXXXXXX
        self.ei2 = self.ei1 / self.f_w_ger(self.roh, self.rohs, 0, self.h, self.d)

    def calc_d(self):
        d = self.h - self.c_nom - self.bw_bg[0] - self.bw[0][0] / 2 #Statische Höhe für Positives Biegemoment
        ds = self.h - self.c_nom - self.bw_bg[0] - self.bw[1][0] / 2 #Statische Höhe für Negatives Biegemoment
        if d <= 0 or ds <= 0:
            print("d of ds<=0. Cross-section is not valid")
        return d, ds

    def calc_mu(self, sign='pos'):
        #in: self
        #out: Biegewiderstand mu [Nm], Druckzonenhöhe x [m], Bewehrungsfläche a_s [m2], Querschnittsklasse qs_klasse []
        b = self.b  #Querschnittsbreite
        fsd = self.rebar_type.fsd
        fcd = self.concrete_type.fcd
        if sign == 'pos':
            [mu, x, a_s, qs_klasse] = self.mu_unsigned(self.bw[0][0], self.bw[0][1], self.d, b, fsd, fcd, self.mr_p)
        elif sign == 'neg':
            [mus, x, a_s, qs_klasse] = self.mu_unsigned(self.bw[1][0], self.bw[1][1], self.ds, b, fsd, fcd, self.mr_n)
            mu = -mus
        else:
            [mu, x, a_s, qs_klasse] = [0, 0, 0, 0]
            print("sign of moment resistance has to be 'neg' or 'pos'")
        return mu, x, a_s, qs_klasse

    @staticmethod
    def mu_unsigned(di, s, d, b, fsd, fcd, mr):
        #in: Bewehrungsdurchmesser di, Bewehrungsabstand s, Statische Höhe d, fsd, fcd, mr
        #out: mu, x, a_s, qs_klasse
        # units input: [m, m, m, m, N/m^2, N/m^2]
        a_s = np.pi * di ** 2 / (4 * s) * b  # [m^2]
        omega = a_s * fsd / (d * b * fcd)  # [-]
        mu = a_s * fsd * d * (1 - omega / 2)  # [Nm]
        x = omega * d / 0.85  # [m]
        if x / d <= 0.35 and mu >= mr:
            return mu, x, a_s, 1
        elif x / d <= 0.5 and mu >= mr:
            return mu, x, a_s, 2
        else:  # zero resistance for x/d>0.5
            epsilon = 1.0e-3
            shift = 0.5
            factor = 1 - 0.5 * (1 + 2 / np.pi * np.arctan((x/d - shift) / epsilon)) #irgendein Faktor, um die Funktion richtig auf 0 gehen zu lassen. Ist keine Formel aus irgendeiner Norm o.Ä., hat auch nichts mit der Statik zu tun#
            return factor*mu, x, a_s, 99  # Querschnitt hat ungenügendes Verformungsvermögen

    def calc_shear_resistance(self, d_installation=0.0):
        # in: self
        # out: Querkraftwiderstand positiv vu_p [N], Querkraftwiderstand negativ vu_n [N], Querkraftbewehrung as_bw [m2]
        #TODO: Anpassung an die SIA 262 (2025)! Ist noch gemäss alter Norm!
        di = self.bw_bg[0]      # diameter
        s = self.bw_bg[1]       # spacing
        n = self.bw_bg[2]       # number of stirrups per spacing
        fck = self.concrete_type.fck        #SIA 262
        fcd = self.concrete_type.fcd        #SIA 262
        tcd = self.concrete_type.tcd        #SIA 262
        dmax = self.concrete_type.dmax      #dmax in mm
        fsk = self.rebar_type.fsk           #SIA 262
        fsd = self.rebar_type.fsd           #SIA 262
        es = self.rebar_type.Es             #SIA 262
        bw = self.b         #Stegbreite
        d = self.d          #Statische Höhe für positives Biegemoment (untere Lagen)
        ds = self.ds        #Statische Höhe für negatives Biegemoment (obere Lagen)
        x_p = self.x_p      #Druckzonenhöhe positives Biegemoment (obere Querschnittsrand)
        x_n = self.x_n      #Druckzonenhöhe negatives Biegemoment (unterer Querschnittsrand)
        as_bw = self.calc_as_bw(di, n, s, d)
        if d_installation < d / 6:
            dv_p = d                    #Wirksame statische Höhe für Querkraft
        else:
            dv_p = d - d_installation   #Wirksame statische Höhe für Querkraft
        if d_installation < ds / 6:
            dv_n = ds                   #Wirksame statische Höhe für Querkraft
        else:
            dv_n = ds - d_installation  #Wirksame statische Höhe für Querkraft
        vu_p = self.vu_unsigned(bw, di, n, s, as_bw, d, dv_p, x_p, fck, fcd, tcd, fsk, fsd, es, dmax)   #Positiver Querkraftwiderstand [N]
        vu_n = self.vu_unsigned(bw, di, n, s, as_bw, ds, dv_n, x_n, fck, fcd, tcd, fsk, fsd, es, dmax)  #Negativer Querkraftwiederstand [N]
        return vu_p, vu_n, as_bw

    @staticmethod
    def calc_as_bw(di, n, s, d):
        #in: Bewehrungsduchmesser di [m], Anzahl Stäbe n [], Bewehrungsabstand s [m], Statische Höhe d [m]
        #out: Bewehrungsquerschnittsfläche Querkraftbewehrung as_bw [mm2]
        as_bw = np.pi * di ** 2 / 4 * n / s * 0.9*d #ToDo: muss die Bügelquerschnittsfläche nicht noch mit der Plattenstärke multipliziert werden?
        return as_bw

    @staticmethod
    def vu_unsigned(bw, di, n, s, as_bw, d, dv, x, fck, fcd, tcd, fsk, fsd, es, dmax=32, alpha=np.pi / 4, kc=0.55):
        rohw = as_bw / min(bw, 0.4)  # SIA 262, Zif. 5.5.2.2
        rohw_min = 0.001 * (fck * 1e-6 / 30) ** 0.5 * 500 / (fsk * 1e-6)  # SIA 262, Zif. 5.5.2.2
        s_max = 25*s  # SIA262, Zif. 5.5.2.2
        if bw < 0.5:  # SIA262, Zif. 5.5.2.3
            n_min = 2
        else:
            n_min = 4
        if rohw < rohw_min or s > s_max or n < n_min:  # cross-section resistance without stirrups
            ev = 1.5 * fsd / es         #SIA 262
            kg = 48 / (16 + dmax)       #SIA 262
            kd = 1 / (1 + ev * d * kg)  #SIA 262
            vrd = kd * tcd * dv
            return vrd #Querkraftwiderstand OHNE Querkraftbewehrung SIA 262
        else:  # cross-section resistance with vertical stirrups
            z = d - 0.85 * x / 2
            vrds = as_bw * z * fsd
            vrdc = bw * z * kc * fcd * np.sin(alpha) * np.cos(alpha)  # unit of alpha: [rad]
            return min(vrds, vrdc) #Querkraftwiderstand MIT Querkraftbewehrung SIA 262

    @staticmethod
    def f_w_ger(roh, rohs, phi, h, d):
        f = (1 - 20 * rohs) / (10 * roh ** 0.7) * (0.75 + 0.1 * phi) * (h / d) ** 3
        #TODO: Prüfen, ob dieser Wert nicht zu konservativ ist! Als Abschätzung für die Vordimensionierung scheint der Wert jedoch schon i.O., ist zumindest nicht komplett willkürlich.
        return f

    @staticmethod
    def fire_resistance(section):
        # fire resistance of 1-D load-bearing plates according to SIA 262, Tab.16
        c_nom = section.c_nom
        h = section.h
        b = section.b
        if c_nom >= 0.04 and h >= 0.15 and b >= 0.4:
            resistance = 180
        elif c_nom >= 0.03 and h >= 0.12 and b >= 0.3:
            resistance = 120
        elif c_nom >= 0.03 and h >= 0.1 and b >= 0.2:
            resistance = 90
        elif c_nom >= 0.02 and h >= 0.08 and b >= 0.15:
            resistance = 60
        elif c_nom >= 0.02 and h >= 0.06 and b >= 0.1:
            resistance = 30
        else:
            resistance = 0
        return resistance


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
class SupStrucRibbedConcrete(Section):
    def __init__(self, section_type, b, b_w, h, h_f, l0, phi=0):
        super().__init__(section_type)
        self.b = b              # flange width [m] (Abstand Rippenachse-Rippenachse)
        self.b_w = b_w          # web width [m]
        self.h = h              # total height [m]
        self.h_f = h_f          # flange height [m]
        self.h_w = h - h_f      # web height [m]
        self.l0 = l0            # Abstand Momentennullpunkte [m]
        self.b_eff = self.calc_beff()               #Effective width [m]
        self.a_brutt = self.calc_area()             #Bruttoquerschnittsfläche [m2]
        self.z_s = self.calc_center_of_gravity()    #center of gravity [m]
        self.iy = self.calc_moment_of_inertia()     #moment of inertia [m4]
        self.w = self.calc_weight()                 #Eigengewicht [n/m]
        self.phi = phi                              #Kriechzahl

    def calc_area(self):
        # in: width b and bw [m], height h and h_f[m]
        # out: area [m2]
        a_brutt = self.b * self.h_f + self.b_w * (self.h - self.h_f)
        return a_brutt

    def calc_beff(self):
        # in: width b [m], bw [m], l_0 [m]
        # out: effective width b_eff
        l_0 = self.l0
        b_eff_i = 0.2 * (self.b - self.b_w) / 2 + 0.1 * l_0  # SIA 262, 4.1.3.3.2 (20)
        if b_eff_i > 0.2 * l_0:
            b_eff_i = 0.2 * l_0
        else:
            pass
        b_eff = 2 * b_eff_i + self.b_w  # SIA 262, 4.1.3.3.2 (19)
        if b_eff > self.b:
            b_eff = self.b
        else:
            pass
        return b_eff

    def calc_center_of_gravity(self):
        # in: Geometry effective width b_eff [m], slab height h_f [m], rib width b_w [m], rib height h_w [m]
        # out: center of gravity z_s [m], z = 0: OK Slab
        z_s = (self.b_eff * self.h_f ** 2 / 2 + self.b_w * self.h_w * (self.h_f +self.h_w/2)) / (
                    self.b_eff * self.h_f + self.b_w * self.h_w)
        return z_s

    def calc_moment_of_inertia(self):
        # in: Geometry effective width b_eff [m], slab height h_f [m], rib width b_w [m], rib height h_w [m], center of gravity z_s [m]
        # out: moment of inertia I_y [m^4]
        i_01 = self.b_eff * self.h_f ** 3 / 12
        as_01 = self.b_eff * self.h_f * abs(self.z_s - self.h_f / 2) ** 2
        i_02 = self.b_w * self.h_w ** 3 / 12
        as_02 = self.b_w * self.h_w * abs(self.z_s - (self.h_f + self.h_w/2)) ** 2
        iy = i_01 + i_02 + as_01 + as_02
        return iy

    def calc_weight(self,
                    spec_weight=25):  #README: Spec-Weight muss automatisch aus Tabelle eingelesen werden können! Ergänzen!
        #  in: specific weight [N/m^3]
        #  out: weight of cross section per m length [N/m]
        w = spec_weight * self.a_brutt
        return w


#.....................................................................................
class RibbedConcrete(SupStrucRibbedConcrete):
    #defines properties of a rectangular, reinforced concrete section
    #di_x_w, n_x_w = diameter and number of longitudinal reinforcement in rib
    def __init__(self, concrete_type, rebar_type, l0, b, b_w, h, h_f, di_xu, s_xu, di_xo, s_xo, di_x_w, n_x_w,
                 di_pb_bw, s_pb_bw, n_pb_bw=2,
                 phi=2.0, c_nom=0.03, xi=0.02, jnt_srch=0.15):
        section_type = "rc_rib"
        super().__init__(section_type, b, b_w, h, h_f, l0, phi)
        self.concrete_type = concrete_type
        self.rebar_type = rebar_type
        self.c_nom = c_nom
        self.bw = [[di_xu, s_xu], [di_xo, s_xo]]  # Slab reinforcement
        self.bw_bg = [0, 0.15, 0]  # Allow for no slab shear reinforcement
        self.bw_r = [di_x_w, n_x_w]  # Longitudinal reinforcement in rib
        self.bw_bg_r = [di_pb_bw, s_pb_bw, n_pb_bw]  # Shear reinforcement in rib
        mr_slab = self.b * self.h ** 2 / 6 * 1.3 * self.concrete_type.fctm  # cracking moment
        mr_pb = self.iy / (self.h - self.z_s) * 1.3 * self.concrete_type.fctm  # cracking moment
        self.mr_p, self.mr_n = mr_slab, -mr_slab
        self.mr_pb_p = mr_pb
        self.mr_pb_n = -mr_pb
        [self.d, self.ds, self.d_PB, self.ds_PB] = self.calc_d()
        [self.mu_max_slab, self.x_p, self.as_p, self.qs_class_p_slab] = self.calc_mu('pos')
        [self.mu_min_slab, self.x_n, self.as_n, self.qs_class_n_slab] = self.calc_mu('neg')
        [self.mu_max, self.x_PB_p, self.as_PB_p, self.qs_class_p] = self.calc_mu_pb('pos')
        [self.mu_min, self.x_PB_n, self.as_PB_n, self.qs_class_n] = self.calc_mu_pb('neg')
        self.roh_slab, self.rohs, self.roh = self.as_p / self.d, self.as_n / self.ds, self.as_PB_p / self.d_PB
        [self.vu_p, self.vu_n, self.as_bw] = self.calc_shear_resistance('Platte')  #Platte "Querrichtung"
        [self.vu_PB_p, self.vu_PB_n, self.as_PB_bw] = self.calc_shear_resistance(
            'Plattenbalken')  #Rippe Plattenbalken "Längsrichtung"
        self.g0k = self.calc_weight(concrete_type.weight)
        a_s_stat = self.as_p + self.as_n + self.as_bw + self.as_PB_p + self.as_PB_n + self.as_PB_bw
        #TODO: Achtung - es fehlt die Spreizbewehrung
        self.joint_surcharge = jnt_srch
        a_s_tot = a_s_stat * (1 + self.joint_surcharge)
        co2_rebar = a_s_tot * self.rebar_type.GWP * self.rebar_type.density  # [kg_CO2_eq/m]
        co2_concrete = (self.a_brutt - a_s_tot) * self.concrete_type.GWP * self.concrete_type.density  # [kg_CO2_eq/m]
        self.ei1 = self.concrete_type.Ecm * self.iy  # elastic stiffness concrete (uncracked behaviour) [Nm^2]
        self.co2 = (co2_rebar + co2_concrete)/self.b
        self.mass = ((self.a_brutt - a_s_tot) * self.concrete_type.density
                     + a_s_tot * self.rebar_type.density) / self.b  # [kg/m²] concrete + rebar
        self.cost = (a_s_tot * self.rebar_type.cost + (self.a_brutt - a_s_tot) * self.concrete_type.cost
                     + self.concrete_type.cost2)
        self.ei_b = self.ei1  #!!!!!!!ANPASSEN AUF PB
        self.xi = xi  # XXXXXXX preset value is an assumption. Has to be verified with literature. XXXXXXX
        self.ei2 = self.ei1 / self.f_w_ger(self.roh, self.rohs, 0, self.h, self.d_PB)  #!!!!!ANPASSEN AUF PB

    def calc_d(self):
        d = self.h_f - self.c_nom - self.bw[0][0] / 2  # Statische Höhe 1. Lage Platte
        ds = self.h_f - self.c_nom - self.bw[1][0] / 2  # Statische Höhe 4. Lage Platte
        d_PB = self.h - self.c_nom - self.bw_bg_r[0] - self.bw_r[
            0] / 2  # Nur eine Lage Längsbewehrung implementiert. ACHTUNG: Check implementieren, ob genug Platz für Längsbewehrung vorhanden!!
        ds_PB = self.h - self.c_nom - self.bw[1][0]  # Mittlere statische Höhe 3./4. Lage Platte
        return d, ds, d_PB, ds_PB

    #Slab = Platte in Querrichtung. ACHTUNG: DURCHLAUFWIRKUNG MUSS NOCH IMPLEMENTIERT WERDEN!
    #Kann man die Berechnung der Platte zusammenführen mit Rectangular Concrete?

    def calc_mu(self, sign='pos'):
        # calculates moment resistence of slab
        b = 1
        fsd = self.rebar_type.fsd
        fcd = self.concrete_type.fcd
        if sign == 'pos':
            [mu, x, a_s, qs_klasse] = self.mu_unsigned(self.bw[0][0], self.bw[0][1], self.d, b, fsd, fcd, self.mr_p)
        elif sign == 'neg':
            [mus, x, a_s, qs_klasse] = self.mu_unsigned(self.bw[1][0], self.bw[1][1], self.ds, b, fsd, fcd, self.mr_n)
            mu = -mus
        else:
            [mu, x, a_s, qs_klasse] = [0, 0, 0, 0]
            print("sign of moment resistance has to be 'neg' or 'pos'")

        return mu, x, a_s, qs_klasse

    def calc_mu_pb(self, sign='pos'):
        # calculates moment resistence of Plattenbalken = PB
        fsd = self.rebar_type.fsd
        fcd = self.concrete_type.fcd
        if sign == 'pos':
            [mu_PB, x, a_s, qs_klasse] = self.mu_unsigned_PB(self.bw_r[0], self.bw_r[1], self.d_PB, self.b_eff,
                                                             self.h_f, fsd, fcd, self.mr_pb_p)
        elif sign == 'neg':
            [mus_PB, x, a_s, qs_klasse] = self.mu_unsigned(self.bw[1][0], self.bw[1][1], self.ds_PB, self.b_w, fsd, fcd,
                                                           self.mr_pb_n)
            mu_PB = - mus_PB
        else:
            [mu_PB, x, a_s, qs_klasse] = [0, 0, 0, 0]
            print("sign of moment resistance has to be 'neg' or 'pos'")

        return mu_PB, x, a_s, qs_klasse

    @staticmethod
    def mu_unsigned(di, s, d, b, fsd, fcd, mr):
        # units input: [m, m, m, m, N/m^2, N/m^2]
        a_s = np.pi * di ** 2 / (4 * s) * b  # [m^2]
        omega = a_s * fsd / (d * b * fcd) # [-]
        mu = a_s * fsd * d * (1 - omega / 2)  # [Nm]
        x = omega * d / 0.85  # [m]
        if x / d <= 0.35 and mu >= mr:
            return mu, x, a_s, 1
        elif x / d <= 0.5 and mu >= mr:
            return mu, x, a_s, 2
        else:
            return mu, x, a_s, 99  # Querschnitt hat ungenügendes Verformungsvermögen

    @staticmethod
    def mu_unsigned_PB(di, n, d, b, h_f, fsd, fcd, mr):
        a_s = np.pi * di ** 2 / 4 * n  # [m^2]
        omega = a_s * fsd / (d * b * fcd)  #[-]
        mu_PB = a_s * fsd * d * (1 - omega / 2)  # [Nm]
        x = omega * d / 0.85
        if x > h_f:
            #print("Druckzonenhöhe > Plattenhöhe")
            mu_PB = 0
            return mu_PB, x, a_s, 99
        else:
            pass

        if x / d <= 0.35 and mu_PB >= mr:
            return mu_PB, x, a_s, 1
        elif x / d <= 0.5 and mu_PB >= mr:
            return mu_PB, x, a_s, 2
        else:
            return mu_PB, x, a_s, 99  # Querschnitt hat ungenügendes Verformungsvermögen

    def calc_shear_resistance(self, bauteil='Platte', d_installation=0.0):
        # calculates shear resistance with d
        di_r = self.bw_bg_r[0]  # diameter
        s_r = self.bw_bg_r[1]  # spacing
        n_r = self.bw_bg_r[2]  # number of stirrups per spacing
        fck = self.concrete_type.fck
        fcd = self.concrete_type.fcd
        tcd = self.concrete_type.tcd
        dmax = self.concrete_type.dmax  # dmax in mm
        fsk = self.rebar_type.fsk
        fsd = self.rebar_type.fsd
        es = self.rebar_type.Es
        bw = self.b
        b_w = self.b_w
        d, d_PB = self.d, self.d_PB
        ds, ds_PB = self.ds, self.ds_PB
        x_p, x_PB_p = self.x_p, self.x_PB_p
        x_n, x_PB_n = self.x_n, self.x_PB_n
        as_bw = 0
        as_PB_bw = np.pi * di_r ** 2 / 4 * n_r / s_r * 0.9 * d

        if bauteil == 'Platte':
            if d_installation < d / 6:  #SIA 262 4.3.3.2.8
                dv_p = d
            else:
                dv_p = d - d_installation
            if d_installation < ds / 6:
                dv_n = ds
            else:
                dv_n = ds - d_installation

            vu_p = self.vu_unsigned(bw, as_bw, d, dv_p, x_p, fck, fcd, tcd, fsk, fsd, es, dmax)
            vu_n = self.vu_unsigned(bw, as_bw, ds, dv_n, x_n, fck, fcd, tcd, fsk, fsd, es, dmax)

            return vu_p, vu_n, as_bw

        else:
            if d_installation < d_PB / 6:  #SIA 262 4.3.3.2.8
                dv_PB_p = d_PB
            else:
                dv_PB_p = d_PB - d_installation
            if d_installation < ds_PB / 6:
                dv_PB_n = ds_PB
            else:
                dv_PB_n = ds_PB - d_installation

            vu_PB_p = self.vu_unsigned(b_w, as_PB_bw, d_PB, dv_PB_p, x_PB_p, fck, fcd, tcd, fsk, fsd, es, dmax)
            vu_PB_n = self.vu_unsigned(b_w, as_PB_bw, ds_PB, dv_PB_n, x_PB_n, fck, fcd, tcd, fsk, fsd, es, dmax)
            return vu_PB_p, vu_PB_n, as_PB_bw

    @staticmethod
    def vu_unsigned(bw, as_bw, d, dv, x, fck, fcd, tcd, fsk, fsd, es, dmax=32, alpha=np.pi / 4, kc=0.55):
        if as_bw == 0:  # cross-section without stirrups
            ev = 1.5 * fsd / es  # SIA 262, 4.3.3.2.2, (39)
            kg = 48 / (16 + dmax)  # SIA 262, 4.3.3.2.1, (37)
            kd = 1 / (1 + ev * d * kg)  # SIA 262, 4.3.3.2.1, (36)
            vrd = kd * tcd * dv  # SIA 262, 4.3.3.2.1, (35)
            return vrd
        else:  # cross-section with vertical stirrups
            z = d - 0.85 * x / 2
            vrds = as_bw * z * fsd  # SIA 262, 4.3.3.4.3, (43)
            vrdc = bw * z * kc * fcd * np.sin(alpha) * np.cos(
                alpha)  # unit of alpha: [rad]    # SIA 262, 4.3.3.4.6, (45)
            rohw = as_bw / bw /(0.9*d)
            rohw_min = 0.001 * (fck * 1e-6 / 30) ** 0.5 * 500 / (fsk * 1e-6)
            if rohw < rohw_min:
                print("minimal reinforcement ratio of stirrups is lower than required according to SIA 262, (110)")
            return min(vrds, vrdc)

    #ÜBERNOMMEN VON RECHTECK-QS, NICHT ANGEPASST
    @staticmethod
    #SIA 262, 4.4.3.2.5: Annahme für den vollständig gerissenen Zustand
    def f_w_ger(roh, rohs, phi, h, d):
        f = (1 - 20 * rohs) / (10 * roh ** 0.7) * (0.75 + 0.1 * phi) * (h / d) ** 3
        return f

    @staticmethod
    def fire_resistance(section):
        # fire resistance of 1-D load-bearing plates according to SIA 262, Tab.16
        c_nom = section.c_nom
        h = section.h
        b = section.b
        if c_nom >= 0.04 and h >= 0.15 and b >= 0.4:
            resistance = 180
        elif c_nom >= 0.03 and h >= 0.12 and b >= 0.3:
            resistance = 120
        elif c_nom >= 0.03 and h >= 0.1 and b >= 0.2:
            resistance = 90
        elif c_nom >= 0.02 and h >= 0.08 and b >= 0.15:
            resistance = 60
        elif c_nom >= 0.02 and h >= 0.06 and b >= 0.1:
            resistance = 30
        else:
            resistance = 0
        return resistance


# .....................................................................................
class SupStrucRibWood(Section):
    def __init__(self, section_type, b, h, a, t2, t3, n, n_inf):
        super().__init__(section_type)
        self.b = b  # rib width [m]
        self.h = h  # rib height [m]
        self.a = a  # spacing between ribs [m]
        self.t2 = t2  # slab height bottom flange [m]
        self.t3 = t3  # slab height top flange [m]
        self.bc_ef = self.calc_bef('comp') + b  # Effective width top flange compression [m]
        self.bt_ef = self.calc_bef('tens') + b  # Effective width bottom flange tension [m]
        self.a_brutt = self.calc_area()
        self.n = n
        self.n_inf = n_inf
        self.z_s = self.calc_center_of_gravity()
        self.iy, self.iy_inf = self.calc_moment_of_inertia()
        self.w = self.calc_weight()

    def calc_area(self):
        # in: width b and bw [m], height h and h_f[m]
        # out: area [m2]
        a_brutt = self.b * self.h / self.a + 1 * self.t2 + 1 * self.t3
        return a_brutt

    def calc_bef(self, sign='comp' ):
        # in: width b and bw [m], Abstand Momentennullpunkte l_0 [m]
        # out: effective width b_eff
        l_0 = self.l0
        if sign == 'comp':
            b_ef_schub = 0.1 * l_0
            b_ef_beulen = 20 * self.t3  # falls Fasern rechtwinklig zu Stegen wären, ist Faktor falsch!
            b_ef = min(b_ef_schub, b_ef_beulen, self.a - self.b)
            return b_ef
        else:
            b_ef_schub = 0.1 * l_0
            b_ef = min(b_ef_schub, self.a - self.b)
            return b_ef

    def calc_center_of_gravity(self):
        # in: Geometry effective width b, h, a, t2, b_ef_t, t3, b_ef_c
        # out: center of gravity z_s [m]
        z_s1 = self.t3 + self.h/2
        z_s2 = self.t3+self.h + self.t2/2
        z_s3 = self. t3/2
        z_s = ((self.b * self.h *z_s1 + self.bt_ef * self.t2 * z_s2 + self.bc_ef * self.t3 * z_s3) /
               (self.b * self.h + self.bt_ef * self.t2 + self.bc_ef * self.t3))
        return z_s

    def calc_moment_of_inertia(self):
        # in: Geometry b, h, t2, bt_ef, t3, bc_ef, zs
        # out: moment of inertia I_y [m^4]

        #z=0: Oberkante obere Beplankung
        z_s1 = self.t3 + self.h/2
        z_s2 = self.t3+self.h + self.t2/2
        z_s3 = self. t3/2

        i_1 = self.n[0] * self.b * self.h ** 3 / 12
        as_1 = self.n[0] * self.b * self.h * abs(self.z_s - z_s1) ** 2
        i_2 = self.n[1] * self.bt_ef * self.t2 ** 3 / 12
        as_2 = self.n[1] * self.bt_ef * self.t2 * abs(self.z_s - z_s2) ** 2
        i_3 = self.n[2] * self.bc_ef * self.t3 ** 3 / 12
        as_3 = self.n[2] * self.bc_ef * self.t3 * abs(self.z_s - z_s3) ** 2
        iy = i_1 + as_1 + i_2 + as_2 + i_3 + as_3
        i_1_inf = self.n_inf[0] * self.b * self.h ** 3 / 12
        as_1_inf = self.n_inf[0] * self.b * self.h * abs(self.z_s - z_s1) ** 2
        i_2_inf = self.n_inf[1] * self.bt_ef * self.t2 ** 3 / 12
        as_2_inf = self.n_inf[1] * self.bt_ef * self.t2 * abs(self.z_s - z_s2) ** 2
        i_3_inf = self.n_inf[2] * self.bc_ef * self.t3 ** 3 / 12
        as_3_inf = self.n_inf[2] * self.bc_ef * self.t3 * abs(self.z_s - z_s3) ** 2
        iy_inf = i_1_inf + as_1_inf + i_2_inf + as_2_inf + i_3_inf + as_3_inf
        return iy, iy_inf

    def calc_weight(self, spec_weight=5):
        #  in: specific weight [N/m^3]
        #  out: weight of cross section per m length [N/m]
        w = spec_weight * self.a_brutt
        return w

#................................................................
class RibWood(SupStrucRibWood):
    # defines properties of ribbed timber slab = "Hohlkastendecke" → box beam floor or "Ripendecke" = → joist floor
    def __init__(self, wood_type_1, wood_type_2, wood_type_3, l0, b, h, a, t2, t3, phi_1=0.6, phi_2=0.6, phi_3=0.6,
                 xi=0.01, ei_b=0.0):  # create a rectangular timber object
        section_type = "wd_rib"
        self.wood_type_1 = wood_type_1
        self.wood_type_2 = wood_type_2
        self.wood_type_3 = wood_type_3

        self.phi_1 = phi_1
        self.phi_2 = phi_2
        self.phi_3 = phi_3
        self.phi = phi_1

        self.l0 = l0

        n, n_inf = self.calc_n()
        super().__init__(section_type, b, h, a, t2, t3, n, n_inf)

        mu1_rand_u, mu1_rand_o, mu2_rand_u, mu2_rand_o, mu3_rand_u, mu3_rand_o = self.calc_mu()
        #print("mu1_rand_u, muq_rand_o, mu2_rand_u, mu2_rand_o, mu3_rand_u, mu3_rand_o =", mu1_rand_u, mu1_rand_o, mu2_rand_u, mu2_rand_o, mu3_rand_u, mu3_rand_o)
        mu_el = max(mu1_rand_u, mu1_rand_o, mu2_rand_u, mu2_rand_o, mu3_rand_u, mu3_rand_o)
        self.mu_max, self.mu_min = [mu_el, -mu_el]
        vu_el = self.calc_vu()
        self.vu_p, self.vu_n = vu_el, vu_el

        self.qs_class_n, self.qs_class_p = [3, 3]  # Required cross-section class: 1:=PP, 2:EP, 3:EE
        self.g0k = self.calc_weight(wood_type_1.weight)
        self.ei1 = self.wood_type_1.Emmean * self.iy  # elastic stiffness [Nm^2], Zeitpunkt t = 0

        self.co2 = (self.b*self.h * self.wood_type_1.GWP * self.wood_type_1.density)/self.a +self.t2 * self.wood_type_2.GWP * self.wood_type_2.density + self.t3 * self.wood_type_3.GWP * self.wood_type_3.density # [kg_CO2_eq/m]
        self.cost = self.b * self.h / self.a * self.wood_type_1.cost + (self.t2 + self.t3)  * self.wood_type_2.cost
        self.ei_b = ei_b  # stiffness perpendicular to direction of span
        self.xi = xi  # damping factor, preset value see: HBT, Page 47 (higher value for some buildups possible)

    def calc_n(self):
        ft0d = 8.5 #C24
        fc0d = 12.4 #C24
        E0mean = 11000 #C24

        factor = 2/3 #Dreischichtplatte 9/9/9 oder 10/10/10

        n1 = self.wood_type_1.Emmean / self.wood_type_1.Emmean          # Wertigkeit Rippe
        n2 = self.wood_type_2.Emmean*factor / self.wood_type_1.Emmean  # Wertigkeit Beplankung unten           #Todo: EMMEAN reduzieren! Stimmt das?
        n3 = self.wood_type_3.Emmean*factor / self.wood_type_1.Emmean  # Wertigkeit Beplankung oben            #Todo: EMMEAN reduzieren!
        n = [n1, n2, n3]
        n1_inf = (self.wood_type_1.Emmean / (1 + self.phi_1)) / (
                self.wood_type_1.Emmean / (1 + self.phi_1))  # Wertigkeit Rippe t=inf
        n2_inf = (self.wood_type_2.Emmean*factor / (1 + self.phi_2)) / (
                self.wood_type_1.Emmean / (1 + self.phi_1))  # Wertigkeit Beplankung unten t=inf    #Todo: EMMEAN reduzieren!
        n3_inf = (self.wood_type_3.Emmean*factor / (1 + self.phi_3)) / (
                self.wood_type_1.Emmean / (1 + self.phi_1))  # Wertigkeit Beplankung oben t=inf     #Todo: EMMEAN reduzieren!
        n_inf = [n1_inf, n2_inf, n3_inf]
        return n, n_inf

    def calc_mu(self):
        #Nachweise nach SIA 5.3.5 Tafelelemente (Biegeelemente)-----PRÜFEN

        fy1 = self.wood_type_1.fmd
        #print("fy1= ", fy1)
        fy2 = 8600000  #self.wood_type_2.fcd      #Festigkeiten für 3S Platten reduzieren
        fy3 = 5900000  #self.wood_type_3.ftd      #Festigkeiten für 3S Platten reduzieren

        mu1_rand_o = min(self.mu_unsigned(fy1, self.iy, (self.z_s - self.t3), self.n[0]),  # z = zs -t3
                       self.mu_unsigned(fy1, self.iy_inf, (self.z_s - self.t3), self.n_inf[0]))
        mu1_rand_u = min(self.mu_unsigned(fy1, self.iy, (self.h + self.t3 - self.z_s), self.n[0]),  # z = h + t3 -zs
                       self.mu_unsigned(fy1, self.iy_inf,(self.h + self.t3 - self.z_s), self.n_inf[0]))


        mu2_rand_o = min(self.mu_unsigned(fy2, self.iy, (self.t3 + self.h - self.z_s ), self.n[1]),  # z = t3 + h - zs
                         self.mu_unsigned(fy2, self.iy_inf, (self.t3 + self.h - self.z_s ), self.n_inf[1]))
        mu2_rand_u = min(self.mu_unsigned(fy2, self.iy, (self.t3 + self.h + self.t2- self.z_s ), self.n[1]),  # z = t3 + h + t2 - zs
                         self.mu_unsigned(fy2, self.iy_inf, (self.t3 + self.h + self.t2- self.z_s ), self.n_inf[1]))

        mu3_rand_o = min(self.mu_unsigned(fy3, self.iy, self.z_s, self.n[2]),  # z = zs
                         self.mu_unsigned(fy3, self.iy_inf, self.z_s, self.n_inf[2]))
        mu3_rand_u = min(self.mu_unsigned(fy3, self.iy, (self.z_s - self.t3), self.n[2]),  # z = zs -t3
                         self.mu_unsigned(fy3, self.iy_inf, (self.z_s - self.t3), self.n_inf[2]))
        return mu1_rand_u, mu1_rand_o, mu2_rand_u, mu2_rand_o, mu3_rand_u, mu3_rand_o

    @staticmethod
    def mu_unsigned(fy, iy, z, n):
        mu = fy * iy / z / n
        return mu

    def calc_vu(self):
        ty1 = self.wood_type_1.fvd
        vu_1 = ty1 * self.b * self.h / 1.5  #nur Rippe angesetzt
        return vu_1

    # FEHLT: Rollschubnachweis!!


#TODO: Aktueller Stand wird kein Abbrand des Hohlkastenquerschnitts berechnet, die Schichtdicken werden gem. Lignum so gewählt, dass der Abbrand nicht Bemessen werden muss.
#TODO: Folgende Zeilen müssen angepasst werden, wenn ein anderes Prinzip gewählt wird.
    @staticmethod
    def fire_resistance(section):
         #bnds = [(0, 240)]
         #t0 = 60
         #max_t = minimize(RectangularWood.fire_minimizer, t0, args=[member], bounds=bnds)
         #t_max = max_t.x[0]
        t2 = section.t2
        t3 = section.t3
        b = section.b
        h = section.h
        if t2 >= 50: #and b > ? and h > ? and t3 >= ?
            resistance = 90
        if t2 >= 26 and b > 60 and h > 180 and t3 >= 27:
            resistance = 60
        if t2 >= 10: #and b > ? and h > ? and t3 >= ?
            resistance = 30
        else:
            resistance = 0
        return resistance

    # @staticmethod
    # def fire_minimizer(t, args):
    #     member = args[0]
    #     rem_sec = RectangularWood.remaining_section(member.section, member.fire, t)
    #     mu_fire = 1.8 * rem_sec.mu_max
    #     vu_fire = 1.8 * rem_sec.vu_p  # SIA 265 (51)
    #     qd_fire = member.psi[2] * member.qk + member.gk
    #     qd_fire_zul = min(mu_fire / (max(member.system.alpha_m) * member.system.l_tot ** 2),
    #                           vu_fire / (max(member.system.alpha_v) * member.system.l_tot))
    #     to_opt = abs(qd_fire - qd_fire_zul)
    #     return to_opt
    #
    #     staticmethod
    # def remaining_section(section, fire, t=60, dred=0.007):
    #     betan = section.wood_type.burn_rate
    #     dcharn = betan * t
    #     d_ef = dcharn + dred
    #     h_fire = max(section.h - d_ef * (fire[0] + fire[2]))
    #     b_fire = max(section.b - d_ef * (fire[1] + fire[3]), 0)
    #     rem_sec = RectangularWood(section.wood_type, b_fire, h_fire)
    #     return rem_sec

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
class MatLayer:  # create a material layer
    def __init__(self, mat_name, h_input, roh_input, database):  # get initial data from database
        self.name = mat_name
        connection = sqlite3.connect(database)
        cursor = connection.cursor()

        # ── auto-detect column names (works with annotated OR plain naming) ───
        cursor.execute("PRAGMA table_info(floor_struc_prop)")
        db_cols = [row[1] for row in cursor.fetchall()]

        def _col(*substrings):
            """Return the first column whose name contains ALL substrings
            (case-insensitive). Falls back to substrings[0] so SQLite raises
            a clear OperationalError if truly absent."""
            for cn in db_cols:
                cl = cn.lower()
                if all(s.lower() in cl for s in substrings):
                    return cn
            return substrings[0]

        c_name    = _col("name")
        c_h_fix   = _col("h_fix")
        c_density = _col("density")
        c_weight  = _col("weight")
        c_gwp     = _col("gwp")
        # E-modulus: annotated "E [float, N/m^2]" or plain "E" / "e"
        c_e = next(
            (cn for cn in db_cols if cn.upper() == "E"
             or cn.lower().startswith("e [float")),
            "E [float, N/m^2]")   # fallback (will OperationalError if absent)

        inquiry = (f'SELECT "{c_h_fix}", "{c_e}", "{c_density}", "{c_weight}", "{c_gwp}" '
                   f'FROM floor_struc_prop WHERE "{c_name}"={mat_name}')
        cursor.execute(inquiry)
        result = cursor.fetchall()
        connection.close()

        if not result:
            raise ValueError(
                f"MatLayer: material {mat_name} not found in 'floor_struc_prop'.\n"
                f"  Searched column: '{c_name}'  |  Available columns: {db_cols}")

        h_fix, e, density, weight, self.GWP = result[0]
        if h_input is False:
            self.h = h_fix
        else:
            self.h = h_input
        if roh_input is False:
            self.density = density
            if weight is not None:
                self.weight = weight
            elif density is not None:
                self.weight = density * 10  # derive N/m³ from kg/m³ (g ≈ 10 m/s²)
        else:
            self.density = roh_input
            self.weight = roh_input * 10
        if e == None:
            self.ei = 0.0
        else:
            i = 1 * self.h ** 3 / 12
            self.ei = e * i
        self.gk = self.weight * self.h  # weight per area in N/m^2
        self.co2 = (self.density or 0.0) * self.h * (self.GWP or 0.0)  # CO2-eq per area in kg-CO2/m^2


class FloorStruc:  # create a floor structure
    def __init__(self, mat_layers, database_name):
        self.layers = []
        self.co2 = 0
        self.gk_area = 0
        self.h = 0
        self.ei = 0
        for mat_name, h_input, roh_input in mat_layers:
            current_layer = MatLayer(mat_name, h_input, roh_input, database_name)
            self.layers.append(current_layer)
            self.co2 += current_layer.co2
            self.gk_area += current_layer.gk
            self.h += current_layer.h
            self.ei = max(self.ei, current_layer.ei)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
class BeamSimpleSup:
    """
    Definiert die statischen Eigenschaften (Faktoren) eines Einfeldträgers
    :M = ql^2/8, 0
    :V = ql/2, -ql/2
    :w = 5/384·ql4/EI
    """
    def __init__(self, length):
        self.l_tot = length
        self.li_max = self.l_tot    # max span (used for calculation of admissible deflections)
        self.alpha_m = [0, 1 / 8]   # Faktor zur Berechung des Momentes unter verteilter Last
        self.alpha_v = [0, 1 / 2]   # Faktor zur Berechung der Querkarft unter verteilter Last
        self.qs_cl_erf = [3, 3]     # Querschnittsklasse: 1 == PP, 2 == EP, 3 == EE
        self.alpha_w = 5 / 384      # Faktor zur Berechung der Durchbiegung unter verteilter Last
        self.kf2 = 1.0              # Hilfsfaktor zur Brücksichtigung der Spannweitenverhältnisse bei Berechnung f1 gem. HBT, S. 46
        self.alpha_w_f_cd = 1 / 48  # Faktor zur Berechnung der Durchbiegung unter Einzellast

class BeamTwoSpan:
    def __init__(self, length):
        self.l_tot = length
        self.li_max = self.l_tot  # max span (used for calculation of admissible deflections)
        self.alpha_m = [-0.125, 0.0703]  # Faktor zur Berechung des Momentes unter verteilter Last
        self.alpha_v = [3/8, 5/8]  # Faktor zur Berechung der Querkarft unter verteilter Last
        self.qs_cl_erf = [3, 3]  # Querschnittsklasse: 1 == PP, 2 == EP, 3 == EE
        self.alpha_w = 2 / 369  # Faktor zur Berechung der Durchbiegung unter verteilter Last
        self.kf2 = 1.0  # Hilfsfaktor zur Brücksichtigung der Spannweitenverhältnisse bei Berechnung f1 gem. HBT, S. 46
        self.alpha_w_f_cd = 1/(48*5**0.5)  # Faktor zur Berechnung der Durchbiegung unter Einzellast

class BeamContinuousSupEl:
    def __init__(self, length):
        self.l_tot = length
        self.li_max = self.l_tot  # max span (used for calculation of admissible deflections)
        self.alpha_m = [-1/12, 1/24]  # Faktor zur Berechung des Momentes unter verteilter Last
        self.alpha_v = [0.5, 0.5]  # Faktor zur Berechung der Querkarft unter verteilter Last
        self.qs_cl_erf = [3, 3]  # Querschnittsklasse: 1 == PP, 2 == EP, 3 == EE
        self.alpha_w = 1 / 384  # Faktor zur Berechung der Durchbiegung unter verteilter Last
        self.kf2 = 1.0  # Hilfsfaktor zur Brücksichtigung der Spannweitenverhältnisse bei Berechnung f1 gem. HBT, S. 46
        self.alpha_w_f_cd = 1/192  # Faktor zur Berechung der Durchbiegung unter Einzellast

class BeamContinuousSupPl:
    def __init__(self, length):
        self.l_tot = length
        self.li_max = self.l_tot  # max span (used for calculation of admissible deflections)
        self.alpha_m = [-3/48, 3/48]  # Faktor zur Berechung des Momentes unter verteilter Last
        self.alpha_v = [0.5, 0.5]  # Faktor zur Berechung der Querkarft unter verteilter Last
        self.qs_cl_erf = [3, 3]  # Querschnittsklasse: 1 == PP, 2 == EP, 3 == EE
        self.alpha_w = 1 / 384  # Faktor zur Berechung der Durchbiegung unter verteilter Last
        self.kf2 = 1.0  # Hilfsfaktor zur Brücksichtigung der Spannweitenverhältnisse bei Berechnung f1 gem. HBT, S. 46
        self.alpha_w_f_cd = 1/192  # Faktor zur Berechung der Durchbiegung unter Einzellast

class Slab:
    """
    Nimmt die Faktoren für die Beanspruchung der Platte aus der Tabelle slab_properties.db, welche mit FE (Cedrus) ermittelt wurden
    Tabelle wird direkt im Skript "create_slab_properties.py" erstellt.
    """

    def __init__(self, length_x, length_y, support):
        self.raender = support
        self.lx = length_x
        self.ly = length_y
        self.li_max = max(length_x, length_y)
        self.l_tot = max(length_x, length_y)
        conn = sqlite3.connect("slab_properties.db")
        cursor = conn.cursor()
        # get mechanical properties from database
        result = cursor.execute(
                    """
                    SELECT NAME, RAENDER, LX, LY, MX_POS, MY_POS, MX_NEG, MY_NEG, V_POS, V_NEG, W, F 
                    FROM slab_properties
                    WHERE RAENDER = ? AND LX = ? AND LY = ? """, (self.raender, self.lx, self.ly)).fetchall()

        self.result = result[0]
        #Faktor alpha_m_x: Bewehrungfür l_max
        #x-Richtung = Richtung mit maximaler Spannweite
        self.alpha_m_x = (float(self.result[4]), float(self.result[5]))
        #Faktor alpha_m_x: Bewehrungfür l_min
        #y-Ritchtun = Richtung mit minimaler Spannweite
        self.alpha_m_y = (float(self.result[6]), float(self.result[7]))
        self.alpha_v = (float(self.result[8]), float(self.result[9]))
        self.qs_cl_erf = [2, 1]
        self.alpha_w = float(self.result[10])
        self.kf2 = 1.0
        self.alpha_w_f_cd = 10000

        #self.factors = [self.alpha_m, self.alpha_v, self.qs_cl_erf, self.alpha_w, self.kf2, self.alpha_w_f_cd]


class Member1D:
    def __init__(self, section, system, floorstruc, requirements, g2k=0.0, qk=2e3, psi0=0.7, psi1=0.5, psi2=0.3,
                 fire_b=True, fire_l=False, fire_t=False, fire_r=False):
        self.section = section
        self.system = system
        self.floorstruc = floorstruc
        self.requirements = requirements
        self.li_max = self.system.li_max
        self.g0k = self.section.g0k
        self.g1k = self.floorstruc.gk_area
        self.g2k = g2k
        self.gk = self.g0k + self.g1k + self.g2k
        self.qk = qk
        self.psi = [psi0, psi1, psi2]
        self.q_rare = self.gk + self.qk  #TODO: Wieso ist hier psi = 1.0?
        self.q_freq = self.gk + self.psi[1] * self.qk
        self.q_per = self.gk + self.psi[2] * self.qk
        self.m = self.q_per / 10
        self.w_install_adm = self.system.li_max / self.requirements.lw_install
        self.w_use_adm = self.system.li_max / self.requirements.lw_use
        self.w_app_adm = self.system.li_max / self.requirements.lw_app
        self.qu = self.calc_qu()
        self.mkd_n = self.system.alpha_m[0] * (self.gk + self.qk) * self.system.l_tot ** 2
        self.mkd_p = self.system.alpha_m[1] * (self.gk + self.qk) * self.system.l_tot ** 2
        self.qk_zul_gzt = float
        self.fire = [0, 0, 0, 0]  # fire from bottom, left, top, right (0: no fire; 1: fire)
        if fire_b is True:
            self.fire[0] = 1
        if fire_l is True:
            self.fire[1] = 1
        if fire_t is True:
            self.fire[2] = 1
        if fire_r is True:
            self.fire[3] = 1
        self.fire_resistance = []

        # calculation of deflections uncracked (plus cracked for concrete sections self.section.section_type[0:2] = rc))
        section_material = self.section.section_type[0:2]
        unit_def = self.system.alpha_w * self.system.l_tot ** 4 / self.section.ei1  # deflection for q = 1, phi = 0

        if self.requirements.install == "ductile":
            self.w_install = unit_def * (self.q_freq + self.q_per * (self.section.phi - 1))
            if section_material == "rc":  # Alternative Durchbiegungsberechnung für Betonquerschnitte gem. SIA262,(102)
                self.w_install_ger = unit_def * (
                        self.q_per * RectangularConcrete.f_w_ger(self.section.roh, self.section.rohs, self.section.phi,
                                                                 self.section.h, self.section.d)
                        + (self.q_freq - self.q_per) * RectangularConcrete.f_w_ger(self.section.roh, self.section.rohs,
                                                                                   0, self.section.h, self.section.d)
                        - self.q_per
                )

        elif self.requirements.install == "brittle":
            self.w_install = unit_def * (self.q_rare + self.q_per * (self.section.phi - 1))
            if section_material == "rc":  # Alternative Durchbiegungsberechnung für Betonquerschnitte gem. SIA262,(102)
                self.w_install_ger = unit_def * (
                        self.q_per * RectangularConcrete.f_w_ger(self.section.roh, self.section.rohs, self.section.phi,
                                                                 self.section.h, self.section.d)
                        + (self.q_rare - self.q_per) * RectangularConcrete.f_w_ger(self.section.roh, self.section.rohs,
                                                                                   0, self.section.h, self.section.d)
                        - self.q_per
                )
        self.w_use = unit_def * (self.q_freq - self.gk)
        if section_material == "rc":  # Alternative Durchbiegungsberechnung für Betonquerschnitte gem. SIA262,(102)
            self.w_use_ger = unit_def * (
                    (self.q_freq - self.q_per) * RectangularConcrete.f_w_ger(self.section.roh, self.section.rohs, 0,
                                                                             self.section.h, self.section.d)
            )
        self.w_app = unit_def * (self.q_per * (1 + self.section.phi))
        if section_material == "rc":  # Alternative Durchbiegungsberechnung für Betonquerschnitte gem. SIA262,(102)
            self.w_app_ger = unit_def * (
                    self.q_per * RectangularConcrete.f_w_ger(self.section.roh, self.section.rohs, self.section.phi,
                                                             self.section.h, self.section.d)
            )
        self.co2 = system.l_tot * (self.floorstruc.co2 + self.section.co2)

        # calculation first frequency (uncracked cross-section, method for cracked cross-section is not implemented jet)
        self.f1 = self.calc_f1()
        # calculation of further vibration criteria for wooden cross-sections
        section_material = self.section.section_type[0:2]
        if section_material == "wd" or section_material == "rc":  # check for material type
            self.ei_b = max(self.section.ei_b,
                            self.floorstruc.ei)  # Berücksichtigung n.t. Bodenaufbau gemäss Beispielsammlung HBT)
            # if floor structure has no stiffness (e.g. raised access floor, E=None),
            # fall back to section stiffness to avoid bm_rech = 0 → division by zero
            if self.ei_b == 0:
                self.ei_b = self.section.ei1
            self.bm_rech = self.system.li_max / 1.1 * (self.ei_b / self.section.ei1) ** 0.25  # HBT Seite 46
            self.a_ed = self.calc_vib1()
            self.wf_ed, self.ve_ed = self.calc_vib2()
            if self.section.xi < 0.015:
                self.r1 = 1.0  # HBT S. 48
            elif self.section.xi < 0.025:
                self.r1 = 1.15  # HBT S. 48
            else:
                self.r1 = 1.25  # HBT S. 48
            # clip exponent to prevent float overflow (100^x overflows for x > ~154)
            # a very large ve_cd simply means the vibration check passes easily
            _exp = min(self.f1 * self.section.xi - 1, 150.0)
            self.ve_cd = self.requirements.alpha_ve_cd * 100.0 ** _exp

    def calc_qu(self):
        # calculates maximal load qu in respect to bearing moment mu_max, mu_min and static system
        alpha_m = self.system.alpha_m
        alpha_v = self.system.alpha_v
        qs_class_erf = self.system.qs_cl_erf  # z.B. [0, 2]
        qs_class_vorh = [self.section.qs_class_n, self.section.qs_class_p]

        if min(alpha_m) >= 0 or abs(alpha_m[1]) >= abs(alpha_m[0]):
            if qs_class_vorh[1] <= qs_class_erf[1]:
                # if cross-section fulfills the ductility criterion (e.g. required: PP, present PP) then assign the full
                # bending strength
                qu_bend = self.section.mu_max / (max(alpha_m) * self.system.l_tot ** 2)
            else:
                # if the cross-section is not fulfilling the ductility criterion (e.g. required: EP, present PP) then
                # assign a value, which drops from the full bending strength fast towards 0 (for concrete sections)
                # or a value of 0 (for all other sections)
                if self.section.section_type[0:2] == "rc":
                    # for reinforced concrete cross-sections: smooth change to 0 load bearing capacity when roh<roh_min
                    # or roh>roh_zul (enables more efficient optimization)
                    epsilon = 1.0e-3
                    if qs_class_vorh[1] == 1:
                        shift = 0.35
                    else:
                        shift = 0.5
                    x_d = self.section.x_p / self.section.d
                    factor = min(0.5 * (1 + 2 / np.pi * np.arctan((self.section.mu_max - self.section.mr_p) / epsilon)),    #README: Wieso wird hier mit diesem factor gearbeitet? und nicht ienfahc mit qu_bend = 0 wie beim Mehrfehldträger?
                                 1 - 0.5 * (1 + 2 / np.pi * np.arctan((x_d - shift) / epsilon)))
                    qu_bend = factor * self.section.mu_max / (max(alpha_m) * self.system.l_tot ** 2)
                else:
                    # for all other cross-sections bending strength = 0
                    qu_bend = 0
            qu_shear = self.section.vu_p / (max(alpha_v) * self.system.l_tot)
        else:
            if qs_class_vorh[1] <= qs_class_erf[1]:
                qu_bend = self.section.mu_min / (min(alpha_m) * self.system.l_tot ** 2)
            else:
                # if the cross-section is not fulfilling the ductility criterion (e.g. required: EP, present PP) then
                # assign a value, which drops from the full bending strength fast towards 0 (for concrete sections)
                # or a value of 0 (for all other sections)
                if self.section.section_type[0:2] == "rc":
                    # for reinforced concrete cross-sections: smooth change to 0 load bearing capacity when roh<roh_min
                    # or roh>roh_zul (enables more efficient optimization)
                    epsilon = 1.0e-3
                    if qs_class_vorh[0] == 1:
                        shift = 0.35
                    else:
                        shift = 0.5
                    x_d = self.section.x_n / self.section.d
                    factor = min(0.5 * (1 + 2 / np.pi * np.arctan((self.section.mu_min - self.section.mr_n) / epsilon)),
                                 1 - 0.5 * (1 + 2 / np.pi * np.arctan((x_d - shift) / epsilon)))
                    qu_bend = factor * self.section.mu_min / (min(alpha_m) * self.system.l_tot ** 2)
                else:
                    # for all other cross-sections bending strength = 0
                    qu_bend = 0
            qu_shear = self.section.vu_n / (min(alpha_v) * self.system.l_tot)
        return min(qu_bend, qu_shear)

    def calc_qk_zul_gzt(self, gamma_g=1.35, gamma_q=1.5):
        self.qk_zul_gzt = (self.qu - gamma_g * self.gk) / gamma_q

    def calc_f1(self):
        # calculates first frequency of system according to HBT, Seite 46
        kf2 = self.system.kf2
        l_rech = self.system.li_max
        section_material = self.section.section_type[0:2]
        if section_material == "rc":  # take cracked stiffness for calculation of concrete sections if section is cracked
            if self.mkd_p < self.section.mr_p and self.mkd_n > self.section.mr_n:
                eil = self.section.ei1
            else:
                eil = self.section.ei2
        else:
            eil = self.section.ei1
        m = self.m
 #       print("m =", m)
        # Guard: eil (stiffness) and m (mass) must be > 0; invalid optimizer
        # candidates can produce negative/zero values — clamp with abs/fallback.
        f1 = kf2 * np.pi / (2 * l_rech ** 2) * (abs(eil / m) ** 0.5 if m else 0.0)  # HBT, Seite 46
        return f1

    def calc_vib1(self, f0=700):
        # calculates a_Ed according to HBT, Seite 47
        f1 = self.f1
        m_gen = self.m * self.system.li_max / 2 * self.bm_rech
        xi = self.section.xi
        if f1 <= 5.1:
            alpha = 0.2
            ff = f1
        elif f1 <= 6.9:
            alpha = 0.06
            ff = f1
        else:
            alpha = 0.06
            ff = 6.9
        a_ed = 0.4 * f0 * alpha / m_gen * 1 / (
                    ((f1 / ff) ** 2 - 1) ** 2 + (2 * xi * f1 / ff) ** 2) ** 0.5  # HBT, Seite 47
        return a_ed

    def calc_vib2(self, f=1000):
        # calculates W_F,ED according to to HBT, Seite 48
        wf_ed = self.system.alpha_w_f_cd * f * self.system.li_max ** 3 / (self.bm_rech * self.section.ei1)
        section_material = self.section.section_type[0:2]
        if section_material == "rc":  # take cracked stiffness for calculation of concrete sections
            eil = self.section.ei2
        else:
            eil = self.section.ei1
        _base = self.m ** 3 * eil * 1e6  # must be > 0 for valid section
        ve_ed = (364 / (self.bm_rech * _base ** 0.25)
                 if _base > 0 else float('inf'))  # HBT; guard against invalid optimizer candidates
        return wf_ed, ve_ed

    def get_fire_resistance(self):
        # evaluate fire resistance
        if self.section.section_type == "rc_rec":
            fire_resistance = RectangularConcrete.fire_resistance(self.section)
        elif self.section.section_type == "wd_rec":
            fire_resistance = RectangularWood.fire_resistance(self)
        elif self.section.section_type == "rc_rib":
            fire_resistance = RibbedConcrete.fire_resistance(self.section)
        elif self.section.section_type == "wd_rib":
            fire_resistance = RibWood.fire_resistance(self.section)
        else:
            #print("fire resistance for is not defined for that cross-section type.")
            fire_resistance = None
        self.fire_resistance = fire_resistance


class Member2D:
    def __init__(self, section, system, floorstruc, requirements, g2k=0.0, qk=2e3, psi0=0.7, psi1=0.5, psi2=0.3,
                     fire_b=True, fire_l=False, fire_t=False, fire_r=False):
        """
        Definiert ein 2-Dimensionales Bauteil (Platte) mit Eigenschaften
        :section:
        :system:
        """
        self.section = section
        self.system = system
        self.floorstruc = floorstruc
        self.requirements = requirements
        self.li_min = min(self.system.lx, self.system.ly)
        self.li_max = self.system.li_max
        self.g0k = self.section.g0k
        self.g1k = self.floorstruc.gk_area
        self.g2k = g2k
        self.gk = self.g0k + self.g1k + self.g2k
        self.qk = qk
        self.psi = [psi0, psi1, psi2]
        self.q_rare = self.gk + self.qk
        self.q_freq = self.gk + self.psi[1] * self.qk
        self.q_per = self.gk + self.psi[2] * self.qk
        self.m = self.q_per / 10
        self.w_install_adm = self.li_min / self.requirements.lw_install
        self.w_use_adm = self.li_min / self.requirements.lw_use
        self.w_app_adm = self.li_min / self.requirements.lw_app
        self.qu = self.calc_qu()
        self.mkd_n = self.system.alpha_m_x[0] * (self.gk + self.qk) * self.li_max ** 2
        self.mkd_p = self.system.alpha_m_x[1] * (self.gk + self.qk) * self.li_max ** 2
        self.mkd_n_y = self.system.alpha_m_y[0] * (self.gk + self.qk) * self.li_min ** 2
        self.mkd_p_y = self.system.alpha_m_y[1] * (self.gk + self.qk) * self.li_min ** 2
        #TODO: Everything for lx and ly!
        self.qk_zul_gzt = float
        self.fire = [0, 0, 0, 0]  # fire from bottom, left, top, right (0: no fire; 1: fire)
        if fire_b is True:
            self.fire[0] = 1
        if fire_l is True:
            self.fire[1] = 1
        if fire_t is True:
            self.fire[2] = 1
        if fire_r is True:
            self.fire[3] = 1
        self.fire_resistance = []

        # calculation of deflections uncracked (plus cracked for concrete sections self.section.section_type[0:2] = rc))
        section_material = self.section.section_type[0:2]
        unit_def = self.system.alpha_w * self.system.l_tot ** 4 / self.section.ei1  # deflection for q = 1, phi = 0


        if self.requirements.install == "ductile":
            self.w_install = unit_def * (self.q_freq + self.q_per * (self.section.phi - 1))
            if section_material == "rc":  # Alternative Durchbiegungsberechnung für Betonquerschnitte gem. SIA262,(102)
                self.w_install_ger = unit_def * (
                        self.q_per * RectangularConcrete.f_w_ger(self.section.roh, self.section.rohs, self.section.phi, self.section.h, self.section.d)
                        + (self.q_freq - self.q_per) * RectangularConcrete.f_w_ger(self.section.roh, self.section.rohs,0, self.section.h,self.section.d)
                        - self.q_per
                        )
            elif self.requirements.install == "brittle":
                self.w_install = unit_def * (self.q_rare + self.q_per * (self.section.phi - 1))
                if section_material == "rc":  # Alternative Durchbiegungsberechnung für Betonquerschnitte gem. SIA262,(102)
                        self.w_install_ger = unit_def * (
                        self.q_per * RectangularConcrete.f_w_ger(self.section.roh, self.section.rohs, self.section.phi, self.section.h, self.section.d)
                        + (self.q_rare - self.q_per) * RectangularConcrete.f_w_ger(self.section.roh, self.section.rohs,0, self.section.h, self.section.d)
                        - self.q_per
                        )
            self.w_use = unit_def * (self.q_freq - self.gk)
            if section_material == "rc":  # Alternative Durchbiegungsberechnung für Betonquerschnitte gem. SIA262,(102)
                self.w_use_ger = unit_def * (
                            (self.q_freq - self.q_per) * RectangularConcrete.f_w_ger(self.section.roh,
                                                                                     self.section.rohs, 0,
                                                                                     self.section.h, self.section.d)
                    )
                self.w_app = unit_def * (self.q_per * (1 + self.section.phi))
                if section_material == "rc":  # Alternative Durchbiegungsberechnung für Betonquerschnitte gem. SIA262,(102)
                    self.w_app_ger = unit_def * (
                            self.q_per * RectangularConcrete.f_w_ger(self.section.roh, self.section.rohs,
                                                                     self.section.phi,
                                                                     self.section.h, self.section.d)
                    )
                self.co2 = system.l_tot * (self.floorstruc.co2 + self.section.co2)

                # calculation first frequency (uncracked cross-section, method for cracked cross-section is not implemented jet)
                """
                self.f1 = self.calc_f1()
                # calculation of further vibration criteria for wooden cross-sections
                section_material = self.section.section_type[0:2]
                if section_material == "wd" or section_material == "rc":  # check for material type
                    self.ei_b = max(self.section.ei_b,
                                    self.floorstruc.ei)  # Berücksichtigung n.t. Bodenaufbau gemäss Beispielsammlung HBT)
                    # if floor structure has no stiffness (e.g. raised access floor, E=None),
                    # fall back to section stiffness to avoid bm_rech = 0 → division by zero
                    if self.ei_b == 0:
                        self.ei_b = self.section.ei1
                    self.bm_rech = self.system.li_max / 1.1 * (self.ei_b / self.section.ei1) ** 0.25  # HBT Seite 46
                    self.a_ed = self.calc_vib1()
                    self.wf_ed, self.ve_ed = self.calc_vib2()
                    if self.section.xi < 0.015:
                        self.r1 = 1.0  # HBT S. 48
                    elif self.section.xi < 0.025:
                        self.r1 = 1.15  # HBT S. 48
                    else:
                        self.r1 = 1.25  # HBT S. 48
                    self.ve_cd = self.requirements.alpha_ve_cd * 100 ** (self.f1 * self.section.xi - 1)
                """
        self.a_ed = 0
        self.wf_ed = 0
        self.ve_cd = 0
        self.ve_ed = 0
        self.r1 = 1.15
        self.f1 = 100

    def calc_qu(self):
        """
        Idea: qu von maximaler Spannweite definiert
        Schauen, welche Nachweise man alles in beide Richtungen machen muss und bei welchen einfach l_max ausreicht!
        """
        # calculates maximal load qu in respect to bearing moment mu_max, mu_min and static system
        alpha_m = self.system.alpha_m_x
        alpha_v = self.system.alpha_v
        qs_class_erf = self.system.qs_cl_erf  # z.B. [0, 2]
        qs_class_vorh = [self.section.qs_class_n, self.section.qs_class_p]

        if min(alpha_m) == 0:
            if qs_class_vorh[1] <= qs_class_erf[1]:
                # if cross-section fulfills the ductility criterion (e.g. required: PP, present PP) then assign the full
                # bending strength
                qu_bend = self.section.mu_max / (max(alpha_m) * self.system.l_tot ** 2)
            else:
                # if the cross-section is not fulfilling the ductility criterion (e.g. required: EP, present PP) then
                # assign a value, which drops from the full bending strength fast towards 0 (for concrete sections)
                # or a value of 0 (for all other sections)
                if self.section.section_type[0:2] == "rc":
                    # for reinforced concrete cross-sections: smooth change to 0 load bearing capacity when roh<roh_min
                    # or roh>roh_zul (enables more efficient optimization)
                    epsilon = 1.0e-3
                    if qs_class_vorh[1] == 1:
                        shift = 0.35
                    else:
                        shift = 0.5
                    x_d = self.section.x_p / self.section.d
                    factor = min(0.5 * (1 + 2 / np.pi * np.arctan((self.section.mu_max - self.section.mr_p) / epsilon)),    #README: Wieso wird hier mit diesem factor gearbeitet? und nicht ienfahc mit qu_bend = 0 wie beim Mehrfehldträger?
                                 1 - 0.5 * (1 + 2 / np.pi * np.arctan((x_d - shift) / epsilon)))
                    qu_bend = factor * self.section.mu_max / (max(alpha_m) * self.system.l_tot ** 2)
                else:
                    # for all other cross-sections bending strength = 0
                    qu_bend = 0
            qu_shear = self.section.vu_p / (max(alpha_v) * self.system.l_tot)
        else:
            if qs_class_vorh[0] <= qs_class_erf[0] & qs_class_vorh[1] <= qs_class_erf[1]:
                qu_bend = min(self.section.mu_max / (max(alpha_m) * self.system.l_tot ** 2), self.section.mu_min /
                              (min(alpha_m) * self.system.l_tot ** 2))
            else:
                qu_bend = 0
            qu_shear = min(self.section.vu_p / (max(alpha_v) * self.system.l_tot),
                           self.section.vu_n / (min(alpha_v) * self.system.l_tot))
        return min(qu_bend, qu_shear)

    def calc_qk_zul_gzt(self, gamma_g=1.35, gamma_q=1.5):
        self.qk_zul_gzt = (self.qu - gamma_g * self.gk) / gamma_q
    '''
    def calc_f1(self):
        # calculates first frequency of system according to HBT, Seite 46
        kf2 = self.system.kf2
        l_rech = self.system.li_max
        section_material = self.section.section_type[0:2]
        if section_material == "rc":  # take cracked stiffness for calculation of concrete sections if section is cracked
            if self.mkd_p < self.section.mr_p and self.mkd_n > self.section.mr_n:
                eil = self.section.ei1
            else:
                eil = self.section.ei2
        else:
            eil = self.section.ei1
        m = self.m
 #       print("m =", m)
        # Guard: eil (stiffness) and m (mass) must be > 0; invalid optimizer
        # candidates can produce negative/zero values — clamp with abs/fallback.
        f1 = kf2 * np.pi / (2 * l_rech ** 2) * (abs(eil / m) ** 0.5 if m else 0.0)  # HBT, Seite 46
        return f1

    def calc_vib1(self, f0=700):
        # calculates a_Ed according to HBT, Seite 47
        f1 = self.f1
        m_gen = self.m * self.system.li_max / 2 * self.bm_rech
        xi = self.section.xi
        if f1 <= 5.1:
            alpha = 0.2
            ff = f1
        elif f1 <= 6.9:
            alpha = 0.06
            ff = f1
        else:
            alpha = 0.06
            ff = 6.9
        a_ed = 0.4 * f0 * alpha / m_gen * 1 / (
                    ((f1 / ff) ** 2 - 1) ** 2 + (2 * xi * f1 / ff) ** 2) ** 0.5  # HBT, Seite 47
        return a_ed

    def calc_vib2(self, f=1000):
        # calculates W_F,ED according to to HBT, Seite 48
        wf_ed = self.system.alpha_w_f_cd * f * self.system.li_max ** 3 / (self.bm_rech * self.section.ei1)
        section_material = self.section.section_type[0:2]
        if section_material == "rc":  # take cracked stiffness for calculation of concrete sections
            eil = self.section.ei2
        else:
            eil = self.section.ei1
        _base = self.m ** 3 * eil * 1e6  # must be > 0 for valid section
        ve_ed = (364 / (self.bm_rech * _base ** 0.25)
                 if _base > 0 else float('inf'))  # HBT; guard against invalid optimizer candidates
        return wf_ed, ve_ed
    '''
    def get_fire_resistance(self):
        # evaluate fire resistance
        if self.section.section_type == "rc_rec":
            fire_resistance = RectangularConcrete.fire_resistance(self.section)
        # elif self.section.section_type == "wd_rec":
        #     fire_resistance = RectangularWood.fire_resistance(self)
        # elif self.section.section_type == "rc_rib":
        #     fire_resistance = RibbedConcrete.fire_resistance(self.section)
        # elif self.section.section_type == "wd_rib":
        #     fire_resistance = RibWood.fire_resistance(self.section)
        else:
            #print("fire resistance for is not defined for that cross-section type.")
            fire_resistance = None
        self.fire_resistance = fire_resistance


class Requirements:
    def __init__(self, install="ductile", lw_install=350, lw_use=350, lw_app=300, f1=8, a_cd=0.1, w_f_cdr1=1.0e-3,
                 alpha_ve_cd=1 / 3, fire='R60', cover_top=20.0):
        self.install = install
        self.lw_install = lw_install  # preset value: SIA 260
        self.lw_use = lw_use  # preset value: SIA 260
        self.lw_app = lw_app  # preset value: SIA 260
        self.f1 = f1  # preset value: HBT, Seite 46
        self.a_cd = a_cd  # preset value: HBT, Seite 46
        self.w_f_cdr1 = w_f_cdr1  # preset value: HBT, Seite 48
        self.alpha_ve_cd = alpha_ve_cd  # preset value: HBT, Seite 49
        self.t_fire    = int(fire[1:])  # required fire resistance period [min]
        self.cover_top = cover_top      # nominal concrete cover, top surface → mesh bar face [mm]
                                        # EN 1992-1-1: c_nom = c_min,dur + 10 mm (XC1 → 20 mm)
                                        # EN 1994-1-2 fire min axis dist: R30→10, R60→20, R90→30, R120→40 mm
                                        # axis distance = cover_top + d_s/2
