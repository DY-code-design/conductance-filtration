# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/lib/L1_materials/materials_real.py
#   sha256(src) : 2b8c97b104306566
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 45 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical. In addition, 9 top-level definitions (and 10 module-level
# statements used only by them) that the figures and verification scripts of
# the paper do not use were omitted (they are listed by name in rb_manifest);
# the remaining definitions are unchanged.
# ---------------------------------------------------------------------------
# Material property table with literature references (filler conductivity,
# interfacial resistances per unit area for each surface treatment, matrix
# conductivity). pair(filler_i, treatment_i, filler_j, treatment_j, matrix)
# returns the MaterialTable used by the engine; each value carries a ref key.
import sys as _sys, pathlib as _pathlib
_LIB = _pathlib.Path(__file__).resolve().parent.parent
for _d in ("L0_engine", "L1_materials", "L2_structures"):
    _p = str(_LIB / _d)
    if _p not in _sys.path:
        _sys.path.insert(0, _p)


from thermal_sheaf_filtration import Material, MaterialTable

REFERENCES = {
    "hBN_Yao2024":
        "Yao et al., 'Thermally Conductive Hexagonal Boron Nitride/Polymer "
        "Composites for Efficient Heat Transport', Adv. Funct. Mater. (2024), "
        "doi:10.1002/adfm.202405235 - h-BN in-plane 300-600 / through-plane 2-30 W/(m K), "
        "representative values in-plane ~400 / through-plane ~30",
    "hBN_review2024":
        "'Thermal Conductivity Enhancement of Polymeric Composites Using "
        "Hexagonal Boron Nitride: Design Strategies and Challenges', "
        "PMC10893155 (2024) - anisotropy of h-BN and effective values in composites",
    "Al2O3_epoxy_ITC":
        "W. Yang, K. Wang, Y. Fu, K. Zheng, Y. Chen, Y. Ma, "
        "'Interfacial Thermal Conductance between Alumina and Epoxy', "
        "J. Phys.: Conf. Ser. 2109, 012018 (2021), "
        "doi:10.1088/1742-6596/2109/1/012018 - TDTR measurement (plasma treatment, 3 min). "
        "interfacial thermal conductance of the alumina/epoxy interface: untreated 9.0 MW/(m^2 K), "
        "treated 26.3 MW/(m^2 K). As thermal resistance: 1.11e-7 / 3.80e-8 m^2 K/W. "
        "WARNING: the ITC is a per-interface quantity, so the implementation treating R_s_matrix as one side is correct.",
    "ZnO_bulk":
        "'Enhanced thermal conductivity of epoxy composites filled with "
        "tetrapod-shaped ZnO', PMC9079263 - ZnO bulk, 60 W/(m K) class",
    "ZnO_poly":
        "'Thermal conductivity of nanoscale polycrystalline ZnO thin films', "
        "J. Cryst. Growth (2010), doi:10.1016/j.jcrysgro.2010.10.038 - "
        "polycrystal about 1/2 of single crystal; nanoscale thin films 1.4-6.5 W/(m K)",
    "TIM_review2022":
        "'Recent Advances in Thermal Interface Materials for Thermal "
        "Management of High-Power Electronics', PMC9565324 (2022) - "
        "effective thermal conductivity of practical TIMs 0.8-4.2 W/(m K); positioning of AlN, BN, etc.",
    "BOB1977":
        "Batchelor & O'Brien, Proc. R. Soc. A 355:313 (1977), "
        "doi:10.1098/rspa.1977.0100 - near-field conduction and intra-particle spreading limit",
    "Nan1997":
        "Nan, Birringer, Clarke & Gleiter, J. Appl. Phys. 81:6692 (1997), "
        "doi:10.1063/1.365209 - Kapitza radius a_K = R_K*k_m",
    "UNVERIFIED_kappa":
        "handbook-level standard values",
    "UNVERIFIED_Rc":
        "no literature value is available for the filler-filler interfacial resistance R_c; "
        "ceramic-ceramic TBR is expected to be of order 1e-9 to 1e-8 m^2 K/W. "
        "Treated as a sensitivity parameter",
}


KAPPA_TABLE = {
    "AlN":     (170.0, "UNVERIFIED_kappa",
                "sintered grade. Around the middle between single crystal 285-320 and commercial sintered bodies 140-200"),
    "Al2O3":   (30.0, "UNVERIFIED_kappa", "99.5% alumina. Lower end of 30-40"),
    "ZnO":     (25.0, "ZnO_poly",
                "filler-grade polycrystal. At most 1/2 of bulk single crystal 60-100"),
    "ZnO_hi":  (60.0, "ZnO_bulk", "upper-bound assumption close to bulk single crystal"),
    "hBN":     (30.0, "hBN_Yao2024",
                "**through-plane** representative value. Lower-side value when platelet particles are used with random orientation"),
    "hBN_ip":  (400.0, "hBN_Yao2024",
                "**in-plane** representative value. Ideal upper bound for perfect alignment (not reached in practice)"),
    "hBN_or":  (100.0, "hBN_review2024",
                "assumed effective value for partial alignment. Between in-plane 400 and through-plane 30"),
    "diamond": (2000.0, "UNVERIFIED_kappa", "single-crystal diamond 2000-2200"),
    "SiO2":    (1.4, "UNVERIFIED_kappa", "fused silica 1.3-1.4. For low-kappa comparison"),
    "Thermalnite": (170.0, "UNVERIFIED_kappa",
                    "AlN single-crystal fibre. Bond branch within the same mol = interfacial resistance removed"),
}
KAPPA = {k: v[0] for k, v in KAPPA_TABLE.items()}

K_MATRIX_TABLE = {
    "epoxy":        (0.20, "TIM_review2022", "unfilled epoxy 0.15-0.25"),
    "silicone":     (0.25, "TIM_review2022", "unfilled silicone 0.15-0.30"),
    "high_k_resin": (1.00, "TIM_review2022",
                     "a hypothetical high-conductivity matrix, used for sensitivity analysis; not a real single resin"),
}
K_MATRIX = {k: v[0] for k, v in K_MATRIX_TABLE.items()}

TREATMENT_TABLE = {
    "untreated": dict(
        R_s_contact=(1.0e-8, "UNVERIFIED_Rc", "assumed TBR between ceramics"),
        R_s_matrix=(1.11e-7, "Al2O3_epoxy_ITC",
                    "reciprocal of alumina/epoxy untreated 9.0 MW/(m^2 K)"),
        t_ip_um=(0.0, "-", ""), k_ip=(None, "-", "")),
    "treated": dict(
        R_s_contact=(2.0e-9, "UNVERIFIED_Rc", "assumed reduction by treatment"),
        R_s_matrix=(3.80e-8, "Al2O3_epoxy_ITC",
                    "reciprocal of alumina/epoxy plasma-treated 26.3 MW/(m^2 K)"),
        t_ip_um=(0.0, "-", ""), k_ip=(None, "-", "")),
    "lowk_coat": dict(
        R_s_contact=(1.0e-6, "UNVERIFIED_Rc", "hypothetical; an extreme case on the unfavourable side"),
        R_s_matrix=(1.0e-6, "UNVERIFIED_Rc", "hypothetical; about ten times the untreated value"),
        t_ip_um=(0.05, "-", "coating layer thickness"),
        k_ip=(0.05, "-", "thermal conductivity of the coating layer")),
}
TREATMENT = {name: {k: v[0] for k, v in d.items()}
             for name, d in TREATMENT_TABLE.items()}
TREATMENT["silane"] = TREATMENT["treated"]
TREATMENT_TABLE["silane"] = TREATMENT_TABLE["treated"]


def material(filler, treatment="untreated", name=None):
    t = TREATMENT[treatment]
    return Material(name or f"{filler}/{treatment}", kappa=KAPPA[filler], **t)


def table(spec, matrix="epoxy", overrides=None):
    mats = {int(t): material(f, tr) for t, (f, tr) in spec.items()}
    return MaterialTable(materials=mats, k_matrix=K_MATRIX[matrix],
                         contact_Rc_override=overrides or {})


def single(filler, treatment="untreated", matrix="epoxy"):
    return table({1: (filler, treatment)}, matrix)


def pair(fi, ti, fj, tj, matrix="epoxy", overrides=None):
    return MaterialTable(materials={1: material(fi, ti), 2: material(fj, tj)},
                         k_matrix=K_MATRIX[matrix],
                         contact_Rc_override=overrides or {})


FILLERS_GRID = ["AlN", "Al2O3", "ZnO", "hBN", "diamond", "SiO2"]
TREATMENTS_GRID = ["untreated", "treated", "lowk_coat"]
MATRICES_GRID = ["epoxy", "silicone", "high_k_resin"]

R_TILDE_GRID = [0.25, 0.5, 1.25, 2.5, 5.0, 12.5, 25.0]


def radii_from_rt(r_tilde):
    return 2.0 * r_tilde, 2.0 * r_tilde


CASES = [
    ("AlN / epoxy (untreated)",     "AlN",     "untreated", "epoxy"),
    ("AlN / epoxy (treated)",       "AlN",     "treated",   "epoxy"),
    ("AlN / epoxy (low-k coat)",    "AlN",     "lowk_coat", "epoxy"),
    ("Al2O3 / epoxy (untreated)",   "Al2O3",   "untreated", "epoxy"),
    ("Al2O3 / epoxy (treated)",     "Al2O3",   "treated",   "epoxy"),
    ("ZnO / epoxy (untreated)",     "ZnO",     "untreated", "epoxy"),
    ("ZnO / epoxy (treated)",       "ZnO",     "treated",   "epoxy"),
    ("ZnO (bulk-like) / epoxy",     "ZnO_hi",  "untreated", "epoxy"),
    ("h-BN through-plane / epoxy",  "hBN",     "untreated", "epoxy"),
    ("h-BN oriented / epoxy",       "hBN_or",  "untreated", "epoxy"),
    ("h-BN in-plane / epoxy",       "hBN_ip",  "untreated", "epoxy"),
    ("diamond / epoxy (untreated)", "diamond", "untreated", "epoxy"),
    ("SiO2 / epoxy (untreated)",    "SiO2",    "untreated", "epoxy"),
    ("AlN / silicone (treated)",    "AlN",     "treated",   "silicone"),
    ("AlN / high-k resin (treated)", "AlN",    "treated",   "high_k_resin"),
]

