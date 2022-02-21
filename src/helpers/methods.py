from PyQt5 import QtWidgets
import numpy

from scipy.optimize import fsolve, newton
from src import constants


def handle_file(mode, parent, suffix, name_filter):
    dialog = QtWidgets.QFileDialog(parent)
    dialog.setDefaultSuffix(suffix)
    dialog.setNameFilter(name_filter)
    dialog.setAcceptMode(mode)
    return_code = dialog.exec_()
    if return_code != 0:
        return str((dialog.selectedFiles())[0])
    else:
        return None


def save_file(parent, suffix="", name_filter=""):
    return handle_file(1, parent, suffix, name_filter)


def load_file(parent, suffix="", name_filter=""):
    return handle_file(0, parent, suffix, name_filter)


def create_bandpass_menu(parent):
    rootmenu = QtWidgets.QMenu(parent)

    # johnson
    johnson = rootmenu.addMenu("Johnson")
    ju = johnson.addAction("U")
    ju.setObjectName("Johnson U")
    jb = johnson.addAction("B")
    jb.setObjectName("Johnson B")
    jv = johnson.addAction("V")
    jv.setObjectName("Johnson V")
    jr = johnson.addAction("R")
    jr.setObjectName("Johnson R")
    ji = johnson.addAction("I")
    ji.setObjectName("Johnson I")

    # stromgren
    stromgren = rootmenu.addMenu("Stromgren")
    su = stromgren.addAction("u")
    su.setObjectName("Stromgren u")
    sv = stromgren.addAction("v")
    sv.setObjectName("Stromgren v")
    sb = stromgren.addAction("b")
    sb.setObjectName("Stromgren b")
    sy = stromgren.addAction("y")
    sy.setObjectName("Stromgren y")

    # cousins
    cousins = rootmenu.addMenu("Cousins")
    crc = cousins.addAction("Rc")
    crc.setObjectName("Cousins Rc")
    cic = cousins.addAction("Ic")
    cic.setObjectName("Cousins Ic")

    # bessel
    bessel = rootmenu.addMenu("Bessel")
    bux = bessel.addAction("UX")
    bux.setObjectName("Bessel UX")
    bbx = bessel.addAction("BX")
    bbx.setObjectName("Bessel BX")
    bb = bessel.addAction("B")
    bb.setObjectName("Bessel B")
    bv = bessel.addAction("B")
    bv.setObjectName("Bessel V")
    br = bessel.addAction("V")
    br.setObjectName("Bessel R")
    bi = bessel.addAction("I")
    bi.setObjectName("Bessel I")

    # tycho
    tycho = rootmenu.addMenu("Tycho")
    tbt = tycho.addAction("Bt")
    tbt.setObjectName("Tycho Bt")
    tvt = tycho.addAction("Vt")
    tvt.setObjectName("Tycho Vt")

    # corot
    corot = rootmenu.addMenu("Corot")
    csis = corot.addAction("SIS")
    csis.setObjectName("COROT SIS")
    cexo = corot.addAction("EXO")
    cexo.setObjectName("COROT EXO")

    # geneva
    geneva = rootmenu.addMenu("Geneva")
    gu = geneva.addAction("U")
    gu.setObjectName("Geneva U")
    gb = geneva.addAction("B")
    gb.setObjectName("Geneva B")
    gb1 = geneva.addAction("B1")
    gb1.setObjectName("Geneva B1")
    gb2 = geneva.addAction("B2")
    gb2.setObjectName("Geneva B2")
    gv = geneva.addAction("V")
    gv.setObjectName("Geneva V")
    gv1 = geneva.addAction("V1")
    gv1.setObjectName("Geneva V1")
    gg = geneva.addAction("G")
    gg.setObjectName("Geneva G")

    # vilnius
    vilnius = rootmenu.addMenu("Vilnius")
    vu = vilnius.addAction("U")
    vu.setObjectName("Vilnius U")
    vp = vilnius.addAction("P")
    vp.setObjectName("Vilnius P")
    vx = vilnius.addAction("X")
    vx.setObjectName("Vilnius X")
    vy = vilnius.addAction("Y")
    vy.setObjectName("Vilnius Y")
    vz = vilnius.addAction("Z")
    vz.setObjectName("Vilnius Z")
    vv = vilnius.addAction("V")
    vv.setObjectName("Vilnius V")
    vs = vilnius.addAction("S")
    vs.setObjectName("Vilnius S")

    # milone
    milone = rootmenu.addMenu("Milone")
    miz = milone.addAction("iz")
    miz.setObjectName("Milone iz")
    mij = milone.addAction("iJ")
    mij.setObjectName("Milone iJ")
    mih = milone.addAction("iH")
    mih.setObjectName("Milone iH")
    mik = milone.addAction("iK")
    mik.setObjectName("Milone iK")
    m230 = milone.addAction("230")
    m230.setObjectName("Milone 230")
    m250 = milone.addAction("250")
    m250.setObjectName("Milone 250")
    m270 = milone.addAction("270")
    m270.setObjectName("Milone 270")
    m290 = milone.addAction("290")
    m290.setObjectName("Milone 290")
    m310 = milone.addAction("310")
    m310.setObjectName("Milone 310")
    m330 = milone.addAction("330")
    m330.setObjectName("Milone 330")

    # yms94
    yms94 = rootmenu.addMenu("YMS94")
    yiz = yms94.addAction("iz")
    yiz.setObjectName("YMS94 iz")
    yij = yms94.addAction("iJ")
    yij.setObjectName("YMS94 iJ")
    yih = yms94.addAction("iH")
    yih.setObjectName("YMS94 iH")
    yik = yms94.addAction("iK")
    yik.setObjectName("YMS94 iK")
    yil = yms94.addAction("iL")
    yil.setObjectName("YMS94 iL")
    yill = yms94.addAction("iL'")
    yill.setObjectName("YMS94 iL'")
    yim = yms94.addAction("iM")
    yim.setObjectName("YMS94 iM")
    yin = yms94.addAction("in")
    yin.setObjectName("YMS94 in")
    yinn = yms94.addAction("iN")
    yinn.setObjectName("YMS94 iN")

    # sloandds
    sloandds = rootmenu.addMenu("Sloan DSS")
    sdu = sloandds.addAction("u'")
    sdu.setObjectName("Sloan DSS u'")
    sdg = sloandds.addAction("g'")
    sdg.setObjectName("Sloan DSS g'")
    sdr = sloandds.addAction("r'")
    sdr.setObjectName("Sloan DSS r'")
    sdi = sloandds.addAction("i'")
    sdi.setObjectName("Sloan DSS i'")
    sdz = sloandds.addAction("z'")
    sdz.setObjectName("Sloan DSS z'")

    # hststis
    hststis = rootmenu.addMenu("HST STIS")
    hly = hststis.addAction("Ly alpha")
    hly.setObjectName("HST STIS Ly alpha")
    hlf = hststis.addAction("Fclear")
    hlf.setObjectName("HST STIS Fclear")
    hlfc = hststis.addAction("Fsrf2")
    hlfc.setObjectName("HST STIS Fsrf2")
    hlfq = hststis.addAction("Fqtz")
    hlfq.setObjectName("HST STIS Fqtz")
    hlc3 = hststis.addAction("C III")
    hlc3.setObjectName("HST STIS C III")
    hlm2 = hststis.addAction("Mg II")
    hlm2.setObjectName("HST STIS Mg II")
    hlnc = hststis.addAction("Nclear")
    hlnc.setObjectName("HST STIS Nclear")
    hlns = hststis.addAction("Nsfr2")
    hlns.setObjectName("HST STIS Nsrf2")
    hlnq = hststis.addAction("Nqtz")
    hlnq.setObjectName("HST STIS Nqtz")
    hlcn = hststis.addAction("cn182")
    hlcn.setObjectName("HST STIS cn182")
    hlcn2 = hststis.addAction("cn270")
    hlcn2.setObjectName("HST STIS cn270")
    hlo = hststis.addAction("Oclear")
    hlo.setObjectName("HST STIS Oclear")
    hloc = hststis.addAction("Oclear-lp")
    hloc.setObjectName("HST STIS Oclear-lp")
    hlo2 = hststis.addAction("[O II]")
    hlo2.setObjectName("HST STIS [O II]")
    hlo3 = hststis.addAction("[O III]")
    hlo3.setObjectName("HST STIS [O III]")

    # 2mass
    twomass = rootmenu.addMenu("2MASS")
    tmj = twomass.addAction("J")
    tmj.setObjectName("2MASS J")
    tmh = twomass.addAction("H")
    tmh.setObjectName("2MASS H")
    tmjks = twomass.addAction("Ks")
    tmjks.setObjectName("2MASS Ks")

    # gaia
    gaia = rootmenu.addMenu("Gaia")
    gg2006 = gaia.addAction("G (2006)")
    gg2006.setObjectName("Gaia G (2006)")
    gg2010 = gaia.addAction("G (2010)")
    gg2010.setObjectName("Gaia G (2010)")
    ggbp = gaia.addAction("Gbp")
    ggbp.setObjectName("Gaia Gbp")
    ggrp = gaia.addAction("Grp")
    ggrp.setObjectName("Gaia Grp")
    ggrvs = gaia.addAction("Grvs")
    ggrvs.setObjectName("Gaia Grvs")

    rootmenu.addSeparator()

    # single menu bands
    hipparchos = rootmenu.addAction("Hipparchos Hp")
    hipparchos.setObjectName("Hipparchos Hp")
    kepler = rootmenu.addAction("KEPLER")
    kepler.setObjectName("KEPLER")
    swasp = rootmenu.addAction("SWASP")
    swasp.setObjectName("SWASP")
    most = rootmenu.addAction("MOST")
    most.setObjectName("MOST")
    triplet = rootmenu.addAction("Ca II triplet")
    triplet.setObjectName("Ca II triplet")
    wire = rootmenu.addAction("WIRE V+R")
    wire.setObjectName("WIRE V+R")
    lut = rootmenu.addAction("Lunar Ultraviolet Telescope")
    lut.setObjectName("Lunar Ultraviolet Telescope")
    tess = rootmenu.addAction("TESS")
    tess.setObjectName("TESS")

    return rootmenu


def compute_temp_from_color(reference, color, error):
    def _compute(cfs, clr):
        return 10 ** (cfs[0] + (cfs[1] * clr) + (cfs[2] * clr ** 2) + (cfs[3] * clr ** 3) +
                      (cfs[4] * clr ** 4) + (cfs[5] * clr ** 5) + (cfs[6] * clr ** 6) + (cfs[7] * clr ** 7))

    param_dict = {
        "gray": [0, 0, 0, 0, 0, 0, 0, 0],  # not actually used
        "gray_cool": [3.981, -0.4728, 0.2434, -0.0620, 0, 0, 0, 0],  # do not reference this
        "gray_hot": [3.981, 0.0142, 16.3618, 81.8910, 161.5075, 0, 0, 0],  # do not reference this
        "flower": [3.979145, -0.654499, 1.740690, -4.608815, 6.792600, -5.396910, 2.192970, -0.359496],
        "drilling_landolt": [3.97758849, -0.63784759, 1.7966422, -4.37949535, 5.05075322, -2.60480351, 0.44186047,
                             0.02579832],
        "popper": [3.97981681, -0.82365352, 2.10679135, -3.47184197, 2.48191975, -0.34341818, -0.4029675,
                   0.14262473],
        "tokunaga_vk": [3.99391909, -3.15370242e-01, 2.33380328e-01, -1.27336537e-01, 3.76914126e-02,
                        -6.01480225e-03, 4.91975597e-04, -1.62340638e-05],
        "tokunaga_jh": [3.99250922, -2.44541306, 16.60593436, -26.76755927, -235.35579444, 1132.65031894,
                        -1810.38727869, 991.27146106],
        "tokunaga_hk": [4.04566110, -8.28316187, 7.97352126e+01, -4.30370358e+02, 7.87866519e+02, 2.05805024e+03,
                        -1.00848335e+04, 1.07400393e+04]
    }

    coeffs = param_dict[reference]

    if reference == "gray" and color < 0.0:
        coeffs = param_dict["gray_hot"]

    elif reference == "gray" and color >= 0.0:
        coeffs = param_dict["gray_cool"]

    temp = _compute(coeffs, color)
    temp_upper = _compute(coeffs, color - error)
    temp_lower = _compute(coeffs, color + error)

    temp_err = ((temp - temp_lower) + (temp_upper - temp)) / 2.0

    return int(numpy.round(temp)), int(numpy.round(temp_err))


def convert_jd_to_ut(jd, add_24=False):
    jd = float(jd) + 0.5

    if add_24 and jd < 2400000.0:
        jd = jd + 2400000.0

    fractional, integral = numpy.modf(jd)
    integral = int(integral)

    A = numpy.trunc((integral - 1867216.25) / 36524.25)

    B = None

    if integral > 2299160:
        B = integral + 1 + A - numpy.trunc(A / 4.)
    else:
        B = integral

    C = B + 1524
    D = numpy.trunc((C - 122.1) / 365.25)
    E = numpy.trunc(365.25 * D)
    G = numpy.trunc((C - E) / 30.6001)

    day = C - E + fractional - numpy.trunc(30.6001 * G)
    day_fract, day = numpy.modf(day)

    month = None
    if G < 13.5:
        month = G - 1
    else:
        month = G - 13

    year = None
    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715

    hour = day_fract * 24.0
    fract_hour, hour = numpy.modf(hour)

    minute = fract_hour * 60.0
    fract_minute, minute = numpy.modf(minute)

    second = fract_minute * 60.0

    return int(year), int(month), int(day), int(hour), int(minute), second


def convert_ut_to_jd(year, month, day, hour, minute, second):

    jd = 367.0 * year - int(7.0 * (year + int((month + 9.0) / 12.0)) / 4.0) + int(275.0 * month / 9.0) - \
         int(3.0 * ((year + (month - 9.0) / 7.0) / 100.0 + 1.0) / 4.0) + 1721028.5 + day + hour / 24.0 + \
         minute / 1440.0 + second / 86400.0

    return jd


def convert_hjd_to_jd(hjd, ra_h, ra_m, ra_s, dec_d, dec_m, dec_s):
    # This calculation is adopted from IRAF setjd package.

    if dec_d < 0:
        dec_m = dec_m * -1.0
        dec_s = dec_s * -1.0

    r = ((ra_h + ra_m / 60.0 + ra_s / 3600.0) * 15.0) * numpy.pi / 180.0
    d = (dec_d + dec_m / 60.0 + dec_s / 3600.0) * numpy.pi / 180.0

    t = (hjd - 2415020.0) / 36525.0

    manom = 358.47583 + t * (35999.04975 - t * (0.000150 + t * 0.000003))
    lperi = 101.22083 + t * (1.7191733 + t * (0.000453 + t * 0.000003))
    oblq = 23.452294 - t * (0.0130125 + t * (0.00000164 - t * 0.000000503))
    eccen = 0.01675104 - t * (0.00004180 + t * 0.000000126)

    manom = numpy.mod(manom, 360.0)
    lperi = numpy.mod(lperi, 360.0)

    manom = numpy.pi * manom / 180.0
    lperi = numpy.pi * lperi / 180.0
    oblq = numpy.pi * oblq / 180.0

    tanom = manom + (2.0 * eccen - 0.25 * eccen ** 3) * numpy.sin(manom) + 1.25 * eccen ** 2 * numpy.sin(
        2 * manom) + 13.0 / 12.0 * eccen ** 3 * numpy.sin(3.0 * manom)

    slong = lperi + tanom + numpy.pi

    ll = (numpy.arctan2((numpy.sin(r) * numpy.cos(oblq) + numpy.tan(d) * numpy.sin(oblq)), numpy.cos(r)))
    b = (numpy.arcsin(numpy.sin(d) * numpy.cos(oblq) - numpy.cos(d) * numpy.sin(oblq) * numpy.sin(r)))

    rsun = (1. - eccen ** 2) / (1. + eccen * numpy.cos(tanom))
    dt = -0.005770 * rsun * numpy.cos(b) * numpy.cos(ll - slong)

    jd = hjd - dt
    return jd


def compute_omega_potential(q, fractradius, f, d):
    omega = (1.0 / fractradius) + (q * ((1.0 / numpy.absolute(d - fractradius)) -
                                        (fractradius / (d ** 2)))) + (((q + 1.0) / 2.0) * (f ** 2) * (fractradius ** 2))

    return omega


def compute_residuals(obs_x, obs_y, model_x, model_y):
    interp_model = numpy.interp(obs_x, model_x, model_y)
    return obs_y - interp_model


def from_jd_to_phase(x, t0, p):
    """
    :param x: jd
    :param t0: ephemeris
    :param p: period
    :return: x in phases
    """

    _x = [((jd - t0) / p) - int((jd - t0) / p) for jd in x]

    return _x


def from_phase_to_jd(x, t0, p):
    """
    :param x: phases
    :param t0: ephemeris
    :param p: period
    :return: x in jd
    """
    _x = [((phase * p) + t0) for phase in x]

    return _x


def alias_phased_obs_with_phase(x, y, start, end):
    """
    :param x: a list containing phases
    :param y: a list containing observations
    :param start: start phase
    :param end: end phase
    :return: aliased phases and observations
    """

    x = [float(n) for n in x]
    y = [float(n) for n in y]

    if start > end:
        raise ValueError("Start phase can't be larger than stop phase.")
    if len(x) != len(y):
        raise ValueError("x and y must be the same size.")

    distance = int(start - min(x))
    if (distance == 0 and min(x) > start) or (distance < 0 < min(x)):
        distance = distance - 1

    x = [phase + distance for phase in x]

    new_x = x[:]
    new_y = y[:]

    i = 1
    while max(new_x) < end:
        x_temp = [phase + i for phase in x]
        new_x = new_x + x_temp
        new_y = new_y + y[:]
        i = i + 1

    _x = []
    _y = []

    for phase, value in zip(new_x, new_y):
        if start <= phase <= end:
            _x.append(phase)
            _y.append(value)

    return _x, _y


def alias_phased_obs_with_jd(x, y, start, end, t0, p):
    """
    :param x: a list containing phases
    :param y: a list containing observations
    :param start: start jd
    :param end: end jd
    :param t0: ephemeris
    :param p: period
    :return: aliased jd and observations
    """

    if start > end:
        raise ValueError("Start phase can't be larger than stop phase.")
    if len(x) != len(y):
        raise ValueError("x and y must be the same size.")

    new_x = x[:]
    new_y = y[:]

    new_x = from_phase_to_jd(new_x, t0, p)

    new_x = numpy.array(new_x)
    new_y = numpy.array(new_y)

    distance = start - min(new_x)
    distance_in_phase = numpy.modf(distance / p)[1] + (1 * numpy.sign(distance))
    new_x = new_x + (distance_in_phase * p)

    xx = new_x.copy()

    while max(new_x) < end:
        xx = xx + p
        new_x = numpy.insert(new_x, new_x.shape[0], xx.copy())
        new_y = numpy.insert(new_y, new_y.shape[0], y.copy())

    _x = []
    _y = []

    for jd, value in zip(new_x, new_y):
        if start <= jd <= end:
            _x.append(jd)
            _y.append(value)

    return _x, _y

    pass


def alias_jd_obs_with_phase(x, y, start, end, t0, p):
    """
    :param x: a list containing jd
    :param y: a list containing observations
    :param start: start phase
    :param end: end phase
    :param t0: ephemeris
    :param p: period
    :return: aliased phases and observations
    """

    if start > end:
        raise ValueError("Start phase can't be larger than stop phase.")
    if len(x) != len(y):
        raise ValueError("x and y must be the same size.")

    new_x = x[:]
    new_y = y[:]

    new_x = from_jd_to_phase(new_x, t0, p)

    return alias_phased_obs_with_phase(new_x, new_y, start, end)


def compute_roche_potentials(w, e, q, phase, phase_shift, plot_elements=None):
    # This snippet is only for f1, f2 = 1 for now
    # For in-depth discussion about calculating Roche potentials, refer to:
    #  - Eclipsing Binary Stars: Modeling and Analysis (Kallrath & Milone, 2009, Springer)

    true_anomaly = (numpy.pi / 2.0) - w
    eccentric_anomaly = 2.0 * numpy.arctan(numpy.sqrt((1.0 - e) / (1.0 + e)) * numpy.tan(true_anomaly / 2.0))
    mean_anomaly = eccentric_anomaly - e * numpy.sin(eccentric_anomaly)
    conjunction = ((mean_anomaly + w) / (2.0 * numpy.pi)) - 0.25 + phase_shift  # superior conjunction phase
    periastron_passage = 1.0 - mean_anomaly / (2.0 * numpy.pi)
    periastron_phase = conjunction + periastron_passage  # phase of periastron

    while periastron_phase > 1.0:
        periastron_phase = periastron_phase - int(periastron_phase)

    M = 2.0 * numpy.pi * (phase - periastron_phase)

    while M < 0.0:
        M = M + 2.0 * numpy.pi

    f_E = lambda E: E - e * numpy.sin(E) - M
    E = fsolve(f_E, M)

    separation_at_phase = 1.0 - e * numpy.cos(E)

    q_inverse = q
    if q >= 1.0:
        q_inverse = 1.0 / q

    f = 1.0

    f_critical_inner = lambda x: (-1.0 / x ** 2) - \
                           (q * ((x - separation_at_phase) / pow(numpy.absolute(separation_at_phase - x), 3))) + \
                           (x * f ** 2 * (q + 1)) - (q / separation_at_phase ** 2)  # Appendix E.12.4

    inner_critical_x = fsolve(f_critical_inner, separation_at_phase / 2.0)

    inner_potential = (1 / inner_critical_x) + (q * ((1 / numpy.absolute(separation_at_phase - inner_critical_x))
                                                     - (inner_critical_x / (separation_at_phase ** 2)))) + \
                      (((q + 1) / 2) * (f ** 2) * (inner_critical_x ** 2))  # Appendix E.12.8

    mu = None
    if q >= 1.0:
        mu = (1.0 / 3.0) * (1.0 / q) / (1.0 + (1.0 / q))
    else:
        mu = (1.0 / 3.0) * q / (1.0 + q)

    outer_critical_estimation = 1.0 + mu ** (1.0 / 3.0) + (1.0 / 3.0) * mu ** \
                                (2.0 / 3.0) + (1.0 / 9.0) * mu ** (3.0 / 3.0)

    f_critical_outer = lambda x: (-1.0 / x ** 2) - \
                                 (q_inverse * ((x - separation_at_phase) /
                                               pow(numpy.absolute(separation_at_phase - x), 3))) \
                                            + (x * f ** 2 * (q_inverse + 1)) - (q_inverse / separation_at_phase ** 2)
    # Appendix E.12.4 with inverse q check  ^^^^^

    f_critical_outer_deriv = lambda x: (2.0 / x ** 3) + \
                                       ((2.0 * q_inverse) / pow(numpy.abs(separation_at_phase - x), 3)) + \
                                       f ** 2.0 * (q_inverse + 1.0)  # Appendix E.12.5

    outer_critical_x = newton(f_critical_outer, outer_critical_estimation, tol=1e-6, fprime=f_critical_outer_deriv)

    if q >= 1.0:
        outer_critical_x = 1.0 - outer_critical_x

    outer_potential = (1.0 / numpy.abs(outer_critical_x)) + (q * (
        (1.0 / numpy.absolute(separation_at_phase - outer_critical_x)) - (
            outer_critical_x / (separation_at_phase ** 2)))) + (
                            ((q + 1.0) / 2.0) * (f ** 2) * (outer_critical_x ** 2))  # Appendix E.12.8

    if plot_elements is not None:
        # only  used for calculating left / right limits
        f_outer_critical = lambda x: ((1.0 / numpy.abs(x)) + (q * (
            (1.0 / numpy.absolute(separation_at_phase - x)) - (x / (separation_at_phase ** 2)))) + (
                                ((q + 1.0) / 2.0) * (f ** 2) * (x ** 2))) - outer_potential
        f_inner_critical = lambda x: ((1.0 / numpy.abs(x)) + (q * (
            (1.0 / numpy.absolute(separation_at_phase - x)) - (x / (separation_at_phase ** 2)))) + (
                                ((q + 1.0) / 2.0) * (f ** 2) * (x ** 2))) - inner_potential

        x_axis_temporary = 0.0
        if q >= 1.0:
            x_axis_temporary = separation_at_phase
        z_axis_temporary = numpy.linspace(0, 5, 5000)
        (X_temp, Z_temp) = numpy.meshgrid(x_axis_temporary, z_axis_temporary)

        all_pots_temp = ((1 / numpy.sqrt(X_temp ** 2 + Z_temp ** 2)) + (q * (
                (1 / numpy.sqrt(
                    (separation_at_phase ** 2) - (2 * X_temp * separation_at_phase) + (
                                numpy.sqrt(X_temp ** 2 + Z_temp ** 2) ** 2))) - (
                        X_temp / (separation_at_phase ** 2)))) + (
                                     0.5 * (f ** 2) * (q + 1) * ((X_temp ** 2) + (Z_temp ** 2))))

        critical_limit = None
        if e == 0.0:
            critical_limit = outer_potential
        else:
            critical_limit = inner_potential

        upper_limit = None
        lower_limit = None

        for y_val in all_pots_temp:
            if y_val <= critical_limit:
                idx = numpy.where(all_pots_temp == y_val)[0]
                upper_limit = (float(idx) * (5.0 / 5000.0)) + 0.01
                lower_limit = upper_limit * -1.0
                break

        left_limit = None
        right_limit = None

        if e == 0.0:
            left_limit = fsolve(f_outer_critical, -0.5)
            right_limit = fsolve(f_outer_critical, 1.5)
        else:
            left_limit = fsolve(f_inner_critical, -0.25)
            right_limit = fsolve(f_inner_critical, 1.5)

        x_axis = numpy.linspace(left_limit, right_limit, 2000)
        z_axis = numpy.linspace(lower_limit, upper_limit, 2000)
        (X, Z) = numpy.meshgrid(x_axis, z_axis)

        all_pots = ((1 / numpy.sqrt(X ** 2 + Z ** 2)) + (q * (
            (1 / numpy.sqrt(
                (separation_at_phase ** 2) - (2 * X * separation_at_phase) + (numpy.sqrt(X ** 2 + Z ** 2) ** 2))) - (
                X / (separation_at_phase ** 2)))) + (0.5 * (f ** 2) * (q + 1) * ((X ** 2) + (Z ** 2))))

        center_of_mass = (separation_at_phase / (1 + (1 / q)))

        x_axis_primary = numpy.linspace(left_limit, inner_critical_x, 2000)
        z_axis_primary = numpy.linspace(-0.6, 0.6, 2000)
        (X_primary, Z_primary) = numpy.meshgrid(x_axis_primary, z_axis_primary)

        x_axis_secondary = numpy.linspace(inner_critical_x, right_limit, 2000)
        z_axis_secondary = numpy.linspace(-0.6, 0.6, 2000)
        (X_secondary, Z_secondary) = numpy.meshgrid(x_axis_secondary, z_axis_secondary)

        all_pots_primary = ((1 / numpy.sqrt(X_primary ** 2 + Z_primary ** 2)) + (q * (
            (1 / numpy.sqrt(
                (separation_at_phase ** 2) - (2 * X_primary * separation_at_phase) +
                (numpy.sqrt(X_primary ** 2 + Z_primary ** 2) ** 2))) - (
                X_primary / (separation_at_phase ** 2)))) +
                            (0.5 * (f ** 2) * (q + 1) * ((X_primary ** 2) + (Z_primary ** 2))))

        all_pots_secondary = ((1 / numpy.sqrt(X_secondary ** 2 + Z_secondary ** 2)) + (q * (
            (1 / numpy.sqrt(
                (separation_at_phase ** 2) - (2 * X_secondary * separation_at_phase) +
                (numpy.sqrt(X_secondary ** 2 + Z_secondary ** 2) ** 2))) - (
                X_secondary / (separation_at_phase ** 2)))) +
                              (0.5 * (f ** 2) * (q + 1) * ((X_secondary ** 2) + (Z_secondary ** 2))))

        axis = plot_elements[0]
        pot1 = plot_elements[1]
        pot2 = plot_elements[2]

        axis.contour(X, Z, all_pots, inner_potential, colors=constants.COLOR_RED)
        if e == 0.0:
            axis.contour(X, Z, all_pots, outer_potential, colors=constants.COLOR_BLUE)
        axis.contourf(X_primary, Z_primary, all_pots_primary, [pot1, pot1 * 1000], cmap='Dark2', vmin=0, vmax=1)
        axis.contourf(X_secondary, Z_secondary, all_pots_secondary, [pot2, pot2 * 1000], cmap='Dark2', vmin=0, vmax=1)
        axis.plot([0, separation_at_phase, center_of_mass], [0, 0, 0], linestyle="",
                  marker="+", markersize=10, color=constants.COLOR_RED)

    return inner_potential, outer_potential


def compute_conjunction_phases(w, e, phase_shift):
    # Calculations below are adopted from JKTEBOP code;
    #  - The equation for phase difference comes from Hilditch (2001) page 238 equation 5.66,
    #  - originally credited to the monograph by Kopal (1959).

    e_fac = numpy.sqrt(1.0 - e**2)
    term1 = 2.0 * numpy.arctan((e * numpy.cos(w)) / e_fac)
    term2 = 2.0 * e * numpy.cos(w) * e_fac / (1.0 - (e * numpy.sin(w)) ** 2)
    phase_diff = (numpy.pi + term1 + term2) / (2.0 * numpy.pi)

    # Calculations below are adopted from;
    #  - Eclipsing Binary Stars: Modeling and Analysis (Kallrath & Milone, 2009, Springer)

    true_anomaly = (numpy.pi / 2.0) - w
    eccentric_anomaly = 2.0 * numpy.arctan(numpy.sqrt((1.0 - e) / (1.0 + e)) * numpy.tan(true_anomaly / 2.0))
    mean_anomaly = eccentric_anomaly - e * numpy.sin(eccentric_anomaly)
    conjunction = ((mean_anomaly + w) / (2.0 * numpy.pi)) - 0.25 + phase_shift
    periastron_passage = 1.0 - mean_anomaly / (2.0 * numpy.pi)
    periastron_phase = conjunction + periastron_passage
    while periastron_phase > 1.0:
        periastron_phase = periastron_phase - int(periastron_phase)

    phase_of_primary_eclipse = conjunction
    phase_of_first_quadrature = conjunction + phase_diff / 2.0
    phase_of_secondary_eclipse = conjunction + phase_diff
    phase_of_second_quadrature = conjunction + (phase_diff / 2.0) + 0.5
    phase_of_periastron = periastron_phase
    phase_of_apastron = periastron_phase + 0.5

    while phase_of_apastron > 1.0:
        phase_of_apastron = phase_of_apastron - int(phase_of_apastron)

    if e == 0.0:
        phase_of_periastron = 0.25
        phase_of_apastron = 0.25

    return phase_of_primary_eclipse, phase_of_first_quadrature, phase_of_secondary_eclipse, \
           phase_of_second_quadrature, phase_of_periastron, phase_of_apastron
