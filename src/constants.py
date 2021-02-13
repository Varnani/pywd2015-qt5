import inspect
import os


# info
MAIN_VERSION = "2015-qt5"
CONFIG_VERSION = "0.2.0"
EXPORT_DELIMITER = "\t"

# paths
MAIN_DIR = os.path.dirname(os.path.abspath(inspect.stack()[-1][1]))  # script dir
PYWD_CONFIG_PATH = os.path.join(MAIN_DIR, "pywd_config.dat")
RESOURCES_DIR = os.path.join(MAIN_DIR, "resources")
MAIN_ICON_PATH = os.path.join(RESOURCES_DIR, "pywd.ico")
PLAY_ICON_PATH = os.path.join(RESOURCES_DIR, "play.png")
PAUSE_ICON_PATH = os.path.join(RESOURCES_DIR, "pause.png")
NEXT_ICON_PATH = os.path.join(RESOURCES_DIR, "next.png")
PREV_ICON_PATH = os.path.join(RESOURCES_DIR, "prev.png")
MONO_FONT_PATH = os.path.join(RESOURCES_DIR, "PTM55FT.ttf")

# hex colors and plot details
MARKER_SIZE = 3
COLOR_BLUE = "#4286f4"
COLOR_RED = "#f5425a"
COLOR_GREEN = "#2aa34b"
COLOR_ORANGE = "#f5a142"
DC_RESULTS_FONTSIZE = 9

# config sections
CONFIG_SECTION_INFO = "Info"
CONFIG_SECTION_MAIN = "Main"
CONFIG_SECTION_SYSTEM = "System"
CONFIG_SECTION_SURFACE = "Surface"
CONFIG_SECTION_3B = "Third_body"
CONFIG_SECTION_LC2015 = "LC2015"
CONFIG_SECTION_DC2015 = "DC2015"
CONFIG_SECTION_KEEPS = "DC Keeps"
CONFIG_SECTION_DELS = "DC Dels"
CONFIG_SECTION_SPOT_PARAMS = "Spot Parameters"
CONFIG_SECTION_ECLIPSE_TIMINGS = "Eclipse Timings"
CONFIG_SECTION_CURVE_COUNTS = "Curve Counts"
CONFIG_SECTION_LIGHT_CURVE_STUB = "Light Curve "
CONFIG_SECTION_VELOCITY_CURVE_STUB = "Velocity Curve"
CONFIG_SECTION_STAR1_SPOTS = "Star 1 Spot "
CONFIG_SECTION_STAR2_SPOTS = "Star 2 Spot "

# dicts
BANDPASS_ID_DICT = {
            "Stromgren u": "1",
            "Stromgren v": "2",
            "Stromgren b": "3",
            "Stromgren y": "4",
            "Johnson U": "5",
            "Johnson B": "6",
            "Johnson V": "7",
            "Johnson R": "8",
            "Johnson I": "9",
            "Johnson J": "10",
            "Johnson K": "11",
            "Johnson L": "12",
            "Johnson M": "13",
            "Johnson N": "14",
            "Cousins Rc": "15",
            "Cousins Ic": "16",
            "Bessel UX": "17",
            "Bessel BX": "18",
            "Bessel B": "19",
            "Bessel V": "20",
            "Bessel R": "21",
            "Bessel I": "22",
            "Tycho Bt": "23",
            "Tycho Vt": "24",
            "Hipparchos Hp": "25",
            "KEPLER": "26",
            "COROT SIS": "27",
            "COROT EXO": "28",
            "Geneva U": "29",
            "Geneva B": "30",
            "Geneva B1": "31",
            "Geneva B2": "32",
            "Geneva V": "33",
            "Geneva V1": "34",
            "Geneva G": "35",
            "Vilnius U": "36",
            "Vilnius P": "37",
            "Vilnius X": "38",
            "Vilnius Y": "39",
            "Vilnius Z": "40",
            "Vilnius V": "41",
            "Vilnius S": "42",
            "Milone iz": "43",
            "Milone iJ": "44",
            "Milone iH": "45",
            "Milone iK": "46",
            "YMS94 iz": "47",
            "YMS94 iJ": "48",
            "YMS94 iH": "49",
            "YMS94 iK": "50",
            "YMS94 iL": "51",
            "YMS94 iL'": "52",
            "YMS94 iM": "53",
            "YMS94 in": "54",
            "YMS94 iN": "55",
            "Sloan DSS u'": "56",
            "Sloan DSS g'": "57",
            "Sloan DSS r'": "58",
            "Sloan DSS i'": "59",
            "Sloan DSS z'": "60",
            "HST STIS Ly alpha": "61",
            "HST STIS Fclear": "62",
            "HST STIS Fsrf2": "63",
            "HST STIS Fqtz": "64",
            "HST STIS C III": "65",
            "HST STIS Mg II": "66",
            "HST STIS Nclear": "67",
            "HST STIS Nsrf2": "68",
            "HST STIS Nqtz": "69",
            "HST STIS cn182": "70",
            "HST STIS cn270": "71",
            "HST STIS Oclear": "72",
            "HST STIS Oclear-lp": "73",
            "HST STIS [O II]": "74",
            "HST STIS [O III]": "75",
            "2MASS J": "76",
            "2MASS H": "77",
            "2MASS Ks": "78",
            "SWASP": "79",
            "MOST": "80",
            "Gaia G (2006)": "81",
            "Gaia G (2010)": "82",
            "Gaia Gbp": "83",
            "Gaia Grp": "84",
            "Gaia Grvs": "85",
            "Milone 230": "86",
            "Milone 250": "87",
            "Milone 270": "88",
            "Milone 290": "89",
            "Milone 310": "90",
            "Milone 330": "91",
            "Ca II triplet": "92",
            "WIRE V+R": "93",
            "Lunar Ultraviolet Telescope": "94"
}
ID_BANDPASS_DICT = {
    value: key for key, value in list(BANDPASS_ID_DICT.items())  # reverse of the above
}
KEEPS_NAME_ID_DICT = {
    "Spot 1 Latitude": 1.0,
    "Spot 1 Longitude": 2.0,
    "Spot 1 Ang. Rad.": 3.0,
    "Spot 1 Temp. Factor": 4.0,

    "Spot 2 Latitude": 5.0,
    "Spot 2 Longitude": 6.0,
    "Spot 2 Ang. Rad.": 7.0,
    "Spot 2 Temp. Factor": 8.0,

    "a": 9.0,
    "e": 10.0,
    "Omega": 11.0,
    "F1": 12.0,
    "F2": 13.0,
    "Phase Shift": 14.0,
    "V Gamma": 15.0,
    "i": 16.0,
    "g1": 17.0,
    "g2": 18.0,
    "T1": 19.0,
    "T2": 20.0,
    "Alb1": 21.0,
    "Alb2": 22.0,
    "Pot1": 23.0,
    "Pot2": 24.0,
    "q (M2/M1)": 25.0,
    "Ephemeris": 26.0,
    "Period": 27.0,
    "dP/dt": 28.0,
    "d(Omega)/dt": 29.0,
    "a (3B)": 30.0,
    "Period (3B)": 31.0,
    "i (3B)": 32.0,
    "e (3B)": 33.0,
    "Omega (3B)": 34.,
    "Ephemeris (3B)": 35.0,
    "Log(d)": 41.0,
    "Desig. Ext.": 42.0,

    "Spot 1 Tstart": 43.0,
    "Spot 1 Tmax1": 44.0,
    "Spot 1 Tmax2": 45.0,
    "Spot 1 Tend": 46.0,

    "Spot 2 Tstart": 47.0,
    "Spot 2 Tmax1": 48.0,
    "Spot 2 Tmax2": 49.0,
    "Spot 2 Tend": 50.0,

    "L1": 56.0,
    "L2": 57.0,
    "X1": 58.0,
    "X2": 59.0,
    "L3": 60.0
}
KEEPS_ID_NAME_DICT = {
    value: key for key, value in list(KEEPS_NAME_ID_DICT.items())
}

ID_LATEX_DICT = {
            "1": "$Co-Latitude_{Spot1}$ $(^{\circ})$",
            "2": "$Longitude_{Spot1}$ $(^{\circ})$",
            "3": "$Radius_{Spot1}$ $(^{\circ})$",
            "4": "$Temperature~Factor_{Spot1}$",
            "5": "$Co-Latitude_{Spot2}$ $(^{\circ})$",
            "6": "$Longitude_{Spot2}$ $(^{\circ})$",
            "7": "$Radius_{Spot2}$ $(^{\circ})$",
            "8": "$Temperature~Factor_{Spot2}$",
            "9": "$a$~({\mbox{$R_{\odot}$}})",
            "10": "$e$",
            "11": "$\omega~(^{\circ})$",
            "12": "$F_{1}$",
            "13": "$F_{2}$",
            "14": "$Phase~Shift$",
            "15": "$V_{\gamma}$",
            "16": "$i~(^{\circ})$",
            "17": "$g_{1}$",
            "18": "$g_{2}$",
            "19": "$T_{1}(K)$",
            "20": "$T_{2}(K)$",
            "21": "$A_{{1}}$",
            "22": "$A_{{2}}$",
            "23": "$\Omega_{1}$",
            "24": "$\Omega_{2}$",
            "25": "$q~(=M{_2}/M{_1})$",
            "26": "$HJD_{Min1}$",
            "27": "$Period$",
            "28": "$dP/dt$",
            "29": "$d\omega/dt$",
            "30": "$a$ ({\mbox{$R_{\odot}$}}) (3B)$",
            "31": "$P~(3B)$",
            "32": "$i~(3B)~(^{\circ})$",
            "33": "$e~(3B)$",
            "34": "$\omega~(3B)$",
            "35": "$Ephemeris~(3B)$",
            "41": "$log(d)$",
            "42": "$Designated~Extinction$",
            "43": "nan",
            "44": "nan",
            "45": "nan",
            "46": "nan",
            "47": "nan",
            "48": "nan",
            "49": "nan",
            "50": "nan",
            "56": "$L_{{1}}$/$(L_{{1}}+L_{{2}})_{{{band}}}$",
            "57": "nan",
            "58": "$x{{_1}}_{{{band}}}$",
            "59": "$x{{_2}}_{{{band}}}$",
            "60": "$L_{{3}}$/$(L_{{1}}+L_{{2}}+L_{{3}})_{{{band}}}$"
        }

JDPHS_DICT = {
    "Time": 1,
    "Phase": 2
}
MAGLITE_DICT = {
    "Flux": 0,
    "Magnitude": 1
}
MODE_DICT = {
     "Mode -1": -1,
     "Mode 0": 0,
     "Mode 1": 1,
     "Mode 2": 2,
     "Mode 3": 3,
     "Mode 4": 4,
     "Mode 5": 5,
     "Mode 6": 6
}
ATM_DICT = {
    "Stellar Atmosphere": 1,
    "Blackbody": 0
}
LIMB_DICT = {
    "Linear Cosine": 1,
    "Logarithmic": 2,
    "Square Root": 3
}
DERIV_DICT = {
    "Symmetrical": 1,
    "Asymmetrical": 0
}
NOMAX_DICT = {
    "Triangular": 1,
    "Trapezoidal": 0
}
NOISE_DICT = {
    "No Scaling": 0,
    "Square Root": 1,
    "Linear": 2
}
