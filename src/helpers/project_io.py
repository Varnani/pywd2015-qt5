"""Parse and validate a complete project before mutating live widgets."""
from configparser import ConfigParser, Error
from pathlib import Path
import math
import locale
import os
import tempfile
from src import constants
from src.helpers.project_fields import FIELDS, CURVE_FIELDS
from src.helpers.observations import read_timings


def atomic_write_parser(parser, path):
    target = Path(path)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=target.parent,
                                         prefix=".pywd-", delete=False) as output:
            temporary = output.name
            parser.write(output)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, target)
    finally:
        if temporary is not None and os.path.exists(temporary):
            os.unlink(temporary)


def read_ini_text(path):
    try:
        return Path(path).read_text(encoding="utf-8-sig")
    except UnicodeDecodeError:
        # Releases before this correction wrote INI with the machine locale.
        return Path(path).read_text(encoding=locale.getpreferredencoding(False))


class ProjectModel:
    @classmethod
    def read(cls, path, window):
        model = cls()
        model.parser = ConfigParser(interpolation=None)
        try:
            model.parser.read_string(read_ini_text(path))
            model.validate(path, window)
        except Error as error:
            raise ValueError("Invalid project: " + str(error)) from error
        return model

    def validate(self, path, window):
        parser = self.parser
        for key, expected in (("version", constants.MAIN_VERSION), ("config", constants.CONFIG_VERSION)):
            if parser.get(constants.CONFIG_SECTION_INFO, key) != expected:
                raise ValueError("Unsupported project " + key)
        for name, getter, section, option in FIELDS:
            widget = window
            for component in name.split("."):
                widget = getattr(widget, component)
            value = getattr(parser, getter)(section, option)
            if getter in ("getfloat", "getint"):
                if not math.isfinite(value):
                    raise ValueError(section + "." + option + " must be finite.")
                if hasattr(widget, "minimum") and not widget.minimum() <= value <= widget.maximum():
                    raise ValueError(section + "." + option + " is outside the supported range.")
                if hasattr(widget, "currentIndex") and not 0 <= value < widget.count():
                    raise ValueError(section + "." + option + " is not a valid selection.")
        # Validate all dynamic records independently, including selected observation columns.
        from src.interfaces.curveproperties_interface import CurveProperty
        from PyQt5 import QtWidgets
        from gui.curveproperties_dialog import Ui_CurvePropertiesDialog
        curve_dialog = QtWidgets.QDialog()
        curve_ui = Ui_CurvePropertiesDialog()
        curve_ui.setupUi(curve_dialog)
        self.light_curves = []
        self.velocity_curves = [None, None]
        count = parser.getint(constants.CONFIG_SECTION_CURVE_COUNTS, "light_curves")
        if count < 0 or count > len(parser.sections()):
            raise ValueError("Invalid number of light curves.")
        def resolve(section):
            original = Path(parser.get(section, "filepath"))
            candidates = [original] if original.is_absolute() else [Path(path).parent / original]
            candidates.append(Path(path).parent / original.name)
            found = next((candidate for candidate in candidates if candidate.is_file()), None)
            if found is None:
                raise ValueError("Observation file not found: " + str(original))
            parser.set(section, "filepath", str(found.resolve()))
            return str(found.resolve())
        def curve(section, kind):
            resolve(section)
            result = CurveProperty(kind, window)
            result.populate_from_parser(parser, section)
            if str(result.band_id) not in constants.ID_BANDPASS_DICT:
                raise ValueError("Invalid observation band.")
            for key, value in vars(result).items():
                if isinstance(value, (int, float)) and not math.isfinite(value):
                    raise ValueError("Non-finite curve property: " + key)
            if result.sigma < 0 or result.ksd < 0 or result.noise not in (0, 1, 2):
                raise ValueError("Invalid curve weight/noise setting.")
            for widget_name, attribute in CURVE_FIELDS:
                if attribute in ("col1_idx", "col2_idx", "col3_idx"):
                    continue  # Data-dependent column bounds are checked by get_data below.
                widget = getattr(curve_ui, widget_name)
                value = getattr(result, attribute)
                if hasattr(widget, "minimum") and not widget.minimum() <= value <= widget.maximum():
                    raise ValueError("Curve property is outside its supported range: " + attribute)
                if hasattr(widget, "currentIndex") and not 0 <= value < widget.count():
                    raise ValueError("Invalid curve selection: " + attribute)
            result.get_data()
            return result
        for index in range(count):
            self.light_curves.append(curve(constants.CONFIG_SECTION_LIGHT_CURVE_STUB + str(index + 1), "light"))
        for index in range(2):
            if parser.getboolean(constants.CONFIG_SECTION_CURVE_COUNTS, "velocity_curve_" + str(index + 1)):
                self.velocity_curves[index] = curve(constants.CONFIG_SECTION_VELOCITY_CURVE_STUB + str(index + 1), "velocity")
        for prefix in (constants.CONFIG_SECTION_STAR1_SPOTS, constants.CONFIG_SECTION_STAR2_SPOTS):
            count = parser.getint(constants.CONFIG_SECTION_SPOT_PARAMS, prefix + "count")
            if count < 0 or count > len(parser.sections()):
                raise ValueError("Invalid spot count.")
            for index in range(count):
                section = prefix + str(index + 1)
                for key in ("a", "b"):
                    parser.getboolean(section, key)
                for key in ("lat", "lon", "radsp", "temsp", "tstart", "tmax1", "tmax2", "tend"):
                    value = parser.getfloat(section, key)
                    if not math.isfinite(value) or not 0 <= value <= 9999999:
                        raise ValueError("Invalid spot value: " + key)
        if parser.get(constants.CONFIG_SECTION_ECLIPSE_TIMINGS, "filepath") != "None":
            read_timings(resolve(constants.CONFIG_SECTION_ECLIPSE_TIMINGS),
                         parser.getboolean(constants.CONFIG_SECTION_ECLIPSE_TIMINGS, "cw"))
