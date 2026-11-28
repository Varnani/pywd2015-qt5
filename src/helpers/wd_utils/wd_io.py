from .wd_containers import _ParameterContainer
import os
import shutil
import tempfile
import threading
import time
import subprocess


class _WDIO:
    def __init__(self, container, wd_path, wd_binary_name):
        self.parameters = container
        self._input = ""
        self._binary_directory = os.path.abspath(wd_path)
        self._binary_path = os.path.abspath(os.path.join(wd_path, wd_binary_name))
        self._job = tempfile.TemporaryDirectory(prefix="pywd-job-")
        self._cwd = self._job.name
        self.command_arguments = []
        self.timeout = 3600.0
        self._cancelled = threading.Event()
        self.succeeded = False
        self._type = ""
        self._wd_binary_name = wd_binary_name

        # TODO implement error checking for common input errors
        self.warning = ""
        self.error = ""
        self.has_warning = False
        self.has_error = False

        self.process = None

    def set_working_directory(self, path):
        self._binary_directory = os.path.abspath(path)
        self._binary_path = os.path.join(self._binary_directory, self._wd_binary_name)

    def _get_input_path(self):
        return os.path.join(self._cwd, self._type + "in.active")

    def _get_output_path(self):
        return os.path.join(self._cwd, self._type + "out.active")

    def save(self):
        self.succeeded = False
        # WD opens auxiliary data relative to cwd, not relative to its executable.
        for name in os.listdir(self._binary_directory):
            if name.lower().endswith((".dat", ".tab")):
                source = os.path.join(self._binary_directory, name)
                if os.path.isfile(source):
                    shutil.copyfile(source, os.path.join(self._cwd, name))
        with open(self._get_input_path(), "w") as output:
            output.write(self._input)
        return self

    def cancel(self):
        # Persistent even when the child has not been spawned yet.
        self._cancelled.set()

    def _run_process(self):
        self.succeeded = False
        if self._cancelled.is_set():
            raise RuntimeError("WD job cancelled.")
        if not os.path.isfile(self._binary_path):
            raise OSError("Cannot find WD binary: " + self._binary_path)
        output = self._get_output_path()
        if os.path.lexists(output):
            os.unlink(output)
        started = time.monotonic()
        try:
            with open(os.path.join(self._cwd, "stdout.log"), "wb") as stdout, \
                    open(os.path.join(self._cwd, "stderr.log"), "wb") as stderr:
                self.process = subprocess.Popen([self._binary_path] + self.command_arguments,
                                                cwd=self._cwd, shell=False, stdin=subprocess.DEVNULL,
                                                stdout=stdout, stderr=stderr,
                                                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
                while self.process.poll() is None:
                    if self._cancelled.is_set():
                        raise RuntimeError("WD job cancelled.")
                    if time.monotonic() - started >= self.timeout:
                        raise RuntimeError("WD job timed out.")
                    try:
                        self.process.wait(timeout=0.05)
                    except subprocess.TimeoutExpired:
                        pass
                if self._cancelled.is_set():
                    raise RuntimeError("WD job cancelled.")
                if self.process.returncode != 0:
                    raise RuntimeError("WD exited with code {0}. Logs: {1}".format(self.process.returncode, self._cwd))
            if not os.path.isfile(output) or os.path.getsize(output) == 0:
                raise RuntimeError("WD produced no output file. Logs: " + self._cwd)
            self.succeeded = True
        finally:
            if self.process is not None:
                if self.process.poll() is None:
                    self.process.kill()
                self.process.wait()
                self.process = None
        return self

    def run(self):
        from PyQt5 import QtCore, QtWidgets
        app = QtWidgets.QApplication.instance()
        try:
            if app is not None and QtCore.QThread.currentThread() == app.thread():
                from src.helpers.wd_utils.qt_job import run_with_progress
                run_with_progress(self)
            else:
                self._run_process()
        finally:
            type(self).last_run = self
        return self

    @staticmethod
    def _validate_numeric_table(table, minimum_columns):
        import math
        if len(table) < minimum_columns or not table[0]:
            raise ValueError("Required WD output table is missing or incomplete.")
        if len({len(column) for column in table}) != 1:
            raise ValueError("WD output table has unequal column lengths.")
        if any(not isinstance(value, (int, float)) or not math.isfinite(value)
               for column in table for value in column):
            raise ValueError("WD output table contains invalid or non-finite numbers.")
        return table

    @staticmethod
    def _format_eccentricity(ipt):
        ipt = float(ipt.get())
        if ipt >= 1.0 or ipt < 0.0:
            raise ValueError("Invalid eccentricity value: " + repr(ipt))
        else:
            output = "{:6.5f}".format(ipt)
            return output[1:]

    def _format_spots(self):
        def _format_spot(spt):
            return spt["xlat"].format(9, 5, "F") + \
                   spt["xlong"].format(9, 5, "F") + \
                   spt["radsp"].format(9, 5, "F") + \
                   spt["temsp"].format(9, 5, "F") + \
                   spt["tstart"].format(14, 5, "F") + \
                   spt["tmax1"].format(14, 5, "F") + \
                   spt["tmax2"].format(14, 5, "F") + \
                   spt["tfinal"].format(14, 5, "F") + "\n"

        star1_spot_lines = ""
        for spot in self.parameters.star1_spots:
            star1_spot_lines = star1_spot_lines + _format_spot(spot)

        star2_spot_lines = ""
        for spot in self.parameters.star2_spots:
            star2_spot_lines = star2_spot_lines + _format_spot(spot)

        return star1_spot_lines, star2_spot_lines

    @staticmethod
    def _slice_with_splitmap(line, splitmap, string=False):
        if splitmap[0] != 0:
            splitmap.insert(0, 0)
        splitted_line = []
        i = 0
        while i < len(splitmap) - 1:
            value = line[splitmap[i]:splitmap[i + 1]]
            value = value.rstrip(" ")
            value = value.strip(" ")
            splitted_line.append(_WDIO._tidy_value(value, string=string))
            i = i + 1

        return splitted_line

    @staticmethod
    def _tidy_value(value, string=False):
        if string:
            return value
        else:
            if "*" in value:
                return float("nan")
            else:
                try:
                    return float(value.replace("D", "e"))
                except ValueError:
                    return value

    @staticmethod
    def _tidy_table(table):
        if len(table) == 0:
            return []
        if any(len(row) != len(table[0]) for row in table):
            raise ValueError("WD output table has inconsistent row widths.")
        columns = [[] for _ in table[0]]
        for line in table:
            for index, data in enumerate(line):
                columns[index].append(data)
        return columns

    @staticmethod
    def _read_table(source, header, offset=1, occurence=1, splitmap=None, tidy=True, string=False):
        table = []
        flag = False
        start = 0
        occured = 0
        from contextlib import nullcontext
        with (nullcontext(source) if hasattr(source, "read") else open(source, "r")) as src:
            for line in src:
                if header in line:
                    occured = occured + 1
                    if occured == occurence:
                        flag = True
                if flag is True:
                    if start < offset:
                        start = start + 1
                    else:
                        if not line.strip():
                            break
                        else:
                            if splitmap is not None:
                                table.append(_WDIO._slice_with_splitmap(line, splitmap, string=string))
                            else:
                                table.append([_WDIO._tidy_value(x, string=string) for x in line.split()])

        if tidy:
            return _WDIO._tidy_table(table)
        else:
            return table

    @staticmethod
    def _read_all_tables(source, header, offset=1, splitmap=None, tidy=True, string=False):
        with open(source, "r") as src:
            splitted_source = src.read().split(header)

        if len(splitted_source) == 1:
            return []

        splitted_source.pop(0)  # we do not care about prior data
        tables = []

        for segment in splitted_source:
            splitted_segment = segment.split("\n")
            current_offset = 0
            while offset > current_offset:
                splitted_segment.pop(0)
                current_offset = current_offset + 1
            table = []
            for line in splitted_segment:
                if not line.split():
                    break
                else:
                    if splitmap is not None:
                        table.append(_WDIO._slice_with_splitmap(line, splitmap, string=string))
                    else:
                        table.append([_WDIO._tidy_value(x, string=string) for x in line.split()])
            if tidy:
                tables.append(_WDIO._tidy_table(table))
            else:
                tables.append(table)
        return tables

    def check_container_type(self):
        expectation = None
        if self._type == "lc":
            expectation = "LC"
        elif self._type == "dc":
            expectation = "DC"

        if self.parameters.name != expectation:
            raise TypeError("Expected container: " + expectation + "\n"
                            "Found container: " + self.parameters.name)

    def __str__(self):
        return self._input


class LCIO(_WDIO):
    def __init__(self, container, wd_path=os.getcwd(), lc_binary_name="LC"):
        _WDIO.__init__(self, container, wd_path=wd_path, wd_binary_name=lc_binary_name)
        self._type = "lc"
        self.check_container_type()

    def _fill_input(self, mpage, ktstep=0):

        self.parameters.check_values()
        from src.helpers.observations import validate_interval
        prefix = ("hjdst", "hjdsp", "hjdin") if self.parameters["jdphs"].get() == 1 else ("phstrt", "phstop", "phin")
        validate_interval(*(self.parameters[key].get() for key in prefix))


        line1 = str(mpage) + " " + \
                self.parameters["nref"].format(1, 0, "") + " " + \
                self.parameters["mref"].format(1, 0, "") + " " + \
                self.parameters["ifsmv1"].format(1, 0, "") + " " + \
                self.parameters["ifsmv2"].format(1, 0, "") + " " + \
                self.parameters["icor1"].format(1, 0, "") + " " + \
                self.parameters["icor2"].format(1, 0, "") + " " + \
                self.parameters["if3b"].format(1, 0, "") + " " + \
                self.parameters["ld1"].format(2, 0, "", signed=True) + " " + \
                self.parameters["ld2"].format(2, 0, "", signed=True) + " " + \
                self.parameters["kspev"].format(1, 0, "") + " " + \
                self.parameters["kspot"].format(1, 0, "") + " " + \
                self.parameters["nomax"].format(1, 0, "") + " " + \
                self.parameters["ifcgs"].format(1, 0, "") + " " + \
                ((" " * (6 - len(str(ktstep)))) + str(ktstep)) + "\n"

        line2 = self.parameters["jdphs"].format(1, 0, "") + \
                self.parameters["hjd0"].format(15, 6, "F") + \
                self.parameters["pzero"].format(17, 10, "D") + \
                self.parameters["dpdt"].format(14, 6, "D") + \
                self.parameters["pshift"].format(10, 6, "D") + \
                self.parameters["delph"].format(8, 5, "D") + \
                self.parameters["nga"].format(3, 0, "") + \
                self.parameters["stdev"].format(11, 4, "D") + \
                self.parameters["noise"].format(2, 0, "") + \
                self.parameters["seed"].format(11, 0, "F") + "\n"

        line3 = self.parameters["hjdst"].format(14, 6, "F") + \
                self.parameters["hjdsp"].format(15, 6, "F") + \
                self.parameters["hjdin"].format(13, 6, "F") + \
                self.parameters["phstrt"].format(12, 6, "F") + \
                self.parameters["phstop"].format(12, 6, "F") + \
                self.parameters["phin"].format(12, 6, "F") + \
                self.parameters["phn"].format(12, 6, "F") + \
                self.parameters["phobs"].format(10, 4, "F") + \
                self.parameters["lsp"].format(2, 0, "") + \
                self.parameters["tobs"].format(8, 4, "F") + "\n"

        line4 = self.parameters["mode"].format(2, 0, "") + \
                self.parameters["ipb"].format(2, 0, "") + \
                self.parameters["ifat1"].format(2, 0, "") + \
                self.parameters["ifat2"].format(2, 0, "") + \
                self.parameters["n1"].format(4, 0, "") + \
                self.parameters["n2"].format(4, 0, "") + \
                self.parameters["perr"].format(13, 6, "F") + \
                self.parameters["dperdt"].format(14, 6, "D") + \
                self.parameters["the"].format(8, 5, "F") + \
                self.parameters["vunit"].format(8, 2, "F") + "\n"

        line5 = self._format_eccentricity(self.parameters["e"]) + \
                self.parameters["a"].format(13, 6, "D") + \
                self.parameters["f1"].format(10, 4, "F") + \
                self.parameters["f2"].format(10, 4, "F") + \
                self.parameters["vga"].format(10, 4, "F") + \
                self.parameters["xincl"].format(9, 3, "F") + \
                self.parameters["gr1"].format(7, 3, "F") + \
                self.parameters["gr2"].format(7, 3, "F") + \
                self.parameters["abunin"].format(7, 2, "F") + \
                self.parameters["fspot1"].format(10, 4, "F") + \
                self.parameters["fspot2"].format(10, 4, "F") + "\n"

        tavh_n = _ParameterContainer.Parameter("tavh_n", float, self.parameters["tavh"].get() / 10000.0)
        tavc_n = _ParameterContainer.Parameter("tavc_n", float, self.parameters["tavc"].get() / 10000.0)

        line6 = tavh_n.format(7, 4, "F") + " " + \
                tavc_n.format(7, 4, "F") + \
                self.parameters["alb1"].format(7, 3, "F") + \
                self.parameters["alb2"].format(7, 3, "F") + \
                self.parameters["phsv"].format(13, 6, "D") + \
                self.parameters["pcsv"].format(13, 6, "D") + \
                self.parameters["rm"].format(13, 6, "D") + \
                self.parameters["xbol1"].format(7, 3, "F") + \
                self.parameters["xbol2"].format(7, 3, "F") + \
                self.parameters["ybol1"].format(7, 3, "F") + \
                self.parameters["ybol2"].format(7, 3, "F") + \
                self.parameters["dpclog"].format(8, 5, "F") + "\n"

        line7 = self.parameters["a3b"].format(12, 6, "D") + \
                self.parameters["p3b"].format(14, 7, "D") + \
                self.parameters["xincl3b"].format(11, 5, "F") + \
                self.parameters["e3b"].format(9, 6, "F") + \
                self.parameters["perr3b"].format(10, 7, "F") + \
                self.parameters["tc3b"].format(17, 8, "F") + "\n"

        line8 = self.parameters.synthetic_curve["iband"].format(3, 0, "") + \
                self.parameters.synthetic_curve["hla"].format(13, 7, "D") + \
                self.parameters.synthetic_curve["cla"].format(13, 7, "D") + \
                self.parameters.synthetic_curve["x1a"].format(7, 3, "F") + \
                self.parameters.synthetic_curve["x2a"].format(7, 3, "F") + \
                self.parameters.synthetic_curve["y1a"].format(7, 3, "F") + \
                self.parameters.synthetic_curve["y2a"].format(7, 3, "F") + \
                self.parameters.synthetic_curve["el3a"].format(12, 4, "D") + \
                self.parameters.synthetic_curve["opsfa"].format(11, 4, "D") + \
                self.parameters.synthetic_curve["zero"].format(8, 3, "F") + \
                self.parameters.synthetic_curve["factor"].format(8, 4, "F") + \
                self.parameters.synthetic_curve["wla"].format(10, 6, "F") + \
                self.parameters.synthetic_curve["aextinc"].format(8, 4, "F") + \
                self.parameters.synthetic_curve["calib"].format(12, 5, "D") + "\n"

        star1_line_profiles = ""
        star2_line_profiles = ""

        if mpage == 3:
            star1_line_profiles = self.parameters["binwm1"].format(11, 5, "D") + \
                                  self.parameters["sc1"].format(9, 4, "F") + \
                                  self.parameters["sl1"].format(9, 2, "F") + \
                                  self.parameters["nf1"].format(3, 0, "") + "\n"

            for line in self.parameters.star1_lines:
                star1_line_profiles = star1_line_profiles + \
                                      line["wll"].format(9, 6, "F") + \
                                      line["ewid"].format(12, 5, "D") + \
                                      line["depth"].format(10, 5, "F") + \
                                      line["kks"].format(5, 0, "") + "\n"

            star1_line_profiles = star1_line_profiles + "-1.\n"

            star2_line_profiles = self.parameters["binwm2"].format(11, 5, "D") + \
                                  self.parameters["sc2"].format(9, 4, "F") + \
                                  self.parameters["sl2"].format(9, 2, "F") + \
                                  self.parameters["nf2"].format(3, 0, "") + "\n"

            for line in self.parameters.star2_lines:
                star2_line_profiles = star2_line_profiles + \
                                      line["wll"].format(9, 6, "F") + \
                                      line["ewid"].format(12, 5, "D") + \
                                      line["depth"].format(10, 5, "F") + \
                                      line["kks"].format(5, 0, "") + "\n"

            star2_line_profiles = star2_line_profiles + "-1.\n"

        star1_spots, star2_spots = self._format_spots()

        eclipse_data = ""
        if mpage == 6 and ktstep == 0:

            if len(self.parameters.data["eclipse_times"]) == 0:
                raise ValueError("Eclipse times must be provided for mpage: 6, ktstep: 0")

            jd_formatter = _ParameterContainer.Parameter("jd", float)
            type_formatter = _ParameterContainer.Parameter("type", int)

            jd_list, type_list = self.parameters.data["eclipse_times"]

            for data in zip(jd_list, type_list):
                jd_formatter.set(data[0])
                type_formatter.set(data[1])

                eclipse_data = eclipse_data + jd_formatter.format(14, 5, "F") + type_formatter.format(6, 0, "") + "\n"

            eclipse_data = eclipse_data + "-10000.\n"

        self._input = line1 + line2 + line3 + line4 + line5 + line6 + line7 + line8 + \
                star1_line_profiles + star2_line_profiles + \
                star1_spots + \
                "300.00000  0.00000  0.00000  0.00000       0.00000       0.00000       0.00000       0.00000\n" + \
                star2_spots + \
                "300.00000  0.00000  0.00000  0.00000       0.00000       0.00000       0.00000       0.00000\n" + \
                "150.\n" + \
                eclipse_data + \
                "9"

        return self

    def fill_for_synthetic_light_curve(self):
        return self._fill_input(1)

    def fill_for_synthetic_velocity_curve(self):
        return self._fill_input(2)

    def fill_for_spectral_lines(self):
        return self._fill_input(3)

    def fill_for_component_dimensions(self):
        return self._fill_input(4)

    def fill_for_star_positions(self):
        return self._fill_input(5)

    def fill_for_etv(self):
        return self._fill_input(6)

    def fill_for_conjunction(self, ktstep):
        return self._fill_input(6, ktstep=ktstep)

    def read_synthetic_light_curve(self):
        lc = self._read_table(self._get_output_path(),
                              "      JD               Phase     light 1       light 2")
        return self._validate_numeric_table(lc, 4)

    def read_cgs_synthetic_light_curve(self):
        lc = self._read_table(self._get_output_path(),
                              "      JD               Phase       cgs1          cgs2          cgstot")
        return self._validate_numeric_table(lc, 4)

    def read_synthetic_velocity_curve(self):
        vc = self._read_table(self._get_output_path(),
                              "      JD              Phase     V Rad 1")
        return self._validate_numeric_table(vc, 4)

    def read_spectral_lines(self):
        star1_spec_lines = self._read_all_tables(self._get_output_path(),
                                                 "                              star 1\n",
                                                 offset=2)
        star2_spec_lines = self._read_all_tables(self._get_output_path(),
                                                 "                              star 2\n",
                                                 offset=2)
        return star1_spec_lines, star2_spec_lines

    def read_component_dimensions(self):
        dimensions = self._read_table(self._get_output_path(),
                                      "      JD             Phase     r1pol      r1pt")
        return self._validate_numeric_table(dimensions, 10)

    def read_star_positions(self):
        positions = self._read_all_tables(self._get_output_path(),
                                          "   Y Sky Coordinate    Z Sky Coordinate\n")
        if not positions:
            raise ValueError("LC returned no star positions.")
        return [self._validate_numeric_table(table, 2) for table in positions]

    def read_etv(self):
        etv = self._read_table(self._get_output_path(),
                               "eclipse timing   type         wt.",
                               offset=2)
        return self._validate_numeric_table(etv, 6)

    def read_conjunction(self):
        conjunction = self._read_table(self._get_output_path(),
                                       "conj. time   type         wt.",
                                       offset=2)
        return self._validate_numeric_table(conjunction, 2)

    def read_abs_params(self):
        abs_params = self._read_table(self._get_output_path(),
                                       " Star         M/Msun   (Mean Radius)/Rsun     M Bol    Log g (cgs)")

        teffs = self._read_table(self._get_output_path(),
                                       "  T1      T2     Alb 1  Alb 2")
        sma = self._read_table(self._get_output_path(),
                                       "  ecc     s-m axis       F1         F2       Vgam")
        lds = self._read_table(self._get_output_path(),
                                       "band      x1        x2        y1        y2")
        lums = self._read_table(self._get_output_path(),
                                       "band         L1           L2         x1      x2      y1      y2")
        return abs_params, teffs, sma, lds, lums

    def read_K1_2_params(self):
        par_set_1 = self._read_table(self._get_output_path(),
                                       "JDPHS     J.D. zero       P zero           dPdt      Ph. shift")
        par_set_2 = self._read_table(self._get_output_path(),
                                       "  ecc     s-m axis       F1         F2       Vgam       Incl")
        par_set_3 = self._read_table(self._get_output_path(),
                                       "  T1      T2     Alb 1  Alb 2    Pot 1        Pot 2           M2/M1")

        p, e, a, i, q = float(par_set_1[2][0]), float(par_set_2[0][0]), float(par_set_2[1][0]), \
                         float(par_set_2[5][0]), float(par_set_3[6][0])
        return p, e, a, i, q


class DCIO(_WDIO):
    def __init__(self, container, wd_path=os.getcwd(), dc_binary_name="DC"):
        _WDIO.__init__(self, container, wd_path=wd_path, wd_binary_name=dc_binary_name)
        self._type = "dc"
        self.check_container_type()

    def fill_for_solution(self):
        def _format_keeps(keep):
            block1 = " " + keep["spot_a_lat"].format(1, 0, "") + \
                     keep["spot_a_long"].format(1, 0, "") + \
                     keep["spot_a_rad"].format(1, 0, "") + \
                     keep["spot_a_tempf"].format(1, 0, "") + " "

            block2 = keep["spot_b_lat"].format(1, 0, "") + \
                     keep["spot_b_long"].format(1, 0, "") + \
                     keep["spot_b_rad"].format(1, 0, "") + \
                     keep["spot_b_tempf"].format(1, 0, "") + " "

            block3 = keep["a"].format(1, 0, "") + \
                     keep["e"].format(1, 0, "") + \
                     keep["perr"].format(1, 0, "") + \
                     keep["f1"].format(1, 0, "") + \
                     keep["f2"].format(1, 0, "") + \
                     keep["pshift"].format(1, 0, "") + \
                     keep["vga"].format(1, 0, "") + " "

            block4 = keep["xincl"].format(1, 0, "") + \
                     keep["g1"].format(1, 0, "") + \
                     keep["g2"].format(1, 0, "") + \
                     keep["tavh"].format(1, 0, "") + \
                     keep["tavc"].format(1, 0, "") + " "

            block5 = keep["alb1"].format(1, 0, "") + \
                     keep["alb2"].format(1, 0, "") + \
                     keep["phsv"].format(1, 0, "") + \
                     keep["pcsv"].format(1, 0, "") + \
                     keep["rm"].format(1, 0, "") + " "

            block6 = keep["hjd0"].format(1, 0, "") + \
                     keep["pzero"].format(1, 0, "") + \
                     keep["dpdt"].format(1, 0, "") + \
                     keep["dperdt"].format(1, 0, "") + \
                     keep["a3b"].format(1, 0, "") + " "

            block7 = keep["p3b"].format(1, 0, "") + \
                     keep["xincl3b"].format(1, 0, "") + \
                     keep["e3b"].format(1, 0, "") + \
                     keep["perr3b"].format(1, 0, "") + \
                     keep["t03b"].format(1, 0, "") + " "

            block8 = "11111 "  # unused block

            block9 = keep["dpclog"].format(1, 0, "") + \
                     keep["desextinc"].format(1, 0, "") + \
                     keep["spot_a_tstart"].format(1, 0, "") + \
                     keep["spot_a_tmax1"].format(1, 0, "") + \
                     keep["spot_a_tmax2"].format(1, 0, "") + " "

            block10 = keep["spot_a_tend"].format(1, 0, "") + \
                      keep["spot_b_tstart"].format(1, 0, "") + \
                      keep["spot_b_tmax1"].format(1, 0, "") + \
                      keep["spot_b_tmax2"].format(1, 0, "") + \
                      keep["spot_b_tend"].format(1, 0, "") + " "

            block11 = "11111 "  # unused block

            block12 = keep["hla"].format(1, 0, "") + \
                      keep["cla"].format(1, 0, "") + \
                      keep["x1a"].format(1, 0, "") + \
                      keep["x2a"].format(1, 0, "") + \
                      keep["el3a"].format(1, 0, "") + " "

            block13 = keep["niter"].format(2, 0, "") + \
                      keep["xlamda"].format(10, 3, "D") + \
                      keep["vlr"].format(6, 3, "F") + "\n"

            return block1 + block2 + block3 + block4 + block5 + \
                   block6 + block7 + block8 + block9 + block10 + \
                   block11 + block12 + block13

        def _format_lc_vc_data(x, y, w):
            data_line = ""

            time_formatter = _ParameterContainer.Parameter("time", float)
            observation_formatter = _ParameterContainer.Parameter("obs", float)
            weight_formatter = _ParameterContainer.Parameter("weight", float)

            for xyw in zip(x, y, w):
                time_formatter.set(xyw[0])
                observation_formatter.set(xyw[1])
                weight_formatter.set(xyw[2])

                data_line = data_line + \
                            time_formatter.format(14, 5, "D") + \
                            observation_formatter.format(11, 6, "D") + \
                            weight_formatter.format(8, 3, "D") + "\n"

            return data_line + "  -10001.00000\n"

        def _format_velocity_curve(vc):
            if vc is None:
                return "", ""

            else:
                vc_info_line = vc["iband"].format(3, 0, "") + \
                               vc["hla"].format(13, 6, "D") + \
                               vc["cla"].format(13, 6, "D") + \
                               vc["x1a"].format(7, 3, "F") + \
                               vc["x2a"].format(7, 3, "F") + \
                               vc["y1a"].format(7, 3, "F") + \
                               vc["y2a"].format(7, 3, "F") + \
                               vc["opsfa"].format(10, 3, "D") + \
                               vc["sigma"].format(12, 5, "D") + \
                               vc["sphas1"].format(8, 5, "F") + \
                               vc["sphas2"].format(8, 5, "F") + \
                               vc["sphas3"].format(8, 5, "F") + \
                               vc["sphas4"].format(8, 5, "F") + \
                               vc["wla"].format(10, 6, "F") + \
                               vc["ksd"].format(2, 0, "") + "\n"

                x, y, w = vc.data["velocity_data"]

                vc_data_line = _format_lc_vc_data(x, y, w)

                return vc_info_line, vc_data_line

        def _format_light_curve(lc):
            if lc is None:
                return "", "", ""

            else:
                lc_info_line = lc["iband"].format(3, 0, "") + \
                               lc["hla"].format(13, 6, "D") + \
                               lc["cla"].format(13, 6, "D") + \
                               lc["x1a"].format(7, 3, "F") + \
                               lc["x2a"].format(7, 3, "F") + \
                               lc["y1a"].format(7, 3, "F") + \
                               lc["y2a"].format(7, 3, "F") + \
                               lc["el3a"].format(12, 4, "D") + \
                               lc["opsfa"].format(10, 3, "D") + \
                               lc["noise"].format(2, 0, "") + \
                               lc["sigma"].format(12, 5, "D") + \
                               lc["sphas1"].format(8, 5, "F") + \
                               lc["sphas2"].format(8, 5, "F") + \
                               lc["sphas3"].format(8, 5, "F") + \
                               lc["sphas4"].format(8, 5, "F") + \
                               lc["ksd"].format(2, 0, "") + "\n"

                lc_extra_line = lc["wla"].format(9, 6, "F") + \
                                lc["aextinc"].format(8, 4, "F") + \
                                lc["xunit"].format(11, 4, "D") + \
                                lc["calib"].format(12, 5, "D") + "\n"

                x, y, w = lc.data["light_data"]

                lc_data_line = _format_lc_vc_data(x, y, w)

                return lc_info_line, lc_extra_line, lc_data_line

        # all del's use same formatting
        del_width = 7
        del_precision = 4
        del_exponent = "d"

        del1 = " " + self.parameters.dels["spot_a_lat"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["spot_a_long"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["spot_a_rad"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["spot_a_tempf"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["spot_b_lat"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["spot_b_long"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["spot_b_rad"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["spot_b_tempf"].format(del_width, del_precision, del_exponent) + "\n"

        del2 = " " + self.parameters.dels["a"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["e"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["perr"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["f1"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["f2"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["pshift"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["xincl"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["g1"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["g2"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["tavh"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["tavc"].format(del_width, del_precision, del_exponent) + " " + "\n"

        del3 = " " + self.parameters.dels["alb1"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["alb2"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["phsv"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["pcsv"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["rm"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["hla"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["cla"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["x1a"].format(del_width, del_precision, del_exponent) + " " + \
                self.parameters.dels["x2a"].format(del_width, del_precision, del_exponent) + "\n"

        keeps = _format_keeps(self.parameters.keeps)

        line5 = self.parameters["kspa"].format(3, 0, "") + \
                self.parameters["nspa"].format(3, 0, "") + \
                self.parameters["kspb"].format(3, 0, "") + \
                self.parameters["nspb"].format(3, 0, "") + "\n"

        line6 = self.parameters["ifvc1"].format(1, 0, "") + " " + \
                self.parameters["ifvc2"].format(1, 0, "") + " " + \
                self.parameters["nlc"].format(2, 0, "") + \
                self.parameters["iftime"].format(2, 0, "") + \
                self.parameters["ko"].format(2, 0, "") + \
                self.parameters["kdisk"].format(2, 0, "") + \
                self.parameters["isym"].format(2, 0, "") + \
                self.parameters["nppl"].format(2, 0, "") + \
                self.parameters["ifder"].format(2, 0, "") + \
                self.parameters["iflcin"].format(2, 0, "") + \
                self.parameters["ifoc"].format(2, 0, "") + "\n"

        line7 = self.parameters["nref"].format(1, 0, "") + " " + \
                self.parameters["mref"].format(1, 0, "") + " " + \
                self.parameters["ifsmv1"].format(1, 0, "") + " " + \
                self.parameters["ifsmv2"].format(1, 0, "") + " " + \
                self.parameters["icor1"].format(1, 0, "") + " " + \
                self.parameters["icor2"].format(1, 0, "") + " " + \
                self.parameters["if3b"].format(1, 0, "") + " " + \
                self.parameters["ld1"].format(2, 0, "", signed=True) + " " +  \
                self.parameters["ld2"].format(2, 0, "", signed=True) + " " +  \
                self.parameters["kspev"].format(1, 0, "") + " " + \
                self.parameters["kspot"].format(1, 0, "") + " " + \
                self.parameters["nomax"].format(1, 0, "") + " " + \
                self.parameters["ifcgs"].format(1, 0, "") + " " + \
                self.parameters["maglite"].format(1, 0, "") + " " + \
                self.parameters["linkext"].format(1, 0, "") + " " + \
                self.parameters["desextinc"].format(7, 4, "F") + "\n"

        line8 = self.parameters["jdphs"].format(1, 0, "") + \
                self.parameters["hjd0"].format(15, 6, "F") + \
                self.parameters["pzero"].format(17, 10, "D") + \
                self.parameters["dpdt"].format(14, 6, "D") + \
                self.parameters["pshift"].format(10, 6, "D") + \
                self.parameters["delph"].format(8, 5, "D") + \
                self.parameters["nga"].format(3, 0, "") + "\n"

        line9 = self.parameters["mode"].format(2, 0, "") + \
                self.parameters["ipb"].format(2, 0, "") + \
                self.parameters["ifat1"].format(2, 0, "") + \
                self.parameters["ifat2"].format(2, 0, "") + \
                self.parameters["n1"].format(4, 0, "") + \
                self.parameters["n2"].format(4, 0, "") + \
                self.parameters["n1l"].format(4, 0, "") + \
                self.parameters["n2l"].format(4, 0, "") + \
                self.parameters["perr"].format(13, 6, "F") + \
                self.parameters["dperdt"].format(13, 5, "D") + \
                self.parameters["the"].format(8, 5, "F") + \
                self.parameters["vunit"].format(9, 3, "F") + "\n"

        line10 = self._format_eccentricity(self.parameters["e"]) + \
                 self.parameters["a"].format(13, 6, "D") + \
                 self.parameters["f1"].format(10, 4, "F") + \
                 self.parameters["f2"].format(10, 4, "F") + \
                 self.parameters["vga"].format(10, 4, "F") + \
                 self.parameters["xincl"].format(9, 3, "F") + \
                 self.parameters["gr1"].format(7, 3, "F") + \
                 self.parameters["gr2"].format(7, 3, "F") + \
                 self.parameters["abunin"].format(7, 2, "F") + \
                 self.parameters["fspot1"].format(10, 4, "F") + \
                 self.parameters["fspot2"].format(10, 4, "F") + "\n"

        tavh_n = _ParameterContainer.Parameter("tavh_n", float, self.parameters["tavh"].get() / 10000.0)
        tavc_n = _ParameterContainer.Parameter("tavc_n", float, self.parameters["tavc"].get() / 10000.0)

        line11 = tavh_n.format(7, 4, "F") + \
                 tavc_n.format(8, 4, "F") + \
                 self.parameters["alb1"].format(7, 3, "F") + \
                 self.parameters["alb2"].format(7, 3, "F") + \
                 self.parameters["phsv"].format(13, 6, "D") + \
                 self.parameters["pcsv"].format(13, 6, "D") + \
                 self.parameters["rm"].format(13, 6, "D") + \
                 self.parameters["xbol1"].format(7, 3, "F") + \
                 self.parameters["xbol2"].format(7, 3, "F") + \
                 self.parameters["ybol1"].format(7, 3, "F") + \
                 self.parameters["ybol2"].format(7, 3, "F") + \
                 self.parameters["dpclog"].format(9, 5, "F") + "\n"

        line12 = self.parameters["a3b"].format(12, 6, "D") + \
                 self.parameters["p3b"].format(14, 7, "D") + \
                 self.parameters["xincl3b"].format(11, 5, "F") + \
                 self.parameters["e3b"].format(9, 6, "F") + \
                 self.parameters["perr3b"].format(10, 7, "F") + \
                 self.parameters["tc3b"].format(17, 8, "F") + "\n"

        star1_spots, star2_spots = self._format_spots()

        vc1_dependent_line, vc1_data = _format_velocity_curve(self.parameters.velocity_curves[0])
        vc2_dependent_line, vc2_data = _format_velocity_curve(self.parameters.velocity_curves[1])

        lc_dependent_lines = ""
        lc_extra_dependent_lines = ""
        lc_data = ""
        for lc_container in self.parameters.light_curves:
            info, extra, data = _format_light_curve(lc_container)
            lc_dependent_lines = lc_dependent_lines + info
            lc_extra_dependent_lines = lc_extra_dependent_lines + extra
            lc_data = lc_data + data

        eclipse_line = ""
        eclipse_data = ""
        if self.parameters.eclipse_timings is not None:
            eclipse_line = (" " * 82) + \
                           self.parameters.eclipse_timings["sigma"].format(10,8,"F") + \
                           (" " * 34) + \
                           self.parameters.eclipse_timings["ksd"].format(1,1,"") + "\n"

            hjd_formatter = _ParameterContainer.Parameter("hjd", float)
            type_formatter = _ParameterContainer.Parameter("type", int)
            weights_formatter = _ParameterContainer.Parameter("weights", float)

            x, y, z = self.parameters.eclipse_timings.data["eclipse_data"][0], \
                      self.parameters.eclipse_timings.data["eclipse_data"][1], \
                      self.parameters.eclipse_timings.data["eclipse_data"][2]
            for xyz in zip(x,y,z):
                hjd_formatter.set(xyz[0])
                type_formatter.set(xyz[1])
                weights_formatter.set(xyz[2])

                eclipse_data = eclipse_data + \
                               hjd_formatter.format(14, 5, "D") + \
                               type_formatter.format(6, 0, "") + \
                               weights_formatter.format(13, 3, "D") + "\n"

            eclipse_data = eclipse_data + "  -10001.00000\n"

        subset_line = ""
        for subset in self.parameters.subsets:
            subset_line = subset_line + _format_keeps(subset)

        self._input = del1 + del2 + del3 + keeps + \
                       line5 + line6 + line7 + line8 + line9 + line10 + line11 + line12 + \
                       vc1_dependent_line + vc2_dependent_line + lc_dependent_lines + \
                       eclipse_line + lc_extra_dependent_lines + \
                       star1_spots + "300.00000\n" + star2_spots + "300.00000\n150.\n" + \
                       vc1_data + vc2_data + lc_data + eclipse_data + subset_line + " 2\n"\

        return self

    def read_extinction_values(self):
        extinction_values = self._read_table(self._get_output_path(),
                                 "band   Wave L    extinc      xunit      calib")
        return extinction_values


    def read_results(self, force_tidy_output=False):
        with open(self._get_output_path()) as source:
            count = source.read().count("Input-Output in F Format")
        if count < self.parameters.keeps["niter"].get():
            raise ValueError("DC output is missing completed iterations.")
        results = self._read_table(self._get_output_path(), "Input-Output in F Format", offset=3,
                                   splitmap=[5, 9, 28, 46, 65, 83], occurence=count, tidy=False)
        if not results:
            raise ValueError("DC output has no fitted parameters.")
        self._validate_numeric_table(self._tidy_table(results), 6)
        return self._tidy_table(results) if force_tidy_output else results

    def _final_result_tail(self):
        with open(self._get_output_path()) as source:
            text = source.read()
        position = text.rfind("Input-Output in F Format")
        if position == -1:
            raise ValueError("DC has no final solution block.")
        return text[position:]

    def read_solution_stats(self):
        from io import StringIO
        header = "   Mean residual for input values"
        return self._validate_numeric_table(
            self._read_table(StringIO(self._final_result_tail()), header), 3)

    def read_curve_stats(self):
        import math
        from io import StringIO
        with open(self._get_output_path()) as source:
            text = source.read()
        result_header = "Input-Output in F Format"
        result_position = text.rfind(result_header)
        previous_result = text.rfind(result_header, 0, result_position)
        header = "    Curve   No. of obs.        Std. dev."
        position = text.rfind(header, 0, result_position)
        if result_position == -1 or position <= previous_result:
            raise ValueError("No curve statistics for the selected DC solution.")
        curve_stats = self._read_table(StringIO(text[position:result_position]), header, tidy=False)
        # Eclipse, velocity and light curves have 3, 5 and 7 fields respectively.
        if not curve_stats or any(
                len(row) not in (3, 5, 7) or
                any(not isinstance(value, (int, float)) or not math.isfinite(value) for value in row)
                for row in curve_stats):
            raise ValueError("Malformed or non-finite DC curve statistics.")
        return curve_stats

    def read_component_dimensions(self):
        from io import StringIO
        text = self._final_result_tail()
        dimensions = []
        for star in (1, 2):
            header = "  {0}   pole".format(star)
            if header not in text:
                raise ValueError("Final DC component dimensions are missing.")
            table = self._read_table(StringIO(text), header, offset=0,
                                     splitmap=[3, 10, 24, 38, 52, 66])
            if len(table) != 6:
                raise ValueError("Malformed DC component dimensions.")
            self._validate_numeric_table([table[0]] + table[2:], 5)
            dimensions.append(table)
        return dimensions

    def read_unweighted_observations(self, split_by_observation=False):
        import math
        results = self.read_results()
        with open(self._get_output_path()) as source:
            text = source.read()
        result_header = "Input-Output in F Format"
        result_position = text.rfind(result_header)
        previous_result = text.rfind(result_header, 0, result_position)
        equation_header = "Unweighted Observational Equations"
        equation_position = text.rfind(equation_header, 0, result_position)
        if equation_position <= previous_result:
            raise ValueError("No observational equations for the selected DC solution.")
        ifoc = self.parameters["ifoc"].get()
        if ifoc not in (0, 1):
            raise ValueError("Unsupported DC ifoc value.")
        base_columns = (2 if self.parameters["jdphs"].get() == 1 else 1) + 2 * ifoc + 1
        expected = len(results) + base_columns
        lines = text[equation_position:].splitlines()[3:]
        table, pending = [], []
        for line in lines:
            if not line.strip():
                break
            # Native formats: F20.10 phase records; F13.5 then F20.10 time records.
            # Continuations can use any wrap width; logical width comes from this solution.
            first_width = 13 if self.parameters["jdphs"].get() == 1 else 20
            fields = [line[:first_width]]
            fields.extend(line[i:i + 20] for i in range(first_width, len(line.rstrip()), 20))
            try:
                values = [float(field.replace("D", "E")) for field in fields if field.strip()]
            except ValueError as error:
                raise ValueError("Malformed DC equation field.") from error
            if not all(math.isfinite(value) for value in values):
                raise ValueError("Non-finite DC equation.")
            pending.extend(values)
            if len(pending) > expected:
                raise ValueError("Too many fields in DC observational equation.")
            if len(pending) == expected:
                table.append(pending)
                pending = []
        if pending or not table:
            raise ValueError("Incomplete DC observational equation.")
        expected_rows = sum(len(curve.data["velocity_data"][0]) for curve in self.parameters.velocity_curves if curve is not None)
        expected_rows += sum(len(curve.data["light_data"][0]) for curve in self.parameters.light_curves)
        if self.parameters.eclipse_timings is not None:
            expected_rows += len(self.parameters.eclipse_timings.data["eclipse_data"][0])
        if expected_rows and len(table) != expected_rows:
            raise ValueError("DC equation count does not match input observations.")
        oc_table = self._tidy_table(table)

        if split_by_observation:
            obs_table = []
            split_table = []
            limit = 0
            if self.parameters.velocity_curves[0] is not None:
                vc1_len = len(self.parameters.velocity_curves[0].data["velocity_data"][0])
                split_table.append([limit, limit + vc1_len])
                limit = limit + vc1_len #+ 1
            if self.parameters.velocity_curves[1] is not None:
                vc2_len = len(self.parameters.velocity_curves[1].data["velocity_data"][0])
                split_table.append([limit, limit + vc2_len])
                limit = limit + vc2_len #+ 1
            for lc in self.parameters.light_curves:
                lc_len = len(lc.data["light_data"][0])
                split_table.append([limit, limit + lc_len])
                limit = limit + lc_len #+ 1
            for split in split_table:
                temp_table = []
                for column in oc_table:
                    temp_table.append(column[split[0]:split[1]])
                obs_table.append(temp_table)
            return obs_table
        else:
            return oc_table

    def update_from_results(self):
        # TODO implement this
        raise NotImplementedError
