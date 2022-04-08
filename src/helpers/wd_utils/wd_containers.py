class _ParameterContainer:

    class Parameter:
        def __init__(self, name, ftype, value=None):
            self.name = name
            self._value = None
            self.ftype = ftype

            if value is not None:
                self.set(value)

        def set(self, value):
            self._value = value
            # is enforcing type here beneficial?
            # if self.ftype == type(value):
            #     self._value = value
            # else:
            #     raise TypeError("\nName          : " + self.name + "\n"
            #                     "Expected type : " + str(self.ftype) + "\n"
            #                     "Given type    : " + str(type(value)))

        def get(self):
            return self._value

        def format(self, width, precision, exponent, signed=False):
            if self._value is None:
                raise ValueError(self.name + ": value is None.")

            output = ""

            if self.ftype == float:

                if self._value == 0.0:
                    return (" " * (width - 3)) + "0.0"
                elif 1.0 > self._value > -1.0:
                    output = "{:{width}.{precision}g}".format(self._value, width=width, precision=precision)
                    if "." not in output and "e" in output:
                        output = output.split("e")[0].strip(" ") + ".0e" + output.split("e")[1].strip(" ")
                elif self._value >= 1.0 or self._value <= -1.0:
                    output = "{:{width}.{precision}f}".format(self._value, width=width, precision=precision)

                output = output.rstrip("0")
                if output[-1] == "." or output[-2] in ("+", "-"):
                    output = output + "0"
                output = " " * (width - len(output)) + output

                if len(output) > width:
                    raise ValueError("Can't format float:" +
                                     "\nName        : " + self.name +
                                     "\nWidth     : " + str(width) +
                                     "\nPrecision : " + str(precision) +
                                     "\nInput     : " + repr(self._value) +
                                     "\nOutput    : " + output)

                else:
                    return output.replace("e", exponent)

            elif self.ftype == int:
                str_val = str(int(self._value))

                if signed and self._value > 0:
                    str_val = "+" + str_val

                output = (" " * (width - len(str_val))) + str_val

                if len(output) > width:
                    raise ValueError("Can't format integer:" +
                                     "\nName   : " + self.name +
                                     "\nWidth  : " + str(width) +
                                     "\nInput  : " + repr(self._value) +
                                     "\nOutput : " + output)

                else:
                    return output

            else:
                raise TypeError(self.name + ": Can't format " + str(self.ftype) + ", must be either float or integer.")

        def __repr__(self):
            return "Parameter({0}, {1}, value={2})".format(self.name, self.ftype, self.get())

        def __str__(self):
            return "Parameter {0}: {1}, {2}".format(self.name, self.get(), self.ftype)

    def __init__(self, name):
        self.name = name
        self.data = {}

        self._parameter_dict = {}

    def _add_parameter(self, name, ftype, value=None):
        if name not in list(self._parameter_dict.keys()):
            param = _ParameterContainer.Parameter(name, ftype, value=value)
            self._parameter_dict[name] = param
        else:
            raise KeyError(name + " already exists.")

    def _add_data(self, name, *columns):
        cols = list(columns)
        self.data[name] = cols

    def check_values(self):
        none_list = []
        for parameter in self:
            if parameter.get() is None:
                none_list.append(parameter)

        if len(none_list) > 0:
            output = "    "
            i = 0
            for parameter in none_list:
                output = output + parameter.name + ", "
                i = i + 1
                if i >= 5:
                    output = output + "\n    "
                    i = 0
            raise ValueError("Following parameters are not filled:\n" + output)

        else:
            return True

    def __getitem__(self, item):
        if type(item) is not str:
            raise TypeError("Expected as string, but found " + str(type(item)))
        else:
            return self._parameter_dict[item]

    def __setitem__(self, key, value):
        if type(key) == str:
            if key in list(self._parameter_dict.keys()):
                self._parameter_dict[key].set(value)
            else:
                raise ValueError(key + " does not exist.")
        else:
            raise TypeError("Expected a string, but found " + str(type(key)) + " for key: " + str(key))

    def __iter__(self):
        for parameter in list(self._parameter_dict.values()):
            yield parameter

    def __repr__(self):
        return "<ParameterContainer for " + self.name + \
               ", parameters: " + str(len(list(self._parameter_dict.keys()))) + \
               ", data: " + str(len(list(self.data.keys()))) + ">"

    def __str__(self):
        output = "ParameterContainer " + self.name + ":\n"
        i = 0
        for parameter in self:
            output = output + parameter.name + ": " + str(parameter.get()) + ", "
            if i == 8:
                output = output + "\n"
                i = 0
            i = i + 1

        output = output + "\n"

        if len(self.data) != 0:
            output = output + "\nAvailable data:\n"
            for key in list(self.data.keys()):
                output = output + key + ": " + str(len(self.data[key])) + " columns"
        return output


class _CommonParameterContainer(_ParameterContainer):

    def __init__(self, name="Common"):
        _ParameterContainer.__init__(self, name)
        self._populate_general_parameters()
        self.star1_spots = []
        self.star2_spots = []

    def _populate_general_parameters(self):
        # general config
        self._add_parameter("jdphs", int)
        self._add_parameter("ifcgs", int)
        self._add_parameter("mode", int)
        self._add_parameter("icor1", int)
        self._add_parameter("icor2", int)

        # system params
        self._add_parameter("hjd0", float)
        self._add_parameter("pzero", float)
        self._add_parameter("dpdt", float)
        self._add_parameter("pshift", float)
        self._add_parameter("delph", float)
        self._add_parameter("nga", int)
        self._add_parameter("e", float)
        self._add_parameter("a", float)
        self._add_parameter("f1", float)
        self._add_parameter("f2", float)
        self._add_parameter("vga", float)
        self._add_parameter("xincl", float)
        self._add_parameter("tavh", float)
        self._add_parameter("tavc", float)
        self._add_parameter("phsv", float)
        self._add_parameter("pcsv", float)
        self._add_parameter("rm", float)
        self._add_parameter("perr", float)
        self._add_parameter("dperdt", float)
        self._add_parameter("the", float)
        self._add_parameter("vunit", float)
        self._add_parameter("abunin", float)
        self._add_parameter("dpclog", float)

        # surface params
        self._add_parameter("ifat1", int)
        self._add_parameter("ifat2", int)
        self._add_parameter("gr1", float)
        self._add_parameter("gr2", float)

        self._add_parameter("ipb", int)
        self._add_parameter("mref", int)
        self._add_parameter("nref", int)
        self._add_parameter("n1", int)
        self._add_parameter("n2", int)
        self._add_parameter("alb1", float)
        self._add_parameter("alb2", float)
        self._add_parameter("ld1", int)
        self._add_parameter("ld2", int)
        self._add_parameter("xbol1", float)
        self._add_parameter("xbol2", float)
        self._add_parameter("ybol1", float)
        self._add_parameter("ybol2", float)

        # spot params
        self._add_parameter("nomax", int)
        self._add_parameter("kspev", int)
        self._add_parameter("kspot", int)
        self._add_parameter("fspot1", float)
        self._add_parameter("fspot2", float)
        self._add_parameter("ifsmv1", int)
        self._add_parameter("ifsmv2", int)

        # third body
        self._add_parameter("if3b", int)
        self._add_parameter("a3b", float)
        self._add_parameter("p3b", float)
        self._add_parameter("xincl3b", float)
        self._add_parameter("e3b", float)
        self._add_parameter("perr3b", float)
        self._add_parameter("tc3b", float)

    def add_spot(self, star, xlat, xlong, radsp, temsp, tstart, tmax1, tmax2, tfinal):
        spot_container = _ParameterContainer("StarSpots")

        spot_container._add_parameter("xlat", float, value=xlat)
        spot_container._add_parameter("xlong", float, value=xlong)
        spot_container._add_parameter("radsp", float, value=radsp)
        spot_container._add_parameter("temsp", float, value=temsp)
        spot_container._add_parameter("tstart", float, value=tstart)
        spot_container._add_parameter("tmax1", float, value=tmax1)
        spot_container._add_parameter("tmax2", float, value=tmax2)
        spot_container._add_parameter("tfinal", float, value=tfinal)

        if star == 1:
            self.star1_spots.append(spot_container)
        elif star == 2:
            self.star2_spots.append(spot_container)

    def remove_spot(self, star, index):
        if star == 1:
            self.star1_spots.pop(index)
        elif star == 2:
            self.star2_spots.pop(index)

    def reset_spots(self):
        self.star1_spots = []
        self.star2_spots = []

    def __str__(self):
        output = _ParameterContainer.__str__(self) + "\n\nStar 1 spots:\n"

        for spot in self.star1_spots:
            output = output + str(spot)

        output = output + "\nStar 2 spots:\n"

        for spot in self.star2_spots:
            output = output + str(spot)

        return output

    def __repr__(self):
        output = _ParameterContainer.__repr__(self)[:-1] + \
                 ", star1 spots: " + str(len(self.star1_spots)) + \
                 ", star2 spots: " + str(len(self.star2_spots)) + ">"

        return output


class LCParameterContainer(_CommonParameterContainer):

    def __init__(self):
        _CommonParameterContainer.__init__(self, name="LC")
        self.synthetic_curve = None
        self._populate_lc_parameters()
        self.star1_lines = []
        self.star2_lines = []

    def _populate_lc_parameters(self):
        # noise
        self._add_parameter("stdev", float)
        self._add_parameter("noise", int)
        self._add_parameter("seed", float, value=138472375)

        # steps
        self._add_parameter("hjdst", float)
        self._add_parameter("hjdsp", float)
        self._add_parameter("hjdin", float)
        self._add_parameter("phstrt", float)
        self._add_parameter("phstop", float)
        self._add_parameter("phin", float)
        self._add_parameter("phn", float, value=0.25)

        # temperature estimation params
        self._add_parameter("phobs", float)
        self._add_parameter("lsp", int)
        self._add_parameter("tobs", float)

        # line profile parameters
        self._add_parameter("binwm1", float, value=0.00001)
        self._add_parameter("sc1", float, value=1.0)
        self._add_parameter("sl1", float, value=0.0)
        self._add_parameter("nf1", int, value=1)
        self._add_parameter("binwm2", float, value=0.00001)
        self._add_parameter("sc2", float, value=1.0)
        self._add_parameter("sl2", float, value=0.0)
        self._add_parameter("nf2", int, value=1)

        # synthetic curve params
        self.synthetic_curve = _ParameterContainer("SyntheticCurve")
        self.synthetic_curve._add_parameter("iband", int)
        self.synthetic_curve._add_parameter("hla", float)
        self.synthetic_curve._add_parameter("cla", float)
        self.synthetic_curve._add_parameter("x1a", float)
        self.synthetic_curve._add_parameter("x2a", float)
        self.synthetic_curve._add_parameter("y1a", float)
        self.synthetic_curve._add_parameter("y2a", float)
        self.synthetic_curve._add_parameter("el3a", float)
        self.synthetic_curve._add_parameter("opsfa", float)
        self.synthetic_curve._add_parameter("zero", float)
        self.synthetic_curve._add_parameter("factor", float)
        self.synthetic_curve._add_parameter("wla", float)
        self.synthetic_curve._add_parameter("aextinc", float)
        self.synthetic_curve._add_parameter("calib", float)

        self.set_dummy_synthetic_curve()

        self.data["eclipse_times"] = []

    def set_synthetic_curve(self, iband, hl, cl, xh, xc, yh, yc, el3, opsf, zero, factor, wl, aextinc, calib):
        self.synthetic_curve["iband"] = iband
        self.synthetic_curve["hla"] = hl
        self.synthetic_curve["cla"] = cl
        self.synthetic_curve["x1a"] = xh
        self.synthetic_curve["x2a"] = xc
        self.synthetic_curve["y1a"] = yh
        self.synthetic_curve["y2a"] = yc
        self.synthetic_curve["el3a"] = el3
        self.synthetic_curve["opsfa"] = opsf
        self.synthetic_curve["zero"] = zero
        self.synthetic_curve["factor"] = factor
        self.synthetic_curve["wla"] = wl
        self.synthetic_curve["aextinc"] = aextinc
        self.synthetic_curve["calib"] = calib

    def set_dummy_synthetic_curve(self):
        self.synthetic_curve["iband"] = 7
        self.synthetic_curve["hla"] = 1.0
        self.synthetic_curve["cla"] = 1.0
        self.synthetic_curve["x1a"] = 0.0
        self.synthetic_curve["x2a"] = 0.0
        self.synthetic_curve["y1a"] = 0.0
        self.synthetic_curve["y2a"] = 0.0
        self.synthetic_curve["el3a"] = 0.0
        self.synthetic_curve["opsfa"] = 0.0
        self.synthetic_curve["zero"] = 8.0
        self.synthetic_curve["factor"] = 1.0
        self.synthetic_curve["wla"] = 0.55
        self.synthetic_curve["aextinc"] = 0.0
        self.synthetic_curve["calib"] = 0.0

    def add_spectral_line(self, star, wll, ewid, depth, kks):
        line = _ParameterContainer("SpectralLine")

        line._add_parameter("wll", float, value=wll)
        line._add_parameter("ewid", float, value=ewid)
        line._add_parameter("depth", float, value=depth)
        line._add_parameter("kks", int, value=kks)

        if star == 1:
            self.star1_lines.append(line)
        elif star == 2:
            self.star2_lines.append(line)

    def remove_spectral_line(self, star, index):
        if star == 1:
            self.star1_lines.pop(index)
        elif star == 2:
            self.star2_lines.pop(index)

    def reset_spectral_lines(self):
        self.star1_lines = []
        self.star2_lines = []

    def add_eclipse_times(self, times, types):
        self._add_data("eclipse_times", times, types)

    def remove_eclipse_times(self):
        self.data = {}

    def __str__(self):
        output = _CommonParameterContainer.__str__(self) + "\nStar1 spectral lines:\n"
        for line in self.star1_lines:
            output = output + str(line)

        output = output + "\nStar2 spectral lines:\n"
        for line in self.star2_lines:
            output = output + str(line)

        return output

    def __repr__(self):
        output = _CommonParameterContainer.__repr__(self)[:-1] + \
                 ", star1 lines: " + str(len(self.star1_lines)) + \
                 ", star2 lines: " + str(len(self.star2_lines)) + ">"
        return output


class DCParameterContainer(_CommonParameterContainer):

    def __init__(self):
        _CommonParameterContainer.__init__(self, name="DC")
        self.light_curves = []
        self.velocity_curves = [None, None]
        self.eclipse_timings = None
        self.dels = None
        self.keeps = None
        self.subsets = []

        self._populate_dc_parameters()
        self.reset_dels()
        self.reset_keeps()

    def _populate_dc_parameters(self):
        # general dc params
        self._add_parameter("ifvc1", int, value=0)
        self._add_parameter("ifvc2", int, value=0)
        self._add_parameter("nlc", int, value=0)
        self._add_parameter("iftime", int, value=0)
        self._add_parameter("ko", int, value=0)
        self._add_parameter("kdisk", int, value=0)
        self._add_parameter("isym", int, value=1)
        self._add_parameter("nppl", int, value=1)
        self._add_parameter("ifder", int, value=1)
        self._add_parameter("iflcin", int, value=0)
        self._add_parameter("ifoc", int, value=1)
        self._add_parameter("maglite", int, value=None)
        self._add_parameter("linkext", int, value=None)
        self._add_parameter("desextinc", float, value=None)
        self._add_parameter("n1l", int, value=None)
        self._add_parameter("n2l", int, value=None)

        # which spots to fit?
        self._add_parameter("kspa", int, value=0)
        self._add_parameter("nspa", int, value=0)
        self._add_parameter("kspb", int, value=0)
        self._add_parameter("nspb", int, value=0)

    def reset_dels(self):
        # DEL's, they are filled with default values from WD manual example
        dels = _ParameterContainer("DELs")

        dels._add_parameter("spot_a_lat", float, value=0.02)
        dels._add_parameter("spot_a_long", float, value=0.02)
        dels._add_parameter("spot_a_rad", float, value=0.001)
        dels._add_parameter("spot_a_tempf", float, value=0.02)

        dels._add_parameter("spot_b_lat", float, value=0.02)
        dels._add_parameter("spot_b_long", float, value=0.02)
        dels._add_parameter("spot_b_rad", float, value=0.001)
        dels._add_parameter("spot_b_tempf", float, value=0.02)

        dels._add_parameter("a", float, value=0.05)
        dels._add_parameter("e", float, value=0.001)
        dels._add_parameter("perr", float, value=0.01)
        dels._add_parameter("f1", float, value=0.01)
        dels._add_parameter("f2", float, value=0.01)
        dels._add_parameter("pshift", float, value=0.002)
        dels._add_parameter("xincl", float, value=0.2)
        dels._add_parameter("g1", float, value=0.01)
        dels._add_parameter("g2", float, value=0.01)
        dels._add_parameter("tavh", float, value=0.02)
        dels._add_parameter("tavc", float, value=0.02)

        dels._add_parameter("alb1", float, value=0.05)
        dels._add_parameter("alb2", float, value=0.05)
        dels._add_parameter("phsv", float, value=0.02)
        dels._add_parameter("pcsv", float, value=0.02)
        dels._add_parameter("rm", float, value=0.003)
        dels._add_parameter("hla", float, value=0.01)
        dels._add_parameter("cla", float, value=0.01)
        dels._add_parameter("x1a", float, value=0.01)
        dels._add_parameter("x2a", float, value=0.01)

        self.dels = dels

    def reset_keeps(self):
        keeps = self.add_subset(_ret=True)
        keeps.name = "MainSet"
        self.keeps = keeps

    def add_subset(self, _ret=False):
        subset = _ParameterContainer("Subset")

        # KEEP's, all of them are off by default
        subset._add_parameter("spot_a_lat", int, value=1)
        subset._add_parameter("spot_a_long", int, value=1)
        subset._add_parameter("spot_a_rad", int, value=1)
        subset._add_parameter("spot_a_tempf", int, value=1)

        subset._add_parameter("spot_b_lat", int, value=1)
        subset._add_parameter("spot_b_long", int, value=1)
        subset._add_parameter("spot_b_rad", int, value=1)
        subset._add_parameter("spot_b_tempf", int, value=1)

        subset._add_parameter("a", int, value=1)
        subset._add_parameter("e", int, value=1)
        subset._add_parameter("perr", int, value=1)
        subset._add_parameter("f1", int, value=1)
        subset._add_parameter("f2", int, value=1)
        subset._add_parameter("pshift", int, value=1)
        subset._add_parameter("vga", int, value=1)
        subset._add_parameter("xincl", int, value=1)
        subset._add_parameter("g1", int, value=1)
        subset._add_parameter("g2", int, value=1)
        subset._add_parameter("tavc", int, value=1)
        subset._add_parameter("tavh", int, value=1)
        subset._add_parameter("alb1", int, value=1)
        subset._add_parameter("alb2", int, value=1)
        subset._add_parameter("phsv", int, value=1)
        subset._add_parameter("pcsv", int, value=1)
        subset._add_parameter("rm", int, value=1)
        subset._add_parameter("hjd0", int, value=1)
        subset._add_parameter("pzero", int, value=1)
        subset._add_parameter("dpdt", int, value=1)
        subset._add_parameter("dperdt", int, value=1)
        subset._add_parameter("a3b", int, value=1)
        subset._add_parameter("p3b", int, value=1)
        subset._add_parameter("xincl3b", int, value=1)
        subset._add_parameter("e3b", int, value=1)
        subset._add_parameter("perr3b", int, value=1)
        subset._add_parameter("t03b", int, value=1)
        subset._add_parameter("dpclog", int, value=1)
        subset._add_parameter("desextinc", int, value=1)

        subset._add_parameter("spot_a_tstart", int, value=1)
        subset._add_parameter("spot_a_tmax1", int, value=1)
        subset._add_parameter("spot_a_tmax2", int, value=1)
        subset._add_parameter("spot_a_tend", int, value=1)

        subset._add_parameter("spot_b_tstart", int, value=1)
        subset._add_parameter("spot_b_tmax1", int, value=1)
        subset._add_parameter("spot_b_tmax2", int, value=1)
        subset._add_parameter("spot_b_tend", int, value=1)

        subset._add_parameter("hla", int, value=1)
        subset._add_parameter("cla", int, value=1)
        subset._add_parameter("x1a", int, value=1)
        subset._add_parameter("x2a", int, value=1)
        subset._add_parameter("el3a", int, value=1)

        # iteration params
        subset._add_parameter("niter", int, value=1)
        subset._add_parameter("xlamda", float, value=0.00001)
        subset._add_parameter("vlr", float, value=1.0)

        if _ret:
            return subset
        else:
            self.subsets.append(subset)

    def remove_subset(self, index):
        self.subsets.pop(index)

    def add_light_curve(self, iband, hla, cla, x1a, x2a, y1a, y2a, opsfa, sigma,
                        ksd, el3a, noise, aextinc, calib, times, observations, weights,
                        wla, xunit, spha1, spha2, spha3, spha4 ):

        lc = _ParameterContainer("LightCurve")

        lc._add_parameter("iband", int, value=iband)
        lc._add_parameter("hla", float, value=hla)
        lc._add_parameter("cla", float, value=cla)
        lc._add_parameter("x1a", float, value=x1a)
        lc._add_parameter("x2a", float, value=x2a)
        lc._add_parameter("y1a", float, value=y1a)
        lc._add_parameter("y2a", float, value=y2a)
        lc._add_parameter("opsfa", float, value=opsfa)
        lc._add_parameter("sigma", float, value=sigma)
        lc._add_parameter("sphas1", float, value=spha1)
        lc._add_parameter("sphas2", float, value=spha2)
        lc._add_parameter("sphas3", float, value=spha3)
        lc._add_parameter("sphas4", float, value=spha4)
        lc._add_parameter("ksd", int, value=ksd)
        lc._add_parameter("wla", float, value=wla)
        lc._add_parameter("el3a", float, value=el3a)
        lc._add_parameter("noise", int, value=noise)
        lc._add_parameter("aextinc", float, value=aextinc)
        lc._add_parameter("xunit", float, value=xunit)
        lc._add_parameter("calib", float, value=calib)

        lc._add_data("light_data", times, observations, weights)

        self.light_curves.append(lc)
        self["nlc"] = len(self.light_curves)

    def remove_light_curve(self, index):
        self.light_curves.pop(index)
        self["nlc"] = len(self.light_curves)

    def reset_light_curves(self):
        self.light_curves = []
        self["nlc"] = 0

    def add_velocity_curve(self, star, sigma, ksd, wla, times, observations, weights,
                           iband=7, hla=1, cla=1, x1a=0, x2a=0, y1a=0, y2a=0, opsfa=0,
                           sphas1=0.05, sphas2=0.45, sphas3=0.55, sphas4=0.95):

        vc = _ParameterContainer("VelocityCurve")

        vc._add_parameter("iband", int, value=iband)
        vc._add_parameter("hla", float, value=hla)
        vc._add_parameter("cla", float, value=cla)
        vc._add_parameter("x1a", float, value=x1a)
        vc._add_parameter("x2a", float, value=x2a)
        vc._add_parameter("y1a", float, value=y1a)
        vc._add_parameter("y2a", float, value=y2a)
        vc._add_parameter("opsfa", float, value=opsfa)
        vc._add_parameter("sigma", float, value=sigma)
        vc._add_parameter("sphas1", float, value=sphas1)
        vc._add_parameter("sphas2", float, value=sphas2)
        vc._add_parameter("sphas3", float, value=sphas3)
        vc._add_parameter("sphas4", float, value=sphas4)
        vc._add_parameter("ksd", int, value=ksd)
        vc._add_parameter("wla", float, value=wla)

        vc._add_data("velocity_data", times, observations, weights)

        if star == 1:
            self["ifvc1"] = 1
            self.velocity_curves[0] = vc
        elif star == 2:
            self["ifvc2"] = 1
            self.velocity_curves[1] = vc

    def remove_velocity_curve(self, star):
        if star == 1:
            self["ifvc1"] = 0
            self.velocity_curves[0] = None
        elif star == 2:
            self["ifvc2"] = 0
            self.velocity_curves[1] = None

    def reset_velocity_curves(self):
        self.velocity_curves = [None, None]
        self["ifvc1"] = 0
        self["ifvc2"] = 0

    def add_eclipse_times(self, sigma, ksd, times, types, weights):

        et = _ParameterContainer("EclipseTimes")

        et._add_parameter("sigma", float, value=sigma)
        et._add_parameter("ksd", int, value=ksd)

        et._add_data("eclipse_data", times, types, weights)

        self["iftime"] = 1
        self.eclipse_timings = et

    def remove_eclipse_times(self):
        self.eclipse_timings = None
        self["iftime"] = 0

    def __str__(self):
        output = _CommonParameterContainer.__str__(self) + "\nLight Curves:\n"
        for curve in self.light_curves:
            output = output + str(curve)

        output = output + "\n\nVelocity Curves:" \
                          "\nStar 1: " + str(self.velocity_curves[0]) + \
                          "\n---\nStar 2: " + str(self.velocity_curves[1]) + \
                          "\n\nEclipse Timings: " + str(self.eclipse_timings) + \
                          "\n\nKEEP's: \n" + str(self.keeps) + \
                          "\nSubsets: \n" + str(self.subsets) + \
                          "\n\nDEL's: \n" + str(self.dels)

        return output

    def __repr__(self):
        ets = 0
        if self.eclipse_timings is not None:
            ets = 1

        output = _CommonParameterContainer.__repr__(self)[:-1] + \
                 ", light curves: " + str(len(self.light_curves)) + \
                 ", velocity curves: " + str(self["ifvc1"].get() + self["ifvc2"].get()) + \
                 ", eclipse timings: " + str(ets) + \
                 ", subsets: " + str(len(self.subsets)) + ">"
        return output
