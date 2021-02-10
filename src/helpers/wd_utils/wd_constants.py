# a dictionary mapping names to ID's of DC keeps
DC_KEEPS_NAME_ID_DICT = {
    "spot_a_lat": 1.0,
    "spot_a_long": 2.0,
    "spot_a_rad": 3.0,
    "spot_a_tempf": 4.0,

    "spot_b_lat": 5.0,
    "spot_b_long": 6.0,
    "spot_b_rad": 7.0,
    "spot_b_tempf": 8.0,

    "a": 9.0,
    "e": 10.0,
    "perr": 11.0,
    "f1": 12.0,
    "f2": 13.0,
    "pshift": 14.0,
    "vga": 15.0,
    "xincl": 16.0,
    "g1": 17.0,
    "g2": 18.0,
    "tavc": 19.0,
    "tavh": 20.0,
    "alb1": 21.0,
    "alb2": 22.0,
    "phsv": 23.0,
    "pcsv": 24.0,
    "rm": 25.0,
    "hjd0": 26.0,
    "pzero": 27.0,
    "dpdt": 28.0,
    "dperdt": 29.0,
    "a3b": 30.0,
    "p3b": 31.0,
    "xincl3b": 32.0,
    "e3b": 33.0,
    "perr3b": 34.,
    "t03b": 35.0,
    "dpclog": 41.0,
    "desextinc": 42.0,

    "spot_a_tstart": 43.0,
    "spot_a_tmax1": 44.0,
    "spot_a_tmax2": 45.0,
    "spot_a_tend": 46.0,

    "spot_b_tstart": 47.0,
    "spot_b_tmax1": 48.0,
    "spot_b_tmax2": 49.0,
    "spot_b_tend": 50.0,

    "hla": 56.0,
    "cla": 57.0,
    "x1a": 58.0,
    "x2a": 59.0,
    "el3a": 60.0
}

# reverse of the dictionary above
DC_KEEPS_ID_NAME_DICT = {value: key for key, value in list(DC_KEEPS_NAME_ID_DICT.items())}
