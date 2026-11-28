# PyWD2015-Qt5
A Qt5 GUI written in Python for Wilson - Devinney eclipsing binary modeling software.

UPDATE: Version 1.0.0 now has some small but useful updates. Download and try the source code or the precompiled Windows executable file available [here](https://github.com/Varnani/pywd2015-qt5/releases/download/v1.0.0/pywd2015_win64bit_1November2025.exe).

IMPORTANT NOTE: TESS bandpass was added to PyWD2015. However, users must get updated atmosphere and limb darkening files from https://faculty.fiu.edu/~vanhamme/wdfiles/lcdc-standard-data.tar.gz in order to run PyWD2015 properly.

PyWD2015 v1.0.0 is released (16th August 2022). Please visit [releases](https://github.com/Varnani/pywd2015-qt5/releases) page and [download](https://github.com/Varnani/pywd2015-qt5/releases/download/v1.0.0/wd-precompiled.zip) all required files, including executable binaries, auxiliary files and source codes. 

A comprehensive user manual is also available in [this link](https://github.com/Varnani/pywd2015-qt5/releases/download/v1.0.0/PyWD2015_qt5_manual.pdf).

Single executable file is [available](https://github.com/Varnani/pywd2015-qt5/releases/download/v1.0.0/pywd2015_win64bit_1November2025.exe) for Windows 7, 8.x, 10 and 11 (all 64 bit) operating systems.

## First Things First
If you happen to use this software in a scientific work, please cite PyWD2015 release proceeding:  
[CoSka (2020), 50, 535](https://ui.adsabs.harvard.edu/abs/2020CoSka..50..535G/abstract)  

You must also properly cite relevant Wilson - Devinney papers:  
[ApJ (1971), 166, 605](https://ui.adsabs.harvard.edu/abs/1971ApJ...166..605W/abstract)  
[ApJ (1979), 234, 1054](https://ui.adsabs.harvard.edu/abs/1979ApJ...234.1054W/abstract)  
[ApJ (1990), 356, 613](https://ui.adsabs.harvard.edu/abs/1990ApJ...356..613W/abstract)  
[ApJ (2007), 661, 1129](https://ui.adsabs.harvard.edu/abs/2007ApJ...661.1129V/abstract)  
[ApJ (2008), 672, 575](https://ui.adsabs.harvard.edu/abs/2008ApJ...672..575W/abstract)  
[AJ (2012), 144, 73](https://ui.adsabs.harvard.edu/abs/2012AJ....144...73W/abstract)  
[ApJ (2014), 780, 151](https://ui.adsabs.harvard.edu/abs/2014ApJ...780..151W/abstract)

## Installing Dependencies

Use CPython 3.12 and a virtual environment:

```sh
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux: source .venv/bin/activate
python -m pip install --require-hashes -r requirements.lock
python pywd2015.py
```

For development, install `requirements-dev.lock` and run `python -m pytest -q`. Tests run offscreen on Windows and Linux. To run the native LC/DC tests on Windows, set `PYWD_ENGINE_DIR` to the extracted `wd-precompiled` directory. Engine versions and fixture hashes are recorded in `tests/fixtures/wd2015-win64/provenance.json`.

Keep the WD `.dat`/`.tab` files beside the selected LC/DC binaries. Each calculation copies them to its own temporary directory. The input/output viewer opens files from the latest job; export any results you want to keep before exiting. Application settings use the user's config directory, with a fallback for existing settings beside the application.

LC and DC calculations can be cancelled and have an adjustable 60-minute time limit. Failed or cancelled calculations do not add results to the history.

Project loading validates all saved fields and observation references before applying state. Project files keep the existing INI format and literal `%` text. Calendar conversion uses the historical Julian/Gregorian transition: Julian through 1582-10-04, Gregorian from 1582-10-15; the intervening dates are rejected. Supported calendar years are 1–9999. No local timezone, JD/HJD or time-scale conversion is implied.

## Getting Started

After installing required dependencies, you can get the source code of the program by either cloning this repository, or downloading a release from the [Releases](https://github.com/Varnani/pywd2015-qt5/releases) page. Extract the source code and issue this command to start the PyWD2015:

```shell
python pywd2015.py
```

To run calculations, you need to provide the paths of compiled DC and LC programs of the Wilson - Devinney code. Precompiled versions are hosted here with the permission of Dr. Robert E. Wilson. These too can be found on the "Releases" page. You can also compile the WD binaries yourself. Source files can be found at https://github.com/rewilsonufl/ebmodel/ (the WD github Homepage).

GUI is mostly self-explanatory, and most labels contain tooltips, but "Releases" page also contains a user manual. While this manual is based on a deprecated, Qt4 version of PyWD2015, nearly all the concepts discussed there also applies to the Qt5 version.

## Contact

If you have any questions, you can reach the authors from:

```
Ozan Güzel:          ozanguzel35@outlook.com  
Dr. Orkun Özdarcan:  orkun.ozdarcan@ege.edu.tr
``` 
