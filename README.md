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

Use CPython **3.12** in a dedicated virtual environment. Python 2.7 and older Python 3 releases are no longer supported by this source tree. PyQt5 is retained; a Qt6 migration is not required.

```sh
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux: source .venv/bin/activate
python -m pip install --require-hashes -r requirements.lock
python pywd2015.py
```

The lock includes the full runtime dependency set, with platform-specific Qt packages. Development checks use `requirements-dev.lock` and `python -m pytest`. The Windows source path is validated with CPython 3.12.11 and offscreen Qt; Windows/Linux CI runs on each push and pull request; Windows additionally executes the hash-verified reference WD engines. Interactive platform acceptance and native Qt security are separate from Python package advisory results.

Each calculation uses a private temporary job directory and copies the WD `.dat`/`.tab` auxiliary files from the executable directory. Keep those data files beside the selected LC/DC binaries. The executable directory is never used for job output. The last LC/DC input and output buttons open the current session's private job files, which are temporary. Application settings are written to the user's config directory; legacy adjacent config is still readable.

LC calculations show a cancellable modal progress window; the Qt event loop remains responsive. Its time limit defaults to 60 minutes and can be changed while running. DC has the same configurable limit in its calculation controls. Failed, timed-out and cancelled jobs do not enter result/history publication.

Project loading validates all saved fields and observation references before applying state. Project files keep the existing INI format and literal `%` text. Calendar conversion uses the historical Julian/Gregorian transition: Julian through 1582-10-04, Gregorian from 1582-10-15; the intervening dates are rejected. Supported calendar years are 1–9999. No local timezone, JD/HJD or time-scale conversion is implied.

See [audit validation](docs/audit-validation.md) for the per-finding evidence, native engine provenance and remaining limits. New releases must include a clean source commit, `tools/release_manifest.py` output and an SBOM; no new application executable is produced by these source corrections.

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
