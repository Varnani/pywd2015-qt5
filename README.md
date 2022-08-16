# PyWD2015-Qt5
A Qt5 GUI written in Python for Wilson - Devinney eclipsing binary modeling software.

IMPORTANT NOTE: TESS bandpass was added to PyWD2015. However, users must get updated atmosphere and limb darkening files from https://faculty.fiu.edu/~vanhamme/wdfiles/lcdc-datafiles-2019.tar.gz in order to run PyWD2015 properly.

PyWD2015 v1.0.0 is released (16th August 2022). Please visit [releases](https://github.com/Varnani/pywd2015-qt5/releases) page and [download](https://github.com/Varnani/pywd2015-qt5/releases/download/v1.0.0/wd-precompiled.zip) all required files, including executable binaries, auxiliary files and source codes. 

A comprehensive user manual is also available in [this link](https://github.com/Varnani/pywd2015-qt5/releases/download/v1.0.0/PyWD2015_qt5_manual.pdf).

Single executable file is [available](https://github.com/Varnani/pywd2015-qt5/releases/download/v1.0.0/pywd2015_win64bit_16August2022.exe) for Windows 7, 8.x and 10 (all 64 bit) operating systems.

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

> **WARNING** | Some Linux distributions use Python 3.x as their default Python interpreter. In that case 'python' command and package name prefix will actually refer to Python 3.x and 'python2' will refer to Python 2.7. Please refer to your distributions package manager database if you are unsure. On Windows, multiple Python installations may conflict with each other, and your 'python' command might run a different version rather than what you intented. Run 'python --version' from command line to check your Python version. 

PyWD2015 supports Python 2.7 and 3.5+. It also relies on Numpy, Scipy, Matplotlib libraries. These libraries can be installed using pip, the Python package installer. Pip should be automatically installed alongside Python on Windows. On Linux, you may need to install pip with your package manager. On Ubuntu and Debian, you can issue:

```shell
(For Python 2.7) sudo apt install python-pip
(For Python 3.x) sudo apt install python3-pip
```   

to install pip. You can then install the dependencies with the following command:  

```shell
Linux: 
(For Python 2.7) sudo python -m pip install numpy scipy matplotlib
(For Python 3.x) sudo python3 -m pip install numpy scipy matplotlib

Windows:
python -m pip install numpy scipy matplotlib
```

On Debian, you may encounter a "backports.functools_lru_cache" and/or "tkinter" error on a fresh matplotlib installation under Python 2.7. To fix this, you can issue:  

```shell
sudo apt install python-backports.functools-lru-cache python-tk
```  

After that, you need to install the PyQt5 library.

### Linux
Installing PyQt5 on Linux depends on your distribution and package manager. In general, you should not use pip to install PyQt5, instead use your distributions package manager. On Ubuntu and Debian you can use:  

```shell
(For Python 2.7) sudo apt install python-pyqt5
(For Python 3.x) sudo apt install python3-pyqt5
```  

### Windows
You can safely use pip to install PyQt5 on Windows 7 and Windows 8.x operating systems:  

```shell
python -m pip install python-qt5
```  

However, for the latest versions of Python (3.7+) Windows 10 users may have to modify 
the command as shown below:

```shell
python -m pip install PyQT5
```  

In some Windows operating systems with older Python versions (especially 3.5 or 3.6), users may experience "dll import error" when trying to run PyWD2015. This issue is usually related to PyQT5 and may not have a unique reason. In many case, finding a solution might be time-consuming. In order to avoid this issue, uninstalling the older version and then installing the latest available Python version is strongly recommended. After that, users should install required libraries as described above, as well as PyQT5. 

PyWD2015 has not been tested on Windows 11 yet.

## Getting Started

After installing required dependencies, you can get the source code of the program by either cloning this repository, or downloading a release from the [Releases](https://github.com/Varnani/pywd2015-qt5/releases) page. Extract the source code and issue this command to start the PyWD2015:

```shell
python pywd2015.py
```

To run calculations, you need to provide the paths of compiled DC and LC programs of the Wilson - Devinney code. Precompiled versions are hosted here with the permission of Dr. Robert E. Wilson. These too can be found on the "Releases" page. You can also compile the WD binaries yourself. Source files can be found at ftp://ftp.astro.ufl.edu/pub/wilson/lcdc2015/ (the WD Homepage).

GUI is mostly self-explanatory, and most labels contain tooltips, but "Releases" page also contains a user manual. While this manual is based on a deprecated, Qt4 version of PyWD2015, nearly all the concepts discussed there also applies to the Qt5 version.

## Contact

If you have any questions, you can reach the authors from:

```
Ozan Güzel:          ozanguzel35@outlook.com  
Dr. Orkun Özdarcan:  orkun.ozdarcan@ege.edu.tr
``` 
