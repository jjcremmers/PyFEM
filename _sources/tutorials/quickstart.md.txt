# Quickstart

To test whether everything is installed properly, the following two simulations can be run.

## Simple Example
In the directory `examples/ch02` the script `PatchTest.py` can be executed from a terminal (or DOS-shell) by typing:
```bash
python PatchTest.py
```
In Windows, this script can also be executed by double-clicking the icon.

## PyFEM Example
The full finite element code PyFEM can be run by typing `pyfem` in the terminal. In directory `examples/ch04`, for example, the input file `ShallowTrussRiks.pro` is processed by typing:
```bash
pyfem ShallowTrussRiks.pro
```
Here, `ShallowTrussRiks.pro` is the input file, which by definition ends with `.pro`. When opened in a text editor, it looks as follows:
```text
input = "ShallowTrussRiks.dat";
TrussElem  = {
  ...
}
SpringElem = {
  ...
}
solver = {
  ...
}
```
(See documentation for further details)
