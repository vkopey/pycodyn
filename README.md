# pycodyn
Component-oriented acausal modeling of the dynamical systems in Python language on the example of the model of the sucker rod string.
For details click the link below.

If you use pycodyn please cite the following reference in your work (books, articles, reports, etc.):

Kopei VB, Onysko OR, Panchuk VG. 2019. Component-oriented acausal modeling of the dynamical systems in Python language on the example of the model of the sucker rod string. PeerJ Preprints 7:e27612v1 https://doi.org/10.7287/peerj.preprints.27612v1

## Package content:
Pycodyn.mo - models of sucker rod string (Modelica language)  
pycodyn.py - components and solver (Euler method)  
pycodynDAE.py - components and solver (DAE)  
main1.py - model of free vibrations of sucker rod string (Euler method)  
trapComponents.py - components (trapezoidal rule)  
main1T.py - model of free vibrations of sucker rod string (trapezoidal rule)  
main1DAE.py - model of free vibrations of sucker rod string (DAE)  
main1Sym.py - model of free vibrations of sucker rod string (analytical)  
main2s.py - single-section model of pumping process (Euler method)  
main2.py - two-section model of pumping process (Euler method)  
main2V.py - two-section model of string breakage (Euler method, events)  
main2sDAE.py - single-section model of pumping process (DAE)  
main2DAE.py - two-section model of pumping process (DAE)  
eliminate.py - equations elimination for SymPy

## Requirements:
Python 3.7 (recommended) or Python 2.7  
numpy-1.16.4+mkl  
scipy-1.2.2 (optionally)  
mpmath-1.1.0  
sympy-1.4  
Assimulo-2.9 (for pycodynDAE)  
matplotlib-2.2.4 (optionally for plotting)  
OpenModelica 1.12 (for Pycodyn.mo)

## Installation:
Just put all the modules in one folder

## Portable installation on Windows:
Download and install WinPython2.7Zero or WinPython3.7Zero from https://sourceforge.net/projects/winpython/files/  
Download requirements from https://www.lfd.uci.edu/~gohlke/pythonlibs/ and put them to WinPython/scripts/  
Run WinPython Command Prompt.exe as admin and install requirements (pip update may be required). Example for WinPython 32bit:

    pip install numpy-1.16.4+mkl-cp37-cp37m-win32.whl
    pip install mpmath-1.1.0-py2.py3-none-any.whl
    pip install sympy-1.4-py2.py3-none-any.whl
    ...

## Usage example (Python):
    d:\WinPython32_37\python-3.7.4\python.exe main1.py

## Usage example (OpenModelica):
    loadModel(Modelica)
    loadFile("Pycodyn.mo")
    simulate(Pycodyn.Oscillator, stopTime=10)
    plot(mass1.s)
    simulate(Pycodyn.Pumping, stopTime=100)
    plotParametric(motion1.flange.s, spring1.flange_a.f)
