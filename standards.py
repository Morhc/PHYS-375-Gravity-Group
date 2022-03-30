#script: standards.py
#purpose: To provide a standard across the project.

import os

class Standards:
    """This is a class that is used to define constants and standardizable constants. Unless
    otherwise stated, the value of all constants is in SI units and taken from Ryden.
    """

    #Solar standard
    Msun = 1.989e30
    Rsun = 6.955e8
    Lsun = 3.839e26

    #Physical constants
    G = 6.673e-11
    e = 1.602e-19 #charge of the electron
    c = 2.998e8
    h = 6.626e-34
    hbar = 1.055e-34
    k = 1.381e-23 #Boltzmann constant
    sb = 5.670e-8 #Stefan-Boltzmann constant
    mp = 1.673-27
    me = 9.109e-31

    #path of standards.py
    here = os.path.dirname(os.path.realpath(__file__))

    #path to plots folder
    plots_folder = os.path.join(here, 'Plots')

    #path to data folder
    data_folder = os.path.join(here, 'Data')

    #standard palette for plot colours [colour-blind friendly]
    colours = ["#d55e00",
               "#cc79a7",
               "#0072b2",
               "#f0e442"
               "#009e73"]
