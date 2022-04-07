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
    a = 4*sb/c #Radiation constant (see dT/dr differential)

    # Project specific constants
    gamma = -2/3 # not sure about this value
    X = 0.73 # need to adjust these values. Ensure that X+Y+Z = 1
    Y = 0.025
    Z = 0.02
    mu = 1/(2*X + 0.75*Y + 0.5*Z)

    #scaling constants
    g_to_kg = 1000
    cm_to_m = 100
    cmgs_to_mkgs = g_to_kg/(cm_to_m**3)
    ergs_to_W = 1.0e-7

    #path of standards.py
    here = os.path.dirname(os.path.realpath(__file__))

    #path to plots folder
    plots_folder = os.path.join(here, 'Plots')

    #path to data folder
    data_folder = os.path.join(here, 'Data')

    #standard palette for plot colours [colour-blind friendly]
    colours = ["black",
               "limegreen",
               "red",
               "blue",
               "purple"]
