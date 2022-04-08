#script: diagrams.py
#purpose: creates an HR diagram, L-M diagram, and R-M diagram.

import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, rc, ticker as tk
rc('text', usetex=True)

from standards import Standards as s


def HR_diag(data, savepath=""):
    """Creates an H-R diagram with the generated stars.
    INPUTS:
        data - A DataFrame [loaded in from something like star_summary.txt] with a L and Tsurf column.
        savepath - The path to save the data to. If none is given then it will display.
    OUTPUTS:
        Plot of a H-R diagram.
    """
    Lscaled = data.L/s.Lsun

    fig, ax = plt.subplots(figsize=(7,10))

    ax.set_xlim(800, 5e4)
    ax.set_ylim(0.005, 200)
    ax.invert_xaxis()
    ax.set_xscale('log')
    ax.yaxis.set_major_locator(tk.MultipleLocator(50))
    ax.yaxis.set_minor_locator(tk.AutoMinorLocator())
    ax.tick_params(direction='in', which='both')
    ax.tick_params(which='major', length=5)
    ax.tick_params(which='minor', length=3)


    ax_t = ax.secondary_xaxis('top')
    ax.set_yscale('log')
    ax_t.tick_params(direction='in', which='both')
    ax_t.tick_params(which='major', length=5)
    ax_t.tick_params(which='minor', length=3)

    ax.scatter(data.Tsurf, Lscaled, color='black', s=5)

    #from p. 342 in Ryden
    spectral_type = ['O', 'B', 'A', 'F', 'G', 'K', 'M', 'L', 'T']
    spectral_temp = [40000, 20000, 9000, 7000, 5500, 4500, 3000, 2000, 1300]

    ax.set_xticks(spectral_temp, spectral_type)

    ax_t.set_xlabel(r'$T$ (K)')
    ax_t.set_xticks([1000, 5000, 1e4], ['1000', '5000', r'10$^4$'])

    ax.set_ylabel(r'$L$/$L_\odot$', rotation=0)
    ax.set_yticks([0.01, 0.1, 1, 10, 100], ['0.01', '0.1', '1', '10', '100'])


    if savepath != "":
        plt.savefig(savepath)
    else:
        plt.show()

    plt.close('all')

def LM_diag(data, savepath=""):
    """Creates a L-M diagram with the generated stars.
    INPUTS:
        data - A DataFrame [loaded in from something like star_summary.txt] with a L and M column.
        savepath - The path to save the data to. If none is given then it will display.
    OUTPUTS:
        Plot of L/Lsun as a function of M/Msun.
    """

    Lscaled = data.L/s.Lsun
    Mscaled = data.M/s.Msun

    fig, ax = plt.subplots(figsize=(7,10))
    ax.set_facecolor('black')

    ax.scatter(Mscaled, Lscaled, color='gold', s=5)

    ax.set_xlabel(r'$M$/$M_\odot$')
    ax.set_ylabel(r'$L$/$L_\odot$')

    if savepath != "":
        plt.savefig(savepath)
    else:
        plt.show()

    plt.close('all')

def RM_diag(data, savepath=""):
    """Creates a R-M diagram with the generated stars.
    INPUTS:
        data - A DataFrame [loaded in from something like star_summary.txt] with a R and M column.
        savepath - The path to save the data to. If none is given then it will display.
    OUTPUTS:
        Plot of R/Rsun as a function of M/Msun.
    """

    Rscaled = data.R/s.Rsun
    Mscaled = data.M/s.Msun

    fig, ax = plt.subplots(figsize=(7,10))
    ax.set_facecolor('black')

    ax.scatter(Mscaled, Rscaled, color='gold', s=5)

    ax.set_xlabel(r'$M$/$M_\odot$')
    ax.set_ylabel(r'$R$/$R_\odot$')

    if savepath != "":
        plt.savefig(savepath)
    else:
        plt.show()

    plt.close('all')



if __name__ == '__main__':

    #test dataset

    M = 5*np.random.rand(100)*s.Msun
    L = 400*np.random.rand(100)*s.Lsun
    R = 8*np.random.rand(100)*s.Rsun
    T = np.power(L/(4*np.pi*s.sb*np.power(R, 2)), 1/4)

    data = pd.DataFrame(data=zip(M, L, R, T), columns=['M', 'L', 'R', 'Tsurf'])

    savepath = os.path.join(s.plots_folder, 'LM_diag_test.png')
    savepath2 = os.path.join(s.plots_folder, 'RM_diag_test.png')
    savepath3 = os.path.join(s.plots_folder, 'HR_diag_test.png')

    LM_diag(data, savepath)
    RM_diag(data, savepath2)
    HR_diag(data, savepath3)
