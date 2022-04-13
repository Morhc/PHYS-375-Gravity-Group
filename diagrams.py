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

    for xx in data:
        fff, lam = xx

        Lscaled = fff.L/s.Lsun * 1e-7

        ax.plot(fff.Tsurf, Lscaled, label=rf'$\lambda$ = {lam}')

    #from p. 342 in Ryden
    spectral_type = ['O', 'B', 'A', 'F', 'G', 'K', 'M', 'L', 'T']
    spectral_temp = [40000, 20000, 9000, 7000, 5500, 4500, 3000, 2000, 1300]

    ax.set_xticks(spectral_temp, spectral_type)

    ax_t.set_xlabel(r'$T$ (K)')
    ax_t.set_xticks([1000, 5000, 1e4], ['1000', '5000', r'10$^4$'])

    ax.set_ylabel(r'$L$/$L_\odot$', rotation=0)
    ax.set_yticks([0.01, 0.1, 1, 10, 100], ['0.01', '0.1', '1', '10', '100'])

    lgd = plt.legend()

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


    fig, ax = plt.subplots(figsize=(7,10))

    for xx in data:
        fff, lam = xx

        Lscaled = fff.L/s.Lsun * 1e-7  * 1e-7
        Mscaled = fff.M/s.Msun / 1000

        ax.plot(Mscaled, Lscaled, label=rf'$\lambda$ = {lam}')

    #from Ryden p. 330
    Mtest = np.linspace(0.1, 10, 100)
    Ltest = []
    for M in Mtest:
        if M < 0.7:
            Ltest.append(0.35*(M**2.62))
        else:
            Ltest.append(1.02*(M**3.92))

    plt.plot(Mtest, Ltest, label='Expectation', ls='--')

    ax.set_xlabel(r'$M$/$M_\odot$')
    ax.set_ylabel(r'$L$/$L_\odot$')

    plt.legend()

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

    fig, ax = plt.subplots(figsize=(7,10))

    for xx in data:
        fff, lam = xx
        Rscaled = fff.R/s.Rsun / 100
        Mscaled = fff.M/s.Msun / 1000
        ax.plot(Mscaled, Rscaled, label=rf'$\lambda$ = {lam}')

    #from Ryden p. 330
    Mtest = np.linspace(0.1, 10, 100)
    Rtest = []
    for M in Mtest:
        if M < 1.66:
            Rtest.append(1.06*(M**0.945))
        else:
            Rtest.append(1.33*(M**0.555))

    plt.plot(Mtest, Rtest, label='Expectation', ls='--')

    ax.set_xlabel(r'$M$/$M_\odot$')
    ax.set_ylabel(r'$R$/$R_\odot$')

    plt.legend()


    if savepath != "":
        plt.savefig(savepath)
    else:
        plt.show()

    plt.close('all')



if __name__ == '__main__':

    both = []
    for lam in [-100, -10, 0, 10, 100, 500]:
        data = pd.read_csv(os.path.join(s.data_folder, f'stars_lam_{lam}', 'generated_stars.csv'), header=0)[:50]
        both.append((data, lam))

        #hr_path = os.path.join(s.plots_folder, f'HR_diag_lam_{lam}')
        #lm_path = os.path.join(s.plots_folder, f'LM_diag_lam_{lam}')
        #rm_path = os.path.join(s.plots_folder, f'RM_diag_lam_{lam}')

        #HR_diag(data, hr_path)
        #LM_diag(data, lm_path)
        #RM_diag(data, rm_path)
    data = pd.read_csv(os.path.join(s.data_folder, 'test', 'generated_stars.csv'), header=0)[:50]

    both = [(data,0)]
    hr_path = os.path.join(s.plots_folder, 'HR_diag_comp')
    lm_path = os.path.join(s.plots_folder, 'LM_diag_comp')
    rm_path = os.path.join(s.plots_folder, 'RM_diag_comp')

    HR_diag(both, hr_path)
    LM_diag(both, lm_path)
    RM_diag(both, rm_path)
