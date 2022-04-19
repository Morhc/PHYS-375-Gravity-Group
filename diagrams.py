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
    fig, ax = plt.subplots(figsize=(5,5))

    #ax.set_xlim(800, 5e4)
    #ax.set_ylim(0.005, 500)
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

        #filter = fff.R < 3.47e12 #means that the code broke

        Lscaled = fff.L/s.Lsun * 1e-7

        #ax.scatter(fff[filter].Tsurf, Lscaled, label=rf'$\lambda$ = {lam}', edgecolor='black', lw=1)
        f1 = fff.Tsurf > 1e4
        f2 = Lscaled > 0.6
        f4 = fff.Tsurf < 1e4
        f3 = Lscaled < 0.6
        ax.scatter(fff[f1].Tsurf, Lscaled[f1], label='Region 1', edgecolor='black', lw=1)
        ax.scatter(fff[f2*f4].Tsurf, Lscaled[f2*f4], label='Region 2', edgecolor='black', lw=1)
        ax.scatter(fff[f3].Tsurf, Lscaled[f3], label='Region 3', edgecolor='black', lw=1)

    #from p. 342 in Ryden
    spectral_type = ['O', 'B', 'A', 'F', 'G', 'K', 'M', 'L', 'T']
    spectral_temp = [40000, 20000, 9000, 7000, 5500, 4500, 3000, 2000, 1300]

    ax.set_xticks(spectral_temp, spectral_type)

    ax_t.set_xlabel(r'$T$ (K)')
    ax_t.set_xticks([1000, 5000, 1e4], ['1000', '5000', r'10$^4$'])

    ax.set_ylabel(r'$L$/$L_\odot$', rotation=0)
    ax.set_yticks([0.01, 0.1, 1, 10, 100], ['0.01', '0.1', '1', '10', '100'])

    #ax.set_xlim(50000, 1000)
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

    fig, ax = plt.subplots(figsize=(5,5))

    for xx in data:
        fff, lam = xx

        filter = fff.R < 3.47e12 #means that the code broke

        Lscaled = fff[filter].L/s.Lsun * 1e-7
        Mscaled = fff[filter].M/s.Msun / 1000

        #ax.scatter(np.log10(Mscaled), np.log10(Lscaled), label=rf'$\lambda$ = {lam}', edgecolor='black', lw=1)
        ax.scatter(np.log10(Mscaled), np.log10(Lscaled), label='Generated', edgecolor='black', lw=1)


    #from Ryden p. 330
    Mtest = np.linspace(0.1, 10, 100)
    Ltest = []
    for M in Mtest:
        if M < 0.7:
            Ltest.append(0.35*(M**2.62))
        else:
            Ltest.append(1.02*(M**3.92))

    plt.plot(np.log10(Mtest), np.log10(Ltest), label='Expectation', ls='--', color='black')

    ax.set_xlabel(r'log$_{10}$($M$/$M_\odot$)')
    ax.set_ylabel(r'log$_{10}$($L$/$L_\odot$)')

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

    fig, ax = plt.subplots(figsize=(5,5))

    for xx in data:
        fff, lam = xx

        filter = fff.R < 3.47e12 #means that the code broke

        Rscaled = fff[filter].R/s.Rsun / 100
        Mscaled = fff[filter].M/s.Msun / 1000

        #ax.scatter(np.log10(Mscaled), np.log10(Rscaled), label=rf'$\lambda$ = {lam}', edgecolor='black', lw=1)
        ax.scatter(np.log10(Mscaled), np.log10(Rscaled), label='Generated', edgecolor='black', lw=1)

    #from Ryden p. 330
    Mtest = np.linspace(0.1, 10, 100)
    Rtest = []
    for M in Mtest:
        if M < 1.66:
            Rtest.append(1.06*(M**0.945))
        else:
            Rtest.append(1.33*(M**0.555))

    plt.plot(np.log10(Mtest), np.log10(Rtest), label='Expectation', ls='--', color='black')

    ax.set_xlabel(r'log$_{10}$($M$/$M_\odot$)')
    ax.set_ylabel(r'log$_{10}$($R$/$R_\odot$)')

    plt.legend()


    if savepath != "":
        plt.savefig(savepath)
    else:
        plt.show()

    plt.close('all')



if __name__ == '__main__':

    both = []
    for lam in ['-1e7', '-1e3', '100', '0', '100', '1e3', '1e7']:
        data = pd.read_csv(os.path.join(s.data_folder, f'stars_lam_{lam}', 'generated_stars.csv'), header=0)[:50]
        both.append((data, lam))

        #hr_path = os.path.join(s.plots_folder, f'HR_diag_lam_{lam}')
        #lm_path = os.path.join(s.plots_folder, f'LM_diag_lam_{lam}')
        #rm_path = os.path.join(s.plots_folder, f'RM_diag_lam_{lam}')

        #HR_diag(data, hr_path)
        #LM_diag(data, lm_path)
        #RM_diag(data, rm_path)
    #data = pd.read_csv(os.path.join(s.data_folder, 'stars_lam_0', 'generated_stars.csv'), header=0)[:50]

    data = pd.read_csv(os.path.join(s.data_folder, 'stars_lam_0', 'generated_stars.csv'), header=0)

    #data = data[round(data.R/100, ) < 30*s.Rsun]

    #data.Tsurf = np.power(data.L * 1e-7 / (4*np.pi*s.sb*(data.R / 100)**2), 1/4)

    both = [(data,0)]
    hr_path = os.path.join(s.plots_folder, 'HR_diag_lam=0')
    lm_path = os.path.join(s.plots_folder, 'LM_diag_lam=0')
    rm_path = os.path.join(s.plots_folder, 'RM_diag_lam=0')

    HR_diag(both, hr_path)
    LM_diag(both, lm_path)
    RM_diag(both, rm_path)
