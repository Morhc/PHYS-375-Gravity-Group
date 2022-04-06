#script: star_info.py
#purpose: Plots the star info.

import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, rc, ticker as tk
rc('text', usetex=True)

from tools import load_data
from standards import Standards as s

def plot_star(data, info = None, savepath=""):
    """Plots the various star information.
    INPUTS:
        data - A DataFrame with the star data
        info - A summary of the star info in a DataFrame.
        savepath - The path to save the plot to. If no path is given, the plot will be displayed.
    """

    def as_si(x, ndp):
        """Returns scientific notation text of a number.
        from: https://stackoverflow.com/questions/31453422/displaying-numbers-with-x-instead-of-e-scientific-notation-in-matplotlib/31453961
        INPUTS:
            x - the number to make into scientific notation
            ndp - the number of decimal points
        """
        s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
        m, e = s.split('e')
        return r'{m:s}\times 10^{{{e:d}}}'.format(m=m, e=int(e))

    fig, axes = plt.subplots(3,2, figsize=(9,9))
    plt.subplots_adjust(hspace=0, wspace=0.4)

    #if we are given star info
    axes[0,1].remove()
    if not isinstance(None, type(info)):

        tsurf = rf'$T_*$ = {round(info.Tsurf,)} K'
        plt.figtext(0.58,0.83, tsurf, fontsize=20)

        star_info = rf'$\rho_c$ = {round(info.rhoc, 2)} g/cm$^3$' + '\n\n' + \
                    r'$T_c$ = ${0:s}$ K'.format(as_si(info.Tc, 2)) + '\n\n' + \
                    rf'$R_*$ = {round(info.R/s.cm_to_m/s.Rsun, 3)} $R_\odot$' + '\n\n' + \
                    rf'$M_*$ = {round(info.M/s.g_to_kg/s.Msun, 3)} $M_\odot$' + '\n\n' + \
                    r'$L_*$ = ${0:s}$ $L_\odot$'.format(as_si(info.L*s.ergs_to_W/s.Lsun, 2))

        plt.figtext(0.58,0.65, star_info, fontsize=12)

    ## FIRST PLOT - density, mass, luminosity, temperature

    pl_rho = data['rho(r)']/info.rhoc
    pl_M = data['M(r)']/info.M
    pl_T = data['T(r)']/info.Tc
    pl_L = data['L(r)']/info.L

    #need to relativize the coordinates of the text
    axes[0,0].plot(data['r/R'], pl_rho, color = s.colours[0])
    axes[0,0].text(0.03, 0.8, r'$\rho$', color = s.colours[0], fontsize=14)
    axes[0,0].plot(data['r/R'], pl_M, color = s.colours[1], ls='dashed')
    axes[0,0].text(0.6, 0.9, r'$M$', color = s.colours[1], fontsize=14)
    axes[0,0].plot(data['r/R'], pl_T, color = s.colours[2], ls=(0, (5, 10)))
    axes[0,0].text(0.7, 0.15, r'$T$', color = s.colours[2], fontsize=14)
    axes[0,0].plot(data['r/R'], pl_L, color = s.colours[3], ls='dotted')
    axes[0,0].text(0.02, 0.58, r'$L$', color = s.colours[3], fontsize=14)

    axes[0,0].set_ylabel(r'$\rho$/$\rho_c$, $T$/$T_c$, $M$/$M_*$, $L$/$L_*$')

    ## SECOND PLOT - pressures

    Pc = data['P'][0]

    #need to relativize the coordinates of the text
    axes[1,0].plot(data['r/R'], data['P']/Pc, color = s.colours[0])
    axes[1,0].text(0.4, 0.1, r'$P_{tot}$', color = s.colours[0], fontsize=14)
    axes[1,0].plot(data['r/R'], data['Pgas']/Pc, color = s.colours[1], ls='dashed')
    axes[1,0].text(0.15, 0.45, r'$P_{gas}$', color = s.colours[1], fontsize=14)
    axes[1,0].plot(data['r/R'], data['Pdeg']/Pc, color = s.colours[2], ls=(0, (5, 10)))
    axes[1,0].text(0.05, 0.1, r'$P_{deg}$', color = s.colours[2], fontsize=14)
    axes[1,0].plot(data['r/R'], data['Ppho']/Pc, color = s.colours[3], ls='dotted')
    axes[1,0].text(0.2, 0.05, r'$P_{\gamma}$', color = s.colours[3], fontsize=14)

    axes[1,0].set_ylabel(r'$P$/$P_c$')


    ## THIRD PLOT - dlogP/dlogT

    axes[2,0].plot(data['r/R'], data['dlogP/dlogT'], color = s.colours[0])

    axes[2,0].set_ylim(0, 2*np.nanmax(data['dlogP/dlogT']))
    axes[2,0].set_ylabel(r'dlog$P$/dlog$T$')
    axes[2,0].set_xlabel(r'$r$/$R_*$')

    ## FOURTH PLOT - log(kappa)

    #need to relativize the coordinates of the text
    axes[1,1].plot(data['r/R'], np.log10(data['kappa']), color = s.colours[0])
    axes[1,1].text(0.5, 1, r'$\kappa$', color = s.colours[0], fontsize=14)
    axes[1,1].plot(data['r/R'], np.log10(data['kappaff']), color = s.colours[1], ls='dashed')
    axes[1,1].text(0.8, 2, r'$\kappa_{ff}$', color = s.colours[1], fontsize=14)
    axes[1,1].plot(data['r/R'], np.log10(data['kappaHm']), color = s.colours[2], ls=(0, (5, 10)))
    axes[1,1].text(0.85, 9, r'$\kappa_{H^-}$', color = s.colours[2], fontsize=14)
    axes[1,1].plot(data['r/R'], np.log10(data['kappaes']), color = s.colours[3], ls='dotted')
    axes[1,1].text(0.3, -1.5, r'$\kappa_{es}$', color = s.colours[3], fontsize=14)

    axes[1,1].set_ylabel(r'log$_{10}$($\kappa$) (cm$^2$/g)')
    axes[1,1].set_ylim(-2, 10)

    #need to relativize the coordinates of the text
    axes[2,1].plot(data['r/R'], data['dL/dr']*info.R/info.L, color = s.colours[0])
    axes[2,1].text(0.6, 10, r'd$L$/d$r$', color = s.colours[0], fontsize=14)
    axes[2,1].plot(data['r/R'], data['dLpp/dr']*info.R/info.L, color = s.colours[2], ls=(0, (5, 10)))
    axes[2,1].text(0.6, 9, r'd$L_{pp}$/d$r$', color = s.colours[2], fontsize=14)
    axes[2,1].plot(data['r/R'], data['dLcno/dr']*info.R/info.L, color = s.colours[3], ls='dotted')
    axes[2,1].text(0.6, 8, r'd$L_{CNO}$/d$r$', color = s.colours[3], fontsize=14)

    axes[2,1].set_ylabel(r'd$L_*$/d$r$ (L$_*$/R$_*$)')
    axes[2,1].set_ylim(0, 2*np.nanmax(data['dL/dr'])*info.R/info.L)
    axes[2,1].set_xlabel(r'$r$/$R_*$')

    #settings aesthetics for all plots
    for i, ax in enumerate(axes.flatten()):
        ax.set_xlim(0,1)
        ax.xaxis.set_major_locator(tk.MultipleLocator(0.2))
        ax.xaxis.set_minor_locator(tk.MultipleLocator(0.1))
        #ax.tick_params(which='major', length=7)
        #ax.tick_params(which='minor', length=4)


    if savepath != "":
        plt.savefig(savepath)
    else:
        plt.show()

    plt.close('all')


if __name__ == '__main__':

    highmass, lowmass, summary = load_data()
    for i, data in enumerate([highmass, lowmass]):
        info = summary.iloc[i,:]
        savepath = os.path.join(s.plots_folder, f'Test_{i}.png')
        plot_star(data, info = info, savepath = savepath)
