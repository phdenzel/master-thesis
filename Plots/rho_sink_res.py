#!/usr/bin/env python
"""
@author: dephil

Program for computing and visualizing initial conditions for (rotating) non-singular isothermal sphere collapse runs in RAMSES(-RT)

"""
from numpy import pi, sqrt, exp, inf, linspace, logspace, log10, array, abs
from scipy import integrate
from matplotlib import pyplot as plt
from matplotlib import rcParams
#rcParams['mathtext.fontset'] = "stix"
rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
rcParams['mathtext.it'] = 'Bitstream Vera Sans'
#rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'

# custom colors
my_blue = '#39b3e6'
my_red = '#fe4365'
my_green = '#91c442'
my_purple = '#9950cb'
my_orange = '#da9605'
my_yellow = '#edde45'
my_turquoise = '#00d1a4'
my_pink = '#fd4d83'
my_neongreen = '#0bf759'
my_grey = '#aab1b7'
my_brown = '#d88c4e'
my_purpleblue = '#603dd0'
my_darkgreen = '#0d8e66'
my_winered = '#7c3658'
my_gold = '#c18d1d'
my_greyishgreen = '#57776f'
my_colors = [my_blue, my_red, my_green, my_purple, my_orange, my_yellow, my_turquoise, my_pink, my_neongreen, my_grey, my_brown, my_purpleblue, my_darkgreen, my_winered, my_gold, my_greyishgreen]


# constants
AU = 1.495978707e13  # cm
M_sol = 1.98855e33   # g
Gcgs = 6.674e-8      # cm^3/(g s^2)

cs = 20000 # cm/s
conversion_to_gcc = 5.939299e-7

def rho_sink(l):
    return cs**2/(16*Gcgs*l**2)

def rho_kappa(l, tau, kappa):
    return tau/(kappa*l)

def schematic_plot():
    r = logspace(log10(0.1), log10(50), 1e4)
    rho_s = [rho_sink(ri*AU) for ri in r]
    rho_k = [rho_kappa(ri*AU, 1000, 0.15) for ri in r]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.loglog(r, rho_s, color=my_purpleblue, lw=2, ls='-', label=r'$\rho_{sink}$')
    ax.loglog(r, rho_k, color=my_red, lw=2, ls='-', label=r'$\rho_{\kappa}$')
    ax.axhline(1e-13, color=my_grey, ls='--', lw=2, alpha=0.5)
    ax.axvline(4, color=my_grey, ls='--', lw=2, alpha=0.5)
    ax.text(0.25, 0.7, r'$\propto$L$^{-2}$', fontsize=18, transform=ax.transAxes)
    ax.text(0.15, 0.5, r'$\propto$L$^{-1}$', fontsize=18, transform=ax.transAxes)
    ax.text(0.5, 1.0075, r'$\sim$4 AU', fontsize=16, transform=ax.transAxes)
    ax.set_xlabel('L [AU]', fontsize=18)
    ax.set_ylabel(r'$\rho$ [g/cc]', fontsize=18)
    plt.legend(fancybox=True, fontsize=20).get_frame().set_alpha(0.5)
    #plt.savefig('Plots/sink_res.pdf')
    plt.show()


if __name__=="__main__":
    schematic_plot()
