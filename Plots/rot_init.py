#!/usr/bin/env python
"""
@author: dephil

Program for computing and visualizing initial conditions for (rotating) singular isothermal sphere collapse runs in RAMSES(-RT)

"""
from numpy import pi, sqrt, exp, inf, linspace, logspace, log10, array, abs
from scipy import integrate
from matplotlib import pyplot as plt
from matplotlib import rcParams
rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
rcParams['mathtext.it'] = 'Bitstream Vera Sans'

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
G_cgs = 6.674e-8      # cm^3/(g s^2)
cms_to_kms = 1e-5


def escape_vel(M, R_trunc, factor=1.):
    """
    Returns escape velocity for a singular isothermal sphere
    """
    return factor*sqrt(M*G_cgs/(4*pi*R_trunc))


def solid_rot(r, M, R_trunc):
    """
    Returns solid body rotation profile
    """
    sigma = M*G_cgs/(4*pi*R_trunc)
    omega = sqrt(36.*pi*sigma)/R_trunc
    return array([omega*ri if (ri**2<=R_trunc**2) else omega*ri*exp(10*(R_trunc**2-ri**2)) for ri in r])


def saturated_rot(r, M, R_trunc, R_vortex, factor):
    """
    Returns a rotation profile
    """
    vesc = sqrt(G_cgs*M/R_trunc)
    omega = vesc*sqrt(1./R_trunc**2+1./R_vortex**2)*0.1
    return array([factor*ri*omega*1./sqrt(1+(ri**2/R_vortex**2)) if (ri**2<=R_trunc**2) else factor*ri*omega*1./sqrt(1+(ri**2/R_vortex**2))*exp(10*(R_trunc**2-ri**2))*1./ri  for ri in r])


def plot_profile(ax, r, profile, v_escape=-1., label=None):
    """
    Plots a radial profile
    """
    if v_escape>0: ax.axhline(v_escape, color=my_red, ls='-', lw=2, alpha=0.5)
    ax.loglog(r, profile, color=my_purpleblue, lw=2, ls='-', label=label)


if __name__=="__main__":

    # parameters
    lowM  =   4.*M_sol
    highM = 100.*M_sol
    lowR_trunc  =   4000.*AU
    highR_trunc = 100000.*AU
    R_vortex = 4000.*AU
    #highResR_min = 10*AU
    #lowResR_min  = 250*AU

    # radii
    lowr = linspace(1.*AU, lowR_trunc*1.1, 1000)
    highr = linspace(1.*AU, highR_trunc*1.1, 100000)

    # profiles
    low_vrot   = solid_rot(lowr, lowM, lowR_trunc)
    low_vrot  *= cms_to_kms
    high_vrot  = saturated_rot(highr, highM, highR_trunc, R_vortex, factor=3*10)
    high_vrot *= cms_to_kms
    print high_vrot


    # escape velocity
    low_vesc  = escape_vel(lowM, lowR_trunc, factor=6*sqrt(pi))*cms_to_kms
    high_vesc = escape_vel(highM, highR_trunc, factor=0.1)*cms_to_kms

    # plotting
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot_profile(ax, lowr/AU, low_vrot, v_escape=low_vesc)
    plot_profile(ax, highr/AU, high_vrot, v_escape=high_vesc)
    plt.show()
