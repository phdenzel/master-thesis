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


def unit_parameters(M, R_min, R_trunc, R_vortex):
    """
    Input number of AU
    Sets the parameters
        R_trunc:    truncation radius of the non-singular isothermal sphere
        R_min:      minimal radius of the non-singular isothermal sphere; rho(R_min) = rho(0)/2
        R_vortex:   radius at which rotation velocity almost doesn't increase anymore
    Returns given parameters in AU
    """
    return [M*M_sol, R_min*AU, R_trunc*AU, R_vortex*AU]


def nsis_profile(r, pars):
    """
    Returns the density at a given radius of a non-singular isothermal truncated sphere profile with given parameters
        rho = M / (4*pi*R_trunc*(r^2 + R_min^2))
    """
    M, R_min, R_trunc, R_vortex = pars[0], pars[1], pars[2], pars[3]
    r2 = r*r
    rho = M/(4*pi*R_trunc*(r2+R_min*R_min))
    if (r2 > R_trunc*R_trunc): rho *= 1e-4*exp(10*(R_trunc-r))
    return rho

def nsis_shell(r, pars):
    """
    Function to integrate non-singular isothermal profile
    """
    return nsis_profile(r, pars)*4*pi*r*r

def get_integrated_mass(pars):
    """
    Get mass from integration of density profile in units of solar masses
    """
    mass_number = pars[0]
    pars = unit_parameters(pars[0], pars[1], pars[2], pars[3])
    mass = integrate.quad(nsis_shell, 0, 8*pars[2], args=pars)
    return mass_number/(mass[0]/M_sol)

def omega_alpha(pars, alpha):
    """
    Simple solid body rotation in cgs
    """
    M, R_min, R_trunc, R_vortex = pars[0], pars[1], pars[2], pars[3]
    return sqrt(9*M*Gcgs*M_sol*alpha/((AU*R_trunc)**3))

def omega_rot(r, pars):
    """
    Returns the rotation frequency for given parameters
    """
    M, R_min, R_trunc, R_vortex = pars[0], pars[1], pars[2], pars[3]
    r2 = r*r
    R_trunc2 = R_trunc*R_trunc
    R_vortex2 = R_vortex*R_vortex
    omega = 0.1*sqrt((R_vortex2+R_trunc2)/(R_vortex2*R_trunc2))*sqrt(M/R_trunc) # Kepler frequency limit
    omega1 = omega * 1./sqrt(1+(r2/R_vortex2))
    if (r2 > R_trunc2):
        omega1 *= exp(10*(R_trunc2-r2))*1./r
        omega2 = omega * exp(10*(R_trunc2-r2))*1./r
    else:
        omega2 = omega * 1./sqrt(1+(r2/R_vortex2))
    return [omega1, omega2]


def v_rot(r, pars):
    """
    Returns the rotation velocity for given parameters
    """
    return [omega_rot(r, pars)[0]*r, omega_rot(r, pars)[1]*r]

def plot_nsis(pars):
    """
    Plots an nsis profile
    """
    M, R_min, R_trunc, R_vortex = pars[0], pars[1], pars[2], pars[3]
    r = logspace(log10(0.1), log10(26*R_trunc), 1e4)
    conversion_to_gcc = 5.939299e-7
    rho = [conversion_to_gcc*nsis_profile(ri, pars) for ri in r]
    rho_scaled = [conversion_to_gcc*nsis_profile(ri, [100., 10., 100000., 4000.]) for ri in r]
    rho_selfsim = [conversion_to_gcc*nsis_profile(ri, [100., 250., 100000., 4000.]) for ri in r]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.loglog(r, rho, color=my_turquoise, lw=2, ls='-', label=r'low mass, high res.')
    ax.loglog(r, rho_scaled, color=my_purpleblue, lw=2, ls='--', label=r'high mass, high res.')
    ax.loglog(r, rho_selfsim, color=my_red, lw=2, ls='--', label=r'high mass, low res.')
    ax.axvline(R_min, color=my_grey, ls='--', lw=2, alpha=0.5)
    ax.axvline(100000., color=my_grey, ls='--', lw=2, alpha=0.5)
    ax.axvline(4000., color=my_grey, ls='--', lw=2, alpha=0.5)
    ax.axvline(250., color=my_grey, ls='--', lw=2, alpha=0.5)
    ax.axhline(1e-13, color=my_grey, ls='--', lw=2, alpha=0.5)
    ax.axhline(1.6e-16, color=my_grey, ls='--', lw=2, alpha=0.5)
    ax.text(0.7, 0.915, r'$\rho_{sink}$ at high res.', fontsize=16, transform=ax.transAxes)
    ax.text(0.7, 0.6, r'$\rho_{sink}$ at low res.', fontsize=16, transform=ax.transAxes)
    #ax.text(0.01, 1.02, r'$\Delta$x$_{max}$', fontsize=18, transform=ax.transAxes)
    ax.text(0.03, 1.02, r'high res. R$_{min}$', fontsize=16, transform=ax.transAxes)
    ax.text(0.31, 1.02, r'low res. R$_{min}$', fontsize=16, transform=ax.transAxes)
    ax.text(0.85, 1.02, r'R$_{trunc}$ at 100 M$_{\odot}$', fontsize=16, transform=ax.transAxes)
    ax.text(0.6, 1.02, r'R$_{trunc}$ at 4 M$_{\odot}$', fontsize=16, transform=ax.transAxes)
    ax.text(0.525, 0.625, r'$\propto$r$^{-2}$', fontsize=16, transform=ax.transAxes)
    ax.set_xlim(xmin=1, xmax=200000)
    ax.set_ylim(ymin=1e-21, ymax=7.5e-13)
    ax.set_xlabel('r [AU]', fontsize=16)
    ax.set_ylabel(r'$\rho$ [g/cc]', fontsize=16)
    plt.legend(loc='lower left', fancybox=True, fontsize=16).get_frame().set_alpha(0.5)
    plt.savefig('Plots/new_nsis.pdf', transparent=True)
    plt.show()

def plot_rot(pars):
    """
    Plots a solid body rotation of the nsis
    """
    M, R_min, R_trunc, R_vortex = pars[0], pars[1], pars[2], pars[3]
    r_old = logspace(log10(0.1), log10(1.0*R_trunc), 1e4)
    r = logspace(log10(0.1), log10(25.0*R_trunc), 1e4)
    v_r = [omega_alpha(pars, 0.5)*ri*AU for ri in r_old[:-1]]
    v_new = [omega_rot(ri*AU, [100.*M_sol, 10.*AU, 100000.*AU, 4000.*AU])[0]*ri*AU for i, ri in enumerate(r)]
    v_r.append(1)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.loglog(r_old, v_r, color=my_turquoise, lw=2, ls='-')
    ax.loglog(r, 0.0005*array(v_new), color=my_purpleblue, lw=2, ls='--')
    #ax.axvline(R_min, color=my_grey, ls='--', lw=2, alpha=0.5)
    ax.axvline(R_trunc, color=my_grey, ls='--', lw=2, alpha=0.5)
    #ax.axvline(3.9, color=my_grey, ls='--', lw=2, alpha=0.5)
    #ax.axhline(sqrt(4.5*Gcgs*M*M_sol/(R_trunc*AU)), color=my_grey, ls='--', lw=2, alpha=0.5)
    ax.axhline(2*sqrt(Gcgs*100.*M_sol/(100000.*AU)), color=my_grey, ls='--', lw=2, alpha=0.5)
    ax.text(0.89, 0.91, r'v$_{E}$', fontsize=16, transform=ax.transAxes)
    ax.text(0.83, 0.73, r'10 % v$_{E}$', fontsize=16, transform=ax.transAxes)
    #ax.text(0.01, 1.02, r'$\Delta$x$_{max}$', fontsize=18, transform=ax.transAxes)
    #ax.text(0.145, 1.02, r'R$_{min}$', fontsize=18, transform=ax.transAxes)
    ax.text(0.72, 1.02, r'R$_{trunc}$ at 4 M$_{\odot}$', fontsize=18, transform=ax.transAxes)
    #ax.set_xlim(xmin=3, xmax=6000)
    ax.set_ylim(ymin=1, ymax=1e6)
    ax.set_xlabel('r [AU]', fontsize=18)
    ax.set_ylabel(r'$v_{rot}$ [cm/s]', fontsize=18)
    plt.savefig('Plots/solid_rot.pdf', transparent=True)
    plt.show()


if __name__=="__main__":

    #params = unit_parameters(4., 1e-4, 4000., 4000.)
    params = [4., 10, 4000., 4000.]

    #plot_nsis(params)

    plot_rot(params)
