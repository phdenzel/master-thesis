##!/usr/bin/env python
"""
Program using an example function (polynomial) to explain the finite volume method
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
#rcParams['mathtext.fontset'] = "stix"
rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
rcParams['mathtext.it'] = 'Bitstream Vera Sans'
rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

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

def func(x):
    """
    Example function
    """
    return 0.2*(x-2)**3+1.4*(x-2)**2+(x-2)

def lin_functional(roots, xmax=1.0):
    """
    Returns characteristics in x,t diagram as tuple of points in x from linear functions for each root upto xmax
    """
    characteristics = []
    pos_slopes = [2**i for i in xrange(-4, -1)]
    neg_slopes = [-2**i for i in xrange(-4, -1)]
    for n in roots:
        for a in pos_slopes:
            characteristics.append([n, xmax/a+n]) # each for x = [0 ,xmax]
        for a in neg_slopes:
            characteristics.append([n, xmax/a+n])
    return characteristics

def piecewise(x, data, bins):
    """
    Returns a binned and averaged array of the same length as data
    """
    bin_len = int(float(len(data))/bins)
    pw = []
    new_x = []
    sum_d = 0.
    sum_x = 0.
    for i, d in enumerate(data):
        sum_d += d
        sum_x += x[i]
        if (i+1)%bin_len==0:
            if bar:
                pw += [sum_d/bin_len]#*bin_len
                new_x += [sum_x/bin_len]#*bin_len
            else:
                pw += [sum_d/bin_len]*bin_len
                new_x += [sum_x/bin_len]*bin_len
            sum_d = 0.
            sum_x = 0.
    return new_x, pw # new_x not needed if bar=False

def plot_lagrange():
    global bar
    # continuous data
    x = np.linspace(-4, 4, 10000)
    y = func(x)
    # piecewise constant function no bar
    bar = False
    pos_nobar, p_nobar = piecewise(x, y, 20)
    # piecewise constant function w/ bar
    bar = True
    pos, p = piecewise(x, y, 20)
    # plotting
    fig = plt.figure(figsize=(14, 6), dpi=100)
    ax = fig.add_subplot(111)
    ax.plot(x, y, ls='-', color=my_red, label="continuous", lw=2)
    ax.plot(x, p_nobar, ls='-', color=my_purpleblue, label="piecewise constant", lw=2)
    ax.bar(pos, p, 8./20, fc="None", ec='gray', align='center')
    plt.legend(fancybox=True).get_frame().set_alpha(0.5)
    ax.yaxis.set_ticks([])
    ax.xaxis.set_ticks(list((np.array(pos[1:])+np.array(pos[:-1]))/2))
    tlabls = [' ']*len(pos[1:])
    tlabls[len(pos[1:])//2-7] = r"x$_{i-\frac{1}{2}}$"
    tlabls[len(pos[1:])//2-6] = r"x$_{i+\frac{1}{2}}$"
    ax.xaxis.set_ticklabels(tlabls)
    ax.tick_params(axis='x', labelsize=22)
    rcParams['xtick.direction'] = 'out'
    rcParams['ytick.direction'] = 'out'
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.set_xlim(xmin=-4, xmax=-1)
    ax.set_ylim(ymin=0, ymax=8)
    ax.set_axisbelow(True)
    ax.annotate(r'u$_{i}$',
                xy=(pos[len(pos[1:])//2-6]-0.1, p[len(pos[1:])//2-6]),
                xytext=(-2.5, 7),
                arrowprops=dict(arrowstyle="simple", fc="0.6", ec="none"),
                fontsize=30,
                size=30)
    #ax.set_axis_off()
    # zoomed inset
    axins = zoomed_inset_axes(ax, 3.5, loc=2, bbox_to_anchor=(-0.005, 0.97), bbox_transform=ax.figure.transFigure)
    axins.plot(x[875:1125], p_nobar[875:1125], ls='-', color=my_purpleblue, lw=2)
    axins.set_ylim(ymin=3.8, ymax=5.08)
    axins.set_xlim(xmin=-3.315, xmax=-3.08)
    axins.yaxis.set_ticks([])
    axins.xaxis.set_ticks([])
    mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.5")
    axins.set_title("Riemann problem", fontsize=20)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(False)
    plt.savefig('Plots/piecewise_u.pdf', transparent=True)
    plt.show()
    #plt.close()

def plot_characteristics():
    global bar
    # x axis
    x = np.linspace(-4, 4, 10000)
    dummy_data = func(x)
    # binning x axis for cell discretization
    bar = True
    pos, p = piecewise(x, dummy_data, 20)
    cs = lin_functional(pos, xmax=.2)
    # start plotting
    fig = plt.figure(figsize=(18, 6), dpi=200)
    ax = fig.add_subplot(111)
    for c in cs:
        ax.plot(c, [0., .2], ls='-', color=my_purpleblue, label="characteristics", lw=2)
    ax.plot([pos[len(pos[1:])//2-7], pos[len(pos[1:])//2-7]], [0., .01], ls='--', color=my_red, label="average solution", lw=3)
    ax.plot([pos[len(pos[1:])//2-7], pos[len(pos[1:])//2-6]], [.00995, .00995], ls='--', color=my_red, label="average solution", lw=3)
    ax.plot([pos[len(pos[1:])//2-6], pos[len(pos[1:])//2-6]], [.01, 0.], ls='--', color=my_red, label="average solution", lw=3)
    ax.yaxis.set_ticks([])
    ax.xaxis.set_ticks(pos)
    tlabls = [' ']*len(pos[1:])
    tlabls[len(pos[1:])//2-7] = r"x$_{i-\frac{1}{2}}$"
    tlabls[len(pos[1:])//2-6] = r"x$_{i+\frac{1}{2}}$"
    ax.xaxis.set_ticklabels(tlabls)
    ax.tick_params(axis='x', labelsize=22)
    rcParams['xtick.direction'] = 'out'
    rcParams['ytick.direction'] = 'out'
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.set_xlim(xmin=-4, xmax=-2)
    ax.set_ylim(ymin=0, ymax=.01)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(False)
    ax.annotate('\nt',
                xy=(-4.0, 0.01),
                xytext=(-4.012, 0.005), # 4.021
                arrowprops=dict(arrowstyle="->", fc="0.2", lw=3.),
                fontsize=25,
                size=25)
    plt.savefig('Plots/characteristics.pdf', transparent=True)
    plt.show()
    plt.close()

def plot_hll():
    t_s = 1.
    x_L = -1.
    x_R = 1.
    sRts = 0.8
    sLts = -0.8
    x = np.linspace(x_L, x_R, 10000)
    # start plotting
    fig = plt.figure(figsize=(18, 6), dpi=100)
    ax = fig.add_subplot(111)
    #ax.axvline(0, color='black')
    # horizontal lines
    ax.plot([x_L, x_R], [t_s, t_s], color=my_grey, lw=1.5, alpha=0.5)
    # vertical lines
    ax.plot([x_L, x_L], [0, t_s], color=my_grey, lw=1.5, alpha=0.5)
    ax.plot([x_R, x_R], [0, t_s], color=my_grey, lw=1.5, alpha=0.5)
    ax.plot([sLts, sLts], [0, t_s], color=my_grey, ls='--', lw=1.5, alpha=0.5)
    ax.plot([sRts, sRts], [0, t_s], color=my_grey, ls='--', lw=1.5, alpha=0.5)
    # characteristics
    ax.plot([0., sRts], [0., t_s], color=my_red, lw=2)
    ax.plot([0., sLts], [0., t_s], color=my_red, lw=2)
    ax.plot([0., 0.3], [0., t_s], color=my_purpleblue, ls='--', lw=2)
    # annotate t axis
    ax.annotate('t',
                xy=(0., 1.02),
                xytext=(-0.015, -0.1), # 4.021
                arrowprops=dict(arrowstyle="->", fc='black', lw=1),
                fontsize=25,
                size=25)
    # axis work
    ax.set_xlim(xmin=-1.02, xmax=1.02)
    ax.set_ylim(ymin=0, ymax=1.02)
    ax.yaxis.set_ticks([])
    ax.xaxis.set_ticks([x_L, sLts, sRts, x_R])
    tlabls = [r'x$_{L}$', r't$_{S}$s$_{L}$', r't$_{S}$s$_{R}$', r'x$_{R}$']
    ax.xaxis.set_ticklabels(tlabls)
    ax.tick_params(axis='x', labelsize=22)
    ax.set_axisbelow(True)
    rcParams['xtick.direction'] = 'out'
    rcParams['ytick.direction'] = 'out'
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(False)
    ax.text(-0.9, .35, r'u$_{L}$', fontsize=20)
    ax.text(0.85, .35, r'u$_{R}$', fontsize=20)
    ax.text(-0.15, .69, r'u$_{HLL}$', fontsize=20)
    # save figure
    plt.savefig('Plots/hll_regions.pdf', transparent=True)
    #plt.show()


def piecewise_linear(pos, pw_const):
    """
    Piecewise linear function around the bin centers
    """
    bin_len = len(pw_const)/len(pos)
    dx = pos[1]-pos[0]
    x = np.linspace(pos[0]-dx/2., pos[-1]+dx/2., len(pw_const))
    pwc_centers = [pw_const[i] for i in range(bin_len/2, len(pw_const)-1, bin_len)]
    print len(pwc_centers)
    pwl = []
    pwl_value = 0
    slope = 3.5
    for i, pwc in enumerate(pw_const):
        pwl_value = pwc + slope*(x[i]-pos[i//bin_len])
        if (i+1)%bin_len==0 and len(pos)-2>(i)//bin_len>0:
            slope = ((pwc_centers[(i+1)//bin_len+1]-pwc_centers[(i-1)//bin_len])*0.5)/(dx)
        elif((i)//bin_len-1==0):
            slope = 4.
        pwl.append(pwl_value)
    return pwl


def plot_pwlinear():
    """
    Plots the piecewise linear data
    """
    global bar
    # continuous data
    x = np.linspace(-4, 4, 10000)
    y = func(x)
    # piecewise constant function no bar
    bar = False
    pos_nobar, p_nobar = piecewise(x, y, 20)
    # piecewise constant function w/ bar
    bar = True
    pos, p = piecewise(x, y, 20)
    # piecewise linear function
    pwl = piecewise_linear(pos, p_nobar)
    # plotting
    fig = plt.figure(figsize=(14, 6), dpi=100)
    ax = fig.add_subplot(111)
    #ax.plot(x, y, ls='--', color=my_red, label="continuous", lw=1)
    ax.plot(x, pwl, ls='-', color=my_purpleblue, label="piecewise linear", lw=2)
    ax.plot(x, p_nobar, ls='--', color=my_red, label="piecewise constant", lw=2)
    ax.bar(pos, p, 8./20, fc="None", ec='gray', align='center')
    plt.legend(fancybox=True).get_frame().set_alpha(0.5)
    ax.yaxis.set_ticks([])
    ax.xaxis.set_ticks(list((np.array(pos[1:])+np.array(pos[:-1]))/2))
    tlabls = [' ']*len(pos[1:])
    tlabls[len(pos[1:])//2-7] = r"x$_{i-\frac{1}{2}}$"
    tlabls[len(pos[1:])//2-6] = r"x$_{i+\frac{1}{2}}$"
    ax.xaxis.set_ticklabels(tlabls)
    ax.tick_params(axis='x', labelsize=22)
    rcParams['xtick.direction'] = 'out'
    rcParams['ytick.direction'] = 'out'
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    #ax.set_xlim(xmin=-3.4, xmax=-1)
    ax.set_xlim(xmin=-3.6, xmax=-1)
    ax.set_ylim(ymin=0, ymax=8)
    ax.set_axisbelow(True)
    #ax.set_axis_off()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(False)
    #plt.savefig('Plots/pw_linear.pdf', transparent=True)
    plt.show()
    #plt.close()

#plot_lagrange()
#plot_characteristics()
plot_hll()
#plot_pwlinear()
