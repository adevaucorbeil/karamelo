import numpy as np
from scipy import math
from numpy.random import random
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import pylab as pylab
import glob 
import matplotlib.ticker as mtick
import matplotlib.colors as mcol
import time
import argparse

parser = argparse.ArgumentParser(description='Plot the evolution of variables.',
                                 add_help=False,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)


parser.add_argument('--quiet', '-q', action='store_const', const=True, default=False,
                    help='do not show plot.')
args = parser.parse_args()

#-------  Setup fonts
plt.rc('font',family='Helvetica',size=14)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('axes',linewidth=1.)
plt.rc('pdf',fonttype=3)
plt.rc('mathtext',fontset='stixsans')
plt.rc("figure", facecolor="none")

markersize = 6

#-------- Setup colours
#bgc = [(252./255.)**0.05, (141./255.)**0.05,(98./255.)**0.05]
#colors = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors


#-------------------------------------------------------------
# Dimer 5: Read Data

fnames = ['Weldox460E_smooth_cylindrical_FEM_with_damage.data', 'Weldox460E_smooth_FEM_with_damage.data', 'Linear/cellsize_0dot2/log.mpm','Cubic_spline/one_particle_per_cell/cellsize_0dot2/log.mpm','Bernstein/one_particle_per_cell/cellsize_0dot2/log.mpm','Linear/cellsize_0dot2/cylinder/log.mpm']

dat = [pylab.genfromtxt(f,skip_header=1) for f in fnames]

x=[]
y=[]
xmax1 = 0.3

#f = plt.subplots()
fig = plt.figure(figsize=(5, 5),facecolor='white')
ax1 = plt.subplot2grid((1, 1), (0, 0))

lw = 1.5
ymax = 800

for i, f in enumerate(fnames):
    if 'FEM' in f:
        x = dat[i][:,0]
        y = dat[i][:,1]
    else:
        x = dat[i][:,5]
        y = dat[i][:,6]*1000
    if 'Linear' in f:
        ax1.plot(x, y, color='red', label='TLMPM - Linear',linewidth=lw)
    elif 'Cubic_spline' in f:
        ax1.plot(x, y, color='blue', label='TLMPM - Cubic splines',linewidth=lw)
    elif 'Bernstein' in f:
        ax1.plot(x, y, color='green', label='TLMPM - Bernstein',linewidth=lw)
    elif 'FEM' in f:
        ax1.plot(x, y, color='black', label='FEM',linewidth=lw)
    else:
        print('ERROR: I don\'t know what to do with {}'.format(f))
    ymax = max(ymax, np.max(y))
    #if '460E' in f:
    #    xmax1 = max(xmax1, np.max(x[i]))
    #elif '700E' in f:
    #    xmax2 = max(xmax2, np.max(x[i]))
    #elif '900E' in f:
    #    xmax3 = max(xmax3, np.max(x[i]))

    
ax1.set_xlabel('Engineering strain')
ax1.set_ylabel('Engineering stress (MPa)')
ax1.set_ylim(0, ymax)
ax1.set_xlim(0, xmax1)

ax1.legend(loc=1, prop={'size': 11})
plt.tight_layout()
plt.savefig('./tensile_with_damage.pdf', format='pdf')
if args.quiet == False:
    plt.show()
