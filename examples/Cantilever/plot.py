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


#-------------------------------------------------------------
# Dimer 5: Read Data

fnames = ['FEM.data','FLIP-099_linear_N-20_2_particles_per_cell-incrementalF.data']#,'SPH_converged.data', 'SPH_HG_0dot1.data', 'SPH_HG_2.data','FLIP-099_linear_N-20_2_particles_per_cell-incrementalF.data']#,'FLIP-099_linear_N-10_1_particles_per_cell.data','FLIP-099_linear_N-20_1_particles_per_cell.data','FLIP-099_linear_N-20_1_particles_per_cell-incrementalF.data','FLIP-099_linear_N-20_2_particles_per_cell-incrementalF.data']#'APIC_linear_N-10_2_particles_per_cell.data','APIC_cubic-spline_N-10_2_particles_per_cell.data','FLIP-099_linear_N-10_1_particles_per_cell.data','FLIP-099_linear_N-20_1_particles_per_cell.data']

dat = [pylab.genfromtxt(f,skip_header=1) for f in fnames]

x=[]
y=[]

fig = plt.figure(figsize=(15, 10),facecolor='white')
ax = plt.subplot2grid((1, 1), (0, 0))

lw = 1.5

for i, f in enumerate(fnames):
    if "SPH" in f:
        x = dat[i][:,0]
        y = dat[i][:,1]
        ax.plot(x, y, label='SPH',linewidth=lw)
    elif f=='FEM.data':
        x = dat[i][:,0]
        y = dat[i][:,1]
        ax.plot(x, y, label='FEM',linewidth=2*lw)
    else:
        x = dat[i][:,2]
        y = dat[i][:,6]+0.014
        ax.plot(x, y, label=f,linewidth=lw)
        
ax.set_xlabel('Time')
ax.set_ylabel('Displacement')
#ax.set_ylim(-0.3,0.1)
ax.set_xlim(0,3)

ax.legend(loc=2, prop={'size': 11})
plt.tight_layout()
plt.grid(True)
plt.savefig('./plot.png', format='png')
if args.quiet == False:
    plt.show()
