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

fnames = ['FLIP_99_linear_N-10_1_particle_per_cell.data','FLIP_99_cubic-spline_N-10_1_particle_per_cell.data','FLIP_99_cubic-spline_N-20_1_particle_per_cell.data','FLIP_99_linear_N-20_1_particle_per_cell.data','FLIP_99_linear_N-30_1_particle_per_cell.data','FLIP_99_cubic-spline_N-30_1_particle_per_cell.data','FLIP_99_linear_N-10_2_particle_per_cell.data','FLIP_99_cubic-spline_N-10_2_particle_per_cell.data','FLIP_99_cubic-spline_N-20_2_particle_per_cell.data']
fnames = ['FLIP_99_linear_N-30_1_particle_per_cell.data','FLIP_99_cubic-spline_N-30_1_particle_per_cell.data','FLIP_99_linear_N-10_2_particle_per_cell.data','FLIP_99_cubic-spline_N-10_2_particle_per_cell.data','FLIP_99_cubic-spline_N-20_2_particle_per_cell.data','APIC_cubic-spline_N-20_2_particle_per_cell.data','APIC_linear_N-20_2_particle_per_cell.data','APIC_linear_N-20_2_particle_per_cell-separated_F.data']
fnames = ['FLIP_99_linear_N-30_1_particle_per_cell.data','APIC_cubic-spline_N-20_2_particle_per_cell.data','APIC_linear_N-20_2_particle_per_cell.data','APIC_linear_N-20_2_particle_per_cell-separated_F.data','APIC_cubic-spline_N-20_2_particle_per_cell-separated_F.data']
fnames = ['FLIP_99_linear_N-30_1_particle_per_cell.data','APIC_cubic-spline_N-20_2_particle_per_cell.data','APIC_cubic-spline_N-20_2_particle_per_cell-separated_F.data', 'FLIP99_linear_N-30_1_particle_per_cell-separated_F.data', 'FLIP99_cubic-spline_N-30_1_particle_per_cell-separated_F.data','FLIP99_cubic-spline_N-30_1_particle_per_cell.data','FLIP99_linear_N-30_1_particle_per_cell_central_point.data','FLIP99_cubic-spline_N-30_1_particle_per_cell_central_point.data']
fnames = ['FLIP_99_linear_N-30_1_particle_per_cell.data','APIC_cubic-spline_N-20_2_particle_per_cell.data','APIC_linear_N-20_2_particle_per_cell.data','APIC_linear_N-20_2_particle_per_cell-separated_F.data','APIC_cubic-spline_N-20_2_particle_per_cell-separated_F.data']
fnames = ['FLIP_99_linear_N-30_1_particle_per_cell.data','FLIP_99_cubic-spline_N-30_1_particle_per_cell.data','FLIP99_linear_N-30_1_particle_per_cell_central_point.data','FLIP99_cubic-spline_N-30_1_particle_per_cell_central_point.data']
fnames = ['FLIP99_linear_N-30_1_particle_per_cell_central_point.data','FLIP99_cubic-spline_N-30_1_particle_per_cell_central_point.data','APIC_cubic-spline_N-30_1_particle_per_cell_central_point-separateF.data','FEM.data','APIC_linear_N-30_2_particle_per_cell_central_point-separated_F.data','APIC_linear_N-30_2_particle_per_cell_corner_point-separated_F.data']
fnames = ['FEM.data','APIC_linear_N-30_2_particle_per_cell_corner_point-separated_F.data','PIC_linear_N-30_2_particle_per_cell_corner_point-separated_F.data','APIC_cubic-spline_N-30_2_particle_per_cell_corner_point-separateF.data']#,'FLIP_99_linear_N-30_2_particle_per_cell_corner_point-separated_F.data','PIC_cubic_spline_N-30_2_particle_per_cell_corner_point-separated_F.data','FLIP_99_cubic-spline_N-30_2_particle_per_cell_corner_point-separated_F.data']
dat = [pylab.genfromtxt(f,skip_header=1) for f in fnames]

x=[]
y=[]

fig = plt.figure(figsize=(15, 10),facecolor='white')
ax = plt.subplot2grid((1, 1), (0, 0))

lw = 1.5

for i, f in enumerate(fnames):
    if f=='SPH.data':
        x = dat[i][:,0]
        y = dat[i][:,1]
        ax.plot(x, y, label='SPH',linewidth=lw)
    elif f=='FEM.data':
        x = dat[i][:,0]*1e3
        y = dat[i][:,1]
        ax.plot(x, y, label='FEM',linewidth=2*lw)
    else:
        x = dat[i][:,2]
        y = dat[i][:,6]
        ax.plot(x, y, label=f,linewidth=lw)
        
ax.set_xlabel('Time')
ax.set_ylabel('Displacement')
ax.set_ylim(-2,0.5)
ax.set_xlim(0,0.25)

ax.legend(loc=2, prop={'size': 11})
plt.tight_layout()
plt.grid(True)
plt.savefig('./plot.png', format='png')
if args.quiet == False:
    plt.show()
