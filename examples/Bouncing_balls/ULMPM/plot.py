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

fnames = ['log.mpm']

dat = [pylab.genfromtxt(f,skip_header=1) for f in fnames]

x=[]
y=[]

#f = plt.subplots()
fig = plt.figure(figsize=(5, 5),facecolor='white')
ax1 = plt.subplot2grid((1, 1), (0, 0))

lw = 1.5

for i, f in enumerate(fnames):
    time = dat[i][:,2]
    ek = dat[i][:,3]
    es = dat[i][:,4]
    et = ek+es
    ax1.plot(time, ek, color='blue', label=r'$E_k$',linewidth=lw)
    ax1.plot(time, es, color='green', label=r'$E_s$',linewidth=lw)
    ax1.plot(time, et, color='red', label=r'$E_{tot}$',linewidth=lw)
    xmax = max(time)
    ymax = max([max(ek),max(es)])*1.4

ymax = 3.5  
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Energy (J)')
ax1.set_ylim(0, ymax)
ax1.set_xlim(0, xmax)

ax1.legend(loc=1, prop={'size': 11})
plt.tight_layout()
plt.savefig('./bouncing_balls.pdf', format='pdf')
if args.quiet == False:
    plt.show()
