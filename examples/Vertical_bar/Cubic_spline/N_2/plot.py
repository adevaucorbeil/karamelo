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
# Read Data

fnames = ['log.mpm']
dat = [pylab.genfromtxt(f,skip_header=1) for f in fnames]

x=[]
y=[]

fig = plt.figure(figsize=(5, 5),facecolor='white')
ax = plt.subplot2grid((1, 1), (0, 0))

lw = 1.5

for i, f in enumerate(fnames):
    x = dat[i][:,2]
    y = dat[i][:,6]
    ax.plot(x, y, label=f,linewidth=lw)
        
ax.set_xlabel('Time (s)')
ax.set_ylabel('Displacement (m)')
ax.set_ylim(-2,0)
ax.set_xlim(0,0.25)

#ax.legend(loc=2, prop={'size': 11})
plt.tight_layout()
plt.grid(True)
plt.savefig('./plot.pdf', format='pdf')
if args.quiet == False:
    plt.show()
