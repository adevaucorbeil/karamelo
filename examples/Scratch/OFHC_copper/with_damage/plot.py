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
from glob import glob
import os

def find_cells_per_unit_length(f):
    posI = f.find('N_')
    posE = f.find('/log.mpm')
    return int(f[posI+2:posE])

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

fnames = []
start_dir = os.getcwd()
pattern   = "log.mpm"

for dir,_,_ in os.walk(start_dir):
    fnames.extend(glob(os.path.join(dir,pattern))) 

dat = [pylab.genfromtxt(f,skip_header=1) for f in fnames]

x=[]
y=[]
xmax = 0#0.3

#f = plt.subplots()
fig = plt.figure(figsize=(5, 5),facecolor='white')
ax1 = plt.subplot2grid((1, 1), (0, 0))

lw = 1.5
ymax = 0#800

for i, f in enumerate(fnames):
    x = -dat[i][:,4]
    y = -dat[i][:,6]
    label = find_cells_per_unit_length(f)
    ax1.plot(x, y, label="N="+str(label), linewidth=lw)
    ymax = max(ymax, np.max(y))
    xmax = max(xmax, np.max(x))
    
ax1.set_xlabel('Indenter vertical displacement (mm)')
ax1.set_ylabel('Force on indenter (kN)')
#ax1.set_ylim(0, ymax)
#ax1.set_xlim(0, xmax)

ax1.legend(loc=2, prop={'size': 11})
plt.tight_layout()
plt.savefig('./plot.pdf', format='pdf')
if args.quiet == False:
    plt.show()
