import numpy
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
import os
from glob import glob
import pdb

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

fnames = []
start_dir = os.getcwd()
pattern   = "log.mpm"


for dir,_,_ in os.walk(start_dir):
    fnames.extend(glob(os.path.join(dir,pattern))) 

list_l0 = []
for i, f in enumerate(fnames):
    label = f[f.find('cellsize_')+9:f.find('/log.mpm')]
    l0 = label.replace('dot', '.')
    list_l0.append(float(l0))
#pdb.set_trace()

fnames = [x for _,x in sorted(zip(list_l0,fnames), reverse=True)]
print(fnames)
dat = [pylab.genfromtxt(f,skip_header=1,skip_footer=1) for f in fnames]

x=[]
y=[]

axs = []
#fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
fig = plt.figure(figsize=(5, 5),facecolor='white')
axs.append(plt.subplot2grid((1, 1), (0, 0)))
#axs.append(plt.subplot2grid((1, 2), (0, 1)))

lw = 1.5

max_f = []
disp_maxf = []
disp_failure = []
list_l0 = []

for i, f in enumerate(fnames):
    if "old" not in f:
        disp = dat[i][:,5]*5 # convert eps to disp
        force = dat[i][:,6]/0.05 #0.05 to correct from stress to force
        label = f[f.find('cellsize_')+9:f.find('/log.mpm')]
        l0 = label.replace('dot', '.')
        list_l0.append(float(l0))
        #force = force/float(l0)
        max_f.append(numpy.max(force))
        td = disp[disp>0.7]
        tf = force[disp>0.7]
        #print(t[numpy.where(force == max_f[-1])[0][0]])
        #pdb.set_trace()
        #disp_maxf.append(disp[numpy.where(force == max_f[-1])[0][0]])
        #try:
        #    disp_failure.append(td[numpy.where(tf <= 0)[0][0]])
        #except:
        #    disp_failure.append(disp_failure[-1])
        #print(disp_maxf[-1])
        #print("{}\t{}".format(list_l0[-1],disp_failure[-1]))
        axs[0].plot(disp, force, label=r"$l_0="+l0+"$ mm", linewidth=lw)

#axs[0].annotate("", xy=(1.25, 4.5), xytext=(0.5, 2), arrowprops=dict(arrowstyle="<-"))

disp_failure = numpy.array(disp_failure)
#axs[1].plot(list_l0, disp_failure, '.')

#ax.set_xlabel('Time (ms)')
axs[0].set_xlabel('Displacement (mm)')
axs[0].set_ylabel('Force (kN)')
axs[0].set_ylim(0,7.5)
axs[0].set_xlim(0,0.6)


#axs[1].set_xlabel(r'Cell size (mm)')
#axs[1].set_ylabel('Displacement at failure (mm)')

#axs[1].set_xscale('log')
#axs[1].set_yscale('log')
#axs[1].set_ylim(1e-2,2)
#axs[1].set_xlim(1e-2,2e-1)

axs[0].legend(loc=1, prop={'size': 11})
plt.tight_layout()
plt.grid(True)
plt.savefig('./force_displacement.pdf', format='pdf')
if args.quiet == False:
    plt.show()
