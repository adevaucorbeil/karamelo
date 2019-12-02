from ovito.data import *
from ovito.io import import_file
from ovito.modifiers import *
import glob
import os

filepattern = "/home/alban/karamelo/examples/Scratch/OFHC_copper/with_damage/N_10/dump_g.proc-*.*.LAMMPS"
#dir = os.getcwd()
#print(dir)
#filepattern = dir + "/dump_p.proc-*.*.LAMMPS"
pos_wildcard = filepattern.find('*');
N_proc = int(sorted(glob.glob(filepattern))[-1].split(".proc-")[1].split(".")[0])+1

pipeline = import_file(filepattern[0:pos_wildcard] + "0" + filepattern[pos_wildcard+1:])
for i in range(1, N_proc):
        modifier = CombineDatasetsModifier()
        modifier.source.load(filepattern[0:pos_wildcard] + str(i) + filepattern[pos_wildcard + 1:])
        pipeline.modifiers.append(modifier)
pipeline.add_to_scene()
