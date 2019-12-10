from ovito.data import *
from ovito.io import import_file
from ovito.modifiers import *
from ovito.vis import *
from ovito import *
import glob
import os
import ovito

filepattern = "/home/alban/karamelo/examples/FSI/dump_p.proc-*.*.LAMMPS"
#filepattern = "/home/alban/temp/dump_p.proc-*.*.LAMMPS.gz"
#dir = os.getcwd()
#print(dir)
#filepattern = dir + "/dump_p.proc-*.*.LAMMPS"
#N_proc = int(sorted(glob.glob(filepattern))[-1].split(".proc-")[1].split(".")[0])+1

pos_wildcard = filepattern.find('*');
N_proc = 0
while glob.glob(filepattern[0:pos_wildcard] + str(N_proc) + filepattern[pos_wildcard+1:]):
        N_proc += 1

pipeline = import_file(filepattern[0:pos_wildcard] + "0" + filepattern[pos_wildcard+1:])
for i in range(1, N_proc):
        modifier = CombineDatasetsModifier()
        modifier.source.load(filepattern[0:pos_wildcard] + str(i) + filepattern[pos_wildcard + 1:])
        pipeline.modifiers.append(modifier)
pipeline.add_to_scene()


vp = ovito.dataset.viewports.active_vp
vp.zoom_all()
