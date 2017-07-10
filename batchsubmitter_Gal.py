from glob import glob
import os, sys
import numpy as np
import time
nside = str(32)

flines = open('submittemplatemedhighmem.sh').readlines()

jname='Galproj'

comps = [[1,1],[28,14],[56,26]]
specs = [2.0, 0.9]
densities = [1.e-5, 1.e-4, 1.e-3]
magfields = [1.e-11, 1.e-10, 1.e-9]

for comp in comps:
    for spec in specs:
        for den in densities:
            for mag in magfields:
                print comp, spec, den, mag
                flist = glob('Outputs/*/ProjectedMaps5/AtEarthMap2_'+str(comp[0])+'_'+str(comp[1])+'_SpecIndex_'+str(spec)+'_Den_' + str(den) + '_EgalMag_' +str(mag)+'_S*_Map.txt')
                print len(flist), 'result files found'
                if  len(flist):
                    continue
                flist = glob('Outputs/GalSummary6_'+str(comp[0])+'_'+str(comp[1])+'_SpecIndex_'+str(spec)+'_Den_'+str(den)+'_EgalMag_'+str(mag)+'.txt')
                print len(flist)
                #continue
                jobname = jname+str(comp[0])+str(comp[1])+str(spec)+str(den)+str(mag)+'job'
                fout = open(jobname+'.slurm', "w")
                jobline = ' python GenerateAtEarthMap.py -a '+str(comp[0])+' -z '+str(comp[1]) + ' -s ' + str(spec) + ' -d ' + str(den) + ' -b ' +str(mag)
                for line in flines:
                    fout.write(line.replace('__NAME__', jobname).replace('__JOBLINE__', jobline))
                fout.close()
                os.system('chmod +x ' + jobname+'.slurm')
                os.system('sbatch -p icecube '+ jobname+'.slurm')
                raw_input('test')
