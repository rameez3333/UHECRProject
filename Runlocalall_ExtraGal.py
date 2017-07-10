from glob import glob
import os, sys
import numpy as np
import time
nside = str(32)

flines = open('submittemplateNormal.sh').readlines()

jname='Egalproj'

comps = [[56,26]]
specs = [0.9]
densities = [1.e-3]
magfields = [1.e-11, 1.e-10, 1.e-9]

for comp in comps:
    for spec in specs:
        for den in densities:
            for mag in magfields:
                flist = glob('ProjectedMaps3/EGal_Injection_'+str(comp[0])+'_'+str(comp[1])+'_SpecIndex_'+str(spec)+'_Den_' + str(den) + '_EgalMag_' +str(mag)+'_S*_Map.txt')
                print len(flist), 'result files found'
                print comp, spec, den, mag
                #continue
                if len(flist) >20:
                    continue
                jobname = jname+str(comp[0])+str(comp[1])+str(spec)+str(den)+'job'
                #fout = open(jobname+'.slurm', "w")
                
                jobline = 'python GenerateExtraGalacticMap.py -a '+str(comp[0])+' -z '+str(comp[1]) + ' -n 300 -s ' + str(spec) + ' -d ' + str(den) + ' -b ' +str(mag)
                os.system(jobline)
                #for line in flines:
                    #fout.write(line.replace('__NAME__', jobname).replace('__JOBLINE__', jobline))
                #fout.close()
                #os.system('chmod +x ' + jobname+'.slurm')
                #os.system('sbatch -p icecube_guest '+ jobname+'.slurm')
                #raw_input('test')
