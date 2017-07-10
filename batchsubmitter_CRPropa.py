from glob import glob
import os, sys
import numpy as np
import time
nside = str(32)

flines = open('submittemplate.sh').readlines()

jname = 'Crpropa'
#nucs = [[9,4], [11,5], [19, 9], [20, 10], [23, 11], [24, 12], [27,13], [28,14], [31, 15], [32, 16], [40,18], [40, 20], [52, 24]]

#nucs = [[35, 17],[40,18],[39, 19],[45,21],[48,22],[51,23],[55,25]]

nucs = [[1,1], [4, 2], [9,4],[12, 6], [14,7], [20, 10], [24, 12], [28,14], [32, 16],[40,18],[39, 19], [45,21],[51,23],[55,25], [56, 26]]

priornucs=[[1,1], [28,14], [56,26]]

for nuc in nucs:
    A = nuc[0]
    Z = nuc[1]
    n=0
    if nuc in priornucs:
        n=7
    print A,Z
    for seq in range(700,703+n):
        jobname = jname+str(seq)+str(A)+str(Z)+'job'
        fout = open(jobname+'.slurm', "w")
        jobline = ' python GalBack.py -n 10 -a '+str(A)+' -z '+str(Z)+' -f ' +str(seq) + ' -p'
        for line in flines:
            fout.write(line.replace('__NAME__', jobname).replace('__JOBLINE__', jobline))
        fout.close()
        os.system('chmod +x ' + jobname+'.slurm')
        os.system('sbatch -p icecube_guest '+ jobname+'.slurm')
        raw_input('test')


