from glob import glob
import os, sys
import numpy as np
import time
nside = str(32)

flines = open('submittemplate.sh').readlines()

jname = 'Crpropa'

for seq in range(600,610):
    jobname = jname+str(seq)+'job'
    fout = open(jobname+'.slurm', "w")
    jobline = ' python GalBack.py -n 10 -f ' +str(seq)
    for line in flines:
        fout.write(line.replace('__NAME__', jobname).replace('__JOBLINE__', jobline))
    fout.close()
    os.system('chmod +x ' + jobname+'.slurm')
    os.system('sbatch -p icecube '+ jobname+'.slurm')
    raw_input('test')


