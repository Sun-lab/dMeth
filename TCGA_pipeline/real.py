#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 09:58:12 2017

@author: mac
"""

import os
import subprocess
#import numpy as np
st = '''#!/bin/bash
#SBATCH --job-name run{0} 
#SBATCH --mail-user=zhanyu@fredhutch.org
#SBATCH --mail-type=ALL
#SBATCH --nodes=1 #SBATCH --output=Rout/par-%J.out
#SBATCH --error=Rout/par-%J.err
#SBATCH --cpus-per-task=1
ml R/3.6.0-foss-2016b-fh1
R CMD BATCH run{0}.R Rout/run{0}.Rout
'''

Cancer = ["COAD","LUAD","LUSC","SKCM"]
for i in Cancer:
    with open('run'+i+'.sbatch','w') as fi:
       fi.write(st.format(i))
    #os.spawnlp(os.P_NOWAIT,'sbatch','sbatch -M beagle','batchsim'+str(i)+'_'+str(j)+'.sbatch')
    subprocess.call(["sbatch",'-M','beagle','run'+i+'.sbatch'])
