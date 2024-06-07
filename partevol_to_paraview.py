# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 15:44:11 2023

@author: cvh
"""

# import numpy as np

import sys
sys.path.insert(1, '/Users/cvh/.ssh/Documents/Research/Python_functions_lib/')
from readRW3D import read_partevol
import hydroinform.DFS as DFS


pathrw = r'C:\Users\cvh\.ssh\Documents\Research\Nret24\benchmark_models\SkernUps_1layer\q_corrected'
pathout = r'C:\Users\cvh\.ssh\Documents\Research\Nret24\benchmark_models\SkernUps_1layer\q_corrected\paraview_'
pathmike = r'C:\Users\cvh\.ssh\Documents\Research\Nret24\MIKE-SHE_outputs\Benchmark1'
dfs3_flow = f'{pathmike}\\SkjernUps_Sim1_LowLevel_Chris_Benchmark1_3DSZflow.dfs3'

#get axis origins
data = DFS.DFS3.from_file(dfs3_flow) 
x0 = data.x_origin
y0 = data.y_origin

file = 'rw3d_part_evol.dat'
filename = f'{pathrw}\{file}'
fileout = 'particles_all_dt'

plume = read_partevol(filename)

for it in range(len(plume.times)):
    plume_it = plume.data[it]
    plume_it.to_csv(f'{pathout}\{fileout}{it}.csv', columns=['xp', 'yp', 'zp'],index=False)
