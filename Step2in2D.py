#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 13:11:37 2023

@author: Joan
"""

import math
import meep as mp
import numpy as np
import matplotlib.pyplot as plt

resolution = 10
s = 10 #size of the cell
w = 2 #width of the wvg
cell_size = mp.Vector3(s,s,0)
dpml = 1
pml_layers = [mp.PML(dpml)]

#Waveguide
waveguide_material = mp.Medium(epsilon=12)
waveguide_size = mp.Vector3(s,w,0)
waveguide = mp.Block(waveguide_size,material=waveguide_material)

qd1_wvlength = 0.93
freq = 1/qd1_wvlength

sources = [mp.Source(mp.ContinuousSource(frequency = freq),
                     component=mp.Ez,
                     center=mp.Vector3(0,0,0),
                     size=mp.Vector3(0,0,0))]

sim = mp.Simulation(cell_size=cell_size,
                    boundary_layers=pml_layers,
                    geometry=[waveguide],
                    sources=sources,
                    resolution=resolution)

direction1 = mp.X
direction2 = mp.X
# =============================================================================
# f0= mp.FluxRegion(center=mp.Vector3(0,0,0), size=mp.Vector3(s,s,s), 
#                    direction=direction0
#                   )
# =============================================================================
# =============================================================================
#                   
# f1= mp.FluxRegion(center=mp.Vector3(s/2-dpml-0.1,0,0), size=mp.Vector3(0,s,s), 
#                    direction=direction1
#                   )
# f2= mp.FluxRegion(center=mp.Vector3(-s/2+dpml+0.1,0,0), size=mp.Vector3(0,s,s), 
#                    direction=direction2
#                   )
# f3= mp.FluxRegion(center=mp.Vector3(0,s/2-dpml-0.1,0), size=mp.Vector3(s,0,s), 
#                   # direction=direction2
#                   )
# f4= mp.FluxRegion(center=mp.Vector3(0,-s/2+dpml+0.1,0), size=mp.Vector3(s,0,s), 
#                   # direction=direction2
#                   )
# f5= mp.FluxRegion(center=mp.Vector3(0,0,s/2-dpml-0.1), size=mp.Vector3(s,s,0), 
#                   # direction=direction2
#                   )
# f6= mp.FluxRegion(center=mp.Vector3(0,0,-s/2+dpml+0.1), size=mp.Vector3(s,s,0), 
#                   # direction=direction2
#                   )
# 
# # f0_mon = sim.add_flux(freq,0,1,f0)
# f1_mon = sim.add_flux(freq,0,1,f1)
# f2_mon = sim.add_flux(freq,0,1,f2)
# f3_mon = sim.add_flux(freq,0,1,f3)
# f4_mon = sim.add_flux(freq,0,1,f4)
# f5_mon = sim.add_flux(freq,0,1,f5)
# f6_mon = sim.add_flux(freq,0,1,f6)
# 
# 
# =============================================================================
# sim.run(mp.at_beginning(mp.output_epsilon),mp.to_appended("ez",mp.at_every(0.6, mp.output_efield_z)),until=20)
# sim.run(mp.at_every(1/(freq*3),mp.output_efield_z), until=20)

sim.run(until=20)


eps_data = sim.get_array(center=mp.Vector3(), size=cell_size, component=mp.Dielectric)
ez_data = sim.get_array(center=mp.Vector3(), size=cell_size, component=mp.Ez)
plt.figure()
plt.imshow(eps_data.transpose(),interpolation='spline36',cmap='binary')
plt.imshow(ez_data.transpose(),interpolation='spline36',cmap='RdBu',alpha=0.9)
plt.axis('off')
plt.show()
# =============================================================================
# # f0_flux = mp.get_fluxes(f0_mon)[0]
# f1_flux = mp.get_fluxes(f1_mon)[0]
# f2_flux = mp.get_fluxes(f2_mon)[0]
# f3_flux = mp.get_fluxes(f3_mon)[0]
# f4_flux = mp.get_fluxes(f4_mon)[0]
# f5_flux = mp.get_fluxes(f5_mon)[0]
# f6_flux = mp.get_fluxes(f6_mon)[0]
# # print('Flux0:', f0_flux)
# print('Flux1:', f1_flux)
# print('Flux2:', f2_flux)
# print('Flux3:', f3_flux)
# print('Flux4:', f4_flux)
# print('Flux5:', f5_flux)
# print('Flux6:', f6_flux)
# =============================================================================

# sim.run(mp.at_beginning(mp.output_epsilon),mp.to_appended("ez",mp.at_every(0.6, mp.output_efield_z)),until=200)