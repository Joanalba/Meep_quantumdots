#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 13:11:37 2023

This is the simulation for the setup varying the Z coordinate of the quantum 
dot from 0 to 200 nm. We obtain the FOM for every 0.02 nm. To get better
results we need a resolution of 50 in order to measure at different pixels.

@author: Joan
"""

import math
import meep as mp
import numpy as np
import matplotlib.pyplot as plt

resolution = 30 #number of pixels per unit of length of the cell
s = 10 #size of the cell
wt = 0.160 #thickness of the wvg
ww = 0.300 #wideness of the waveguide
cell_size = mp.Vector3(s,s,s)
dpml = 1 #Thickness of the absorptive layers
n = 10 #number of points at which we measure the FOM
time = 30 #time before taking a measurement
pml_layers = [mp.PML(dpml)] #defining the absorptive layers
qd1_wvlength = 0.93 #Wavelength of the ligth emmited by the quantum dot and the laser
freq = 1/qd1_wvlength #Frequency of the sources

def aver(lis):
    """
    
    Parameters
    ----------
    lis : List that we want to compute
        its average
        
    Returns
    -------
    The average of the list

    """
    return sum(lis)/len(lis)

# =============================================================================
# Step 1: adding a quantum dot inside the material to compute flux emited 
# =============================================================================

#Quantum Dot

def QD_no_wvg(detun,time):
    
    """
    
    Function for the first part of the FOM.
    
    Parameters
    ----------
    detun : determines the displacement of the quantum dot
        with respect the center of the cell in the Z direction
    time : until what time the simulation runs before computing
        the outcomes
        
    Returns
    -------
    The total power flux from the quantum dot that is captured
    by a cubic cell of size box_side around it

    """
    
    print('Doing QD no wvg')
    
    qd1_wvlength = 0.93 #Wavelength of the ligth emmited by the quantum dot
    freq = 1/qd1_wvlength
    
    center_qd = mp.Vector3(0,0,detun) #Point where the qd is placed
    
    #Symmetries to consider, for every symmetry the simulation is a factor or 2 faster
    symmetries = [mp.Mirror(direction=mp.X), mp.Mirror(direction=mp.Y)] 
    
    sources = [mp.Source(mp.ContinuousSource(frequency = freq,fwidth = 0.1*freq), 
                         component=mp.Ez, #component of the source
                         center=center_qd,
                         size=mp.Vector3())] #size of the source, we want a point for the qd
    
    sim = mp.Simulation(cell_size=cell_size,
                        boundary_layers=pml_layers,
                        sources=sources,
                        resolution=resolution,
                        progress_interval=99999, #if we increase it, the feedback is less
                        force_complex_fields=True,
                        symmetries=symmetries
                        )
    
    box_side = s-2*dpml
    
    sim.run(until=time)
    
    #Capture the flux in the following surfaces, one for each side of the box
    fluxbox1 = sim.flux_in_box(d = mp.Y, box = mp.Volume(center=mp.Vector3(0,box_side/2,0), size=mp.Vector3(box_side,0,box_side)))
    fluxbox2 = sim.flux_in_box(d = mp.Y, box = mp.Volume(center=mp.Vector3(0,-box_side/2,0), size=mp.Vector3(box_side,0,box_side)))
    fluxbox3 = sim.flux_in_box(d = mp.X, box = mp.Volume(center=mp.Vector3(box_side/2,0,0), size=mp.Vector3(0,box_side,box_side)))
    fluxbox4 = sim.flux_in_box(d = mp.X, box = mp.Volume(center=mp.Vector3(-box_side/2,0,0), size=mp.Vector3(0,box_side,box_side)))
    fluxbox5 = sim.flux_in_box(d = mp.Z, box = mp.Volume(center=mp.Vector3(0,0,box_side/2), size=mp.Vector3(box_side,box_side,0)))
    fluxbox6 = sim.flux_in_box(d = mp.Z, box = mp.Volume(center=mp.Vector3(0,0,-box_side/2), size=mp.Vector3(box_side,box_side,0)))
    total_fluxbox_small = fluxbox1-fluxbox2+fluxbox3-fluxbox4+fluxbox5-fluxbox6
        
    return total_fluxbox_small

# =============================================================================
# Step 2: adding a quantum dot and a wvg to compute beta1 and gamma_tot
# =============================================================================

def QD_wvg(detun,time):
    
    """
    
    Function for the second part of the FOM.
    
    Parameters
    ----------
    detun : determines the displacement of the quantum dot
        with respect the center of the cell in the Z direction
    time : until what time the simulation runs before computing
        the outcomes
        
    Returns
    -------
    The ratio of total power flux from the quantum dot that is captured
    at the ends of the waveguide with respect to the power captured by a 
    cubic cell of size box_side around it

    """
    
    print('Doing QD wvg')
    
    #We introduce the Waveguide
    waveguide_material = mp.Medium(epsilon=12)
    waveguide_size = mp.Vector3(s,wt,ww)
    waveguide = mp.Block(waveguide_size,material=waveguide_material)
    
    #Quantum Dot
    qd1_wvlength = 0.93 #Wavelength of the ligth emmited by the quantum dot
    freq = 1/qd1_wvlength
    
    center_qd = mp.Vector3(0,0,detun) #Point where the qd is placed
    
    symmetries = [mp.Mirror(direction=mp.X), mp.Mirror(direction=mp.Y)]
    
    sources = [mp.Source(mp.ContinuousSource(frequency = freq), 
                         component=mp.Ez,
                         center=center_qd,
                         size=mp.Vector3())]
    
    sim = mp.Simulation(cell_size=cell_size,
                        boundary_layers=pml_layers,
                        geometry=[waveguide],
                        sources=sources,
                        resolution=resolution,
                        progress_interval=99999,
                        force_complex_fields=True,
                        symmetries=symmetries
                        )
    
    box_side = s-2*dpml
    
    sim.run(until=time)
    
    fluxbox1 = sim.flux_in_box(d = mp.Y, box = mp.Volume(center=mp.Vector3(0,box_side/2,0), size=mp.Vector3(box_side,0,box_side)))
    fluxbox2 = sim.flux_in_box(d = mp.Y, box = mp.Volume(center=mp.Vector3(0,-box_side/2,0), size=mp.Vector3(box_side,0,box_side)))
    fluxbox3 = sim.flux_in_box(d = mp.X, box = mp.Volume(center=mp.Vector3(box_side/2,0,0), size=mp.Vector3(0,box_side,box_side)))
    fluxbox4 = sim.flux_in_box(d = mp.X, box = mp.Volume(center=mp.Vector3(-box_side/2,0,0), size=mp.Vector3(0,box_side,box_side)))
    fluxbox5 = sim.flux_in_box(d = mp.Z, box = mp.Volume(center=mp.Vector3(0,0,box_side/2), size=mp.Vector3(box_side,box_side,0)))
    fluxbox6 = sim.flux_in_box(d = mp.Z, box = mp.Volume(center=mp.Vector3(0,0,-box_side/2), size=mp.Vector3(box_side,box_side,0)))
    fluxwvg1 = sim.flux_in_box(d = mp.X, box = mp.Volume(center=mp.Vector3(-s/2+dpml,0,0), size=mp.Vector3(0,wt,ww)))
    fluxwvg2 = sim.flux_in_box(d = mp.X, box = mp.Volume(center=mp.Vector3(s/2-dpml,0,0), size=mp.Vector3(0,wt,ww)))
    #Flux for the surfaces of the cell
    total_fluxbox_big = fluxbox1-fluxbox2+fluxbox3-fluxbox4+fluxbox5-fluxbox6 
    #Flux at the ends of the waveguide
    fluxoutwvg = fluxwvg2-fluxwvg1
    beta1 = fluxoutwvg/total_fluxbox_big
        
    return beta1

# =============================================================================
# Step 3: Reseting meep and computing the variables for the laser without wvg
# =============================================================================

def LaserBeam_no_wvg(detun,time,angle):
    
    """
    
    Function for the third part of the FOM. It computes the total power emitted
    by the laser.
    
    Parameters
    ----------
    detun : determines the displacement of the quantum dot
        with respect the center of the cell in the Z direction
    time : until what time the simulation runs before computing
        the outcomes
    angle : The angle that the laser beam has with respect to the 
        Y axis
        
    Returns
    -------
    The total power flux from the laser.

    """
    
    print('Doing laser no wvg')
    
    #LaserBeam
    laser1_wvlength = 0.93
    fcen = 1/laser1_wvlength
    laser_center = -s/2+dpml+0.2 #Distance in the Y direction of the laser center
    laser_size = s-2*dpml #Size of the laser center
    w0 = 1.0 #Size of the laser waist radius, its profile at the focus point
    rot_angle = angle# CCW rotation angle about z axis (0: +y axis) in degrees
    dx = math.tan(rot_angle*math.pi/180)*laser_center #Position of the laser center for every angle
    beam_x0 = mp.Vector3(dx,-laser_center,0) #beam focus
    beam_kdir = mp.Vector3(0,1,0).rotate(mp.Vector3(0, 0, 1), math.radians(rot_angle)) #beam propagating direction
    beam_w0 = w0#beam waist radius
    beam_E0 = mp.Vector3(0, 0, 1) #Polarization
    
    #Defining symmetries depending on the angle
    if angle == 0:
        symmetries = [mp.Mirror(direction=mp.X), mp.Mirror(direction=mp.Z,phase=-1.0)]
    else:
        symmetries = [mp.Mirror(direction=mp.Z,phase=-1.0)]
    
    sources = [mp.GaussianBeamSource(src=mp.ContinuousSource(fcen), #Source as a laser
                                     center = mp.Vector3(-dx,laser_center,0),
                                     size = mp.Vector3(laser_size,0,laser_size),
                                     beam_x0=beam_x0,
                                     beam_kdir=beam_kdir,
                                     beam_w0=beam_w0,
                                     beam_E0=beam_E0)]
    
    #SimulationParameters
    sim = mp.Simulation(cell_size=cell_size,
                        boundary_layers=pml_layers,
                        sources=sources,
                        resolution=resolution,
                        progress_interval=9999,
                        force_complex_fields=True,
                        symmetries=symmetries
                        )
    
    box_side = s-2*dpml
    
    sim.run(until=time)
    
    #We compute the flux for different surfaces to average it and 
    # get a better approximation of the laser power
    surface= np.linspace(0,box_side/2,n)
    flux_surface = np.zeros(len(surface))
    
    for j in range(len(surface)):
        fluxbox1 = sim.flux_in_box(d = mp.Y, box = mp.Volume(center=mp.Vector3(0,surface[j],0), size=mp.Vector3(box_side,0,box_side)))
        flux_surface[j] = fluxbox1
    
    average = aver(flux_surface)
    
    return average

# =============================================================================
# Step 4: Reseting meep and computing the variables for the laser with wvg
# =============================================================================

def LaserBeam_wvg(detun,time,angle):
    
    """
    
    Parameters
    ----------
    detun : determines the displacement of the quantum dot
        with respect the center of the cell in the Z direction
    time : until what time the simulation runs before computing
        the outcomes
    angle : The angle that the laser beam has with respect to the 
        Y axis
    
    Returns
    -------
    The total power flux from the laser captured at the ends of the waveguide
    and the electric field at the quantum dot spot.

    """
    
    print('Doing laser wvg')
    
    center_qd = mp.Vector3(0,0,detun) #Point where the qd is placed
    
    #Waveguide
    waveguide_material = mp.Medium(epsilon=12)
    waveguide_size = mp.Vector3(s,wt,ww)
    waveguide = mp.Block(waveguide_size,material=waveguide_material)
    
    #LaserBeam
    laser1_wvlength = 0.93
    fcen = 1/laser1_wvlength
    laser_center = -s/2+dpml+0.2 #Distance in the Y direction of the laser center
    laser_size = s-2*dpml #Size of the laser center
    w0 = 1.0 #Size of the laser waist radius, its profile at the focus point
    rot_angle = angle# CCW rotation angle about z axis (0: +y axis) in degrees
    dx = math.tan(rot_angle*math.pi/180)*laser_center #Position of the laser center for every angle
    beam_x0 = mp.Vector3(dx,-laser_center,0) #beam focus
    beam_kdir = mp.Vector3(0,1,0).rotate(mp.Vector3(0, 0, 1), math.radians(rot_angle)) #beam propagating direction
    beam_w0 = w0#beam waist radius
    beam_E0 = mp.Vector3(0, 0, 1) #Polarization
    
    #Defining symmetries depending on the angle
    if angle == 0:
        symmetries = [mp.Mirror(direction=mp.X), mp.Mirror(direction=mp.Z,phase=-1.0)]
    else:
        symmetries = [mp.Mirror(direction=mp.Z,phase=-1.0)]
    
    sources = [mp.GaussianBeamSource(src=mp.ContinuousSource(fcen), #Source as a laser
                                     center = mp.Vector3(-dx,laser_center,0),
                                     size = mp.Vector3(laser_size,0,laser_size),
                                     beam_x0=beam_x0,
                                     beam_kdir=beam_kdir,
                                     beam_w0=beam_w0,
                                     beam_E0=beam_E0)]
    
    #SimulationParameters
    sim = mp.Simulation(cell_size=cell_size,
                        boundary_layers=pml_layers,
                        geometry=[waveguide],
                        sources=sources,
                        resolution=resolution,
                        progress_interval=9999,
                        force_complex_fields=True,
                        symmetries=symmetries
                        )
    
    box_side = s-2*dpml
    
    sim.run(until=time)
    
    #Electric field at a specific point
    ez_point1 = sim.get_array(center=center_qd, size=mp.Vector3(), component = mp.Ez)
    ez_point = abs(ez_point1)**2
    
    #Flux in the ends of the waveguide
    fluxwvg1 = sim.flux_in_box(d = mp.X, box = mp.Volume(center=mp.Vector3(-box_side/2,0,0), size=mp.Vector3(0,wt,ww)))
    fluxwvg2 = sim.flux_in_box(d = mp.X, box = mp.Volume(center=mp.Vector3(box_side/2,0,0), size=mp.Vector3(0,wt,ww)))
    fluxoutwvg = fluxwvg2-fluxwvg1
    
    #We ask for the electric field already squared and the flux comming from the ends
    return ez_point, fluxoutwvg

#Representation of functions
    
dipole_det= np.linspace(0,0.200,n) #Positions for the quantum dot
laser_ang = 0 #Incident angle for the laser

#Different variables
eta = np.zeros(len(dipole_det))
gamma_tot = np.zeros(len(dipole_det))
beta1 = np.zeros(len(dipole_det))
beta2 = np.zeros(len(dipole_det))
FOM = np.zeros(len(dipole_det))
beta2_mod = np.zeros(len(dipole_det))
FOM_mod = np.zeros(len(dipole_det))
list_flux_laser_wvg = np.zeros(len(dipole_det))
laser_no_wvg = LaserBeam_no_wvg(0,time,laser_ang)

for j in range(len(dipole_det)):
    gamma_tot[j] = QD_no_wvg(dipole_det[j],time)/(2*math.pi*freq)
    beta1[j] = QD_wvg(dipole_det[j],time)
    laser_wvg = LaserBeam_wvg(dipole_det[j],time,laser_ang)
    eta[j] = laser_wvg[1]/laser_no_wvg
    list_flux_laser_wvg[j] = laser_wvg[1]
    beta2[j] = (laser_wvg[0]*(2*math.pi))/(gamma_tot[j])
    if j != 0:
        beta2_mod[j] = (laser_wvg[0]*(2*math.pi))/(gamma_tot[j]*beta2_mod[0])
    else:
        beta2_mod[j] = (laser_wvg[0]*(2*math.pi))/(gamma_tot[j])
        
    FOM[j] = eta[j]/(beta1[j]*beta2[j])
    FOM_mod[j] = eta[j]/(beta1[j]*beta2_mod[j])
    print('')
    print('###')
    print(j)
    print('###')
    print('')
    
FOM_mod[0] = FOM_mod[0]*beta2_mod[0]
beta2_mod[0] = 1

#Plotting
plt.figure()
plt.plot(dipole_det, eta, 'go-')
plt.xlabel('Displacement of the QD in Z (um)')
plt.ylabel('Eta')
plt.savefig('Eta')

plt.figure()
plt.plot(dipole_det, gamma_tot, 'bo-')
plt.xlabel('Displacement of the QD in Z (um)')
plt.ylabel('Gamma_tot')
plt.savefig('Beta2')

plt.figure()
plt.plot(dipole_det, beta1, 'co-')
plt.xlabel('Displacement of the QD in Z (um)')
plt.ylabel('Beta1')
plt.savefig('Beta2')

plt.figure()
plt.plot(dipole_det, beta2, 'ro-')
plt.xlabel('Displacement of the QD in Z (um)')
plt.ylabel('Beta2')
plt.savefig('Beta2')

plt.figure()
plt.plot(dipole_det, FOM, 'mo-')
plt.xlabel('Displacement of the QD in Z (um)')
plt.ylabel('FOM')
plt.savefig('FOM')

plt.figure()
plt.plot(dipole_det, beta2_mod, 'ro-')
plt.xlabel('Displacement of the QD in Z (um)')
plt.ylabel('Beta2')
plt.savefig('Beta2_mod')

plt.figure()
plt.plot(dipole_det, FOM_mod, 'mo-')
plt.xlabel('Displacement of the QD in Z (um)')
plt.ylabel('FOM')
plt.savefig('FOM_mod')


with open('FOM.txt','w') as f:
    for j in range(len(dipole_det)):
        f.write(str(dipole_det[j]))
        f.write('\t')
        f.write(str(FOM[j]))
        f.write('\n')
        
f.close()