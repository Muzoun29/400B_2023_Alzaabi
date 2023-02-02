#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HW2
ParticleProperties
Muzoun Alzaabi

This python code will create a new function called ParticleInfo
, which takes as inputs: filename, particle type and particle number.
"""

# Import Modules
import numpy as np 
from astropy import units as u
from ReadFile import Read

filename = "MW_000.txt"

#Define a function takes in a filename, particle type, and particle number as input
def ParticleInfo(filename, particle_type, particle_number):   
    data = Read(filename)
    
    # Get the index of the particle with the specified particle type and number
    index = np.where((data['type'] == particle_type) & (data['m, x, y, z, vx, vy, vz'] == particle_number))[3]
    #magnitude of the distance in kpc
    distance = np.around(np.sqrt(data['x'][index]**2 + data['y'][index]**2 + data['z'][index]**2),3) * u.kpc
   
    #velocit in km/s , rounded into 3 decimal places 
    velocity = np.around(np.sqrt(data['vx'][index]**2 + data['vy'][index]**2 + data['vz'][index]**2),3) * u.km/u.s
    
    #mass in of M⊙
    mass = data['m'][index] * 1e10 * u.Msun
    
    #Convert the Dist. of particle to Ly and rounded in 3 decimal places
    distance_ly = np.around(distance.to(u.lyr),3)
    
    return distance_ly, velocity, mass

distance_ly, velocity, mass = ParticleInfo(filename, 2, 100)

print("Distance: ",distance_ly)
print("Velocity: ",velocity)
print("Mass: ",mass)


