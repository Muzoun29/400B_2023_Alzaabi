

# Homework 6 Template
# G. Besla & R. Li
"""
HW6 

ASTR 400B 
Muzoun Alzaabi

This code will store the separation and relative velocities 
of the center of mass of the simulated MW, M33 and M31 
over the entire simulation and plot the corre- sponding orbits.


This code using Template resource :
 G. Besla & R. Li + 2023
 

Note : This is unfinished code 

"""



# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from ReadFile import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterOfMass2 import CenterOfMass




def OrbitCOM (galaxy, start, end, n):
    """
    

    function that loops over all the desired snapshots to compute the COM pos and vel as a function of time.
    inputs:
    ----------
    galaxy : str
    the name of the galaxy, e.g. “MW”
        
    start : int
    the number of the first snapshot to be read in.
        
    end : int
    the number of the last snapshot to be read in.
        
    n : int
    an integer indicating the intervals over which you will return the COM.
        

    Returns :
        orbit : Array that contains the time, position, and velocity of the COM.
        
    -------
    """


    
    # compose the filename for output
    
    fileout = 'Orbit_'+(galaxy)+'.txt'
    
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    if galaxy == "M33":  
        
        volDec = 4
        
    else:     
        volDec = 2
        
    delta = 0.1
    
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
    snap_id_seq = np.arange(start, end+1 , n )
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros([snap_id_seq.size, 7])
    
    # a for loop 
    for  i, snap_id in enumerate(snap_id_seq):# loop over files
        
        # compose the data filename (be careful about the folder)
        
        """Same as HW5 """
        
        ilbl = '000' + str(snap_id)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        # create filenames
        filename='%s_'%(galaxy)+ilbl+'.txt'
        
        
        # Initialize an instance of CenterOfMass class, using disk particles
        COM = CenterOfMass(filename , 2)
        # Store the COM pos and vel. Remember that now COM_P required VolDec
        COM_pos = COM.COM_P(delta , volDec)
        COM_vel = COM.COM_V(COM_pos[0], COM_pos[1], COM_pos[2])
    
        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        # Store results in orbit array

        # note that you can store 
        # a[i] = var1, *tuple(array1)
        
        orbit[i] = COM.time.value/1000, *tuple(COM_pos.value), *tuple(COM_vel.value)
        



        
        # print snap_id to see the progress
        #print(snap_id)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
    return orbit 


# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first! 




# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt




# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  




# Determine the magnitude of the relative position and velocities 

# of MW and M31

# of M33 and M31




# Plot the Orbit of the galaxies 
#################################




# Plot the orbital velocities of the galaxies 
#################################

