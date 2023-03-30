"""
Muzoun Alzaabi 
ASTR400B

Research Assignment 3 
Mar 30 , 2023 



*This is initial code set 

Topic : Tidal Debris from M33 : Stellar Streams of M33


Question : 
    how/when do M33’s streams form? 
    what is the Jacobi radius and how does this change in time?
    
    
This code will analyze the formation of M33's stellar streams 
and the change in the Jacobi radius over time

Using observational data for M33's Center of Mass position and velocity
and then calculate the Jacobi radius at some points in time.

"""



# import modules
import numpy as np
import astropy.units as u
import astropy.constants as const


# import plotting modules
import matplotlib.pyplot as plt
import matplotlib


from ReadFile import Read
#This to define COM position anf velocity of given galaxy
from CenterOfMass2 import CenterOfMass
from GalaxyMass import ComponentMass




#Resource : G. Besla & H. Foote + 2023   #    HW7 & Lab 3 


class TidalStreams : 
    """This class to perform simulations of the merger and 
    calculate the positions and properties of M33's tidal streams"""
    
    def __init__(self, filename):
        """
        Parameters
        ----------
        filename : `str`
            Name of the file in which to store the orbit


        """
        # get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
        
        # store the output file name
        self.filename = filename
        
        
        # Calculate the COM Position of MW, M31 and M33 
        com_M31 = CenterOfMass('M31_000.txt', 2) # COM object for M31
        com_M33 = CenterOfMass('M33_000.txt', 2) # COM object for M33
        com_MW = CenterOfMass('MW_000.txt', 2) # COM object for MW
        
        # COM Pos and Velocity for M31 , MW , M33
        # For MW and M31 delta = 0.1 and VelDec = 2
        # For M33 dalta = 0.1 and VelDec = 4
        
        # COM Pos and Velocity for M31
        com_p_M31 = com_M31.COM_P(volDec=2, delta=0.1)
        com_v_M31 = com_M31.COM_V(com_p_M31[0], com_p_M31[1], com_p_M31[2])
        # COM Pos and Velocity for for M33
        com_p_M33 = com_M33.COM_P(volDec=4, delta=0.1)
        com_v_M33 = com_M33.COM_V(com_p_M33[0], com_p_M33[1], com_p_M33[2])
        # COM Pos and Velocity for for MW
        com_p_MW = com_MW.COM_P(volDec=2, delta=0.1)
        com_v_MW = com_MW.COM_V(com_p_MW[0], com_p_MW[1], com_p_MW[2])
        

#INCLUDE MW AS Stellar Streams of M33 (during and after the MW-M31 merger)
#Mtwo systems M31-M33 and MW-M33 system
# set up initial conditions (MW pos and velocity relative to M31)
        #self.r0 = (com_p_MW - com_p_M31).value # separation vector in kpc
        #self.v0 = (com_v_MW - com_v_M31).value # relative velocity in km/s """
        
        
        
        
        
        # Need to Define Function to Calculate total Mass of each Galaxy 
        
        
        
        
        
        
        # Need Function to calculate the enclosed mass 
        
        
        
        
        
        
        # Define function to Calculate the Jacobi Radius of M31-M33 and MW-M33 system
        #Recall  # Rj=r*(Msat/2/Mhost)**(1/3)
                    #Where Msat is M33
                    #Mhost enclose mass of MW and M31




        
        
        
        
 
        
        

        