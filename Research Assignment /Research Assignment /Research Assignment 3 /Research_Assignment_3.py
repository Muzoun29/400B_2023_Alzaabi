"""


Research Assignment 3 
Mar 30 , 2023 

*This is initial code set 

Topic : Tidal Debris from M33 : Stellar Streams of M33


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







# Define class "TidalStreams" 
#which performs simulations of the merger 
#and calculates the positions and properties of M33's tidal streams.

class TidalStreams : 
    """This class to perform simulations of the merger """
    
    def __init__(self, filename):
        """
        Initialize the TidalStreams class.
        ----------
        filename : `str`
            Name of the file in which to store the orbit


        """
        # Constants and units    
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
#CONFUSE 
 
        # Set up initial conditions for 
        r0_M33_MW = (com_p_M33 - com_p_MW).value  # separation vector in kpc
        v0_M33_MW = (com_v_M33 - com_v_MW).value  # relative velocity in km/s
        
        r0_M31_M33 = (com_p_M31 - com_p_M33).value  # separation vector in kpc
        v0_M31_M33 = (com_v_M31 - com_v_M33).value  # relative velocity in km/s

        # Define time array
        self.t = np.linspace("***")* u.Gyr

            
        
            
        
        # Need to Define Function to Calculate total Mass of each Galaxy 
        
        
        
        
        
        
        # Need Function to calculate the enclosed mass 
        
        
        
        
        
        
# Recall  Rj = r * (Msat / (2 * Mhost))**(1/3)
# Where Msat is M33 and Mhost is the enclosed mass of MW and M31


#Jacobi radius imp. to identify the tidal debris from M33 
#that has been stripped by the gravitational interaction with MW,M33.
    def calc_Jacobi_radius(self, r_M31_M33, r_MW_M33):"""
    Calculate the Jacobi radius of M31-M33 and MW-M33 system at each time step.

"""
 

    # Calculate the Jacobi radius of M31-M33 and MW-M33 system at each time step
        #M31_M33 =
        #MW_M33 = 
        





        
        
        
        
 
        
        

        