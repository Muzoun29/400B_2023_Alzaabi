"""
Muzoun Alzaabi
ASTR400B 
Research Assignment 5 
April 18 , 2023 


Topic : Tidal Debris from M33 : Stellar Streams of M33
Q:
how/when do M33’s streams form? 
what is the Jacobi radius and how does this change in time?


This Python code that examines tidal debris and stellar stream formation in M33, 
focusing on the Jacobi radius and its change over time . 


Extention plot in analyzing mass change and its impact on tidal stream formation.
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
        
        #initial conditions for M33-MW and M31-M33 systems
        
        # COM Pos and Velocity for M31, MW, M33
        # For MW and M31 delta = 0.1 and VelDec = 2
        # For M33 delta = 0.1 and VelDec = 4
        
        # COM Pos and Velocity for M31
        com_p_M31 = com_M31.COM_P(volDec=2, delta=0.1)
        com_v_M31 = com_M31.COM_V(com_p_M31[0], com_p_M31[1], com_p_M31[2])
        
        # COM Pos and Velocity for M33
        com_p_M33 = com_M33.COM_P(volDec=4, delta=0.1)
        com_v_M33 = com_M33.COM_V(com_p_M33[0], com_p_M33[1], com_p_M33[2])
        
        # COM Pos and Velocity for MW
        com_p_MW = com_MW.COM_P(volDec=2, delta=0.1)
        com_v_MW = com_MW.COM_V(com_p_MW[0], com_p_MW[1], com_p_MW[2])
    
        # Set up initial conditions for M33-MW and M31-M33 systems
        r0_M33_MW = (com_p_M33 - com_p_MW)  # separation vector in kpc
        v0_M33_MW = (com_v_M33 - com_v_MW)  # relative velocity in km/s
        
        r0_M31_M33 = (com_p_M31 - com_p_M33) # separation vector in kpc
        v0_M31_M33 = (com_v_M31 - com_v_M33) # relative velocity in km/s
    
        # Define time array
        self.t = np.linspace(0, 6, 200) * u.Gyr
        
        #Using the ComponentMass function
        # Calculate the total mass of each galaxy
        M_M31 = ComponentMass('M31_000.txt', 2)
        M_M33 = ComponentMass('M33_000.txt', 2)
        M_MW = ComponentMass('MW_000.txt', 2)
    
        # Calculate the enclosed mass
        M_enclosed_M31_M33 = M_M31 + M_M33
        M_enclosed_MW_M33 = M_MW + M_M33
        
        # Calculate the Jacobi radius at each time step for M31-M33 and MW-M33 systems
        self.M31_M33_jacobi_radius = [self.calc_Jacobi_radius(r0_M31_M33, M_enclosed_M31_M33, M_M33) for _ in self.t]
        self.MW_M33_jacobi_radius = [self.calc_Jacobi_radius(r0_M33_MW, M_enclosed_MW_M33, M_M33) for _ in self.t]
    
    def calc_Jacobi_radius(self, r, M_host, M_satellite):
        """
        Calculate the Jacobi radius of the system at each time step.
        
    
        Parameters
        ----------
        r : array
            Separation vector between two galaxies in kpc
        M_host : float
            Mass of the host galaxy ( MW and M31 )
        M_satellite : float
            Mass of the satellite galaxy in (M33)
            
    
        Returns
        -------
        R_jacobi : array
            Jacobi radius of the system at each time step
        """
    
        #Msat is M33, Mhost is enclosed mass of MW and M31
        
        return r * (M_satellite / (2 * M_host))**(1/3)
    
    # Plot the Jacobi radius as a function of time for the M31-M33 and MW-M33 systems.    
    def plot_jacobi_radius(self):
        """
        Plot the Jacobi radius as a function of time for the M31-M33 and MW-M33 systems.
        Reason : To give us the insights into the behavior of tidal streams.
        understanding of the strength of tidal forces and the the impact on M33's structure.
        
        Note of Result: 
            Decrease : Indicates that the tidal forces are becoming stronger
                        tidal disruption of M33 increases
            Increase or Constant : 
                    Tidal forces are weaker
                    M33 will likely experience less tidal disruption.

        """
        plt.plot(self.t, self.M31_M33_jacobi_radius, label='M31-M33')
        plt.plot(self.t, self.MW_M33_jacobi_radius, label='MW-M33')
    
        plt.xlabel('Time (Gyr)')
        plt.ylabel('Jacobi Radius (kpc)')
        plt.title('Jacobi Radius vs. Time')
        
        plt.legend()
        plt.show()
        
        
        

        
#Create instance of a class        
tidal_streams = TidalStreams("TidalStreams.txt")
tidal_streams.plot_jacobi_radius() #calling plot_jacobi_radius() method 


#RESULT : 
    #PLOT 1 : observing how the Jacobi radius changes with time
#Result : The tidal forces are weaker, and M33 will likely experience less tidal disruption.



#2nd plot 


# plan is to investigate how the mass of M33 changes over time,
# providing insight into how the tidal streams form.

#How ? 
#implement the mass change calculation and 
#analyze the impact of mass change on the formation of tidal streams in M33.
   #Help to understand Answering Q1 in how and when M33's streams form

#Code initiate 

#def calculate_mass_change(self, snapshots):
#Calculate the mass change of M33 over time.
#This include an array of snapshot numbers for the simulation
def calculate_mass_change(self, snapshots):
    """
    Calculate the mass change of M33 over time.

    Parameters
    ----------
    snapshots : array
        An array of snapshot numbers for the simulation

    Returns
    -------
    mass_change : array
        Mass change of M33 at each snapshot
    """

    mass_change = []

    # Initial mass of M33
    initial_mass = ComponentMass('M33_000.txt', 2)

    for snapshot in snapshots:
        # Create a string of the form "M33_###.txt"
        filename = f'M33_{snapshot:03d}.txt'

        # Calculate the mass of M33 at the current snapshot
        current_mass = ComponentMass(filename, 2)

        # Calculate the mass change
        mass_diff = current_mass - initial_mass

        mass_change.append(mass_diff)

    return mass_change

def plot_mass_change(self, snapshots):
    """
    Plot the mass change of M33 as a function of time.

    Parameters
    ----------
    snapshots : array
        An array of snapshot numbers for the simulation
    """

    # Calculate the mass change of M33
    mass_change = self.calculate_mass_change(snapshots)

    # Calculate the time corresponding to each snapshot
    time = snapshots * 0.01 * u.Gyr

    # Plot the mass change as a function of time
    plt.plot(time, mass_change, label='M33 Mass Change')

    plt.xlabel('Time (Gyr)')
    plt.ylabel('Mass Change (1e10 M_sun)')
    plt.title('Mass Change vs. Time')

    plt.legend()
    plt.show()

#Create an instance of the TidalStreams class
tidal_streams = TidalStreams("TidalStreams.txt")

#Define the snapshots to analyze
snapshots = np.arange(0, 801, 10)


#Plot the mass change of M33 as a function of time
tidal_streams.plot_mass_change(snapshots)

#RESULT:
#PLOT 2: observing how the mass of M33 changes with time
#Result: The plot shows the mass loss of M33 over time. 
#The mass loss can be related to the formation of tidal streams.

    
          

    
    




        
        
        
        
 
        
        

        