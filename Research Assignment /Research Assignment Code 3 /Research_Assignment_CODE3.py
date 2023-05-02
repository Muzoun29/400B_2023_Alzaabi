"""

HACK DAY 3/3 
May 02 , 23 


Topic : Tidal Debris from M33 : Stellar Streams of M33

Q:
how/when do M33’s streams form? 
what is the Jacobi radius and how does this change in time?

This Python code that examines tidal debris and stellar stream formation in M33, 
focusing on the Jacobi radius and its change over time . 
Extention plot in analyzing mass change and its impact on tidal stream formation.

Two plots are generated to show the Jacobi radius as a function of time 
for the M31-M33 and MW-M33 systems and the mass change of M33 over time . 
"""

# Import modules
import numpy as np
import astropy.units as u
import astropy.constants as const

# Import plotting modules
import matplotlib.pyplot as plt
import matplotlib

from ReadFile import Read
# This to define COM position and velocity of given galaxy
from CenterOfMass2 import CenterOfMass
from GalaxyMass import ComponentMass


# Define class "TidalStreams"
# which performs simulations of the merger
# and calculates the positions and properties of M33's tidal streams.
class TidalStreams:
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
        com_M31 = CenterOfMass('M31_000.txt', 2)  # COM object for M31
        com_M33 = CenterOfMass('M33_000.txt', 2)  # COM object for M33
        com_MW = CenterOfMass('MW_000.txt', 2)  # COM object for MW

        # initial conditions for M33-MW and M31-M33 systems

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

        # Set up initial conditions for
        r0_M33_MW = (com_p_M33 - com_p_MW)  # separation vector in kpc
        v0_M33_MW = (com_v_M33 - com_v_MW)  # relative velocity in km/s

        r0_M31_M33 = (com_p_M31 - com_p_M33)  # separation vector in kpc
        v0_M31_M33 = (com_v_M31 - com_v_M33)  # relative velocity in km/s

        # Define time array
        self.t = np.linspace(0, 6, 200) * u.Gyr

        # Using the ComponentMass function
        # Calculate the total mass of each galaxy
        M_M31 = ComponentMass('M31_000.txt', 2)
        M_M33 = ComponentMass('M33_000.txt', 2)
        M_MW  = ComponentMass('MW_000.txt', 2)

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

        # Msat is M33, Mhost is enclosed mass of MW and M31

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
        plt.plot(self.t.value, np.array([x.value for x in self.M31_M33_jacobi_radius]), label='M31-M33')
        plt.plot(self.t.value, np.array([x.value for x in self.MW_M33_jacobi_radius]), label='MW-M33')
    
        plt.xlabel('Time (Gyr)')
        plt.ylabel('Jacobi Radius (kpc)')
        plt.title('Jacobi Radius vs. Time')
        
        plt.savefig('plot_jacobi_radius.png')
        
        plt.legend()
        plt.show()

    # For the second plot
    # NEW FUNCTION METHOD ADD TO THE CLASS
    # 2nd Plot PURPOSE :
    #   To investigate how the mass of M33 changes over time,
    #   providing insight into how the tidal streams form.



#reads data from three simulation snapshots corresponding to MW, M31, and M33, 
#and calculates the total mass of M33 at each snapshot using the ComponentMass function. 

    def calculate_mass_change(self):
        """
        Calculate the mass change of M33 over time.
        This includes an array of snapshot numbers for the simulation.
        """

        mass_history = []
        time_history = []

        # Define the list of snapshot numbers: 0 for MW, 1 for M31, and 2 for M33
        snapshots = [1, 2, 3]

        for snapshot in snapshots:
            # Map snapshot numbers to the corresponding galaxy file names
            file_map = {1: 'MW_000.txt', 2: 'M31_000.txt', 3: 'M33_000.txt'}

            # Read data from snapshot file
            time, total, data = Read(file_map[snapshot])

            # Calculate the total mass of M33 at the current snapshot
            mass_M33 = ComponentMass(file_map[snapshot], 2)

            # Store mass and time in lists
            mass_history.append(mass_M33)
            time_history.append(time)

        return time_history, mass_history
    
    

    def plot_mass_change(self):

        #Plot the mass change of M33 over time.

        time_history, mass_history = self.calculate_mass_change()
    
        plt.plot([t.value for t in time_history], mass_history, label='M33 Mass Change')
    
        plt.xlabel('Time (Gyr)')
        plt.ylabel('Mass (1e10 M_sun)')
        plt.title('Mass Change of M33 vs. Time')
        plt.savefig('Mass_Change_M33_VS_Time.png')
        
        plt.legend()
        plt.show()
      
        
#Create a TidalStreams object
tidal_streams = TidalStreams("TidalStreams.txt")

#Plot the Jacobi radius as a function of time for the M31-M33 and MW-M33 systems
tidal_streams.plot_jacobi_radius()



#Plot the mass change of M33 over time
tidal_streams.plot_mass_change()        



#RESULT : 
#PLOT 1 : observing how the Jacobi radius changes with time
#Result : The tidal forces are weaker, and M33 will likely experience less tidal disruption.

#PLOT 2 : 
    # Decreasing ? 
#not taking into account the time evolution.
#calculating the mass of M33 at the initial snapshot (000) 
#for each of the three galaxies (MW, M31, and M33). 




#Ideas : 
#make a rough estimate based on the initial properties and the Jacobi radius.
#By creating 
#function to Estimate the time when M33's streams may start forming based 
#on the initial Jacobi radius.






