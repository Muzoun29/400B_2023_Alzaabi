"""
HW7 
ASTR400B 
Muzoun Alzaabi 

This code will predict the future trajectory of M33 in the frame of M31

Using resource template code by Rixin Li & G . Besla + 2023
+ HW 4 solutions provided by  Rixin Li & G . Besla + 2023


Note: This code include the HW Analysis at the end of code 



^Technical code issue UnitConversionError 
This result for undeo plot and output data results
"""


# # Homework 7 Template
# 
# Rixin Li & G . Besla
# 
# Make edits where instructed - look for "****", which indicates where you need to 
# add code. 




# import necessary modules
# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
# astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const
# import Latex module so we can display the results with symbols
from IPython.display import Latex

# **** import CenterOfMass to determine the COM pos/vel of M33
from CenterOfMass2 import CenterOfMass


# **** import the GalaxyMass to determine the mass of M31 for each component
from GalaxyMass import ComponentMass

#M33AnalyticOrbit



class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """
    
    def __init__(self, filename):
        """
        

        Initialize the class
        
        INPUT : 
        filename : string 
            file  store the integrated orbit 

        """
        # **** add inputs
    
        ### get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
        
        ### **** store the output file name
        self.filename = filename
        
        ### get the current pos/vel of M33 
        # **** create an instance of the  CenterOfMass class for M33 
        M33_COM = CenterOfMass('M33_000.txt',2)

        # **** store the position VECTOR of the M33 COM (.value to get rid of units)
        #Refered from HW4 G . Besla + 2023
        #vlDec=0.4
        M33_COM_p = M33_COM.COM_P(0.1,0.4)
        

        # **** store the velocity VECTOR of the M33 COM (.value to get rid of units)
        M33_COM_v = M33_COM.COM_V(M33_COM_p[0],M33_COM_p[1],M33_COM_p[2])
        
        ### get the current pos/vel of M31 
        # **** create an instance of the  CenterOfMass class for M31 
        M31_COM = CenterOfMass('M31_000.txt',2)
        # **** store the position VECTOR of the M31 COM (.value to get rid of units)
        #Refered from HW4 G . Besla + 2023
        M31_COM_p = M31_COM.COM_P(0.1,2)
        # **** store the velocity VECTOR of the M31 COM (.value to get rid of units)
        M31_COM_v = M31_COM.COM_V(M31_COM_p[0],M31_COM_p[1],M31_COM_p[2])
        
        ### store the DIFFERENCE between the vectors posM33 - posM31
        # **** create two VECTORs self.r0 and self.v0 and have them be the
        # relative position and velocity VECTORS of M33
        self.r0 = M33_COM_p - M31_COM_p
        self.v0 = M33_COM_v - M31_COM_v
        
        ### get the mass of each component in M31 
        ### disk
        # **** self.rdisk = scale length (no units)
        #HW instruction to set self.rdisk=5 kpc
        self.rdisk = 5

        # **** self.Mdisk set with ComponentMass function. 
        #Remember to *1e12 to get the right units. Use the right ptype
        self.Mdisk = ComponentMass('M31_000.txt',2)*1e12
        
        ### bulge
        #HW instruction : 
            #self.rbulge=1 kpc
        # **** self.rbulge = set scale length (no units)
        self.rbulge = 1

        # **** self.Mbulge  set with ComponentMass function. Remember to *1e12 to get the right units Use the right ptype
        self.Mbulge = ComponentMass('M31_000.txt',1)*1e12
        # Halo
        # **** self.rhalo = set scale length from HW5 (no units)
        self.rhalo = 62

        # **** self.Mhalo set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mhalo = ComponentMass('M31_000.txt',1)*1e12
    
    
    def HernquistAccel(self , M , ra , r ): 
        """
        

        This fuction will calculate the gravitational acceleration induced by a Hernquist profile
        ----------
        INPUT : 
            
        M : float
            Mass of the system
        ra : float
             Scale length 
        r : array
            Position vector
            
        OUTPUT : 
        Hern : returns the acceleration vector from a Hernquist potential
        -------

        """
        # it is easiest if you take as an input the position VECTOR 
        #Input as HW instruction to use 

        
        ### **** Store the magnitude of the position vector
        rmag = np.sqrt(r[0]**2+r[1]**2+r[2]**2)
        
        ### *** Store the Acceleration
        # use  -G*M/(rmag *(ra + rmag)**2) * r --> where the last r is a VECTOR 
        Hern =   -self.G * M/(rmag*(ra+rmag)**2)*r
        #follow the formula in the HW instructions
        # NOTE: we want an acceleration VECTOR so you need to make sure that in the Hernquist equation you 
        
        
        return Hern
    
    
    
    #rd is self.rdisk , 
    def MiyamotoNagaiAccel(self,M,rd,r):
        """
        

        This function will calculate the acceleration vector from a Miyamoto-Nagai profile
        ----------
        INPUT : 
            
        M : float
            Mass of the system
        rd : float
            Disck Scale length
        r : array
            Position vector
            
        OUTPUT : 
        Returns acceleration vector from a Miyamoto-Nagai profile
        

        """
        # it is easiest if you take as an input a position VECTOR  r 

        
        ### Acceleration **** follow the formula in the HW instructions
        # AGAIN note that we want a VECTOR to be returned  (see Hernquist instructions)
        # this can be tricky given that the z component is different than in the x or y directions. 
        # we can deal with this by multiplying the whole thing by an extra array that accounts for the 
        # differences in the z direction:
        #  multiply the whle thing by :   np.array([1,1,ZSTUFF]) 
        # where ZSTUFF are the terms associated with the z direction
        R = np.sqrt(r[0]**2+r[1]**2)
        zd = rd/5.0
        B = rd + np.sqrt(r[2]**2+zd**2)
        acce_V = (-self.G * M) / ((R**2+B**2)**1.5 * r)
        
        
        
        
       
        return acce_V * np.array([1,1,B/np.sqrt(r[2]**2+zd**2)])
        # the np.array allows for a different value for the z component of the acceleration
     
    
    def M31Accel(self,r):
        """
        

        Thiss function will sums all acceleration vectors from each galaxy component 
        ----------
        INPUT : 
        r : array
            Position vector
            
        OUTPUT : 
        Returns sum acceleration vector


        """
        # input should include the position vector, r


        ### Call the previous functions for the halo, bulge and disk
        # **** these functions will take as inputs variable we defined in the initialization of the class like 
        # self.rdisk etc. 
            
        
         
        halo_acc = self.HernquistAccel(self.Mhalo,self.rhalo,r)
        bulge_acc = self.HernquistAccel(self.Mbulge,self.rbulge,r)
        disk_acc = self.MiyamotoNagaiAccel(self.Mdisk,self.rdisk,r)
        
        # return the SUM of the output of the acceleration functions - this will return a VECTOR
        return halo_acc+bulge_acc+disk_acc
    
    
    
    def LeapFrog(self,dt,r,v):
        """
        

        Leap Frog integration scheme 
        
        ----------
        INPUT :
        dt : float
            time interval for integration
        r : array
            position vector for the M33 COM position relative to the M31
        v : array 
            velocity vector v for the M33 relative to M31
            
        OUTPUT : 
        Returns : new position and velocity vectors


        """
        # take as input r and v, which are VECTORS. Assume it is ONE vector at a time
                               
        
        # predict the position at the next half timestep
        rhalf =  r + v*dt/2.0
        
        # predict the final velocity at the next timestep using the acceleration field at the rhalf position 
        vnew = v + self.M31Accel(rhalf)*dt
        
        # predict the final position using the average of the current velocity and the final velocity
        # this accounts for the fact that we don't know how the speed changes from the current timestep to the 
        # next, so we approximate it using the average expected speed over the time interval dt. 
        rnew = rhalf + vnew * dt/2.0
        
        return rnew, vnew # **** return the new position and velcoity vectors
    
    
    
    def OrbitIntegration(self, t0, dt, tmax):
        """
        

        This function will loop over the LeapFrog integrator to solve the equations of motion and compute the future orbit of M33 for 10 Gyr into the future
        ----------
        INPUT :
        t0 : float
            Starting time
        dt : float
            Time interval
        tmax : float
            Final time
        
        OUTPUT :
            
        Returns
        -------
        None.

        """


        # initialize the time to the input starting time
        t = t0
        
        # initialize an empty array of size :  rows int(tmax/dt)+2  , columns 7
        orbit = np.zeros([int(tmax/dt)+2, 7])
        
        # initialize the first row of the orbit
        orbit[0] = t0, *tuple(self.r0), *tuple(self.v0)
        # this above is equivalent to 
        # orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]
        
        
        # initialize a counter for the orbit.  
        i = 1 # since we already set the 0th values, we start the counter at 1
        
        # start the integration (advancing in time steps and computing LeapFrog at each step)
        while (t <= tmax):  # as long as t has not exceeded the maximal time 
            
            # **** advance the time by one timestep, dt
            t = t + dt
            # **** store the new time in the first column of the ith row
            orbit[i,0] = t
            
            # ***** advance the position and velocity using the LeapFrog scheme
            # remember that LeapFrog returns a position vector and a velocity vector  
            # as an example, if a function returns three vectors you would call the function and store 
            # the variable like:     a,b,c = function(input)
            p,v = self.LeapFrog(dt,orbit[i,1:4],orbit[i,4:7])
            
         
    
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            # TIP:  if you want columns 5-7 of the Nth row of an array called A, you would write : 
            # A[n, 5:8] 
            # where the syntax is row n, start at column 5 and end BEFORE column 8
            orbit[i,1:4] = p
            
            
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            orbit[i,4:7] = v
            
            # **** update counter i , where i is keeping track of the number of rows (i.e. the number of time steps)
            i = i +1
        
        
        # write the data to a file
        np.savetxt(self.filename, orbit, fmt = "%11.3f"*7, comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
        
        # there is no return function




#Store the array into a file
M33AnalyticOrbit = M33AnalyticOrbit('M33_Analytic_Orbit.txt')
M33AnalyticOrbit.OrbitIntegration(0, 0.1, 10)

#From Hw6 solution Function that computes the magnitude of the difference: 
def relative_mag(a, b): 
    """
    Function that computes the magnitude of the difference between two vectors.
    Inputs with shape (3, n) will return n outputs

    PARAMETERS
    ----------
    a : `np.ndarray'
        first vector
    b : 'np.ndarray'
        second vector

    RETURNS
    -------
    mag : `float or np.ndarray`
        |a-b|
    """
    
    # compute the difference vector
    x = a[0] - b[0] 
    y = a[1] - b[1]
    z = a[2] - b[2]

    # return its magnitude
    return np.sqrt(x**2 + y**2 + z**2)

#Out-put data into an array to plot 
M31_orbit = np.genfromtxt('Orbit_M31.txt',dtype=None,names=True) 
M33_orbit = np.genfromtxt('Orbit_M33.txt',dtype=None,names=True) 
M33_output = np.genfromtxt('M33_Analytic_Orbit.txt',dtype=None,names=True)





#compute the seperation of position and velocity
sep_M31_M33, vel_M31_M33 = relative_mag(M31_orbit, M33_orbit)


pos = np.sqrt(M33_output['x']**2 + M33_output['y']**2 + M33_output['z']**2)
vel = np.sqrt(M33_output['vx']**2 + M33_output['vy']**2 + M33_output['vz']**2)




#plotting the difference in seperation

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111) 

#Plot the seperation 
ax.plot(M31_orbit['t'] , color = 'Black' ,linewidth=5, label='Seperation')
ax.plot(t,r,color = 'red' ,linewidth=5, label='Analytic Seperation')

plt.xlabel('Time (Gyr)', fontsize=22)
plt.ylabel('Separation (kpc)', fontsize=22)

#Legend
legend = ax.legend(loc='upper right',fontsize='x-large')

plt.savefig('Seperation.png')

#Plot the velocity
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111) 

#Plot the seperation 
ax.plot(M31_orbit['t'] , vel_M31_M33,color = 'Black' ,linewidth=5, label='Velocity Diffrence')
ax.plot(t,vel,color = 'red' ,linewidth=5, label='Analytic Velocity Diffrence ')

plt.xlabel('Time (Gyr)', fontsize=22)
plt.ylabel('Velocity (km/s)', fontsize=22)

#Legend
legend = ax.legend(loc='upper right',fontsize='x-large')

plt.savefig('Velocity.png')








#ANSWER TO Q1 : 
    
"""
In Q1 , 
In resource of using HW 6 in solutions By using a function called "relative_mag," 
give's the ability to get both the relative position 
and relative velocity for two galaxies throughout their entire orbit. 
This helps to compute the required separation between position and velocity.

"""

#ANSWER TO Q2 : 
    
"""Unfortunatily plot not resulted """
    
#ANSWER TO Q3 : What missing physics could make the difference?



#ANSWER TO Q4 : The MW is missing in these calculations. How might you include its effects?
"""
By Including MW seperation between on the galaxies as did in Analysis Q1 either to  M31 or M33
"""