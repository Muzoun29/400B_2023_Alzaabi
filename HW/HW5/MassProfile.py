
"""
HW5
ASTR 400B 

Muzoun Alzaabi
Files “MW 000.txt”, “M31 000.txt” , “M33 000.txt” and CenterofMass.py
Resource :  G Besla+2023

This code determine the mass distribution of each galaxy "MW , M31 , M33"
@SnapNumber 0 and use this to determine each galaxy’s rotation curve.

"""


# import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl
from ReadFile import Read
from CenterOfMass import CenterOfMass

import matplotlib #plot 
import matplotlib.pyplot as plt

from astropy.constants import G #Gravitational constant. from astropy.constants import G.



#P1 : Mass Profile

class MassProfile:
    """
    Class to calculate the Mass Profile with Galaxy 

    """

    def __init__(self,galaxy,snap):

       # From HW instructions 
       #Add a string of the filenumber to the value “000” ilbl = ‘000’ + str(Snap)
       ilbl = '000' + str(snap)
       #Remove all but the last 3 digits
       ilbl = ilbl[-3:]
       
       self.filename = "%s_"%(galaxy) + ilbl + '.txt' #Note:if need the file “MW 010.txt”, you would input “MW” and 10.
       #Read file 
       self.time,self.total,self.data = Read(self.filename)
       
       #Read in the data for the mass and x,y,z positions
       self.m = self.data['m']
       self.x = self.data['x']*u.kpc
       self.y = self.data['y']*u.kpc
       self.z = self.data['z']*u.kpc
       
       #Store the name of the galaxy as a global property self.gname
       self.gname = galaxy
       
#P2 : Mass Enclosed


    def MassEnclosed(self,ptype,radii ) :
        """
        This function will compute the mass enclosed within a given radius of the COM position 
        for a specified galaxy and a specified component of that galaxy. 
        

        Parameters
        ----------
        INPUT : 
            
        ptype : `int; 1, 2, or 3`
        particle type to use for COM calculations
        
        radii : array
        1D array of radii (kpc)
            
        OUTPUT : 
        Returns an array of masses (in units of M⊙)
        
        """
        
        #Create an index of particles
        index = np.where(self.data['type'] == ptype)

        #Determine the COM position
        COM = CenterOfMass(self.filename,2) #class to determine the COM position of a galaxy
        COM_position = COM.COM_P(1.0) #1.0 specifies that the position should be calculated using all particles within a radius of 1 unit 
        
        ##Define the position and mass using above index and store the components of each variable
        m_COM = self.m[index]
        x_COM = self.x[index] - COM_position[0]
        y_COM = self.y[index] - COM_position[1]
        z_COM = self.z[index] - COM_position[2]
        
        
        #Define the radius using the distance equation
        radius = np.sqrt(x_COM**2. + y_COM**2. + z_COM**2. )

        #Array to result  an array of masses
        mass = np.zeros(np.size(radii))
      
        #Loop over the Radius array
        for i in range(len(radii)):
            #New index that find the mass within the radius
            index_radd = np.where(radius <= radii[i]*u.kpc)
            mass[i] = np.sum(m_COM[index_radd])
            
        return mass*1e10*u.Msun #Return an array in unit of M_sun
        
        
        
#P3 : Total Mass Enclosed

    def MassEnclosedTotal(self,radii):
            """
            

            Parameters
            ----------
            INPUT : 
            radii : array
            1D array of radii (kpc)
            
            OUTPUT : 
            
            Returns an array of masses (in units of M⊙) 
            representing the total enclosed mass (bulge+disk+halo) at each radius of the input array.
            -------


            """
            
            #Since M33 Exeption case with no bulge only disk and halo :
                #Instruction : using self.gname
                #Using if statement 
                #Calculates the total enclosed mass (bulge+disk+halo)
            if (self.gname == 'M33'):
                M_Halo = self.MassEnclosed(1,radii)     #ptype = 1 halo
                M_Disk = self.MassEnclosed(2,radii)     #ptype = 2 disk
                M_Bulge = np.zeros(len(radii))  
                
                Mass_Total = M_Halo+M_Disk
                
            else: #MW , M31
                M_Halo = self.MassEnclosed(1,radii)     #ptype = 1 halo
                M_Disk = self.MassEnclosed(2,radii)     #ptype = 2 disk
                M_Bulge = self.MassEnclosed(3,radii)    #ptype = 3 bulge
       
                Mass_Total = M_Halo+M_Disk+M_Bulge
            
            return Mass_Total
       
#P4 : Hernquist Mass Profile
        
    def HernquistMass(self,radii,a,M_Halo):
           """
           This Function will compute the mass enclosed within a given radius using the theoretical profile:

           Parameters
           ----------
           INPUTS : 
           radii : array
            1D array of radii (kpc)
           a : float
               scale factor
           M_Halo : TYPE
               halo mass
              
           OUTPUTS :      
           Returns the halo mass in units of M⊙.
           -------
           Hernquist Mass = M_Halo*r**2 / (a+r)**2 
           """
           Hernequist_Mass= (M_Halo*radii**2)/(a+radii)**2
           
           return Hernequist_Mass*u.Msun
           
       
#P5 : Circular Velocity
        
    def CircularVelocity(self,ptype,radii):    
            """
            
            Parameters
            ----------
            INPUT : 
                
            ptype : `int; 1, 2, or 3`
            particle type 
            
            radii : array
                DESCRIPTION.
            
            OUTPUT : 
            Returns an array of circular speeds in units of km/s, rounded to two decimal places.

            -------
       

            """
            # Circular speed is computed using the Mass enclosed 
            # @ each radius, assuming spherical symmetry
            MassEnclosed = self.MassEnclosed(ptype,radii)
            G_const = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
            
            return np.round(np.sqrt(G_const*MassEnclosed/radii),2)*u.km/u.s # v = sqrt(GM/r) Rounded to two decimal places
        
#P6 : Total CircularVelocity

    def CircularVelocityTotal(self,radii):
            """
            
            Parameters
            ----------
            INPUT : 
                
            radii : array

            OUTPUT : 
                
            Returns an array of circular velocity (in units of km/s) 
            representing the total Vcirc created by all the galaxy components 
            (bulge+disk+halo) at each radius of the input array.
            -------

            """
            MassEnclosed_Tot = self.MassEnclosedTotal(radii)
            G_const = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
            
            
            return np.round(np.sqrt(G_const*MassEnclosed_Tot/radii),2)*u.km/u.s # v = sqrt(GM/r) Rounded to two decimal places

#P7 : Hernquist Circular Speed
        
    def HernquistVCirc(self,radii,M_Halo,a):
            
           """This function computes the circular speed using the Hernquist mass profile."""
           
           mass = self.HernquistMass(radii,M_Halo,a)
       
           G_const = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
             
            
           return np.round(np.sqrt((G_const*mass)/radii),2)*u.km/u.s
            
#P8 : Plot the Mass Profile for each Galaxy

#Create a plot of the mass profile (mass enclosed as a function of radius) 
#of each com- ponent of the MW to a radius of 30 kpc

#From HW instructions : 
    
MW  =  MassProfile("MW",0) # initialize the MassProfile class for MW
r = np.arange(0.25, 30.5, 1); print(r) # create an array of radii as the input



MW.MassEnclosed(1, r) ; print(MW.MassEnclosed(1, r))# get the enclosed halo masses at each element in 'r'

#Plot of the mass profile for MW  
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

#MW Halo result found in HW3
MW_Halo = 1.97e12 


# Mass profile for ptype 1 , 2 , 3 
plt.semilogy(r, MW.MassEnclosed(1,r), color='green', label='Halo Mass' )
plt.semilogy(r, MW.MassEnclosed(2,r), color='blue', label='Disk Mass')
plt.semilogy(r, MW.MassEnclosed(3,r), color='black', label='Bulge Mass')
plt.semilogy(r, MW.MassEnclosedTotal(r), color='red', label='Total Mass' )
plt.semilogy(r, MW.HernquistMass(r,55,MW_Halo),linestyle='dashed', color='purple', label='Hernquist Mass')
    
                 



       
# Add axis labels
plt.xlabel('Radius [kpc]', fontsize=22)
plt.ylabel('Enclosed Mass ', fontsize=22)

# Add title 
plt.title('Mass Profile for MW ',fontsize = 22)


#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')



plt.savefig('MassProfileforMWGalaxy.png')

#For M33 
#M33 Halo result found in HW3
M33_Halo = 0.187e12 

    
M33  =  MassProfile("M33",0) # initialize the MassProfile class for MW
r = np.arange(0.25, 30.5, 1.5); print(r) # create an array of radii as the input

M33.MassEnclosed(1, r) ; print(MW.MassEnclosed(1, r))# get the enclosed halo masses at each element in 'r'

#Plot of the mass profile for MW  
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)




# Mass profile for ptype 1 , 2 , 3 
plt.semilogy(r, M33.MassEnclosed(1,r), color='green', label='Halo Mass' )
plt.semilogy(r, M33.MassEnclosed(2,r), color='blue', label='Disk Mass')
#plt.semilogy(r, M33.MassEnclosed(3,r), color='black', label='Bulge Mass') #No Bulge
plt.semilogy(r, M33.MassEnclosedTotal(r), color='red', label='Total Mass' )
plt.semilogy(r, M33.HernquistMass(r,55,M33_Halo),linestyle='dashed',color='purple', label='Hernquist Mass')
    
                 



       
# Add axis labels
plt.xlabel('Radius [kpc]', fontsize=22)
plt.ylabel('Enclosed Mass ', fontsize=22)

# Add title 
plt.title('Mass Profile for M33 ',fontsize = 22)


#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')



plt.savefig('MassProfileforM33Galaxy.png')




#For M31
#M31 Halo result found in HW3
M31_Halo = 1.921e12 

    
M31  =  MassProfile("M31",0) # initialize the MassProfile class for MW
r = np.arange(0.25, 30.5, 1.5); print(r) # create an array of radii as the input

M31.MassEnclosed(1, r) ; print(MW.MassEnclosed(1, r))# get the enclosed halo masses at each element in 'r'

#Plot of the mass profile for MW  
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)




# Mass profile for ptype 1 , 2 , 3 
plt.semilogy(r, M31.MassEnclosed(1,r), color='green', label='Halo Mass' )
plt.semilogy(r, M31.MassEnclosed(2,r), color='blue', label='Disk Mass')
plt.semilogy(r, M31.MassEnclosed(3,r), color='black', label='Bulge Mass')
plt.semilogy(r, M31.MassEnclosedTotal(r), color='red', label='Total Mass' )
plt.semilogy(r, M31.HernquistMass(r,55,M31_Halo), color='purple',linestyle='dashed', label='Hernquist Mass')
    
                 



       
# Add axis labels
plt.xlabel('Radius [kpc]', fontsize=22)
plt.ylabel('Enclosed Mass ', fontsize=22)

# Add title 
plt.title('Mass Profile for M31 ',fontsize = 22)


#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')


plt.savefig('MassProfileforM31Galaxy.png')

#P9 : Plot the Rotation Curve for each Galaxy
#circular speed of each galaxy component to a radius of 30 kpc.
#Used same r value in above P8 

#For MW
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
    
plt.semilogy(r, MW.CircularVelocity(1,r), color='red', label='Halo Mass')
plt.semilogy(r, MW.CircularVelocity(2,r), color='blue', label='Disk Mass')
plt.semilogy(r, MW.CircularVelocity(3,r), color='green', label='Bulge Mass')
plt.semilogy(r, MW.CircularVelocityTotal(r), color='black', label='Total Mass') #Overplot the total circular speed for all galaxy component
plt.semilogy(r, MW.HernquistVCirc(r,55,MW_Halo), color='Purple',linestyle='dashed', label='Hernquist Mass')

plt.xlabel('Radius [kpc]', fontsize=22)
plt.ylabel('Circular speed [km/s]', fontsize=22)
# Add title 
plt.title('Circular Speen For MW ',fontsize = 22)
#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')


plt.savefig('CircularSpeenforMWGalaxy.png')

#For M33

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

    
plt.semilogy(r, M33.CircularVelocity(1,r), color='red', label='Halo Mass')
plt.semilogy(r, M33.CircularVelocity(2,r), color='blue', label='Disk Mass')
#plt.semilogy(r, M33.CircularVelocity(3,r), color='green', label='Bulge Mass') #NO Bulge
plt.semilogy(r, M33.CircularVelocityTotal(r), color='black', label='Total Mass')
plt.semilogy(r, M33.HernquistVCirc(r,23,M33_Halo), color='Purple',linestyle='dashed', label='Hernquist Mass')

plt.xlabel('Radius [kpc]', fontsize=22)
plt.ylabel('Circular speed [km/s]', fontsize=22)
# Add title 
plt.title('Circular Speed For M33 ',fontsize = 22)
#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')


plt.savefig('CircularSpeenforM33Galaxy.png')

#For M31

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

    
plt.semilogy(r, M31.CircularVelocity(1,r), color='red', label='Halo Mass')
plt.semilogy(r, M31.CircularVelocity(2,r), color='blue', label='Disk Mass')
plt.semilogy(r, M31.CircularVelocity(3,r), color='green', label='Bulge Mass') 
plt.semilogy(r, M31.CircularVelocityTotal(r), color='black', label='Total Mass')
plt.semilogy(r, M31.HernquistVCirc(r,45,M31_Halo), color='Purple',linestyle='dashed', label='Hernquist Mass')

plt.xlabel('Radius [kpc]', fontsize=22)
plt.ylabel('Circular speed [km/s]', fontsize=22)
# Add title 
plt.title('Circular Speed For M31 ',fontsize = 22)
#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')


plt.savefig('CircularSpeedforM31Galaxy.png')



       
