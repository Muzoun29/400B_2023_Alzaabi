"""
ASTR400B
HW3 
Muzoun Alzaabi 

1.File Data in this program are stored in the file in astr400b/ directory on nimoy using sftp .
These files are for this assignment: MW 000.txt, M31 000.txt, M33 000.txt

2.This python code under the name of "GalaxyMass" will 
            1.Return the total mass of any desired galaxy component .
            2.Compute the total mass of each galaxy (three galaxies (MW, M31, M33))
            3.Compute the total mass of the Local Group
            4.Compute the baryon fraction fbar for each galaxy and the whole Local Group.
            
Recource : ReadFile.py , MW 000.txt, M31 000.txt, M33 000.txt 


"""

# Import Modules 
import numpy as np # import numpy
import astropy.units as u # import astropy units
from ReadFile import Read #import read file 

#Now 1. Write a function called ComponentMass This function aims to return the total mass  M_tot
#Define the function that will take as input: a filename, the particle type.
def ComponentMass (filename, particletype) :
    """
    This function created to find the total mass of the galaxy
    INPUT :
        filename: str
        particletype : float
            Halo type (1), Disk type (2), Bulge type (3)
    OUTPUT : 
        The total mass , 
        units of solar masses ,rounded to three decimal places

    """
    time,total,data = Read(filename) #Those data are what returned in the ReadFile / 
    #The ReadFile has output of  time:Time of snapshot in Myr , total: Total number of particles ,data : array with the particle data
    index = np.where(data['type']==particletype) #Array created to result that each array contains indices of elements in the 'type' column of the data array that match the particletype value
    #Now , All the Mass result to returned in units of 10**12 MSun rounded to three decimal places
    #The file stores 3 diffrent types  Halo type (1), Disk type (2), Bulge type (3)
    #Mass assigned 
    
    #So , Now we calculate the mass of an object in units of solar masses
    #Use 1e10 to convert the mass units to units of solar masses
    #then applies the unit of u.solMass to the result ,giving the final mass in units of solar masses
    massresult = data['m'][index]*u.Msun*1e10
    
    #Sum of masses 
    mass = np.sum(massresult)
    
    #Sum over all masses 
    #The question asked in Q2/Part 4 that The Mass should be returned in units of 10**12 M⊙ rounded to three decimal places
    the_total_mass = np.around(mass,3)/1e12
    
    return the_total_mass #return The Total Mass as in Q2

#Q3 : Now we use the returned define function 
# to compute the total mass of each component of each galaxy (MW, M31, M33)
#(MW_000.txt) --> MilkeyWay  , (M31_000.txt) , (M33_000.txt)


# Recall that the particle types are: Halo type (1), Disk type (2), Bulge type (3).
#And we have Three galaxies (MW, M31, M33) 
#So , define the mass for three galaxies (MW, M31, M33) 


#Data located in MW_000.txt and particle types (Halo , Disk , Bulge )
Halo_MW = np.around(ComponentMass('MW_000.txt',particletype=1),3)
Disk_MW = np.around(ComponentMass('MW_000.txt',particletype=2),3)
Bulge_MW = np.around(ComponentMass('MW_000.txt',particletype=3),3)

#Data located in M31_000.txt and particle types (Halo , Disk , Bulge )
Halo_M31 = np.around(ComponentMass('M31_000.txt',particletype=1),3)
Disk_M31 = np.around(ComponentMass('M31_000.txt',particletype=2),3)
Bulge_M31 = np.around(ComponentMass('M31_000.txt',particletype=3),3)

#Data located in M33_000.txt and particle types (Halo , Disk , No Bulge )
#M33 does not possess a bulge 
Halo_M33 = np.around(ComponentMass('M33_000.txt',particletype=1),3)
Disk_M33 = np.around(ComponentMass('M33_000.txt',particletype=2),3)


#Compute the total mass of each galaxy 
MW_comp = np.around((Halo_MW + Disk_MW + Bulge_MW) ,3)
M31_comp = np.around((Halo_M31 + Disk_M31 + Bulge_M31),3)
M33_comp = np.around((Halo_M33 + Disk_M33),3)

# Compute the baryon fraction 
# fbar = total stellar mass / total mass (dark+stellar)
# for each galaxy (MW , M31 , M33)
# and the whole Local Group

fbar_MW = np.around(((Disk_MW + Bulge_MW) /MW_comp),3)
fbar_M31 = np.around(((Disk_M31 + Bulge_M31) /M31_comp),3)
fbar_M33 = np.around(((Disk_M33)/M33_comp),3)
fbar_Whole_Local_Group = np.around(((Disk_MW + Bulge_MW + Disk_M31 + Bulge_M31 + Disk_M33 + Halo_M33)/(MW_comp + M31_comp + M33_comp)),3)
    
print ("Mass MW Halo = ", Halo_MW)
print ("Mass MW Disk = ", Disk_MW)
print ("Mass MW Bulge = ", Bulge_MW)
print ("Mass MW TOTAL = ", MW_comp)
                        
print ("Mass M31 Halo = ", Halo_M31)
print ("Mass M31 Disk = ", Disk_M31)
print ("Mass M31 Bulge = ", Bulge_M31)
print ("Mass M31 TOTAL = ", M31_comp)
                        
print ("Mass M33 Halo = ", Halo_M33)
print ("Mass M33 Disk = ", Disk_M33)
print ("Mass M33 TOTAL = ", M33_comp)

print ("fbar for MW :", fbar_MW)
print ("fbar for M31 :", fbar_M31)
print ("fbar for M33 :", fbar_M33)
print ("fbar for Local Group :", fbar_Whole_Local_Group)