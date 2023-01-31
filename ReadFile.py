
"""
HW2
Muzoun Alzaabi 

This python code will 
1. open and read the MW 000.txt data file.
2. returns: the (first two lines of the file).
3. returns: particle type, mass, x,y,z, vx,vy,vz columns as a data array 

"""



# Import Modules  by NumPy and AstroPy
import numpy as np
import astropy.units as u


# open the file  
#f = open('MW_000.txt','r')

# Define filename 
filename = "MW_000.txt"


#Defining READ function that will read MW_000.txt file 
def Read (filename) :
    # open the file "MW_000.txt " that defined in prev. steps 
    file = open(filename , 'r' )
    
    #Then following the recpmmended steps in the HW 
    #Read the first line and store the time in units of Myr
    line1 = file.readline()
    label, value = line1.split()
    time = float(value)*u.Myr
    
    #Read 2nd line for the total number of particles
    line2 = file.readline()
    label, value = line2.split()
    num_total_particles = int(value)
    
    
    #close file   
    file.close()
   
    #store the remainder of the file using the NumPy function
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    
    #return the time, total number of particles and the full data array as separate quantities
    return time, num_total_particles, data
#Call function
#time, num_total_particles, data = Read(filename)
    

#Test the code 
#print(data[’type’][1] ) 











