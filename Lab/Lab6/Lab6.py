"""
Lab 6 
Muzoun Alzaabi
ASTR 400B

"""

# In Class Lab 6
# Surface Brightness Profiles




# Load Modules
import numpy as np
import astropy.units as u

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib
#get_ipython().run_line_magic('matplotlib', 'inline')

# my modules
from ReadFile import Read
from CenterOfMass import CenterOfMass
from MassProfile import MassProfile
from GalaxyMass import ComponentMass


# # Lab 6: Sersic Profiles
# 
# In this lab we will use Homework 5 solutions to compute the mass profile of the Milky Way's bulge. 
# We will turn the mass profile into a density profile and see if we can fit it reasonably well 
# with a sersic profile. 

# # Part A : 
# 
# Create a function called `sersicE` that returns the Sersic Profile in terms of the effective radius $R_e$ (i.e. the half light radius).
# 
# $I(r) = I_e exp^{-7.67 ( (r/R_e)^{1/n} - 1)}$
# 
# Where 
# 
# $ L = 7.2 I_e \pi R_e^2$
# 
# and  $R_e$ is the half light radius.  We will assume a mass to light ratio for the stellar bulge of 1, so this is also the half mass radius.
# 
# The function should take as input: the radius, $R_e$, $n$ and the total stellar mass of the system.
# 
def sersicE(r , re , n , mtot ) :
    """
    

    finction that compute the sersic profie for an elliptical galxy 
    assuming M/L ~ 1 
    ----------
    INPUTS : 
    r : float 
        distance from the center of the galaxy (kpc) 
    re : float 
        effictive radius (2D-radius) contain 1/2 the light in [kpc]
    n : float
        sersic index
    mtot : float
        the total stellar mass in Msun
        
    OUTPUTS : 

    Returns 
    I : array of floats 
    the surface btigtness profile of the elliptical Lsun/kpc^2
    -------
    None.

    """
    #Assumming Mlight = 1 
    lum = mtot 
    Ie = lum/7.2/np.pi/re**2 #Effictive surface brigtness 
    
    
    a=(r/re)**(1/2)
    b=-7.67*(a-1)
    I= Ie*np.exp(b)
    
    return I 




  


# # Part B
# 
# a) Create an instance of the MassProfile Class for M31. Store it as a variable `M31`. 
# 
M31 = MassProfile('M31',0 )




# b) Create an array of radii from 0.1 kpc to 30 kpc in increments of 0.1
# 
r = np.arange(0.1,30,0.1)



# c) Define a new array called `bulge_mass`, that uses the function `MassEnclosed` within MassProfile to compute the mass profile of the bulge.  Get rid of astropy units in `bulge_mass` by adding `.value` 
# 
bulge_mass = M31.massEnclosed(3,r).value




# d) Compute the surface mass density profile for the simulated bulge and store it as an array called `bulge_I`. Assuming M/L ~ 1 this is also the surface brightness profile in Lsun/kpc^2

bulge_I = bulge_mass/4/np.pi/r**2




# # Part C
# 
# Compute $R_e$, the half mass radius, for the bulge

#total mass of bulge
bulge_total = ComponentMass("M31_000.txt",3)
print (f"{bulge_total:.2e}")

#half

b_half=bulge_total/2
print(b_half)

#where the bulge mass > 1/2 total mass
index=np.where(bulge_mass>b_half)
print(bulge_mass[index][0])

re_bulge=r[index][0]
print(re_bulge)

# # Part D
# 
# a) Plot the surface density profile of the simulated bulge
# 
# b) Plot the Sersic profile, assuming a de Vaucouleurs Profile.
# 
# c) If the profiles don't match, try changing either $R_e$ or $n$

# In[ ]:


# Plot the Bulge density profile vs 
# the Sersic profile
####################################


fig = plt.figure(figsize=(8,8))
ax = plt.subplot(111)


# plot the bulge mass density as a proxy for surface brighntess
plt.semilogy(r,bulge_I, color='black',linewidth=3, label='Simulated Bulge')

# YOU ADD HERE: Sersic fit to the surface brightness Sersic fit
# Sersic 
plt.semilogy(r,sersicE(r,re_bulge, 4, bulge_total), color='red',linestyle="-." ,linewidth=3, label='Sersic=4')
plt.semilogy(r,sersicE(r,re_bulge, 5.3 , bulge_total), color='orange',linestyle="-." ,linewidth=3, label='Sersic= 5.3')
plt.semilogy(r,sersicE(r,re_bulge, 6 , bulge_total), color='green',linestyle="-." ,linewidth=3, label='Sersic= 6')




#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size


# Add axis labels
plt.xlabel('Radius (kpc)', fontsize=22)
plt.ylabel('Log(I)  $L_\odot/kpc^2$', fontsize=22)



# add a legend with some customizations.
legend = ax.legend(loc='upper right',fontsize='x-large')


plt.savefig('Lab6.png')






