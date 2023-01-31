"""
Class lab1 
Muzoun Alzaabi
01/26/23
"""
# # In Class Lab 1
# Must be uploaded to your Github repository under a "Labs/Lab1" folder by 5 PM on Jan 31st 2023

# ## Part A:  The Local Standard of Rest
# Proper motion of Sgr A* from Reid & Brunthaler 2004
# $\mu = 6.379$ mas/yr 
# 
# Peculiar motion of the sun, $v_\odot$ = 12.24 km/s  (Schonrich 2010)
# 
# 
# $v_{tan} = 4.74 \frac{\mu}{\rm mas/yr} \frac{R_o}{\rm kpc} = V_{LSR} + v_\odot$
# 
# 
# ### a)
# 
# Create a function called VLSR to compute the local standard of res (V$_{LSR}$).
# 
# The function should take as input: the solar radius (R$_o$), the proper motion (mu)
# and the peculiar motion of the sun in the $v_\odot$ direction.
# 
# Compute V$_{LSR}$ using three different values R$_o$: 
# 1. Water Maser Distance for the Sun :  R$_o$ = 8.34 kpc   (Reid 2014 ApJ 783) 
# 2. GRAVITY Collaboration Distance for the Sun:  R$_o$ = 8.178 kpc   (Abuter+2019 A&A 625)
# 3. Value for Distance to Sun listed in Sparke & Gallagher : R$_o$ = 7.9 kpc 
# 


#A:


# Import Modules 
import numpy as np # import numpy
import astropy.units as u # import astropy units
from astropy import constants as const # import astropy constants

def VLSR (Ro , mu=6.379 , vsun=12.24*u.km/u.s) :
    """
    This function will compute the velocity at the local standard of rest 
    
    VLSR=4.74*mu*Ro-vsun 
    
    INPUTS : 
    
    # Ro:'astropy quantity '
    --> unit so , distance from the sun to the glactic center in kpc
    # mu : 'Float'
        The proper motion of Sag A* in mas/yr .
        Default is from Reid & Brunthaler 2004 
    # vsun : 'astropy quantity '
        The peculiar motion of the sun in the v-dir. (km/s)
        Defult is from schonoric +2010 
        
     OUTPUTS : 
        VLSR : 'astropy quantity '
            The velocity of the local standard of rest (km/s)
"""
    return 4.74*mu*(Ro/u.kpc)*u.km/u.s-vsun 


#Define our distances (Dist.)
RoReid = 8.34*u.kpc #Dist. from Reid et al.2004 in kpc 
RoGravity = 8.178*u.kpc #Dist. from the Gravity collab Abuter+2019 in kpc 
RoSG = 7.9*u.kpc #Dist. from TB Sparke & Gallagher 

#compute VLSR Using Ro from Reid 2014 
VLSR_Reid = VLSR(RoReid)
print (VLSR_Reid)

#Compute VLSR using Ro from Gravity Collab 
VLSR_Gravity = VLSR(RoReid)
print (np.round(VLSR_Gravity))

#Compute VLSR using Ro from Sparke & Gallagher 
VLSR_SG = VLSR(RoSG)
print (np.round(VLSR_SG))



# ### b)
# 
# compute the orbital period of the sun in Gyr using R$_o$ from the GRAVITY Collaboration (assume circular orbit)
# 
# Note that 1 km/s $\sim$ 1kpc/Gyr
#A :

def TorbSun (R , V) :
    """ 
    This Function will compute the orbital period of the sun 
        T=2*pi*R / V
    INPUTS : 
        # R : 'astropy quantity ' --> Dist. in Kpc (Dist.To Galactic center)
        # V : 'astropy quantity ' --> Velocity in km/s (Velocity of the sun in v-Direction)
        
    OUTPUTS : 
        'astropy quantity ' --> Orbital period in Gyr
         
    """
    VkpcGyr = V.to(u.kpc/u.Gyr) # convert v from km/s to kpc/Gyr
    T= 2*np.pi*R/VkpcGyr #Orbital period 
    
    return T 

#Velocity of Sun = VLSR + Peculiar Motion 
VsunPeculiar = 12.24*u.km/u.s
VSun = VLSR_Gravity+VsunPeculiar

T_Grav= TorbSun(RoGravity , VSun)
print(T_Grav)






# ### c)
# 
# Compute the number of rotations about the GC over the age of the universe (13.8 Gyr)

#A :
#Age of Universe / Orbital Period 
Age = 13.8*u.Gyr  #Age of the Universe 
print(Age/T_Grav)



# ## Part B  Dark Matter Density Profiles
# 
# ### a)
# Try out Fitting Rotation Curves 
# [here](http://wittman.physics.ucdavis.edu/Animations/RotationCurve/GalacticRotation.html)
# 
# 
# ### b)
# 
# 
# In the Isothermal Sphere model, what is the mass enclosed within the solar radius (R$_o$) in units of M$_\odot$? 
# 
# Recall that for the Isothermal sphere :
# $\rho(r) = \frac{V_{LSR}^2}{4\pi G r^2}$
# 
# Where $G$ = 4.4985e-6 kpc$^3$/Gyr$^2$/M$_\odot$, r is in kpc and $V_{LSR}$ is in km/s
# 
# What about at 260 kpc (in units of  M$_\odot$) ? 


#Gravitational constant  
Grav = const.G.to(u.kpc**3/u.Gyr**2/u.Msun)
"""
Gravitational constant 
Density profile rho = VLSR^2/(4*pi*G*R^2)
MASS = Integrate rho dv 
    = rho *4*pi*4**2 dr 
    =VLSR **2 / G / (4*pr*r**2) *(4*pi*r**2) dr
    = VLSR**2 / G * r 

"""
def MassIso(r,VLSR) :
    """
    this function will compute the darkmatter mass enclosed with a given Dist.
    Assumming an Isothermal Sphere Modle for the dark matter 
    M= VLSR**2 / G*r 
    
    INPUTS : 
        # r : 'astropy quantity ' --> Dist. to the Glactic center (Kpc)
        # VLSR : 'astropy quantity ' --> Velocity of the local standard of REst (km/s)
        
    OUTPUTS : 
        M: Mass enclosed within r in units of Msun
    
    """
    VLSRkpcGyr = VLSR.to(u.kpc/u.Gyr) #conv. km/s --> kpc/Gyr
    M=VLSRkpcGyr**2 / Grav * r  # Mass for isothermal sphere 
    
    return M 


MIsoSolar = MassIso (RoGravity , VLSR_Gravity)
print (MIsoSolar)

#print(f"{MIsoSolar:.2e }")

#compute mass  within 260 kpc 
MIso260= MassIso(260*u.kpc ,VLSR_Gravity)
#print(f" {MIsoSolar:.2e } ")

# ## c) 
# 
# The Leo I satellite is one of the fastest moving satellite galaxies we know. 
# 
# 
# It is moving with 3D velocity of magnitude: Vtot = 196 km/s at a distance of 260 kpc (Sohn 2013 ApJ 768)
# 
# If we assume that Leo I is moving at the escape speed:
# 
# $v_{esc}^2 = 2|\Phi| = 2 \int G \frac{\rho(r)}{r}dV $ 
# 
# and assuming the Milky Way is well modeled by a Hernquist Sphere with a scale radius of $a$= 30 kpc, what is the minimum mass of the Milky Way (in units of M$_\odot$) ?  
# 
# How does this compare to estimates of the mass assuming the Isothermal Sphere model at 260 kpc (from your answer above)





