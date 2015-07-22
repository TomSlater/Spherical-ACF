# -*- coding: utf-8 -*-
"""
sph_abs_function.py

A script to plot the absorption correction factor across a spherical nanoparticle.

@author: Tom Slater
"""

import numpy as np
import scipy as sp
import pylab as pl

def sph_absorption(r, p, mass_atten, ang=0, u=0, v=0):
    r_m = ((r**2)-(v**2))**0.5
    t = 2*(r_m**2-u**2)**0.5
    
    A_1 = -(p*mass_atten)
    
    i_1 = lambda z: np.exp(A_1*((r_m**2-(z-0.5*t+u*np.tan(ang))**2*np.cos(ang)**2)**0.5-u*np.cos(ang)-(0.5*t-z)*np.sin(ang)))
    
    result, error = sp.integrate.quad(i_1, 0, t)
    I_1 = result
    
    return(I_1)
    
def sph_ACF(r, p, mass_atten1, mass_atten2, ang=0, u=0, v=0, det1=True, det2=False, det3=False, det4=False):
    det1fac = 0
    det2fac = 0
    det3fac = 0
    det4fac = 0
        
    if det1 == True:
        det1fac = 1
    if det2 == True:
        det2fac = 1
    if det3 == True:
        det3fac = 1
    if det4 == True:
        det4fac = 1
    
    ACF = (det1fac*sph_absorption(r, p, mass_atten1, ang, u, v)+det2fac*sph_absorption(r, p, mass_atten1, ang, -v, -u)+det3fac*sph_absorption(r, p, mass_atten1, ang, -u, -v)+det4fac*sph_absorption(r, p, mass_atten1, ang, v, u))/(det1fac*sph_absorption(r, p, mass_atten2, ang, u, v)+det2fac*sph_absorption(r, p, mass_atten2, ang, -v, -u)+det3fac*sph_absorption(r, p, mass_atten2, ang, -u, -v)+det4fac*sph_absorption(r, p, mass_atten2, ang, v, u))
    
    return(ACF)

#Main function
if __name__ == "__main__":
    thetaE = 0
    dens = 6.67
    mass_atten_Ma = 4650 #pow(0.98*10,3)
    mass_atten_La = 47.9 #pow(1.3*10,2)
    
    radius = 100*10**-7
    
    img = np.zeros([256,256])
    
    thetaE=0*np.pi/180
    
    for y in range(0, 255):
        for x in range(0, 255):
            if (x-128)**2 + (y-128)**2 <= (radius*10**7)**2:
                img[y,x] = sph_ACF(radius, dens, mass_atten_Ma, mass_atten_La, thetaE, (x-128)*10**-7, (y-128)*10**-7,det1=True,det2=True,det3=True,det4=True)
    
    plot = pl.imshow(img)
    #plot.set_cmap('Greys_r')
    plot.set_clim(0.5,1)
    pl.colorbar()