# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 10:11:04 2015

A script to calculate the absorption correction in Au nanoparticles using EDX detectors at 4 and 18 degree elevation angles.

@author: Tom Slater
"""


from scipy.integrate import quad
import numpy as np
import pylab as pl

thetaE = 0
dens = 19.32
mass_atten_Ma = 1.02*10**3
mass_atten_La = 1.28*10**2

radius = np.linspace(1, 300, 300, endpoint=True)
radius = radius*10**-7

def sph_absorption(r, ang, p, mass_atten1, mass_atten2):
    A_2 = -(p*mass_atten2*r)/np.cos(ang)
    A_1 = -(p*mass_atten1*r)/np.cos(ang)

    I_1=[]
    for A in A_1:
        result, error = quad(lambda x:np.sin(x)*np.exp(A*np.sin(x)), 0, np.pi - 2*ang)
        I_1.append(result)
    I_1 = np.array(I_1)

    I_2=[]
    for A in A_2:
        result, error = quad(lambda x:np.sin(x)*np.exp(A*np.sin(x)), 0, np.pi - 2*ang)
        I_2.append(result)
    I_2 = np.array(I_2)

    Full_1 = r*I_1+(np.sin(ang)/(p*mass_atten1))*(1-np.exp(-2*p*mass_atten1*np.sin(ang)*r))
    Full_2 = r*I_2+(np.sin(ang)/(p*mass_atten2))*(1-np.exp(-2*p*mass_atten2*np.sin(ang)*r))

    Abs = Full_1/Full_2
    return(Abs)

sph18_pointsx = np.array([11.9,29,53,86,97,217,250,280,380])
sph18_pointsy = np.array([1.00,0.98,0.96,0.93,0.94,0.89,0.85,0.84,0.8])
sph18_errory = np.array([0.04,0.01,0.004,0.01,0.02,0.01,0.02,0.02,0.02])
sph18_errorx = np.array([1,4,7,6,5,10,10,10,10])
pl.errorbar(sph18_pointsx,sph18_pointsy,sph18_errory,sph18_errorx,'k.')

sph4_pointsx = np.array([10,100,200,320,350,380,])
sph4_pointsy = np.array([0.99,0.96,0.89,0.83,0.83,0.81])
sph4_errory = np.array([0.02,0.01,0.01,0.02,0.01,0.01])
sph4_errorx = np.array([1,4,7,10,10,10])
pl.errorbar(sph4_pointsx,sph4_pointsy,sph4_errory,sph4_errorx,'b.')

thetaE=4*np.pi/180
angle4=sph_absorption(radius, thetaE, dens, mass_atten_Ma, mass_atten_La)
pl.plot(radius*2*10**7,angle4,'b',label=r'$\theta_E$= 4$^\circ$ Au')

Ptangle4=sph_absorption(radius, thetaE, 20.39, 137, 125)
pl.plot(radius*2*10**7,Ptangle4,'b--',label=r'$\theta_E$= 4$^\circ$ PtAu')

NiOangle4=sph_absorption(radius, thetaE, 6.67, 4650, 47.9)
pl.plot(radius*2*10**7,NiOangle4,'b-.',label=r'$\theta_E$= 4$^\circ$ NiO')

thetaE=18*np.pi/180
angle18=sph_absorption(radius, thetaE, dens, mass_atten_Ma, mass_atten_La)
pl.plot(radius*2*10**7,angle18,'k',label=r'$\theta_E$= 18$^\circ$ Au')

Ptangle18=sph_absorption(radius, thetaE, 20.39, 137, 125)
pl.plot(radius*2*10**7,Ptangle18,'k--',label=r'$\theta_E$= 18$^\circ$ PtAu')

NiOangle18=sph_absorption(radius, thetaE, 6.67, 4650, 47.9)
pl.plot(radius*2*10**7,NiOangle18,'k-.',label=r'$\theta_E$= 18$^\circ$ NiO')

pl.xlabel('Particle Diameter (nm)')
pl.ylabel(r'Normalized Count Ratio')
pl.xlim(0,405)
pl.ylim(0.75,1.1)
pl.legend(loc=1,prop={'size':10})
pl.show()