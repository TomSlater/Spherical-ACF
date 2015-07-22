# -*- coding: utf-8 -*-
"""
Created on Sun Feb 22 14:48:24 2015

A script to calculate the absorption correction factor in thin films and spherical nanoparticles.

@author: Tom Slater
"""
from scipy.integrate import quad
import numpy as np
import pylab as pl

thetaE = 0
dens = 19.32
mass_atten_Ma = 1.02*10**3
mass_atten_La = 1.28*10**2

diam = np.linspace(1, 410, 410, endpoint=True)
diam = diam*10**-7

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

thetaE=18*np.pi/180
pl.plot(sph_absorption(diam/2, thetaE, dens, mass_atten_Ma, mass_atten_La),'k', label="Spherical nanoparticle Au")
pl.plot(sph_absorption(diam/2, thetaE, 20.39, 137, 125),'k--', label="Spherical nanoparticle PtAu")
pl.plot(sph_absorption(diam/2, thetaE, 6.67, 4650, 47.9),'k-.', label="Spherical nanoparticle NiO")

radius = np.linspace(1, 410, 410, endpoint=True)
radius = radius*10**-7

def tf_absorption(t, ang, p, mass_atten1, mass_atten2):

    Full_1 = mass_atten2*(1-np.exp(-mass_atten1*dens*t/np.sin(ang)))
    Full_2 = mass_atten1*(1-np.exp(-mass_atten2*dens*t/np.sin(ang)))

    Abs = Full_1/Full_2
    return(Abs)

tf_pointsx = np.array([23,27,41,64,82])
tf_pointsy = np.array([0.931,0.93,0.911,0.861,0.822])
tf_errory = np.array([0.007,0.004,0.003,0.006,0.009])
tf_errorx = np.array([2,3,4,6,8])
pl.errorbar(tf_pointsx,tf_pointsy,yerr=tf_errory,xerr=tf_errorx,fmt='r.')

sph_pointsx = np.array([11.9,29,53,86,97,217,250,280,380])
sph_pointsy = np.array([1.00,0.98,0.96,0.93,0.94,0.89,0.85,0.84,0.8])
sph_errory = np.array([0.04,0.01,0.004,0.01,0.02,0.01,0.02,0.02,0.02])
sph_errorx = np.array([1,4,7,6,5,10,10,10,10])
pl.errorbar(sph_pointsx,sph_pointsy,yerr=sph_errory,xerr=sph_errorx,fmt='k.')

pl.plot(tf_absorption(radius, thetaE, dens, mass_atten_Ma, mass_atten_La),'r', label="Thin film Au")
pl.plot(tf_absorption(radius, thetaE, 20.39, 137, 125),'r--', label="Thin film PtAu")
pl.plot(tf_absorption(radius, thetaE, 6.67, 4650, 47.9),'r-.', label="Thin film NiO")
pl.xlim(0, 405)
pl.ylim(0.2,1.1)
pl.xlabel('Particle Diameter / Film Thickness (nm)')
pl.ylabel(r'Normalized Count Ratio')
pl.legend(loc=4,prop={'size':10})
pl.show()