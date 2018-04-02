#!/usr/bin/env python
# Module holding the functions needed to compute the fields
import warnings
import numpy as np
from geometry import geometry as geom
from params import *

if not geom.aset:
    warnings.warn('Geometry not yet defined!',category=ImportWarning)


#------------------------------------------------------------
# Basis functions and their partial derivatives
#------------------------------------------------------------

# For the atmosphere
def Fi(i,x,y):
    w=geom.aftable[i]
    Nx=w['Nx']
    Ny=w['Ny']
    if w['typ']=='A':
        return np.sqrt(2)*np.cos(Ny*y)
    elif w['typ']=='K':
        return 2*np.cos(model.nr*Nx*x)*np.sin(Ny*y)
    else:
        return 2*np.sin(model.nr*Nx*x)*np.sin(Ny*y)

# Partial derivatives
def dxFi(i,x,y):  # in x
    w=geom.aftable[i]
    Nx=w['Nx']
    Ny=w['Ny']
    if w['typ']=='A':
        return 0 #np.sqrt(2)*np.cos(Ny*y)
    elif w['typ']=='K':
        return -2*model.nr*Nx*np.sin(model.nr*Nx*x)*np.sin(Ny*y) #2*np.cos(model.nr*Nx*x)*np.sin(Ny*y)
    else:
        return 2*model.nr*Nx*np.cos(model.nr*Nx*x)*np.sin(Ny*y) #2*np.sin(model.nr*Nx*x)*np.sin(Ny*y)

def dyFi(i,x,y):  # in y
    w=geom.aftable[i]
    Nx=w['Nx']
    Ny=w['Ny']
    if w['typ']=='A':
        return -Ny*np.sqrt(2)*np.sin(Ny*y) #np.sqrt(2)*np.cos(Ny*y)
    elif w['typ']=='K':
        return 2*Ny*np.cos(model.nr*Nx*x)*np.cos(Ny*y) #2*np.cos(model.nr*Nx*x)*np.sin(Ny*y)
    else:
        return 2*Ny*np.sin(model.nr*Nx*x)*np.cos(Ny*y) #2*np.sin(model.nr*Nx*x)*np.sin(Ny*y)

# For the ocean
def phi(i,x,y):
    w=geom.oftable[i]
    Nx=w['Nx']
    Ny=w['Ny']
    return 2*np.sin(model.nr*Nx*x)*np.sin(Ny*y)

# Partial derivatives
def dxphi(i,x,y): # in x
    w=geom.oftable[i]
    Nx=w['Nx']
    Ny=w['Ny']
    return 2*model.nr*Nx*np.cos(model.nr*Nx*x)*np.sin(Ny*y) #2*np.exp(-al*x)*np.sin(model.nr*Nx*x)*np.sin(Ny*y)

def dyphi(i,x,y):
    w=geom.oftable[i]
    Nx=w['Nx']
    Ny=w['Ny']
    return 2*Ny*np.sin(model.nr*Nx*x)*np.cos(Ny*y) #2*np.exp(-al*x)*np.sin(model.nr*Nx*x)*np.sin(Ny*y)

# Function returning the fields based on the coefficients
# ---------------------------------------------------------

# For the atmospheric streamfunction
def astream(x,y,pt):
    Z=pt[0]*Fi(1,x,y)
    for i in range(1,geom.amod):
        Z+=pt[i]*Fi(i+1,x,y)
    return Z

# For the atmospheric wind fields
def avec(x,y,pt):
    U=-pt[0]*dyFi(1,x,y)
    V=pt[0]*dxFi(1,x,y)
    for i in range(1,geom.amod):
        U-=pt[i]*dyFi(i+1,x,y)
        V+=pt[i]*dxFi(i+1,x,y)
    return U,V

# Average of the oceanic basis functions (used to maintain mass conservation)
def Ci(i):
    w=geom.oftable[i]
    Nx=w['Nx']
    Ny=w['Ny']

    return 2*(-1 + (-1)**(Nx*2))*(-1 + (-1)**Ny)/((Nx*2)*Ny*np.pi**2)
    
# Oceanic conserved fields
def ostream_cons(x,y,a):
    Z=a[0]*phi(1,x,y)-a[0]*Ci(1)
    for i in range(1,geom.omod):
        Z+=a[i]*phi(i+1,x,y)-a[i]*Ci(i+1)
    return Z

# Oceanic non-conserved fields (used to compute the temperature fields)
def ostream(x,y,a):
    Z=a[0]*phi(1,x,y)
    for i in range(1,geom.omod):
        Z+=a[i]*phi(i+1,x,y)
    return Z

# Oceanic current vector field
def ovec(x,y,a):
    U=-a[0]*dyphi(1,x,y)
    V=+a[0]*dxphi(1,x,y)
    for i in range(1,geom.omod):
        U-=a[i]*dyphi(i+1,x,y)
        V+=a[i]*dxphi(i+1,x,y)
    return U,V

# Function that return the geopotential height difference
# between North (pi/n,3pi/4) and South (pi/n,pi/4)
def geodiff(pt):
    Zn=pt[0]*Fi(1,np.pi/model.nr,3*np.pi/4)
    Zs=pt[0]*Fi(1,np.pi/model.nr,np.pi/4)
    for i in range(1,geom.amod):
        Zn+=pt[i]*Fi(i+1,np.pi/model.nr,3*np.pi/4)
        Zs+=pt[i]*Fi(i+1,np.pi/model.nr,np.pi/4)
    return Zs-Zn


def compute_frame(line,X,Y):
    y=line.split()
    y=np.array(y,dtype=np.float64)
    psi=y[general.sdi['psi']:general.sdi['theta']]
    theta=y[general.sdi['theta']:general.sdi['A']]
    aa=y[general.sdi['A']:general.sdi['T']]
    tt=y[general.sdi['T']:geom.ndim+1]

    psi=psi*dimension.dimdv['psi']
    theta=theta*dimension.dimdv['theta']
    aa=aa*dimension.dimdv['A']
    tt=tt*dimension.dimdv['T']

    Z=[None,None,None,None]

    if 'op' in view.Zsel:
        Z[view.Zsel.index('op')]=ostream_cons(X,Y,aa)
    if 'ot' in view.Zsel:
        Z[view.Zsel.index('ot')]=ostream(X,Y,tt)
    if 'at' in view.Zsel:
        Z[view.Zsel.index('at')]=astream(X,Y,theta)
    if 'ap' in view.Zsel:
        Z[view.Zsel.index('ap')]=astream(X,Y,psi)

    if 'uo' in view.Zsel or 'vo' in view.Zsel:
        U,V=ovec(X,Y,aa/dimension.dimd['length'])
        if 'uo' in view.Zsel:
            Z[view.Zsel.index('uo')]=U
        if 'vo' in view.Zsel:
            Z[view.Zsel.index('vo')]=V

    if 'ua' in view.Zsel or 'va' in view.Zsel:
        U,V=avec(X,Y,psi*(dimension.dimd['strfunc']/(dimension.dimdv['psi']*dimension.dimd['length'])))
        if 'ua' in view.Zsel:
            Z[view.Zsel.index('ua')]=U
        if 'va' in view.Zsel:
            Z[view.Zsel.index('va')]=V

    if 'ua1' in view.Zsel or 'va1' in view.Zsel or 'ua3' in view.Zsel or 'va3' in view.Zsel:
        U,V=avec(X,Y,psi*(dimension.dimd['strfunc']/(dimension.dimdv['psi']*dimension.dimd['length'])))
        U1,V1=avec(X,Y,theta*(dimension.dimd['strfunc']/(dimension.dimdv['theta']*dimension.dimd['length'])))
        if 'ua1' in view.Zsel:
            Z[view.Zsel.index('ua1')]=U+U1
        if 'va1' in view.Zsel:
            Z[view.Zsel.index('va1')]=V+V1
        if 'ua3' in view.Zsel:
            Z[view.Zsel.index('ua3')]=U-U1
        if 'va3' in view.Zsel:
            Z[view.Zsel.index('va3')]=V-V1

    if 'p1' in view.Zsel or 'p3' in view.Zsel:
        pp=astream(X,Y,psi*(dimension.dimd['strfunc']/dimension.dimdv['psi']))
        pt=astream(X,Y,theta*(dimension.dimd['strfunc']/dimension.dimdv['theta']))
        if 'p1' in view.Zsel:
            Z[view.Zsel.index('p1')]=pp+pt
        if 'p3' in view.Zsel:
            Z[view.Zsel.index('p3')]=pp-pt

    if 'dt' in view.Zsel:
        Z[view.Zsel.index('dt')]=ostream(X,Y,x[3][:,i])-astream(X,Y,x[1][:,i])

    if 'ap1' in view.Zsel or 'ap3' in view.Zsel:
        pp=astream(X,Y,psi)
        pt=astream(X,Y,theta*(dimension.dimdv['psi']/dimension.dimdv['theta']))
        if 'ap1' in view.Zsel:
            Z[view.Zsel.index('ap1')]=pp+pt
        if 'ap3' in view.Zsel:
            Z[view.Zsel.index('ap3')]=pp-pt
    return Z

def compute_quant(line):
    y=line.split()
    y=np.array(y,dtype=np.float64)
    psi=y[general.sdi['psi']:general.sdi['theta']]
    theta=y[general.sdi['theta']:general.sdi['A']]
    aa=y[general.sdi['A']:general.sdi['T']]
    tt=y[general.sdi['T']:geom.ndim+1]

    psi=psi*dimension.dimdv['psi']
    theta=theta*dimension.dimdv['theta']
    aa=aa*dimension.dimdv['A']
    tt=tt*dimension.dimdv['T']

    comp=[psi,theta,aa,tt]

    if 'diff' in view.Isel:
        if 'geo' in view.IIsel[view.Isel.index('diff')]:
            geoap=geodiff(psi)
        else:
            geoap=None
    else:
        geoap=None
    return geoap,comp
