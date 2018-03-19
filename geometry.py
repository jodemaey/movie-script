#!/usr/bin/env python
# Module containing the geometry

class geometry(object):
    ass=None
    oss=None

    ams=None
    oms=None

    amod=None
    omod=None

    ndim=None

    aftable={}
    oftable={}

    aset=False


def set_geometry(ageom,ogeom):

    ass=map(int,ageom.split('x'))
    oss=map(int,ogeom.split('x'))

    ams=[[i,j] for i in range(1,ass[0]+1) for j in range(1,ass[1]+1)]
    oms=[[i,j] for i in range(1,oss[0]+1) for j in range(1,oss[1]+1)]

    amod=2*ass[0]*ass[1]+ass[1]
    omod=oss[0]*oss[1]

    ndim=amod*2+omod*2

    geometry.ass=ass
    geometry.oss=oss

    geometry.ams=ams
    geometry.oms=oms

    geometry.amod=amod
    geometry.omod=omod

    geometry.ndim=ndim

    # Compute the relation table functions index -> functions wavenumbers and type
    ii=0
    aftable={}
    for w in ams:
        if w[0]==1:
            ii+=1
            x={'typ':'A','Nx':0,'Ny':w[1],'Nxi':0,'Nyi':w[1]}
            aftable[ii]=x
            ii+=1
            x={'typ':'K','Nx':w[0],'Ny':w[1],'Nxi':w[0],'Nyi':w[1]}
            aftable[ii]=x
            ii+=1
            x={'typ':'L','Nx':w[0],'Ny':w[1],'Nxi':w[0],'Nyi':w[1]}
            aftable[ii]=x
        else:
            ii+=1
            x={'typ':'K','Nx':w[0],'Ny':w[1],'Nxi':w[0],'Nyi':w[1]}
            aftable[ii]=x
            ii+=1
            x={'typ':'L','Nx':w[0],'Ny':w[1],'Nxi':w[0],'Nyi':w[1]}
            aftable[ii]=x
    geometry.aftable=aftable

    ii=0
    oftable={}
    for w in oms:
        ii+=1
        x={'typ':'L','Nx':w[0]/2.,'Ny':w[1],'Nxi':w[0],'Nyi':w[1]}
        oftable[ii]=x
    geometry.oftable=oftable
    geometry.aset=True

