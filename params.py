#!/usr/bin/env python
# Module containing all the parameters

import numpy as np

class view(object):
#-----------------------------------    
# Organization of the figures layout
#-----------------------------------

# -----------------------------------
# |  info1   |   2      |      4    |
# -----------------------------------
# |  info2   |   1      |      3    |
# -----------------------------------
# info1-2 : information views
# 1--4 : spatial fields representations

# Selection of the spatial fields
#-----------------------------------

# Zsel : List holding the displayed fields in the order of the layout
# 
# Available labels:   - at : atmospheric temperature at 500 mb
#                     - ap : atmospheric pressure at 500 mb
#                     - ot : oceanic temperature
#                     - op : oceanic streamfunction
#                     - p3 : atmospheric lower layer streamfunction psi^3
#                     - p1 : atmospheric upper layer streamfunction psi^1
#                     - ap3 : atmospheric pressure at 750 mb
#                     - ap1 : atmospheric pressure at 250 mb
#                     - dt : ocean-atmosphere temperature difference
#                     - ua : atmospheric U wind component at 500 mb
#                     - va : atmospheric V wind component at 500 mb
#                     - uo : oceanic U current component
#                     - vo : oceanic V current component
#                     - ua1: atmospheric upper U wind component (250 mb)
#                     - va1: atmospheric upper V wind component (250 mb)
#                     - ua3: atmospheric lower U wind component (750 mb)
#                     - va3: atmospheric lower V wind component (750 mb)


    Zsel=['ap','at','op','ot']

# Selection of the 1-2 infoviews modes:
#---------------------------------

    Isel=["diff","3D"] # first relates to info1, second to info2

#Isel components can be : - diff : Difference plot of various quantities (i.e. geopot. height diff.)
#                         - yprof : Profile of various quantities along the spatial direction y
#                         - xprof : Profile of various quantities along the spatial direction x
#                         - mode : Instaneous spectral modes contents
#                         - 3D : 3D projection of the attractor, with locator (warning: both 
#                                infoviews cannot be simultaneously in 3D mode)

    IIsel=[['geo'],[]]

# IIsel : List holding the content to be shown in the infoviews
#         Again, first list relates to info1, second to info2
#
# Available content:
#     For the "diff" mode:
#               - geo : Time evolution of the North-South geopotential height difference at 500mb 
#               - sp(n): Difference between n-th displayed spatial field maximum and minimum
#     For the "yprof" mode:
#               - sp(n): Profile of the (n)-th displayed spatial field in the y spatial direction 
#                         and at the middle of domain
#               - spa(n): Profile of the (n)-th displayed spatial field in the y spatial direction
#                         and averaged in the x direction
#     For the "xprof" mode:
#               - sp(n): Profile of the (n)-th displayed spatial field in the x spatial direction 
#                         and at the middle of domain
#               - spa(n): Profile of the (n)-th displayed spatial field in the x spatial direction
#                         and averaged in the y direction
#     For the "mode" mode, only one of the following:
#               - op: Ocean streamfunction modes
#               - ot: Ocean temperature modes
#               - ap: Atmosphere streamfunction modes
#               - at: Atmosphere temperature modes

class model(object):
#--------------------------------------------------------------------------
# Setting of some general model parameters (those in general do not change)
#--------------------------------------------------------------------------
    nr=1.5                   # aspect ratio
    f0=0.0001032             # Coriolis parameter
    L=5000000./np.pi         # characteristic spatial scale
    rpr=L**2*f0              # streamfunction scaling
    RR=287.                  # Gas constant of dry air
    RK=rpr*f0/RR             # Temperature scaling
    at=365.25                # Duration of a year in days
    ct=(1/(f0*24*3600))/at   # Time scaling from timeunits to years
    geo=f0/9.81              # Geopotential scaling in meters

class mail(object):
#------------------------------------
# Mailserver configuration (optional)
#------------------------------------

# Defining mail address from where and to which send mails
# Not used if no addresses are defined

    fromaddr = ""
    toaddr = ""

# Server to contact
    servername='localhost'

class dimension(object):
    dim=True # dimensionalize the fields

class general(object):
    # Index of the components in the data
    sdi={'psi':1,'theta':amod+1,'A':2*amod+1,'T':2*amod+omod+1,'time':0}


class labels(object):
#-----------------------------
# Defining and ordering labels
#-----------------------------

    Zlabel={'at':r"Atm. Temperature $\theta_a$",'ot':r'Ocean Temperature $\theta_o$','dt':'Oc.-Atm. Temperature diff.','ap':r'Atmospheric $\psi_a$','op':r'Ocean $\psi_o$','p3':r'Atm. low. layer $\psi_a^3$','p1':r'Atm. up. layer $\psi_a^1$','uo':'Ocean U current','vo':'Ocean V current','ua':'Atm. 500mb U wind','va':'Atm. 500mb V wind','ua3':'Atm. 750mb U wind','va3':'Atm. 750mb V wind','ua1':'Atm. 250mb U wind','va1':'Atm. 250mb V wind','ap3':r'Atm. low. layer $\psi_a^3$','ap1':r'Atm. up. layer $\psi_a^1$'}

    strm=r" (m$^2$s$^{-1}$)"
    strg=" (m)"
    strt=r"($^\circ\!$C)"
    strw=r" (ms$^{-1}$)"

    Zlabelunit={'at':strt,'ot':strt,'dt':"years",'ap':strg,'op':strm,'p3':strm,'p1':strm,'ua':strw,'va':strw,'uo':strw,'vo':strw,'ua3':strw,'va3':strw,'ua1':strw,'va1':strw,'ap3':strg,'ap1':strg}

    Zlabelmini={'at':"Atm. T$^\circ$",'ot':'Oc. T$^\circ$','dt':'Oc.-Atm. T$^\circ$ diff.','ap':r'Atm. $\psi_a$','op':r'Oc. $\psi_o$','p3':r'Atm. $\psi_a^3$','p1':r'Atm. $\psi_a^1$','ua':'Atm. U wind','va':'Atm. V wind','uo':'Ocean U current','vo':'Ocean V current','ua3':'Atm. low. U wind','va3':'Atm. low. V wind','ua1':'Atm. up. U wind','va1':'Atm. up. V wind','ap3':r'Atm. $\psi_a^3$','ap1':r'Atm. $\psi_a^1$'}

    #Defining some labels to be used later
    sdd={'ap':'psi','at':'theta','op':'A','ot':'T'}
    sd={'psi':0,'theta':1,'A':2,'T':3,'time':4}
    vl={'psi':r'\psi_{a,','theta':r'\theta_{a,','A':r'\psi_{o,','T':r'\theta_{o,','time':r't'}

    # Infoview labels

    Ivtit={'diff':"Differences plot",'yprof':"Meridional profile",'xprof':"Zonal profile","mode":r"% Modes distribution","3D":'3-D phase space projection'}


def params_initialize():
    labels.Zlab=[]
    labels.Zlabmini=[]
    for x in view.Zsel:
        labels.Zlab.append(labels.Zlabel[x])
        if dimension.dim:
            labels.Zlab[-1]+=labels.Zlabelunit[x]
        labels.Zlabmini.append(labels.Zlabelmini[x])

    # Defining the dico of the dimensionalization
    if dimension.dim:
        dimension.dimd={'geo':model.geo,'strfunc':model.rpr,'temp':model.RK,'timey':model.ct,'times':1/model.f0,'length':model.L}
    else:
        dimension.dimd={'geo':1,'strfunc':1,'temp':1,'timey':1,'times':1,'length':1}

    # Defining the dico of the variables dimensionalization    
    dimension.dimdv={'psi':dimension.dimd['strfunc']*dimension.dimd['geo'],'theta':2*dimension.dimd['temp'],'A':dimension.dimd['strfunc'],'T':dimension.dimd['temp'],'time':dimension.dimd['timey']}
