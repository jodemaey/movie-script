#!/usr/bin/env python
# Code to compute videos of the output of the MAOOAM model

# Copyright :
# 2016-2018 Jonathan Demaeyer.
# See LICENSE.txt for license information.  

# Usage : ./movie-script.py <data-filename> <ageom> <ogeom>
# Example : ./movie-script.py test.dat 2x4 2x4

# This code needs: 
# - mencoder
# - matplotlib >= 1.5
# - numpy

# TODO : - Move the parameters at the beginning of the code
#        - Generate frames "on the fly"

# Loading of the libraries

import numpy as np
import matplotlib
# matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
# print plt.get_backend()
import matplotlib.animation as anim
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import rc
rc('font',**{'family':'serif','sans-serif':['Times'],'size':14})

import time
import smtplib
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText
import gzip
import sys
import subprocess

#---------------------------------------------
# Parameter section
#---------------------------------------------

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


# Mailserver configuration (optional)
#------------------------------------

# Defining mail address from where and to which send mails
# Not used if no addresses are defined

fromaddr = ""
toaddr = ""

# Server to contact
servername='localhost'

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


#--------------------------
# Preparation
#--------------------------

# Asking for adimensionalization

dim=raw_input('Adimensionalize? (y/N)')
if dim in ['y','Y']:
    dim=False
else:
    dim=True

# Defining and ordering labels
#-----------------------------

Zlabel={'at':r"Atm. Temperature $\theta_a$",'ot':r'Ocean Temperature $\theta_o$','dt':'Oc.-Atm. Temperature diff.','ap':r'Atmospheric $\psi_a$','op':r'Ocean $\psi_o$','p3':r'Atm. low. layer $\psi_a^3$','p1':r'Atm. up. layer $\psi_a^1$','uo':'Ocean U current','vo':'Ocean V current','ua':'Atm. 500mb U wind','va':'Atm. 500mb V wind','ua3':'Atm. 750mb U wind','va3':'Atm. 750mb V wind','ua1':'Atm. 250mb U wind','va1':'Atm. 250mb V wind','ap3':r'Atm. low. layer $\psi_a^3$','ap1':r'Atm. up. layer $\psi_a^1$'}

strm=r" (m$^2$s$^{-1}$)"
strg=" (m)"
strt=r"($^\circ\!$C)"
strw=r" (ms$^{-1}$)"

Zlabelunit={'at':strt,'ot':strt,'dt':"years",'ap':strg,'op':strm,'p3':strm,'p1':strm,'ua':strw,'va':strw,'uo':strw,'vo':strw,'ua3':strw,'va3':strw,'ua1':strw,'va1':strw,'ap3':strg,'ap1':strg}

Zlabelmini={'at':"Atm. T$^\circ$",'ot':'Oc. T$^\circ$','dt':'Oc.-Atm. T$^\circ$ diff.','ap':r'Atm. $\psi_a$','op':r'Oc. $\psi_o$','p3':r'Atm. $\psi_a^3$','p1':r'Atm. $\psi_a^1$','ua':'Atm. U wind','va':'Atm. V wind','uo':'Ocean U current','vo':'Ocean V current','ua3':'Atm. low. U wind','va3':'Atm. low. V wind','ua1':'Atm. up. U wind','va1':'Atm. up. V wind','ap3':r'Atm. $\psi_a^3$','ap1':r'Atm. $\psi_a^1$'}

Zlab=[]
Zlabmini=[]
for x in Zsel:
    Zlab.append(Zlabel[x])
    if dim:
        Zlab[-1]+=Zlabelunit[x]
    Zlabmini.append(Zlabelmini[x])

#Defining some labels to be used later
sdd={'ap':'psi','at':'theta','op':'A','ot':'T'}
sd={'psi':0,'theta':1,'A':2,'T':3,'time':4}
vl={'psi':r'\psi_{a,','theta':r'\theta_{a,','A':r'\psi_{o,','T':r'\theta_{o,','time':r't'}

# Defining the dico of the dimensionalization
if dim:
    dimd={'geo':geo,'strfunc':rpr,'temp':RK,'timey':ct,'times':1/f0,'length':L}
else:
    dimd={'geo':1,'strfunc':1,'temp':1,'timey':1,'times':1,'length':1}

# Defining the dico of the variables dimensionalization    
dimdv={'psi':dimd['strfunc']*dimd['geo'],'theta':2*dimd['temp'],'A':dimd['strfunc'],'T':dimd['temp'],'time':dimd['timey']}

# Infoview labels

Ivtit={'diff':"Differences plot",'yprof':"Meridional profile",'xprof':"Zonal profile","mode":r"% Modes distribution","3D":'3-D phase space projection'}


# Utility functions
#------------------

#Count the number of line of a file
def linecount(filename):
    if filename[-3:]=='.gz':
        lines = 0
        with gzip.open(filename, "r+") as f:
            for x in f:
                lines += 1
    else:
        lines=int(subprocess.Popen("cat "+filename+" | wc -l", shell=True, stdout=subprocess.PIPE).stdout.read())
    return lines

# Gives the order of a number
def order(n):
    if n==0.:
	return 0
    h=np.abs(n)
    if h<1.:
        return int(np.log10(h))-1
    else:
        return int(np.log10(h))+1

# Formating function to be used by matplotlib

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    if x==0.:
	return u'0'
    elif x<0.:
        return unicode(r'{}$\times 10^{{{}}}$'.format(a, b)).replace(u'-',u'\u2212',1)
    else:
        return unicode(r'{}$\times 10^{{{}}}$'.format(a, b))

strf=ticker.FuncFormatter(fmt)

if dim:
    Zformat={'at':None,'ot':None,'dt':None,'ap':None,'op':strf,'p3':strf,'p1':strf,'ua':None,'va':None,'uo':None,'vo':None,'ua3':None,'va3':None,'ua1':None,'va1':None,'ap3':None,'ap1':None}
else:
    Zformat={'at':None,'ot':None,'dt':None,'ap':None,'op':None,'p3':None,'p1':None,'ua':None,'va':None,'uo':None,'vo':None,'ua3':None,'va3':None,'ua1':None,'va1':None,'ap3':None,'ap1':None}

Zcm={'at':cm.coolwarm,'ot':cm.coolwarm,'dt':cm.coolwarm,'ap':cm.gist_rainbow_r,'op':cm.gist_rainbow_r,'p3':cm.jet,'p1':cm.jet,'ua':cm.hsv_r,'va':cm.hsv_r,'uo':cm.hsv_r,'vo':cm.hsv_r,'ua3':cm.hsv_r,'va3':cm.hsv_r,'ua1':cm.hsv_r,'va1':cm.hsv_r,'ap3':cm.jet,'ap1':cm.jet}

# Update function for the bar3D plot

def update_Poly3D(p, x, y, z, dx, dy, dz):
    minx, miny, minz = 1e20, 1e20, 1e20
    maxx, maxy, maxz = -1e20, -1e20, -1e20

    polys = []
    for xi, yi, zi, dxi, dyi, dzi in zip(x, y, z, dx, dy, dz):
        minx = min(xi, minx)
        maxx = max(xi + dxi, maxx)
        miny = min(yi, miny)
        maxy = max(yi + dyi, maxy)
        minz = min(zi, minz)
        maxz = max(zi + dzi, maxz)

        polys.extend([
            ((xi, yi, zi), (xi + dxi, yi, zi),
                (xi + dxi, yi + dyi, zi), (xi, yi + dyi, zi)),
            ((xi, yi, zi + dzi), (xi + dxi, yi, zi + dzi),
                (xi + dxi, yi + dyi, zi + dzi), (xi, yi + dyi, zi + dzi)),

            ((xi, yi, zi), (xi + dxi, yi, zi),
                (xi + dxi, yi, zi + dzi), (xi, yi, zi + dzi)),
            ((xi, yi + dyi, zi), (xi + dxi, yi + dyi, zi),
                (xi + dxi, yi + dyi, zi + dzi), (xi, yi + dyi, zi + dzi)),

            ((xi, yi, zi), (xi, yi + dyi, zi),
                (xi, yi + dyi, zi + dzi), (xi, yi, zi + dzi)),
            ((xi + dxi, yi, zi), (xi + dxi, yi + dyi, zi),
                (xi + dxi, yi + dyi, zi + dzi), (xi + dxi, yi, zi + dzi)),
        ])
    p.set_verts(polys)
 

# Loading the geometries given as arguments or by the users
#-----------------------------------------------------------

if len(sys.argv)==4:
    ageom=sys.argv[2][:]
    ogeom=sys.argv[3][:]
else:
    ageom=raw_input('Atm. geometry ? (default 2x4)')
    ogeom=raw_input('Oc. geometry ? (default 2x4)')

if not ageom:
    ageom='2x4'
if not ogeom:
    ogeom='2x4'

ass=map(int,ageom.split('x'))
oss=map(int,ogeom.split('x'))

ams=[[i,j] for i in range(1,ass[0]+1) for j in range(1,ass[1]+1)]
oms=[[i,j] for i in range(1,oss[0]+1) for j in range(1,oss[1]+1)]

amod=2*ass[0]*ass[1]+ass[1]
omod=oss[0]*oss[1]
# print ams
# print oms
# print amod,omod

ndim=amod*2+omod*2

# Compute the relation table functions index -> functions wavenumbers and type
aftable={}
ii=0
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


oftable={}
ii=0
for w in oms:
    ii+=1
    x={'typ':'L','Nx':w[0]/2.,'Ny':w[1],'Nxi':w[0],'Nyi':w[1]}
    oftable[ii]=x


# Defining the basis functions and their partial derivatives
#------------------------------------------------------------

# For the atmosphere
def Fi(i,x,y):
    w=aftable[i]
    Nx=w['Nx']
    Ny=w['Ny']
    if w['typ']=='A':
        return np.sqrt(2)*np.cos(Ny*y)
    elif w['typ']=='K':
        return 2*np.cos(nr*Nx*x)*np.sin(Ny*y)
    else:
        return 2*np.sin(nr*Nx*x)*np.sin(Ny*y)

# Partial derivatives
def dxFi(i,x,y):  # in x
    w=aftable[i]
    Nx=w['Nx']
    Ny=w['Ny']
    if w['typ']=='A':
        return 0 #np.sqrt(2)*np.cos(Ny*y)
    elif w['typ']=='K':
        return -2*nr*Nx*np.sin(nr*Nx*x)*np.sin(Ny*y) #2*np.cos(nr*Nx*x)*np.sin(Ny*y)
    else:
        return 2*nr*Nx*np.cos(nr*Nx*x)*np.sin(Ny*y) #2*np.sin(nr*Nx*x)*np.sin(Ny*y)

def dyFi(i,x,y):  # in y
    w=aftable[i]
    Nx=w['Nx']
    Ny=w['Ny']
    if w['typ']=='A':
        return -Ny*np.sqrt(2)*np.sin(Ny*y) #np.sqrt(2)*np.cos(Ny*y)
    elif w['typ']=='K':
        return 2*Ny*np.cos(nr*Nx*x)*np.cos(Ny*y) #2*np.cos(nr*Nx*x)*np.sin(Ny*y)
    else:
        return 2*Ny*np.sin(nr*Nx*x)*np.cos(Ny*y) #2*np.sin(nr*Nx*x)*np.sin(Ny*y)

# For the ocean
def phi(i,x,y):
    w=oftable[i]
    Nx=w['Nx']
    Ny=w['Ny']
    return 2*np.sin(nr*Nx*x)*np.sin(Ny*y)

# Partial derivatives
def dxphi(i,x,y): # in x
    w=oftable[i]
    Nx=w['Nx']
    Ny=w['Ny']
    return 2*nr*Nx*np.cos(nr*Nx*x)*np.sin(Ny*y) #2*np.exp(-al*x)*np.sin(nr*Nx*x)*np.sin(Ny*y)

def dyphi(i,x,y):
    w=oftable[i]
    Nx=w['Nx']
    Ny=w['Ny']
    return 2*Ny*np.sin(nr*Nx*x)*np.cos(Ny*y) #2*np.exp(-al*x)*np.sin(nr*Nx*x)*np.sin(Ny*y)

# Function returning the fields based on the coefficients
# ---------------------------------------------------------

# For the atmospheric streamfunction
def astream(x,y,pt):
    Z=pt[0]*Fi(1,x,y)
    for i in range(1,amod):
        Z+=pt[i]*Fi(i+1,x,y)
    return Z

# For the atmospheric wind fields
def avec(x,y,pt):
    U=-pt[0]*dyFi(1,x,y)
    V=pt[0]*dxFi(1,x,y)
    for i in range(1,amod):
        U-=pt[i]*dyFi(i+1,x,y)
        V+=pt[i]*dxFi(i+1,x,y)
    return U,V

# Average of the oceanic basis functions (used to maintain mass conservation)
def Ci(i):
    w=oftable[i]
    Nx=w['Nx']
    Ny=w['Ny']

    return 2*(-1 + (-1)**(Nx*2))*(-1 + (-1)**Ny)/((Nx*2)*Ny*np.pi**2)
    
# Oceanic conserved fields
def ostream_cons(x,y,a):
    Z=a[0]*phi(1,x,y)-a[0]*Ci(1)
    for i in range(1,omod):
        Z+=a[i]*phi(i+1,x,y)-a[i]*Ci(i+1)
    return Z

# Oceanic non-conserved fields (used to compute the temperature fields)
def ostream(x,y,a):
    Z=a[0]*phi(1,x,y)
    for i in range(1,omod):
        Z+=a[i]*phi(i+1,x,y)
    return Z

# Oceanic current vector field
def ovec(x,y,a):
    U=-a[0]*dyphi(1,x,y)
    V=+a[0]*dxphi(1,x,y)
    for i in range(1,omod):
        U-=a[i]*dyphi(i+1,x,y)
        V+=a[i]*dxphi(i+1,x,y)
    return U,V

# Function that return the geopotential height difference
# between North (pi/n,3pi/4) and South (pi/n,pi/4)
def geodiff(pt):
    Zn=pt[0]*Fi(1,np.pi/nr,3*np.pi/4)
    Zs=pt[0]*Fi(1,np.pi/nr,np.pi/4)
    for i in range(1,amod):
        Zn+=pt[i]*Fi(i+1,np.pi/nr,3*np.pi/4)
        Zs+=pt[i]*Fi(i+1,np.pi/nr,np.pi/4)
    return Zs-Zn

# User specification of the data file
if len(sys.argv)==2 or len(sys.argv)==4:
    s=sys.argv[1][:]
    print "Loading from file "+s
else:
    print "No data filename specified as argument."
    s=raw_input('Filename of the data ?')


#-----------------------------------------
# Loading the data
#-----------------------------------------

# Opening files and determing data range
#-----------------------------------------

# Possible legend for the data (not used in general)
#leg=raw_input('Legende?')
leg=''
sl=[s] # Nom du fichier de donnee
fl=[leg] # Legende
nlf=linecount(s) # Counting the number of line of data

# Opening the file
evol=[]
for s in sl:
    if s[-3:]=='.gz':
        evol.append(gzip.open(s,'r'))
    else:
        evol.append(open(s,'r'))

# Computing the frequency of sampling
s=sl[0]
if s[-3:]=='.gz':
    f = gzip.open(s,'r')
else:
    f = open(s,'r')

x1=f.readline()
x2=f.readline()

f.close()

x1=dimd['timey']*float(x1.split()[0])
x2=dimd['timey']*float(x2.split()[0])

# Asking the user from which line to which line he want to read
# providing a gross (experimental) estimation of what it will take in the memory
# Warning : Estimation not accurate for the moment
print 'There is '+str(nlf)+' lines of data in the file'
print 'representing '+str(nlf*(ndim+1)*8/1.e6)+' Mbytes ('+str(nlf*(ndim+1)*8/1.e9)+' GB)'
print 'and '+str((x2-x1)*nlf)+' years of data.'
while True:
    sti=raw_input('Where do you want to start reading ? (default=first)')
    ste=raw_input('And where do you want to stop ? (default=last)')
    itr=raw_input('What interval should we use to read lines ? (default=1)')
    if not sti:
        sti=0
    else:
        sti=int(sti)-1
    if not ste:
        ste=nlf-1
    else:
        ste=int(ste)-1
    if not itr:
        itr=1
    else:
        itr=int(itr)
    print 'It will represent '+str((ste-sti)*(ndim+1)*8/1.e6/itr)+' Mbytes of data in the memory. ('+str((ste-sti)*(ndim+1)*8/1.e9/itr)+' GB)'
    print 'and '+str((x2-x1)*(ste-sti))+' years of data.'
    x=raw_input('Do you agree with this ? (y/N)')
    if x in ['y','Y']:
        break
    else:
        print 'No? Ok...'

# Defining the variables that will be shown in the "attractor" view
#-------------------------------------------------------------------

ls=[]
ms=[]
showp=[]
showp2=[]
showp3=[]
n1=[]
n2=[]
n3=[]
s=''
for x in sd.keys():
    s+=x+', '
s=s[:-2]
for i in range(len(sl)):
    print 'Parameter for the file ',sl[i]
    x='A' #raw_input('Var ('+s+') ?')
    showp.append(x)
    if x=='time':
        n1.append(0)
    else:
        x=2 #input('Number?')
        n1.append(x-1)
    x='T' #raw_input('Var2 ('+s[:-6]+') ?')
    showp2.append(x)
    x=2 #input('Number?')
    n2.append(x-1)
    x='psi' #raw_input('Var3 ('+s[:-6]+') ?')
    showp3.append(x)
    if not x:
        n3.append(-1)
    else:
        x=1 #input('Number ?')
        n3.append(x-1)

    x='' #raw_input('Line style ? (default=None)')
    ls.append(x)
    x='' #raw_input('Marker ? (default=pixel)')
    if not x:
        ms.append(',')
    else:
        ms.append(x)

    
# Retrieving the data from the files
#---------------------------------------

sete=[]
for j in range(len(evol)):
    e=evol[j]
    tu=[]
    psi=[]
    theta=[]
    for i in range(amod):
        psi.append([])
        theta.append([])
    aa=[]
    tt=[]
    for i in range(omod):
        aa.append([])
    for i in range(omod):
        tt.append([])

    ef=[]
    ii=0
    for x in e:
        if ii>=sti and ii<=ste:
            if np.mod(ii-sti,itr)==0:
                ef.append(x)
        if ii>ste:
            break
        ii+=1

    for x in ef:
        y=x.split()
        tu.append(float(y[0]))
        for i in range(1,amod+1):
            psi[i-1].append(float(y[i]))
            theta[i-1].append(float(y[i+amod]))
        for i in range(omod):
            aa[i].append(float(y[1+2*amod+i]))
        for i in range(omod):
            tt[i].append(float(y[1+2*amod+omod+i]))

    tu=np.array(tu)
    tup=tu.copy()
    tup.shape=1,len(tu)
    sete.append([np.array(psi)*dimdv['psi'],np.array(theta)*dimdv['theta'],np.array(aa)*dimdv['A'],np.array(tt)*dimdv['T'],tup*dimdv['time']])
    tu=tu*dimdv['time']

ival=1
tl=[]
for i in range(len(sete)):
    tl.append(len(sete[i][0][0]))



# Sending a mail to alert that the first stage is completed (loading the data
# and preparing the plots)

if fromaddr and toaddr:
    msg = MIMEMultipart()
    msg['From'] = fromaddr
    msg['To'] = toaddr
    msg['Subject'] = "Movie run info"

    body = "The run to generate the video for the geometry:\n\n"
    body += "    atm. "+ageom+" -  oc."+ogeom+"\n\n"
    body += "has finished loading data! "

    msg.attach(MIMEText(body, 'plain'))

    server = smtplib.SMTP(servername)
    text = msg.as_string()
    server.sendmail(fromaddr, toaddr, text)
    server.quit()

#------------------------------
# Spatial fields computation
#------------------------------

# Setting the grids
delta=raw_input('Space between points on the grid (default = 0.025) ?') 
if not delta:
    delta=0.025
else:
    delta=float(delta)
x = np.arange(0., 2*np.pi/nr, delta)
y = np.arange(0., np.pi, delta)
X, Y = np.meshgrid(x, y)
sh=X.shape

# Setting the number of frames and the time of the first and the last one
while True:
    print 'Total number of frames:',tl[0]
    if dim:
        print 'Time between each frame is '+str((tu[1]-tu[0])*at)+' days'
    else:
        print 'Time between each frame is '+str(tu[1]-tu[0])+' timeunit'
    sti=raw_input('Start at frame (default = first) ?')
    ste=raw_input('End at frame (default = last) ?')
    ite=raw_input('Interval (default = 1) ?')
    if not sti:
        sti=0
    else:
        sti=int(sti)
    if not ste:
        ste=tl[0]
    else:
        ste=int(ste)
    if not ite:
        ite=1
    else:
        ite=int(ite)
    if dim:
        print 'Time between each frame is now '+str(ite*(tu[1]-tu[0])*at)+' days'
    else:
        print 'Time between each frame is '+str(ite*(tu[1]-tu[0]))+' timeunits'
    print 'Number of frames that will effectively be computed :'+str((ste-sti)/ite)
    if dim:
        print 'for a total time of '+str((ste-sti)*(tu[1]-tu[0]))+' years'
    else:
        print 'for a total time of '+str((ste-sti)*(tu[1]-tu[0]))+' timeunits'
    print 'It will take '+str((ste-sti)*(4*sh[0]*sh[1]+4*3+1)*8/(1.e6*ite))+' Mbytes in the memory! ('+str((ste-sti)*(4*sh[0]*sh[1]+4*3+1)*8/(1.e9*ite))+' GB)'
    x=raw_input('Do you agree (y/N) ?')
    if x in ['y','Y']:
        break
    else:
        print "Let's start again..."

# Loop generating the frame (computing the fields)
#-------------------------------------------------

# Preparing space to store the fields

# spatial field + minmax
Z=[]
Zm=[]
ZM=[]

for i in range(4):
    Z.append([])
    Zm.append([])
    ZM.append([])

# Meridional profile
yprof=[]
yprofmid=[]
yprofmidave=[]
yprofave=[]

for i in range(4):
    yprof.append([])
    yprofmid.append([])
    yprofmidave.append([])
    yprofave.append([])

# Zonal profile
xprof=[]
xprofmid=[]
xprofmidave=[]
xprofave=[]

for i in range(4):
    xprof.append([])
    xprofmid.append([])
    xprofmidave.append([])
    xprofave.append([])

# Geopotential height difference at 500mb
geoap=[]

# Modes distribution

za=[[],[]]
zk=[[],[]]
zl=[[],[]]

if 'mode' in Isel:
    if 'a' in IIsel[Isel.index('mode')][0]:
        ss=ass
    else:
        ss=oss

    z0=np.zeros(ss).flatten()
    nx=np.arange(ss[0])-0.5
    ny=np.arange(ss[1])
    NXa, NYa = np.meshgrid(nx-1.,ny)
    NXa=NXa.flatten()
    NYa=NYa.flatten()
    NXk, NYk = np.meshgrid(nx,ny)
    NXk=NXk.flatten()
    NYk=NYk.flatten()
    NXl, NYl = np.meshgrid(nx+0.5,ny+0.5)
    NXl=NXl.flatten()
    NYl=NYl.flatten()

#overall fields max and min
mmin=np.zeros((4))
mmax=np.zeros((4))

startt=time.time()
x=sete[0]
ifsmax=[None,None]
ifsmin=[None,None]
for i in range(sti,ste,ite):
    if np.mod(i-sti,100*ite)==0:
        print 'Generating the fields in the frame ',i,'('+str((i-sti)/ite)+')'
        if dim:
            print 'At time t=',x[4][0,i],'years'
        else:
            print 'At time t=',x[4][0,i],'timeunits'
    if 'op' in Zsel:
        Z[Zsel.index('op')].append(ostream_cons(X,Y,x[2][:,i]))
    if 'ot' in Zsel:
        Z[Zsel.index('ot')].append(ostream(X,Y,x[3][:,i]))
    if 'at' in Zsel:
        Z[Zsel.index('at')].append(astream(X,Y,x[1][:,i]))
    if 'ap' in Zsel:
        Z[Zsel.index('ap')].append(astream(X,Y,x[0][:,i]))

    if 'uo' in Zsel or 'vo' in Zsel:
        U,V=ovec(X,Y,x[2][:,i]/dimd['length'])
        if 'uo' in Zsel:
            Z[Zsel.index('uo')].append(U)
        if 'vo' in Zsel:
            Z[Zsel.index('vo')].append(V)

    if 'ua' in Zsel or 'va' in Zsel:
        U,V=avec(X,Y,x[0][:,i]*(dimd['strfunc']/(dimdv['psi']*dimd['length'])))
        if 'ua' in Zsel:
            Z[Zsel.index('ua')].append(U)
        if 'va' in Zsel:
            Z[Zsel.index('va')].append(V)

    if 'ua1' in Zsel or 'va1' in Zsel:
        U,V=avec(X,Y,x[0][:,i]*(dimd['strfunc']/(dimdv['psi']*dimd['length'])))
        U1,V1=avec(X,Y,x[1][:,i]*(dimd['strfunc']/(dimdv['theta']*dimd['length'])))
        if 'ua1' in Zsel:
            Z[Zsel.index('ua1')].append(U+U1)
        if 'va1' in Zsel:
            Z[Zsel.index('va1')].append(V+V1)

    if 'ua3' in Zsel or 'va3' in Zsel:
        U,V=avec(X,Y,x[0][:,i]*(dimd['strfunc']/(dimdv['psi']*dimd['length'])))
        U1,V1=avec(X,Y,x[1][:,i]*(dimd['strfunc']/(dimdv['theta']*dimd['length'])))
        if 'ua3' in Zsel:
            Z[Zsel.index('ua3')].append(U-U1)
        if 'va3' in Zsel:
            Z[Zsel.index('va3')].append(V-V1)

    if 'p1' in Zsel:
        Z[Zsel.index('p1')].append(astream(X,Y,x[0][:,i]*(dimd['strfunc']/dimdv['psi']))+astream(X,Y,x[1][:,i]*(dimd['strfunc']/dimdv['theta'])))
    if 'p3' in Zsel:
        Z[Zsel.index('p3')].append(astream(X,Y,x[0][:,i]*(dimd['strfunc']/dimdv['psi']))-astream(X,Y,x[1][:,i]*(dimd['strfunc']/dimdv['theta'])))
    if 'dt' in Zsel:
        Z[Zsel.index('dt')].append(ostream(X,Y,x[3][:,i])-astream(X,Y,x[1][:,i]))

    if 'ap1' in Zsel:
        Z[Zsel.index('ap1')].append(astream(X,Y,x[0][:,i])+astream(X,Y,x[1][:,i]*(dimdv['psi']/dimdv['theta'])))
    if 'ap3' in Zsel:
        Z[Zsel.index('ap3')].append(astream(X,Y,x[0][:,i])-astream(X,Y,x[1][:,i]*(dimdv['psi']/dimdv['theta'])))

    if 'diff' in Isel:
        if 'geo' in IIsel[Isel.index('diff')]:
            geoap.append(geodiff(x[0][:,i]))

    for j in range(4):
        ZM[j].append(np.amax(Z[j][-1]))
        Zm[j].append(np.amin(Z[j][-1]))

    for j in range(4):
        mmax[j]=max(mmax[j],ZM[j][-1])
        mmin[j]=min(mmin[j],Zm[j][-1])
    
    if 'yprof' in Isel:
        for j in range(4):
            yprof[j].append(np.mean(Z[j][-1],axis=1))
            yprofmid[j].append(Z[j][-1][:,sh[1]/2])
            if i==sti:
                yprofmidave[j].append(yprofmid[j][-1])
                yprofave[j].append(yprof[j][-1])
            else:
                y=yprofmidave[j][-1]+(yprofmid[j][-1]-yprofmidave[j][-1])/(i-sti)
                yprofmidave[j].append(y)
                y=yprofave[j][-1]+(yprof[j][-1]-yprofave[j][-1])/(i-sti)
                yprofave[j].append(y)
    if 'xprof' in Isel:
        for j in range(4):
            xprof[j].append(np.mean(Z[j][-1],axis=0))
            xprofmid[j].append(Z[j][-1][sh[0]/2,:])
            if i==sti:
                xprofmidave[j].append(xprofmid[j][-1])
                xprofave[j].append(xprof[j][-1])
            else:
                y=xprofmidave[j][-1]+(xprofmid[j][-1]-xprofmidave[j][-1])/(i-sti)
                xprofmidave[j].append(y)
                y=xprofave[j][-1]+(xprof[j][-1]-xprofave[j][-1])/(i-sti)
                xprofave[j].append(y)
    if "mode" in Isel:
        iii=0
        for z in Isel:
            if z=='mode':
                za[iii].append(np.zeros(ss))
                zk[iii].append(np.zeros(ss))
                zl[iii].append(np.zeros(ss))
                y=np.absolute(x[sd[sdd[IIsel[iii][0]]]][:,i])
                y=100*y/y.sum()
                for ii in range(1,len(y)+1):
                    if 'a' in IIsel[iii][0]:
                        af=aftable[ii]
                    else:
                        af=oftable[ii]
                    if af['typ']=='A':
                        za[iii][-1][af['Nxi'],af['Nyi']-1]=y[ii-1]
                    elif af['typ']=='K':
                        zk[iii][-1][af['Nxi']-1,af['Nyi']-1]=y[ii-1]
                    elif af['typ']=='L':
                        zl[iii][-1][af['Nxi']-1,af['Nyi']-1]=y[ii-1]
                za[iii][-1] = za[iii][-1].T
                za[iii][-1] += 1e-10
                zk[iii][-1] = zk[iii][-1].T
                zl[iii][-1] = zl[iii][-1].T
                za[iii][-1]=za[iii][-1].flatten()
                zk[iii][-1]=zk[iii][-1].flatten()
                zl[iii][-1]=zl[iii][-1].flatten()
                ifsmax[iii]=max(ifsmax[iii],np.amax(za[iii][-1]),np.amax(zk[iii][-1]),np.amax(zl[iii][-1]))
            iii+=1

    
ZM=np.array(ZM)
Zm=np.array(Zm)
if 'diff' in Isel:
    if 'geo' in IIsel[Isel.index('diff')]:
        geoap=np.array(geoap)

if "diff" in Isel:
    diff=[[],[]]
    difflab=[[],[]]

if "yprof" in Isel:
    prof=[[],[]]
    profave=[[],[]]
    proflab=[[],[]]

if "xprof" in Isel:
    prof2=[[],[]]
    prof2ave=[[],[]]
    prof2lab=[[],[]]

ii=0
for z in Isel:
    if z=="diff":
        for x in IIsel[ii]:
            if 'sp' in x:
                i=int(x[2])-1
                diff[ii].append(ZM[i]-Zm[i])
                difflab[ii].append(Zlabmini[i])
                if dim:
                    difflab[ii][-1]+=Zlabelunit[Zsel[i]]
            if x=='geo':
                diff[ii].append(geoap)
                difflab[ii].append('Geop. H.')
                if dim:
                    difflab[ii][-1]+=strg

        ifsmax[ii]=max(map(np.amax,diff[ii]))
        ifsmin[ii]=min(map(np.amin,diff[ii]))
    if z=="yprof":
        for x in IIsel[ii]:
            if 'spa' in x:
                i=int(x[3])-1
                prof[ii].append(yprof[i])
                profave[ii].append(yprofave[i])
                proflab[ii].append('Z.A. '+Zlabmini[i])
                if dim:
                    proflab[ii][-1]+=Zlabelunit[Zsel[i]]
            elif 'sp' in x:
                i=int(x[2])-1
                prof[ii].append(yprofmid[i])
                profave[ii].append(yprofmidave[i])
                proflab[ii].append(Zlabmini[i])
                if dim:
                    proflab[ii][-1]+=Zlabelunit[Zsel[i]]
        ifsmax[ii]=[]
        ifsmin[ii]=[]
        for x in prof[ii]:
            ifsmax[ii].append(max(map(np.amax,x)))
            ifsmin[ii].append(min(map(np.amin,x)))
        ifsmax[ii]=max(ifsmax[ii])
        ifsmin[ii]=min(ifsmin[ii])
    if z=="xprof":
        for x in IIsel[ii]:
            if 'spa' in x:
                i=int(x[3])-1
                prof2[ii].append(xprof[i])
                prof2ave[ii].append(xprofave[i])
                prof2lab[ii].append('Z.A. '+Zlabmini[i])
                if dim:
                    prof2lab[ii][-1]+=Zlabelunit[Zsel[i]]
            elif 'sp' in x:
                i=int(x[2])-1
                prof2[ii].append(xprofmid[i])
                prof2ave[ii].append(xprofmidave[i])
                prof2lab[ii].append(Zlabmini[i])
                if dim:
                    prof2lab[ii][-1]+=Zlabelunit[Zsel[i]]
        ifsmax[ii]=[]
        ifsmin[ii]=[]
        for x in prof2[ii]:
            ifsmax[ii].append(max(map(np.amax,x)))
            ifsmin[ii].append(min(map(np.amin,x)))
        ifsmax[ii]=max(ifsmax[ii])
        ifsmin[ii]=min(ifsmin[ii])
    ii+=1

#--------------------------------------
# Setup of the plots
#--------------------------------------

fig=plt.figure(num=10,figsize=(16,9))

# General title

fig.suptitle('Spatial fields in the model MAOOAM')

# General subtitle

suptit=r'atmosphere $'+ageom.replace('x','x$-$')+r'y$  ocean $'+ogeom.replace('x','x$-$')+r'y$'
fig.text(0.42,0.92,'Resolution : '+suptit)

# Setting the six views

if Isel[0] in ["3D","mode"]:
    ax1=fig.add_subplot(2,3,1,projection='3d')
else:
    ax1=fig.add_subplot(2,3,1)
if Isel[1] in ["3D","mode"]:
    ax2=fig.add_subplot(2,3,4,projection='3d')
else:
    ax2=fig.add_subplot(2,3,4)
ax3=fig.add_subplot(2,3,2)
ax4=fig.add_subplot(2,3,3)
ax5=fig.add_subplot(2,3,5)
ax6=fig.add_subplot(2,3,6)

# Views title

ax1.set_title(Ivtit[Isel[0]])
if Isel[0]=="mode":
    ax1.text2D(0.5, 0.9,Zlabelmini[IIsel[0][0]], horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes,fontdict={'size':12})
if Isel[0]=="3D":
    ax1.text2D(0.5, 0.9,'(Non-dimensional units)', horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes,fontdict={'size':12})
ax2.set_title(Ivtit[Isel[1]])
if Isel[1]=="mode":
    ax2.text2D(0.5, 0.9,Zlabelmini[IIsel[1][0]], horizontalalignment='center',verticalalignment='center',transform=ax2.transAxes,fontdict={'size':12})
if Isel[1]=="3D":
    ax2.text2D(0.5, 0.9,'(Non-dimensional units)', horizontalalignment='center',verticalalignment='center',transform=ax2.transAxes,fontdict={'size':12})

ax3.set_title(Zlab[1])
ax4.set_title(Zlab[3]) 
ax5.set_title(Zlab[0])
ax6.set_title(Zlab[2])

axm=[ax1,ax2]

# 3D view axis ranges and labels
#----------------------------------

# Range
if "3D" in Isel:
    axs=axm[Isel.index('3D')]
    mv=0.25 # Overflow factor
    smax=3*[-30000.]
    x=sete[0]
    smax[0]=max(smax[0],np.amax(x[sd[showp[0]]][n1[0],:]))
    smax[1]=max(smax[1],np.amax(x[sd[showp2[0]]][n2[0],:]))
    smax[2]=max(smax[2],np.amax(x[sd[showp3[0]]][n3[0],:]))
    smin=smax[:]
    smin[0]=min(smin[0],np.amin(x[sd[showp[0]]][n1[0],:]))
    smin[1]=min(smin[1],np.amin(x[sd[showp2[0]]][n2[0],:]))
    smin[2]=min(smin[2],np.amin(x[sd[showp3[0]]][n3[0],:]))
    dmm=[]
    for i in range(3):
        dmm.append(smax[i]-smin[i])

    # Rescaling
    smax[0]=smax[0]+dmm[0]*mv
    smax[1]=smax[1]+dmm[1]*mv
    smax[2]=smax[2]+dmm[2]*mv
    smin[0]=smin[0]-dmm[0]*mv
    smin[1]=smin[1]-dmm[1]*mv
    smin[2]=smin[2]-dmm[2]*mv

    # Setting the limits
    axs.set_xlim(smin[0],smax[0])
    axs.set_ylim(smin[1],smax[1])
    axs.set_zlim(smin[2],smax[2])
    fig.canvas.draw()

    # Setting the ticks and the labels on the 3D view

    # x ticks and axis label

    labels = [item.get_text() for item in axs.get_xticklabels()]
    ii=len(labels)-1
    labto=[]
    #	print labels
    jk=0
    for x in labels:
        if x:
            jk+=1
            y=float(x.replace(u'\u2212',u'-'))
            y=y/dimdv[showp[0]]
            if jk==1:
                n=order(y)
            if abs(n)>2:
                y=y*10**(-n)
                y=round(y,2)
                # labto.append(unicode(y).replace(u'-',u'\u2212')+u'e'+unicode(n).replace(u'-',u'\u2212'))
                labto.append(unicode(y).replace(u'-',u'\u2212')+unicode(r'$\times 10^{'+str(n)+r'}$'))
            else:
                y=round(y,2)
                labto.append(unicode(y).replace(u'-',u'\u2212'))
        else:
            ii-=1
            labto.append(x)
    axs.set_xticklabels(labto)
    til1=axs.xaxis.get_major_ticks()
    for j in range(0,len(til1),1):
        til1[j].label1.set_visible(False)
    til1[ii].label1.set_visible(True)
    til1[0].label1.set_visible(True)

    if showp[0]=='time':
        axs.set_xlabel('\n\n'+r'$'+vl[showp[0]]+r'$',fontdict={'size':18})
    else:
    #    if abs(n)>2:
    #        axs.set_xlabel('\n\n'+r'$'+vl[showp[0]]+str(n1[0]+1)+r'}$ ${(\times 10^{'+str(-n)+r'})}$',fontdict={'size':18})
    #    else:
        axs.set_xlabel('\n\n\n'+r'$'+vl[showp[0]]+str(n1[0]+1)+r'}$',fontdict={'size':18})

    # y ticks and axis label

    labels = [item.get_text() for item in axs.get_yticklabels()]
    ii=len(labels)-1
    labto=[]
    jk=0
    for x in labels:
        if x:
            jk+=1
            y=float(x.replace(u'\u2212',u'-'))
            y=y/dimdv[showp2[0]]
            if jk==1:
                n=order(y)
            if abs(n)>2:
                y=y*10**(-n)
                y=round(y,2)
                # labto.append(unicode(y).replace(u'-',u'\u2212')+u'e'+unicode(n).replace(u'-',u'\u2212'))
                labto.append(unicode(y).replace(u'-',u'\u2212')+unicode(r'$\times 10^{'+str(n)+r'}$'))
            else:
                y=round(y,2)
                labto.append(unicode(y).replace(u'-',u'\u2212'))
        else:
            ii-=1
            labto.append(x)
    axs.set_yticklabels(labto)
    til2=axs.yaxis.get_major_ticks()
    for j in range(0,len(til2),1):
        til2[j].label1.set_visible(False)
    til2[ii].label1.set_visible(True)
    til2[0].label1.set_visible(True)

    #if abs(n)>2:
    #    axs.set_ylabel('\n\n'+r'$'+vl[showp2[0]]+str(n2[0]+1)+r'}$ ${(\times 10^{'+str(-n)+r'})}$',fontdict={'size':18})
    #else:
    axs.set_ylabel('\n\n'+r'$'+vl[showp2[0]]+str(n2[0]+1)+r'}$',fontdict={'size':18})

    # z ticks

    labels = [item.get_text() for item in axs.get_zticklabels()]
    ii=len(labels)-1
    labto=[]
    jk=0
    for x in labels:
        if x:
            jk+=1
            y=float(x.replace(u'\u2212',u'-'))
            y=y/dimdv[showp3[0]]
            if jk==1:
                n=order(y)
            if abs(n)>2:
                y=y*10**(-n)
                y=round(y,2)
                # labto.append(unicode(y).replace(u'-',u'\u2212')+u'e'+unicode(n).replace(u'-',u'\u2212'))
                labto.append(unicode(y).replace(u'-',u'\u2212')+unicode(r'$\times 10^{'+str(n)+r'}$'))
            else:
                y=round(y,2)
                labto.append(unicode(y).replace(u'-',u'\u2212'))
        else:
            ii-=1
            labto.append(x)
    axs.set_zticklabels(labto)
    til3=axs.zaxis.get_major_ticks()
    for j in range(0,len(til3),1):
        til3[j].label1.set_visible(False)
    til3[ii-1].label1.set_visible(True)
    til3[0].label1.set_visible(True)

    #if abs(n)>2:
    #    axs.set_zlabel(''+r'$'+vl[showp3[0]]+str(n3[0]+1)+r'}$ ${(\times 10^{'+str(-n)+r'})}$',fontdict={'size':18})    
    #else:
    axs.set_zlabel(''+r'$'+vl[showp3[0]]+str(n3[0]+1)+r'}$',fontdict={'size':18})


    # ticks alignement

    [t.set_va('center') for t in axs.get_yticklabels()]
    [t.set_ha('left') for t in axs.get_yticklabels()]
    [t.set_va('center') for t in axs.get_xticklabels()]
    [t.set_ha('right') for t in axs.get_xticklabels()]
    [t.set_va('center') for t in axs.get_zticklabels()]
    [t.set_ha('left') for t in axs.get_zticklabels()]

# Other views labels
#--------------------

ax3.set_xlabel('$x^\prime$')
ax3.set_ylabel('$y^\prime$')

ax4.set_xlabel('$x^\prime$')
ax4.set_ylabel('$y^\prime$')

ax5.set_xlabel('$x^\prime$')
ax5.set_ylabel('$y^\prime$')

ax6.set_xlabel('$x^\prime$')
ax6.set_ylabel('$y^\prime$')

# Setting the limit of the infoview
#-------------------------------------------------------
i=0
for x in Isel:
    axs=axm[i]
    if x=="diff":
        axs.set_xlim(0.,sete[0][4][0,ste-1]-sete[0][4][0,sti])
        if ifsmin>0.:
            axs.set_ylim(0.9*ifsmin[i],1.4*ifsmax[i])
        else:
            axs.set_ylim(1.1*ifsmin[i],1.4*ifsmax[i])
        if dim:
            axs.set_xlabel('time (years)')
        else:
            axs.set_xlabel('time (timeunits)')
    if x=="yprof":
        axs.set_xlim(Y[0,0],Y[-1,0])
        if ifsmin>0.:
            axs.set_ylim(0.9*ifsmin[i],1.4*ifsmax[i])
        else:
            axs.set_ylim(1.1*ifsmin[i],1.4*ifsmax[i])
        axs.set_xlabel(r'$y^\prime$')
    if x=="xprof":
        axs.set_xlim(X[0,0],X[0,-1])
        if ifsmin>0.:
            axs.set_ylim(0.9*ifsmin[i],1.4*ifsmax[i])
        else:
            axs.set_ylim(1.1*ifsmin[i],1.4*ifsmax[i])
        axs.set_xlabel(r'$x^\prime$')
    if x=="mode":
        axs.set_zlim(0.,ifsmax[i]+10)
    i+=1

#------------------------------
# Initialization of the animation
#------------------------------

# Infoview plot
ii=0
xrlines=[[],[]]
xralines=[[],[]]
pa=[None,None]
pk=[None,None]
pl=[None,None]
for x in Isel:
    axs=axm[ii]
    if x=="diff":
        for i in range(len(diff[ii])):
            xrlines[ii].append(None)
        for i in range(len(diff[ii])):
            xrlines[ii][i],=axs.plot([],[],label=difflab[ii][i])
        axs.legend(fontsize=12)
    if x=="yprof":
        for i in range(len(prof[ii])):
            xrlines[ii].append(None)
        for i in range(len(prof[ii])):
            xrlines[ii][i],=axs.plot([],[],label=proflab[ii][i])
        for i in range(len(prof[ii])):
            xrlines[ii][i].set_xdata(Y[:,0])

        for i in range(len(profave[ii])):
            xralines[ii].append(None)
        for i in range(len(profave[ii])):
            xralines[ii][i],=axs.plot([],[],color=xrlines[ii][i].get_color(),ls=':')
        for i in range(len(profave[ii])):
            xralines[ii][i].set_xdata(Y[:,0])
        axs.plot([],[],color='k',ls=':',label='Time ave.')
        axs.legend(fontsize=12)
    if x=="xprof":
        for i in range(len(prof2[ii])):
            xrlines[ii].append(None)
        for i in range(len(prof2[ii])):
            xrlines[ii][i],=axs.plot([],[],label=prof2lab[ii][i])
        for i in range(len(prof2[ii])):
            xrlines[ii][i].set_xdata(X[0,:])

        for i in range(len(prof2ave[ii])):
            xralines[ii].append(None)
        for i in range(len(prof2ave[ii])):
            xralines[ii][i],=axs.plot([],[],color=xrlines[ii][i].get_color(),ls=':')
        for i in range(len(prof2ave[ii])):
            xralines[ii][i].set_xdata(X[0,:])
        axs.plot([],[],color='k',ls=':',label='Time ave.')
        axs.legend(fontsize=12)
    if x=='mode':
         if 'a' in IIsel[ii][0]:
             axs.set_xlabel('\n\n'+r'$M,\, H$',fontdict={'size':12})
             axs.set_ylabel('\n\n'+r'$P$',fontdict={'size':12})
             lnx=np.arange(ss[0]+1)-1.25
             lny=ny+0.25
             nxl=[' nd']+map(str,range(1,ss[0]+1))
             nyl=map(str,range(1,ss[1]+1))

             dx=np.ones_like(z0)*0.3
             dy=dx.copy()
             pk[ii]=axs.bar3d(NXk,NYk,z0,dx,dy,z0+1.e-10,color='b',zsort='max',alpha=0.4,linewidth=0,edgecolor="none")
             pl[ii]=axs.bar3d(NXl,NYl,z0,dx,dy,z0+1.e-10,color='r',zsort='max',alpha=0.6,linewidth=0,edgecolor="none")
             pa[ii]=axs.bar3d(NXa,NYa,z0,dx,dy,z0+1.e-10,color='c',zsort='average',alpha=0.95,linewidth=0,edgecolor="none")
         else:
             axs.set_xlabel('\n\n'+r'$H_{\rm{o}}$',fontdict={'size':12})
             axs.set_ylabel('\n\n'+r'$P_{\rm{o}}$',fontdict={'size':12})
             lnx=nx+0.75
             lny=ny+0.75
             nxl=map(str,range(1,ss[0]+1))
             nyl=map(str,range(1,ss[1]+1))

             dx=np.ones_like(z0)*0.8
             dy=dx.copy()

             pl[ii]=axs.bar3d(NXl,NYl,z0,dx,dy,z0+1.e-10,color='b',zsort='max',alpha=0.7,linewidth=0,edgecolor="none")
         axs.set_xticks(lnx)
         axs.set_xticklabels(nxl)
         axs.set_yticks(lny)
         axs.set_yticklabels(nyl)
         axs.tick_params(axis='both',labelsize=10)
         if 'a' in IIsel[ii][0]:
             axl=fig.add_axes([0.12,0.77-ii*0.45,0.085,0.075],projection='3d')
             axl.set_axis_off()
             axl.view_init(elev=12,azim=98)

             axl.bar3d(1., 0.5, .75,  1., 0.5, 2., color='r',alpha=0.6,linewidth=0,edgecolor="none")
             axl.bar3d(1., 0.5, 6.5,  1., 0.5, 2., color='b',alpha=0.4,linewidth=0,edgecolor="none")
             axl.bar3d(1., 0.5, 12.,  1., 0.5, 2., color='c',alpha=0.95,linewidth=0,edgecolor="none")

             axl.set_ylim(0.,1.0)
             axl.set_xlim(-1.,2.)
             axl.set_zlim(0.,12.5)


             axl.text(0.25,0.6,-.1,'L type',fontdict={'size':12})
             axl.text(0.5,0.,4.4,'K type',fontdict={'size':12})
             axl.text(0.5,0.,10.5,'A type',fontdict={'size':12})

             axb = fig.add_axes([0.12,0.77-ii*0.45,0.085,0.075])
             axb.xaxis.set_visible(False)
             axb.yaxis.set_visible(False)
             axb.set_zorder(1000)
             # axb.patch.set_alpha(0.)
             axb.patch.set_fill(False)
             axb.patch.set_color('k')
    ii+=1


# Attractor plot
if "3D" in Isel:
    axs=axm[Isel.index('3D')]
    x=sete[0]
    axs.plot(x[sd[showp[0]]][n1[0],::ival],x[sd[showp2[0]]][n2[0],::ival],zs=x[sd[showp3[0]]][n3[0],::ival],marker=ms[0],linestyle=ls[0])#,label=fl[i])
    xpoint,=axs.plot([],[],zs=[],marker=',',linestyle='',color='r')#,label=fl[i])

# Spatial plots

im2=ax3.imshow(Z[1][0],interpolation='bilinear', cmap=Zcm[Zsel[1]], origin='lower', extent=[0,2*np.pi/nr,0,np.pi],vmin=mmin[1],vmax=mmax[1]) 
cl2=fig.colorbar(im2,ax=ax3,format=Zformat[Zsel[1]])

im1=ax4.imshow(Z[3][0],interpolation='bilinear', cmap=Zcm[Zsel[3]], origin='lower', extent=[0,2*np.pi/nr,0,np.pi],vmin=mmin[3],vmax=mmax[3])
cl1=fig.colorbar(im1,ax=ax4,format=Zformat[Zsel[3]])

im3=ax5.imshow(Z[0][0],interpolation='bilinear', cmap=Zcm[Zsel[0]], origin='lower', extent=[0,2*np.pi/nr,0,np.pi],vmin=mmin[0],vmax=mmax[0])
cl3=fig.colorbar(im3,ax=ax5,format=Zformat[Zsel[0]])

im0=ax6.imshow(Z[2][0],interpolation='bilinear', cmap=Zcm[Zsel[2]], origin='lower', extent=[0,2*np.pi/nr,0,np.pi],vmin=mmin[2],vmax=mmax[2]) # ,label='year '+str(ny))
cl0=fig.colorbar(im0,ax=ax6,format=Zformat[Zsel[2]])

# im0=ax6.streamplot(X,Y,Uop[0],Vop[0],color=np.sqrt(Uop[0]**2+Vop[0]**2),linewidth=2,cmap=cm.Reds)

# im0=ax6.quiver(X,Y,Uop[0],Vop[0])

# Pruning data not needed.

for i in range(len(sete[0])):
    sete[0][i]=sete[0][i][:,sti:ste:ite]

# Shifting the time

sete[0][4]=sete[0][4]-sete[0][4][0,0]

# Setting tick locator for the fields plots

ax3.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax4.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax5.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax6.xaxis.set_major_locator(ticker.MultipleLocator(1.0))

ax3.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax4.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax5.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax6.yaxis.set_major_locator(ticker.MultipleLocator(1.0))


# Defining the animation update function

def animate(i):
    if i==1 and "3D" in Isel:
        xpoint.set_marker('o')
    if np.mod(i,100)==0:
        print i

    l=i*ival

    if "3D" in Isel or 'diff' in Isel: 
        x=sete[0]  

    if "3D" in Isel: #update of the attractor plot locator
        xpoint.set_data(x[sd[showp[0]]][n1[0],l:l+1],x[sd[showp2[0]]][n2[0],l:l+1])
        xpoint.set_3d_properties(x[sd[showp3[0]]][n3[0],l:l+1])
    
    # Update of the info view
    ii=0
    for z in Isel:
        if z=='diff':
            for j in range(len(diff[ii])):
                xrlines[ii][j].set_xdata(x[4][0,:l+1:ival])
                xrlines[ii][j].set_ydata(diff[ii][j][:l+1:ival])
        if z=='yprof':
            for j in range(len(prof[ii])):
                xrlines[ii][j].set_ydata(prof[ii][j][l])
                xralines[ii][j].set_ydata(profave[ii][j][l])
        if z=='xprof':
            for j in range(len(prof2[ii])):
                xrlines[ii][j].set_ydata(prof2[ii][j][l])
                xralines[ii][j].set_ydata(prof2ave[ii][j][l])
        if z=='mode':
            if pk[ii]:
                update_Poly3D(pk[ii],NXk,NYk,z0,dx,dy,zk[ii][l])
            if pl[ii]:
                update_Poly3D(pl[ii],NXl,NYl,z0,dx,dy,zl[ii][l])
            if pa[ii]:
                update_Poly3D(pa[ii],NXa,NYa,z0,dx,dy,za[ii][l])
        ii+=1
    

    # Actualising the fields plot
    # ax6.cla()
    # # im0=ax6.streamplot(X,Y,Uop[l],Vop[l],color=np.sqrt(Uop[l]**2+Vop[l]**2),linewidth=2,cmap=cm.Reds)
    # im0=ax6.quiver(X,Y,Uop[l],Vop[l])
    im2.set_data(Z[1][l])
    cl2.on_mappable_changed(im2)
    im1.set_data(Z[3][l])
    cl1.on_mappable_changed(im1)
    im3.set_data(Z[0][l])
    cl3.on_mappable_changed(im3)
    im0.set_data(Z[2][l])
    cl0.on_mappable_changed(im0)

    r=[im0,im1,im2,im3]
    if '3D' in Isel:
        r.insert(0,xpoint)
    if 'xprof' in Isel or 'yprof' in Isel:
        r.insert(0,xralines)
    if 'diff'in Isel or 'xprof' in Isel or 'yprof' in Isel:
        r.insert(0,xrlines)
    if 'mode' in Isel:
        r.insert(0,pl)
        r.insert(0,pk)
        r.insert(0,pa)
    
    return r


# Computing the animation

showb=raw_input('Do you want to plot it or to make a movie (p/M) ?')
if not showb:
    showb='M'

lengf=len(ZM[0])

matplotlib.verbose.set_level('debug')

if showb in ['P','p']:
    ani = anim.FuncAnimation(fig, animate, frames=(lengf)/ival, interval=10, blit=False)
    plt.show()
else:
    while True:
        mival=raw_input('Interval between frames in milliseconds ? (default=40)')
        if not mival:
            mival=40
        else:
            mival=int(mival)
        print 'The FPS will be '+str(int(1./(mival/1000.)))+' frames per second.'
        print 'The movie will last '+str(mival*lengf/ival/1.e3)+' seconds.'
        x=raw_input('Do you agree with this ? (y/N)')
        if x in ['y','Y']:
            break
        else:
            print 'No? Ok...'

    # Setting the metadata of the video
    print ' Setting metadata of the movie, please answer the questions:'
    tit=raw_input('Title of the movie?')
    if not tit:
        tit='test' # Title
    auth=raw_input('Author(s)?')
    if not auth:
        auth='test' # Author
    lic=raw_input('License?')
    if not lic:
        lic='MIT' # License
    year=raw_input('Year?')
    if not year:
        year='1982' # Year of production
    comment=raw_input('Comment?')
    if not comment:
        comment='test' # Comments
    meta={'title':tit,'artist':auth,'copyright':lic,'comment':comment,'year':year}

    # Sending a mail to alert that the video encoding has begun
    if fromaddr and toaddr:
        msg = MIMEMultipart()
        msg['From'] = fromaddr
        msg['To'] = toaddr
        msg['Subject'] = "Movie run info"

        body = "The run to generate the video for the geometry:\n\n"
        body += "    atm. "+ageom+" -  oc."+ogeom+"\n\n"
        body += "starts the video encoding! "
        
        msg.attach(MIMEText(body, 'plain'))

        server = smtplib.SMTP(servername)
        text = msg.as_string()
        server.sendmail(fromaddr, toaddr, text)
        server.quit()

# Actual video generation
    ani = anim.FuncAnimation(fig, animate, frames=(lengf)/ival, interval=mival, blit=False)
    
    ssav=raw_input('Output filename ? (default : "out.mp4")')
    if not ssav:
        ssav='out.mp4'
    ani.save(ssav,writer='mencoder',bitrate=None,codec='mpeg4:vbitrate=3000',metadata=meta) #extra_args=['-vf','scale=1024:576'],
        
endt=time.time()



# Sending a final mail to alert that the video generation has ended.

if fromaddr and toaddr:
    msg = MIMEMultipart()
    msg['From'] = fromaddr
    msg['To'] = toaddr
    msg['Subject'] = "End of the movie run"
 
    body = "The run to generate the video for the geometry:\n\n"
    body += "    atm. "+ageom+" -  oc."+ogeom+"\n\n"
    body += "is finished. "

    m, s = divmod(endt-startt,60)
    h, m = divmod(m,60)

    body += "Time elapsed: "+ "%d:%02d:%02d" % (h,m,s) +' .'

    msg.attach(MIMEText(body, 'plain'))
 
    server = smtplib.SMTP(servername)
    text = msg.as_string()
    server.sendmail(fromaddr, toaddr, text)
    server.quit()
