#!/usr/bin/env python
# Code to compute videos of the output of the MAOOAM model

# Copyright :
# 2016 Jonathan Demaeyer.
# See LICENSE.txt for license information.  


# This code needs: 
# - mencoder
# - matplotlib >= 1.3 
# - numpy

# TODO : - Move the parameters at the beginning of the code
#        - Generate frame on the "fly"

# Loading of the libraries

import numpy as np
#import scipy as scp
import matplotlib
# matplotlib.use('Agg')
matplotlib.verbose.set_level('debug')
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
print plt.get_backend()
import matplotlib.animation as anim
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc
rc('font',**{'family':'serif','sans-serif':['Times'],'size':14})

import time
import smtplib
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText
import gzip
import sys

# Defining mail address from where and to which send mails

fromaddr = "jodemaey@gloin.oma.be"
toaddr = "jodemaey@meteo.be"

# Utility functions

#Count the number of line of a file
def linecount(filename):
    lines = 0
    if filename[-3:]=='.gz':
        with gzip.open(filename, "r+") as f:
            for x in f:
                lines += 1
    else:
        with open(filename, "r+") as f:
            for x in f:
                lines += 1
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

# Load the geometries given as arguments or by the users

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
        x={'typ':'A','Nx':0,'Ny':w[1]}
        aftable[ii]=x
        ii+=1
        x={'typ':'K','Nx':w[0],'Ny':w[1]}
        aftable[ii]=x
        ii+=1
        x={'typ':'L','Nx':w[0],'Ny':w[1]}
        aftable[ii]=x
    else:
        ii+=1
        x={'typ':'K','Nx':w[0],'Ny':w[1]}
        aftable[ii]=x
        ii+=1
        x={'typ':'L','Nx':w[0],'Ny':w[1]}
        aftable[ii]=x


oftable={}
ii=0
for w in oms:
    ii+=1
    x={'Nx':w[0]/2.,'Ny':w[1]}
    oftable[ii]=x

# Setting of some general model parameters (those in general do not change)
nr=1.5
al=0.
f0=0.0001032
L=5000000./np.pi
rpr=L**2*f0
RR=287.
RK=rpr*f0/RR
dim=True
at=365.25
ct=(1/(f0*24*3600))/at 
geo=f0/9.81

#Defining some dico to be used later
sd={'psi':0,'theta':1,'A':2,'T':3,'time':4}
vl={'psi':r'\psi_{a,','theta':r'\theta_{a,','A':r'\psi_{o,','T':r'\theta_{o,','time':r't'}
dimd={'psi':rpr*geo,'theta':2*RK,'A':rpr,'T':RK,'time':ct}


# Defining the basis functions and their partial derivatives

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
    return 2*nr*Nx*np.cos(nr*Nx*x)*np.sin(Ny*y)/L #2*np.exp(-al*x)*np.sin(nr*Nx*x)*np.sin(Ny*y)

def dyphi(i,x,y):
    w=oftable[i]
    Nx=w['Nx']
    Ny=w['Ny']
    return 2*Ny*np.sin(nr*Nx*x)*np.cos(Ny*y)/L #2*np.exp(-al*x)*np.sin(nr*Nx*x)*np.sin(Ny*y)

# Function returning the (vector) fields based on the coefficients

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

# Oceanic non-conserved fields (used to compute the temperature field)
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

# Possible legend for the data (not used in general)
#leg=raw_input('Legende?')
leg=''
sl=[s] # Nom du fichier de donnee
fl=[leg]#,'averaged','determinist'] # Legende
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

x1=ct*float(x1.split()[0])
x2=ct*float(x2.split()[0])

# Asking the user from which line to which line he want to read
# providing a gross (experimental) estimation of what it will take in the memory
# Warning : Estimation not accurate for the moment
print 'There is '+str(nlf)+' lines of data in the file'
print 'representing '+str(nlf*(ndim+1)*8/1.e6)+' Mbytes ('+str(nlf*(ndim+1)*8/1.e9)+' GB)'
print 'and '+str((x2-x1)*nlf)+' years of data.'
while True:
    sti='' #raw_input('Where do you want to start reading ? (default=first)')
    ste='2000000' #raw_input('And where do you want to stop ? (default=last)')
    itr='' #raw_input('What interval should we use to read lines ? (default=1)')
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
    print 'It will represent '+str((ste-sti)*(ndim+1)*8*8/1.e6/itr)+' Mbytes of data in the memory. ('+str((ste-sti)*(ndim+1)*8*8/1.e9/itr)+' GB)'
    print 'and '+str((x2-x1)*(ste-sti))+' years of data.'
    x='y'  #raw_input('Do you agree with this ? (y/n)')
    if x=='y':
        break
    else:
        print 'No? Ok...'



# Defining the variables that will be shown in the "attractor" view
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
    if dim:
        sete.append([np.array(psi)*dimd['psi'],np.array(theta)*dimd['theta'],np.array(aa)*dimd['A'],np.array(tt)*dimd['T'],tup*dimd['time'],(np.array(psi)-np.array(theta))*dimd['psi']])
        tu=tu*ct
    else:
        sete.append([np.array(psi),np.array(theta),np.array(aa),np.array(tt),tup,(np.array(psi)-np.array(theta))])

ival=1
tl=[]
for i in range(len(sete)):
    tl.append(len(sete[i][0][0]))

# Setting the plots

fig=plt.figure(num=10,figsize=(16,9))

# General title
fig.suptitle('Spatial fields in the model MAOOAM')
#fig.text(0.42,0.92,r'Parameters : '+r'$\delta = d/f_0 = %.3e $ , $C_o = ' % dd +str(CO)+r'$')

# General subtitle
suptit=r'atmosphere $'+ageom.replace('x','x$-$')+r'y$  ocean $'+ogeom.replace('x','x$-$')+r'y$'
fig.text(0.42,0.92,'Resolution : '+suptit)

# Setting the six views
ax1=fig.add_subplot(2,3,1)
ax2=fig.add_subplot(2,3,4,projection='3d')
ax3=fig.add_subplot(2,3,2)
ax4=fig.add_subplot(2,3,3)
ax5=fig.add_subplot(2,3,5)
ax6=fig.add_subplot(2,3,6)

# Views title

ax1.set_title("Geopotential height difference (m)")
ax2.set_title('3-D phase space projection')
ax2.text2D(0.5, 0.9,'(Non-dimensional units)', horizontalalignment='center',verticalalignment='center',transform=ax2.transAxes,fontdict={'size':12})
ax3.set_title("Atm. Temperature ($^\circ\!$C)")
ax4.set_title('Ocean Temperature ($^\circ\!$C)') # You can choose between ocean temp. or ocean-atm. temp. difference
# ax4.set_title('Oc.-Atms Temperature diff.')
ax5.set_title(r'Geopotential height $\psi_a f_0/g$ (m)',fontdict={'size':14})
ax6.set_title(r'Ocean $\psi_o$ (m$^2$s$^{-1}$)')#,fontdict={'size':14})

# Attractor axis ranges and labels

# Range
#ax.view_init(20,24)
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
#print smin,smax
ax2.set_xlim(smin[0],smax[0])
ax2.set_ylim(smin[1],smax[1])
ax2.set_zlim(smin[2],smax[2])
fig.canvas.draw()

# Setting the ticks and the labels on the attractor view

# x ticks and axis label

labels = [item.get_text() for item in ax2.get_xticklabels()]
ii=len(labels)-1
labto=[]
#	print labels
jk=0
for x in labels:
    if x:
        jk+=1
        y=float(x.replace(u'\u2212',u'-'))
        y=y/dimd[showp[0]]
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
ax2.set_xticklabels(labto)
til1=ax2.xaxis.get_major_ticks()
for j in range(0,len(til1),1):
    til1[j].label1.set_visible(False)
til1[ii].label1.set_visible(True)
til1[0].label1.set_visible(True)

if showp[0]=='time':
    ax2.set_xlabel('\n\n'+r'$'+vl[showp[0]]+r'$',fontdict={'size':18})
else:
#    if abs(n)>2:
#        ax2.set_xlabel('\n\n'+r'$'+vl[showp[0]]+str(n1[0]+1)+r'}$ ${(\times 10^{'+str(-n)+r'})}$',fontdict={'size':18})
#    else:
    ax2.set_xlabel('\n\n\n'+r'$'+vl[showp[0]]+str(n1[0]+1)+r'}$',fontdict={'size':18})

# y ticks and axis label

labels = [item.get_text() for item in ax2.get_yticklabels()]
ii=len(labels)-1
labto=[]
jk=0
for x in labels:
    if x:
        jk+=1
        y=float(x.replace(u'\u2212',u'-'))
        y=y/dimd[showp2[0]]
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
ax2.set_yticklabels(labto)
til2=ax2.yaxis.get_major_ticks()
for j in range(0,len(til2),1):
    til2[j].label1.set_visible(False)
til2[ii].label1.set_visible(True)
til2[0].label1.set_visible(True)

#if abs(n)>2:
#    ax2.set_ylabel('\n\n'+r'$'+vl[showp2[0]]+str(n2[0]+1)+r'}$ ${(\times 10^{'+str(-n)+r'})}$',fontdict={'size':18})
#else:
ax2.set_ylabel('\n\n'+r'$'+vl[showp2[0]]+str(n2[0]+1)+r'}$',fontdict={'size':18})

# z ticks

labels = [item.get_text() for item in ax2.get_zticklabels()]
ii=len(labels)-1
labto=[]
jk=0
for x in labels:
    if x:
        jk+=1
        y=float(x.replace(u'\u2212',u'-'))
        y=y/dimd[showp3[0]]
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
ax2.set_zticklabels(labto)
til3=ax2.zaxis.get_major_ticks()
for j in range(0,len(til3),1):
    til3[j].label1.set_visible(False)
til3[ii-1].label1.set_visible(True)
til3[0].label1.set_visible(True)

#if abs(n)>2:
#    ax2.set_zlabel(''+r'$'+vl[showp3[0]]+str(n3[0]+1)+r'}$ ${(\times 10^{'+str(-n)+r'})}$',fontdict={'size':18})    
#else:
ax2.set_zlabel(''+r'$'+vl[showp3[0]]+str(n3[0]+1)+r'}$',fontdict={'size':18})


# ticks alignement

[t.set_va('center') for t in ax2.get_yticklabels()]
[t.set_ha('left') for t in ax2.get_yticklabels()]
[t.set_va('center') for t in ax2.get_xticklabels()]
[t.set_ha('right') for t in ax2.get_xticklabels()]
[t.set_va('center') for t in ax2.get_zticklabels()]
[t.set_ha('left') for t in ax2.get_zticklabels()]

# other views labels

ax3.set_xlabel('$x^\prime$')
ax3.set_ylabel('$y^\prime$')

ax4.set_xlabel('$x^\prime$')
ax4.set_ylabel('$y^\prime$')

ax5.set_xlabel('$x^\prime$')
ax5.set_ylabel('$y^\prime$')

ax6.set_xlabel('$x^\prime$')
ax6.set_ylabel('$y^\prime$')

# Sending a mail to alert that the first stage is completed (loading the data
# and preparing the plots)

msg = MIMEMultipart()
msg['From'] = fromaddr
msg['To'] = toaddr
msg['Subject'] = "Movie run info"

body = "The run to generate the video for the geometry:\n\n"
body += "    atm. "+ageom+" -  oc."+ogeom+"\n\n"
body += "has finished loading data! "

msg.attach(MIMEText(body, 'plain'))

server = smtplib.SMTP('localhost')
text = msg.as_string()
server.sendmail(fromaddr, toaddr, text)
server.quit()


# Actually computing the fields

# Setting the grids
delta='' #raw_input('Space between points on the grid (default = 0.025) ?') # 
if not delta:
    delta=0.025
else:
    delta=float(delta)
x = np.arange(0., 2*np.pi/nr, delta)
y = np.arange(0., np.pi, delta)
X, Y = np.meshgrid(x, y)
sh=X.shape

# List holding the fields and their minima and maxima
Zot=[]
Uop=[]
Vop=[]
Zop=[]
Zat=[]
Zap=[]
Zp3=[]
Zdt=[]
mZot=[]
mZop=[]
mZat=[]
mZap=[]
mZp3=[]
mZdt=[]
MZot=[]
MZop=[]
MZat=[]
MZap=[]
MZp3=[]
MZdt=[]
geoap=[]


mmin=[0.,0.,0.,0.,0.,0.]
mmax=[0.,0.,0.,0.,0.,0.]

# Setting the number of frames and the time of the first and the last
#tl[0]=tl[0]/20
#tl[0]=int(tl[0]*0.9)
while True:
    print 'Total number of frames:',tl[0]
    if dim:
        print 'Time between each frame is '+str((tu[1]-tu[0])*at)+' days'
    else:
        print 'Time between each frame is '+str(tu[1]-tu[0])+' timeunit'
    sti=str(nlf/2) #raw_input('Start at frame (default = first) ?')
    mti=min(nlf/2-1,15000)
    ste=str(mti+nlf/2) #raw_input('End at frame (default = last) ?')
    ite='3' #raw_input('Interval (default = 1) ?')
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
    print 'It will take '+str((ste-sti)*len(mmin)*sh[0]*sh[1]*8/(1.e6*ite))+' Mbytes in the memory! ('+str((ste-sti)*len(mmin)*sh[0]*sh[1]*8/(1.e9*ite))+' GB)'
    x='y' #raw_input('Do you agree (y/n) ?')
    if x=='y':
        break
    else:
        print "Let's start again..."

# Loop generating the frame (computing the fields)
startt=time.time()
x=sete[0]
for i in range(sti,ste,ite):
    if np.mod(i-sti,100*ite)==0:
        print 'Generating the fields in the frame ',i,'('+str((i-sti)/ite)+')'
        print 'At time t=',x[4][0,i],'years'
    Zop.append(ostream_cons(X,Y,x[2][:,i]))
    # U,V=ovec(X,Y,x[2][:,i])
    # Uop.append(U)
    # Vop.append(V)
    Zot.append(ostream(X,Y,x[3][:,i]))
    Zat.append(astream(X,Y,x[1][:,i]))
    Zap.append(astream(X,Y,x[0][:,i]))

    # Zp3.append(astream(X,Y,x[5][:,i]))
    # Zdt.append(Zot[-1]-Zat[-1])

    MZop.append(np.amax(Zop[-1]))
    mZop.append(np.amin(Zop[-1]))

    # MZop.append(np.amax(np.sqrt(Uop[-1]**2+Vop[-1]**2)))
    # mZop.append(np.amin(np.sqrt(Uop[-1]**2+Vop[-1]**2)))

    MZot.append(np.amax(Zot[-1]))
    mZot.append(np.amin(Zot[-1]))
    MZap.append(np.amax(Zap[-1]))
    mZap.append(np.amin(Zap[-1]))
    MZat.append(np.amax(Zat[-1]))
    mZat.append(np.amin(Zat[-1]))
    geoap.append(geodiff(x[0][:,i]))
    # MZp3.append(np.amax(Zp3[-1]))
    # mZp3.append(np.amin(Zp3[-1]))
    # MZdt.append(np.amax(Zdt[-1]))
    # mZdt.append(np.amin(Zdt[-1]))

    
    mmax[0]=max(mmax[0],MZop[-1])
    mmin[0]=min(mmin[0],mZop[-1])
    mmax[1]=max(mmax[1],MZot[-1])
    mmin[1]=min(mmin[1],mZot[-1])
    mmax[2]=max(mmax[2],MZat[-1])
    mmin[2]=min(mmin[2],mZat[-1])
    mmax[3]=max(mmax[3],MZap[-1])
    mmin[3]=min(mmin[3],mZap[-1])
    # mmax[4]=max(mmax[4],MZp3[-1])
    # mmin[4]=min(mmin[4],mZp3[-1])
    # mmax[5]=max(mmax[5],MZdt[-1])
    # mmin[5]=min(mmin[5],mZdt[-1])
    
    

mZot=np.array(mZot)
mZop=np.array(mZop)
mZat=np.array(mZat)
mZap=np.array(mZap)
# mZp3=np.array(mZp3)
# mZdt=np.array(mZdt)
MZot=np.array(MZot)
MZop=np.array(MZop)
MZat=np.array(MZat)
MZap=np.array(MZap)
geoap=np.array(geoap)
# MZp3=np.array(MZp3)
# MZdt=np.array(MZdt)

dZot=MZot-mZot
dZop=MZop-mZop
dZap=MZap-mZap
dZat=MZat-mZat
# dZp3=MZp3-mZp3

# smax=max(np.amax(dZp3),np.amax(dZap))
# smin=min(np.amin(dZp3),np.amin(dZap))
# smax=np.amax(dZap)
# smin=np.amin(dZap)
smax=np.amax(geoap)
smin=np.amin(geoap)

# Now that the geopential height difference have been computed
 # we can set the limit of the view showing it
ax1.set_xlim(0.,sete[0][4][0,ste-1]-sete[0][4][0,sti])
# ax1.set_xticks([0.,40000.,80000.])
if smin>0.:
    ax1.set_ylim(0.9*smin,1.1*smax)
else:
    ax1.set_ylim(1.1*smin,1.1*smax)
ax1.set_xlabel('time (years)')
#ax1.set_ylabel('Geopotential height (m)')

diff=[dZot,dZat,dZop,geoap,dZap]#,dZp3]
difflab=['Ocean temp.','Atm. temp.',r'Ocean $\Psi$',r'Atm. $\Psi$',r'$\Psi_3$']

# Initializing the animation


# Geopotential height difference plot

lx=5
xrlines=[]
for i in range(lx):
    xrlines.append(None)
for i in [3]: #,4
    xrlines[i],=ax1.plot([],[],label=difflab[i])

# ax1.legend()

# Attractor plot

x=sete[0]
ax2.plot(x[sd[showp[0]]][n1[0],::ival],x[sd[showp2[0]]][n2[0],::ival],zs=x[sd[showp3[0]]][n3[0],::ival],marker=ms[0],linestyle=ls[0])#,label=fl[i])
xpoint,=ax2.plot([],[],zs=[],marker=',',linestyle='',color='r')#,label=fl[i])

# Atm. Temperature plot

im2=ax3.imshow(Zat[0],interpolation='bilinear', cmap=cm.coolwarm, origin='lower', extent=[0,2*np.pi/nr,0,np.pi],vmin=mmin[2],vmax=mmax[2]) # ,label='year '+str(ny))
cl2=fig.colorbar(im2,ax=ax3)#,format=ticker.FuncFormatter(fmt))

# Ocean Temperature plot

im1=ax4.imshow(Zot[0],interpolation='bilinear', cmap=cm.coolwarm, origin='lower', extent=[0,2*np.pi/nr,0,np.pi],vmin=mmin[1],vmax=mmax[1]) # ,label='year '+str(ny))
# im1=ax4.imshow(Zdt[0],interpolation='bilinear', cmap=cm.coolwarm, origin='lower', extent=[0,2*np.pi/nr,0,np.pi],vmin=mmin[5],vmax=mmax[5]) # ,label='year '+str(ny))
cl1=fig.colorbar(im1,ax=ax4)#,format=ticker.FuncFormatter(fmt))

# Atmospheric stream plot

im3=ax5.imshow(Zap[0],interpolation='bilinear', cmap=cm.gist_rainbow_r, origin='lower', extent=[0,2*np.pi/nr,0,np.pi],vmin=mmin[3],vmax=mmax[3]) # ,label='year '+str(ny))
cl3=fig.colorbar(im3,ax=ax5)#,format=ticker.FuncFormatter(fmt))

# Ocean stream plot

im0=ax6.imshow(Zop[0],interpolation='bilinear', cmap=cm.gist_rainbow_r, origin='lower', extent=[0,2*np.pi/nr,0,np.pi],vmin=mmin[0],vmax=mmax[0]) # ,label='year '+str(ny))
cl0=fig.colorbar(im0,ax=ax6,format=ticker.FuncFormatter(fmt))

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
    if i==1:
        xpoint.set_marker('o')
    if np.mod(i,100)==0:
        print i
    l=i*ival
    x=sete[0]  #update of the attractor plot locator
    xpoint.set_data(x[sd[showp[0]]][n1[0],l:l+1],x[sd[showp2[0]]][n2[0],l:l+1])
    xpoint.set_3d_properties(x[sd[showp3[0]]][n3[0],l:l+1])
    
    # Update of the geopotential plot
    for j in [3]: #,4
        xrlines[j].set_xdata(x[4][0,:l+1:ival])
        xrlines[j].set_ydata(diff[j][:l+1:ival])

    # Actualising the fields plot
    # ax6.cla()
    # # im0=ax6.streamplot(X,Y,Uop[l],Vop[l],color=np.sqrt(Uop[l]**2+Vop[l]**2),linewidth=2,cmap=cm.Reds)
    # im0=ax6.quiver(X,Y,Uop[l],Vop[l])
    im0.set_data(Zop[l])
    cl0.on_mappable_changed(im0)
    im1.set_data(Zot[l])
    cl1.on_mappable_changed(im1)
    im2.set_data(Zat[l])
    cl2.on_mappable_changed(im2)
    im3.set_data(Zap[l])
    cl3.on_mappable_changed(im3)
       
    return xrlines,xpoint,im0,im1,im2,im3#,xilines

# Computing the animation

# showb=raw_input('Do you want to plot it or to make a movie (P/M) ?')
showb='M'

lengf=len(MZot)


if showb=='P' or showb=='p':
    ani = anim.FuncAnimation(fig, animate, frames=(lengf)/ival, interval=10, blit=False)
    plt.show()
else:
    while True:
        # mival=raw_input('Interval between frames in milliseconds ? (default=40)')
        mival=40
        if not mival:
            mival=40
        else:
            mival=int(mival)
        print 'The FPS will be '+str(int(1./(mival/1000.)))+' frames per second.'
        print 'The movie will last '+str(mival*lengf/ival/1.e3)+' seconds.'
        # x=raw_input('Do you agree with this ? (y/n)')
        x='y'
        if x=='y':
            break
        else:
            print 'No? Ok...'

    # Setting the metadata of the video
    print ' Setting metadata of the movie, please answer the questions:'
    #tit=raw_input('Title of the movie?')
    tit='test' # Title
    #auth=raw_input('Author(s)?')
    auth='test' # Author
    lic='test' # License
    #year=raw_input('Year?')
    year='1980' # Year of production
    # comment=raw_input('Comment?')
    comment='test' # Comments
    meta={'title':tit,'artist':auth,'copyright':lic,'comment':comment,'year':year}

    # Sending a mail to alert that the video encoding has begun

    msg = MIMEMultipart()
    msg['From'] = fromaddr
    msg['To'] = toaddr
    msg['Subject'] = "Movie run info"

    body = "The run to generate the video for the geometry:\n\n"
    body += "    atm. "+ageom+" -  oc."+ogeom+"\n\n"
    body += "starts the video encoding! "

    msg.attach(MIMEText(body, 'plain'))

    server = smtplib.SMTP('localhost')
    text = msg.as_string()
    server.sendmail(fromaddr, toaddr, text)
    server.quit()

# Actual video generation
    ani = anim.FuncAnimation(fig, animate, frames=(lengf)/ival, interval=mival, blit=False)
    
    # ssav=raw_input('Output filename ? (default : "out_d'+str(delt)+'_Co'+str(CO)+'.mp4")')
    ssav='out.mp4'
    # if not ssav:
    #     ssav='out_d'+str(delt)+'_Co'+str(CO)+'.mp4'
    ani.save(ssav,writer='mencoder',bitrate=None,codec='mpeg4:vbitrate=3000',metadata=meta) #extra_args=['-vf','scale=1024:576'],
        
endt=time.time()



# Sending a final mail to alert that the video generation has ended.

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
 
server = smtplib.SMTP('localhost')
text = msg.as_string()
server.sendmail(fromaddr, toaddr, text)
server.quit()
