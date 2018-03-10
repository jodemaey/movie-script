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

import geometry as geom
from geofunctions import compute_frame,geodiff
from util import *
from params import *
params_initialize()
#-----------------------------------------------------------
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

geom.set_geometry(ageom,ogeom)


#----------------------------------------
# Opening files
#-----------------------------------------

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
sl=[s] # Data filename
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
e=evol[0]
    
x1=e.readline()
x2=e.readline()

x1=float(x1.split()[0])
x2=float(x2.split()[0])

dt=x2-x1
dty=(x2-x1)*dimd['timey']

#------------------------------------------------
# Data for the 3D mode (if selected)
#------------------------------------------------

if "3D" in view.Isel:

    # Asking the user from which line to which line he want to read
    # providing a gross (experimental) estimation of what it will take in the memory
    # Warning : Estimation not accurate for the moment

    print 'Selecting the data for the 3D view'
    print "----------------------------------"
    print ""
    print 'There are '+str(nlf)+' lines of data in the file'
    print 'representing '+str(dty*nlf)+' years of data.'
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
        print 'It will represent '+str((ste-sti)*(3+1)*8/1.e6/itr)+' Mbytes of data in the memory. ('+str((ste-sti)*(3+1)*8/1.e9/itr)+' GB)'
        print 'and '+str(dty*(ste-sti))+' years of data.'
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

    dserie=[]
    for j in range(len(evol)):
        e=evol[j]
        e.seek(0)
        tu=[]
        aa=[]
        for i in range(3):
            aa.append([])

        ii=0
        for x in e:
            if ii>=sti and ii<=ste and np.mod(ii-sti,itr)==0:
                y=x.split()
                tu.append(float(y[sdi['time']]))
                aa[0].append(float(y[n1[j]+sdi[showp[j]]]))
                aa[1].append(float(y[n2[j]+sdi[showp2[j]]]))
                aa[2].append(float(y[n3[j]+sdi[showp3[j]]]))
            if ii>ste:
                break
            ii+=1

        tu=np.array(tu)
        tup=tu.copy()
        tup.shape=1,len(tu)
        dserie.append([np.array(aa),tup*dimdv['time']])

    # Sending a mail to alert that the data are loaded

    if fromaddr and toaddr:
        msg = MIMEMultipart()
        msg['From'] = fromaddr
        msg['To'] = toaddr
        msg['Subject'] = "Movie run info"

        body = "The run to generate the video for the geometry:\n\n"
        body += "    atm. "+ageom+" -  oc."+ogeom+"\n\n"
        body += "has finished loading data for the 3D view! "

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
    print 'Total number of possible frames:',nlf
    if dim:
        print 'Time between each frame is '+str(dt*at)+' days'
    else:
        print 'Time between each frame is '+str(dt)+' timeunit'
    sti=raw_input('Start at frame (default = first) ?')
    ste=raw_input('End at frame (default = last) ?')
    ite=raw_input('Interval (default = 1) ?')
    if not sti:
        sti=0
    else:
        sti=int(sti)
    if not ste:
        ste=nlf-1
    else:
        ste=int(ste)
    if not ite:
        ite=1
    else:
        ite=int(ite)
    if dim:
        print 'Time between each frame is now '+str(ite*dt*at)+' days'
    else:
        print 'Time between each frame is '+str(ite*dt)+' timeunits'
    print 'Number of frames that will effectively be computed :'+str((ste-sti)/ite)
    if dim:
        print 'for a total time of '+str((ste-sti)*dt)+' years'
    else:
        print 'for a total time of '+str((ste-sti)*dt)+' timeunits'

    x=raw_input('Do you agree (y/N) ?')
    if x in ['y','Y']:
        break
    else:
        print "Let's start again..."

#------------------------------------------------------------
# First loop : Computing the fields extremums and derivatives
#------------------------------------------------------------

startt=time.time()

print "Computing fields extermums..."
print "This may take a while..."

# Preparing space to store the results

# spatial field + minmax
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

NXa=[None,None]
NYa=[None,None]

NXk=[None,None]
NYk=[None,None]

NXl=[None,None]
NYl=[None,None]

z0=[None,None]

ss=[None,None]

if 'mode' in view.Isel:
    ii=0
    for z in view.Isel:
        if z=='mode':
            if 'a' in view.IIsel[ii][0]:
                ss[ii]=ass
            else:
                ss[ii]=oss
            z0[ii]=np.zeros(ss[ii]).flatten()
            nx=np.arange(ss[ii][0])-0.5
            ny=np.arange(ss[ii][1])
            NXa[ii], NYa[ii] = np.meshgrid(nx-1.,ny)
            NXa[ii]=NXa[ii].flatten()
            NYa[ii]=NYa[ii].flatten()
            NXk[ii], NYk[ii] = np.meshgrid(nx,ny)
            NXk[ii]=NXk[ii].flatten()
            NYk[ii]=NYk[ii].flatten()
            NXl[ii], NYl[ii] = np.meshgrid(nx+0.5,ny+0.5)
            NXl[ii]=NXl[ii].flatten()
            NYl[ii]=NYl[ii].flatten()
        ii+=1

#overall fields max and min
mmin=np.zeros((4))
mmax=np.zeros((4))

#infoviews max and min
ifsmax=[None,None]
ifsmin=[None,None]

e=evol[0]
e.seek(0)

for i,line in enumerate(e):
    if i>=sti and np.mod(i-sti,ite)==0:

        if np.mod(i-sti,100*ite)==0:
            print 'Probing the fields in the frame ',i,'('+str((i-sti)/ite)+')'
            z=float(line.split()[0])
            if dim:
                print 'At time t=',z*dimdv['time'],'years'
            else:
                print 'At time t=',z,'timeunits'

        Z=compute_frame(line)

        if 'diff' in view.Isel:
            if 'geo' in view.IIsel[view.Isel.index('diff')]:
                geoap.append(geodiff(psi))

        for j in range(4):
            ZM[j].append(np.amax(Z[j]))
            Zm[j].append(np.amin(Z[j]))

        for j in range(4):
            mmax[j]=max(mmax[j],ZM[j])
            mmin[j]=min(mmin[j],Zm[j])

        if 'yprof' in view.Isel:
            for j in range(4):
                yprof[j].append(np.mean(Z[j],axis=1))
                yprofmid[j].append(Z[j][:,sh[1]/2])
                if i==sti:
                    yprofmidave[j].append(yprofmid[j][-1])
                    yprofave[j].append(yprof[j][-1])
                else:
                    y=yprofmidave[j][-1]+(yprofmid[j][-1]-yprofmidave[j][-1])/(i-sti)
                    yprofmidave[j].append(y)
                    y=yprofave[j][-1]+(yprof[j][-1]-yprofave[j][-1])/(i-sti)
                    yprofave[j].append(y)
        if 'xprof' in view.Isel:
            for j in range(4):
                xprof[j].append(np.mean(Z[j],axis=0))
                xprofmid[j].append(Z[j][sh[0]/2,:])
                if i==sti:
                    xprofmidave[j].append(xprofmid[j][-1])
                    xprofave[j].append(xprof[j][-1])
                else:
                    y=xprofmidave[j][-1]+(xprofmid[j][-1]-xprofmidave[j][-1])/(i-sti)
                    xprofmidave[j].append(y)
                    y=xprofave[j][-1]+(xprof[j][-1]-xprofave[j][-1])/(i-sti)
                    xprofave[j].append(y)
        if "mode" in view.Isel:
            iii=0
            for z in view.Isel:
                if z=='mode':
                    za[iii].append(np.zeros(ss[iii]))
                    zk[iii].append(np.zeros(ss[iii]))
                    zl[iii].append(np.zeros(ss[iii]))
                    y=np.absolute(x[sd[sdd[view.IIsel[iii][0]]]][:,i])
                    y=100*y/y.sum()
                    for ii in range(1,len(y)+1):
                        if 'a' in view.IIsel[iii][0]:
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
    if i>=ste:
        break

    
ZM=np.array(ZM)
Zm=np.array(Zm)
if 'diff' in view.Isel:
    if 'geo' in view.IIsel[view.Isel.index('diff')]:
        geoap=np.array(geoap)

if "diff" in view.Isel:
    diff=[[],[]]
    difflab=[[],[]]

if "yprof" in view.Isel:
    prof=[[],[]]
    profave=[[],[]]
    proflab=[[],[]]

if "xprof" in view.Isel:
    prof2=[[],[]]
    prof2ave=[[],[]]
    prof2lab=[[],[]]

ii=0
for z in view.Isel:
    if z=="diff":
        for x in view.IIsel[ii]:
            if 'sp' in x:
                i=int(x[2])-1
                diff[ii].append(ZM[i]-Zm[i])
                difflab[ii].append(Zlabmini[i])
                if dim:
                    difflab[ii][-1]+=Zlabelunit[view.Zsel[i]]
            if x=='geo':
                diff[ii].append(geoap)
                difflab[ii].append('Geop. H.')
                if dim:
                    difflab[ii][-1]+=strg

        ifsmax[ii]=max(map(np.amax,diff[ii]))
        ifsmin[ii]=min(map(np.amin,diff[ii]))
    if z=="yprof":
        for x in view.IIsel[ii]:
            if 'spa' in x:
                i=int(x[3])-1
                prof[ii].append(yprof[i])
                profave[ii].append(yprofave[i])
                proflab[ii].append('Z.A. '+Zlabmini[i])
                if dim:
                    proflab[ii][-1]+=Zlabelunit[view.Zsel[i]]
            elif 'sp' in x:
                i=int(x[2])-1
                prof[ii].append(yprofmid[i])
                profave[ii].append(yprofmidave[i])
                proflab[ii].append(Zlabmini[i])
                if dim:
                    proflab[ii][-1]+=Zlabelunit[view.Zsel[i]]
        ifsmax[ii]=[]
        ifsmin[ii]=[]
        for x in prof[ii]:
            ifsmax[ii].append(max(map(np.amax,x)))
            ifsmin[ii].append(min(map(np.amin,x)))
        ifsmax[ii]=max(ifsmax[ii])
        ifsmin[ii]=min(ifsmin[ii])
    if z=="xprof":
        for x in view.IIsel[ii]:
            if 'spa' in x:
                i=int(x[3])-1
                prof2[ii].append(xprof[i])
                prof2ave[ii].append(xprofave[i])
                prof2lab[ii].append('Z.A. '+Zlabmini[i])
                if dim:
                    prof2lab[ii][-1]+=Zlabelunit[view.Zsel[i]]
            elif 'sp' in x:
                i=int(x[2])-1
                prof2[ii].append(xprofmid[i])
                prof2ave[ii].append(xprofmidave[i])
                prof2lab[ii].append(Zlabmini[i])
                if dim:
                    prof2lab[ii][-1]+=Zlabelunit[view.Zsel[i]]
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

if view.Isel[0] in ["3D","mode"]:
    ax1=fig.add_subplot(2,3,1,projection='3d')
else:
    ax1=fig.add_subplot(2,3,1)
if view.Isel[1] in ["3D","mode"]:
    ax2=fig.add_subplot(2,3,4,projection='3d')
else:
    ax2=fig.add_subplot(2,3,4)
ax3=fig.add_subplot(2,3,2)
ax4=fig.add_subplot(2,3,3)
ax5=fig.add_subplot(2,3,5)
ax6=fig.add_subplot(2,3,6)

# Views title

ax1.set_title(Ivtit[view.Isel[0]])
if view.Isel[0]=="mode":
    ax1.text2D(0.5, 0.9,Zlabelmini[view.IIsel[0][0]], horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes,fontdict={'size':12})
if view.Isel[0]=="3D":
    ax1.text2D(0.5, 0.9,'(Non-dimensional units)', horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes,fontdict={'size':12})
ax2.set_title(Ivtit[view.Isel[1]])
if view.Isel[1]=="mode":
    ax2.text2D(0.5, 0.9,Zlabelmini[view.IIsel[1][0]], horizontalalignment='center',verticalalignment='center',transform=ax2.transAxes,fontdict={'size':12})
if view.Isel[1]=="3D":
    ax2.text2D(0.5, 0.9,'(Non-dimensional units)', horizontalalignment='center',verticalalignment='center',transform=ax2.transAxes,fontdict={'size':12})

ax3.set_title(Zlab[1])
ax4.set_title(Zlab[3]) 
ax5.set_title(Zlab[0])
ax6.set_title(Zlab[2])

axm=[ax1,ax2]

# 3D view axis ranges and labels
#----------------------------------

# Range
if "3D" in view.Isel:
    axs=axm[view.Isel.index('3D')]
    mv=0.25 # Overflow factor
    smax=3*[-30000.]
    x=dserie[0]
    smax[0]=max(smax[0],np.amax(x[0][0,:]))
    smax[1]=max(smax[1],np.amax(x[0][1,:]))
    smax[2]=max(smax[2],np.amax(x[0][2,:]))
    smin=smax[:]
    smin[0]=min(smin[0],np.amin(x[0][0,:]))
    smin[1]=min(smin[1],np.amin(x[0][1,:]))
    smin[2]=min(smin[2],np.amin(x[0][2,:]))
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
for x in view.Isel:
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
dx=[None,None]
dy=[None,None]
for x in view.Isel:
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
         if 'a' in view.IIsel[ii][0]:
             axs.set_xlabel('\n\n'+r'$M,\, H$',fontdict={'size':12})
             axs.set_ylabel('\n\n'+r'$P$',fontdict={'size':12})
             lnx=np.arange(ss[ii][0]+1)-1.25
             lny=ny+0.25
             nxl=[' nd']+map(str,range(1,ss[ii][0]+1))
             nyl=map(str,range(1,ss[ii][1]+1))

             dx[ii]=np.ones_like(z0[ii])*0.3
             dy[ii]=dx[ii].copy()
             pk[ii]=axs.bar3d(NXk[ii],NYk[ii],z0[ii],dx[ii],dy[ii],z0[ii]+1.e-10,color='b',zsort='max',alpha=0.4,linewidth=0,edgecolor="none")
             pl[ii]=axs.bar3d(NXl[ii],NYl[ii],z0[ii],dx[ii],dy[ii],z0[ii]+1.e-10,color='r',zsort='max',alpha=0.6,linewidth=0,edgecolor="none")
             pa[ii]=axs.bar3d(NXa[ii],NYa[ii],z0[ii],dx[ii],dy[ii],z0[ii]+1.e-10,color='c',zsort='average',alpha=0.95,linewidth=0,edgecolor="none")
         else:
             axs.set_xlabel('\n\n'+r'$H_{\rm{o}}$',fontdict={'size':12})
             axs.set_ylabel('\n\n'+r'$P_{\rm{o}}$',fontdict={'size':12})
             lnx=nx+0.75
             lny=ny+0.75
             nxl=map(str,range(1,ss[ii][0]+1))
             nyl=map(str,range(1,ss[ii][1]+1))

             dx[ii]=np.ones_like(z0[ii])*0.8
             dy[ii]=dx[ii].copy()

             pl[ii]=axs.bar3d(NXl[ii],NYl[ii],z0[ii],dx[ii],dy[ii],z0[ii]+1.e-10,color='b',zsort='max',alpha=0.7,linewidth=0,edgecolor="none")
         axs.set_xticks(lnx)
         axs.set_xticklabels(nxl)
         axs.set_yticks(lny)
         axs.set_yticklabels(nyl)
         axs.tick_params(axis='both',labelsize=10)
         if 'a' in view.IIsel[ii][0]:
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
if "3D" in view.Isel:
    axs=axm[view.Isel.index('3D')]
    x=dserie[0]
    axs.plot(x[0][0,:],x[0][1,:],zs=x[0][2,:],marker=ms[0],linestyle=ls[0])#,label=fl[i])
    xpoint,=axs.plot([],[],zs=[],marker=',',linestyle='',color='r')#,label=fl[i])


# Spatial plots

e=evol[0]
e.seek(0)

for i, line in enumerate(e):
    if i==sti:
        break

Z=compute_frame(line)

im2=ax3.imshow(Z[1],interpolation='bilinear', cmap=Zcm[view.Zsel[1]], origin='lower', extent=[0,2*np.pi/nr,0,np.pi],vmin=mmin[1],vmax=mmax[1]) 
cl2=fig.colorbar(im2,ax=ax3,format=Zformat[view.Zsel[1]])

im1=ax4.imshow(Z[3],interpolation='bilinear', cmap=Zcm[view.Zsel[3]], origin='lower', extent=[0,2*np.pi/nr,0,np.pi],vmin=mmin[3],vmax=mmax[3])
cl1=fig.colorbar(im1,ax=ax4,format=Zformat[view.Zsel[3]])

im3=ax5.imshow(Z[0],interpolation='bilinear', cmap=Zcm[view.Zsel[0]], origin='lower', extent=[0,2*np.pi/nr,0,np.pi],vmin=mmin[0],vmax=mmax[0])
cl3=fig.colorbar(im3,ax=ax5,format=Zformat[view.Zsel[0]])

im0=ax6.imshow(Z[2],interpolation='bilinear', cmap=Zcm[view.Zsel[2]], origin='lower', extent=[0,2*np.pi/nr,0,np.pi],vmin=mmin[2],vmax=mmax[2]) # ,label='year '+str(ny))
cl0=fig.colorbar(im0,ax=ax6,format=Zformat[view.Zsel[2]])

# im0=ax6.streamplot(X,Y,Uop[0],Vop[0],color=np.sqrt(Uop[0]**2+Vop[0]**2),linewidth=2,cmap=cm.Reds)

# im0=ax6.quiver(X,Y,Uop[0],Vop[0])

# Setting tick locator for the fields plots

ax3.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax4.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax5.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax6.xaxis.set_major_locator(ticker.MultipleLocator(1.0))

ax3.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax4.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax5.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax6.yaxis.set_major_locator(ticker.MultipleLocator(1.0))

#Pruning unneeded data

for i in range(len(dserie[0])):
    dserie[0][i]=dserie[0][i][:,sti:ste+1:ite]

print len(dserie[0][i][1,:])
print len(ZM[0])

# Defining the animation update function

def animate(l):
    if l==1 and "3D" in view.Isel:
        xpoint.set_marker('o')
    if np.mod(l,100)==0:
        print l

    if "3D" in view.Isel or 'diff' in view.Isel: 
        x=dserie[0]  

    if "3D" in view.Isel: #update of the attractor plot locator
        xpoint.set_data(x[0][0,l:l+1],x[0][1,l:l+1])
        xpoint.set_3d_properties(x[0][2,l:l+1])
    
    # Update of the info view
    ii=0
    for z in view.Isel:
        if z=='diff':
            for j in range(len(diff[ii])):
                xrlines[ii][j].set_xdata(x[1][0,:l+1:ival])
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
                update_Poly3D(pk[ii],NXk[ii],NYk[ii],z0[ii],dx[ii],dy[ii],zk[ii][l])
            if pl[ii]:
                update_Poly3D(pl[ii],NXl[ii],NYl[ii],z0[ii],dx[ii],dy[ii],zl[ii][l])
            if pa[ii]:
                update_Poly3D(pa[ii],NXa[ii],NYa[ii],z0[ii],dx[ii],dy[ii],za[ii][l])
        ii+=1
    

    # Actualising the fields plot
    # ax6.cla()
    # # im0=ax6.streamplot(X,Y,Uop[l],Vop[l],color=np.sqrt(Uop[l]**2+Vop[l]**2),linewidth=2,cmap=cm.Reds)
    # im0=ax6.quiver(X,Y,Uop[l],Vop[l])
    line=e.readline()
    Z=compute_frame(line)

    im2.set_data(Z[1])
    cl2.on_mappable_changed(im2)
    im1.set_data(Z[3])
    cl1.on_mappable_changed(im1)
    im3.set_data(Z[0])
    cl3.on_mappable_changed(im3)
    im0.set_data(Z[2])
    cl0.on_mappable_changed(im0)

    r=[im0,im1,im2,im3]
    if '3D' in view.Isel:
        r.insert(0,xpoint)
    if 'xprof' in view.Isel or 'yprof' in view.Isel:
        r.insert(0,xralines)
    if 'diff'in view.Isel or 'xprof' in view.Isel or 'yprof' in view.Isel:
        r.insert(0,xrlines)
    if 'mode' in view.Isel:
        r.insert(0,pl)
        r.insert(0,pk)
        r.insert(0,pa)
    
    return tuple(r)


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
    print ' (Leave blank if not needed!)'
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
