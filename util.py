#!/usr/bin/env python
# Module with utility functions

import numpy as np
import matplotlib.ticker as ticker
import matplotlib.cm as cm
import gzip
import subprocess
from params import dimension as dim

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

if dim.dim:
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

