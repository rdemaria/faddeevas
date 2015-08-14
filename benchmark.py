from numpy import *
import scipy

from wwerf import ccperrfr as wf90
wf90=vectorize(wf90)

from wwerf2 import errf as wf
wf=vectorize(wf)

from mywwerf import wwerf as wpy

from scipy.special import erfc, wofz
def wsci(x,y):
    z=x+1j*y
    w=wofz(z)
    return w.real,w.imag

#https://github.com/PyCOMPLETE/PyHEADTAIL/master/spacecharge/spacecharge.py

#http://www.bigdft.org/devel-doc/d0/dfa/wofz_8f90_source.html
#http://www.bigdft.org/devel-doc/d0/dfa/wofz_8f90.html
#http://dl.acm.org/citation.cfm?id=77629
# latest
#www.ccsenet.org/journal/index.php/jmr/article/viewFile/41877/24151 
# http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package
# arxiv.org/pdf/1407.0748 In root with source code


r=15
xx=10**linspace(-r,r,301)
yy=10**linspace(-r,r,301)
x,y=meshgrid(xx,yy)

x1,y1=wsci(x,y)
x2,y2=wf90(x,y)
x3,y3=wf(x,y)

clf();
imshow(log10(abs(x1-x2)),origin="bottom",extent=[-r,r,-r,r],aspect='auto')
colorbar()


figure()
clf();imshow(x1);colorbar()
figure()
clf();imshow(x2);colorbar()



ref90=[ (x,y)+wf90(x,y) for x in xx for y in yy]
repy =[ (x,y)+wpy(x,y)  for x in xx for y in yy]
resci=[ (x,y)+wsci(x,y) for x in xx for y in yy]

(array(ref90)-array(resci)).std()



