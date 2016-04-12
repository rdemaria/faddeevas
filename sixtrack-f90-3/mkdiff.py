
import mpmath

import numpy as np

from matplotlib.pyplot import *

data=open('sixtrack-3').read().replace('D','E')
a,b,x,y,u,v=np.fromstring(data,sep='\n').reshape(-1,6).T

mpmath.mp.dps = 50

def wofz(x, y):
    z = mpmath.mpc(x, y)
    w = mpmath.exp(-z**2) * mpmath.erfc(z * -1j)
    return w.real, w.imag

wofz = np.vectorize(wofz)

ur,vr=wofz(x,y)

dr=np.array(ur-u,dtype=float).reshape(101,101)
di=np.array(vr-v,dtype=float).reshape(101,101)

figure();
imshow(np.log10(abs(dr)),aspect='auto',origin='bottom')
clim(-15,0)
colorbar()
savefig('1.png')

figure();
imshow(np.log10(abs(di)),aspect='auto',origin='bottom')
clim(-15,0)
colorbar()
savefig('2.png')



