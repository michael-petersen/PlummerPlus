
import numpy as np

# local imports
from .positions import unitv



#-------- Python3 compatibility: to be deprecated
try:
    xrange
except NameError:
    xrange = range
#----------------------


def generate_einstein_sphere(args,w,r,rbit,theta):
    r2 = np.power(r,2)

    # circular velocity in a Plummer potential
    v = np.sqrt( np.divide(r2, np.power( np.sqrt( r2+1.0 ),3.0 ) ))

    for i in xrange(args.n):
    	# create unit vectors e1,e2 in plane normal to r
    	# assuming ri[0], ri[1] non-zero
    	vr,vt = unitv(w[i,1:4],rbit,theta)
    	w[i,4:] = vt*v[i]

    return w
