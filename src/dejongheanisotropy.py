"""
Anisotropy from Dejonghe (1986)

The hypergeometrics need to be checked for very negative q here

"""


import numpy as np
from math import pi, sqrt, cos, sin, acos
from scipy.interpolate import interp1d
import scipy.special as sci

# local imports
from .positions import unitv


#-------- Python3 compatibility: to be deprecated
try:
    xrange
except NameError:
    xrange = range
#----------------------


def make_dejonghe_anistropy(args,w,r,rbit,theta):
    assert args.q <= +2, " q value needs to be in range (-inf,+2) "
    #
    # Define distribution function as in Dejonghe (1987)
    #
    # sci.gamma(a+b) replaced with 1.0 as cancels with sci.gamma(0.5*q) in Fq
    # for x <= 1, sci.gamma(a+d) = 1
    # cost1 coeff for x=<1, cost2 coeff for x>1
    cost1 = 1.0/( sci.gamma(4.5-args.q))
    cost2 = 1.0/( sci.gamma(1.0-0.5*args.q) * sci.gamma(4.5-0.5*args.q))

    def H(a, b, c, d, x ):
    	if x  <= 1:
    		return cost1*(x**a)*sci.hyp2f1(a+b, 1+a-c, a+d, x)
    	else:
    		y = 1.0 / x
    		return cost2*(y**b)*sci.hyp2f1(a+b, 1+b-d, b+c, y)

    #remove factor sci.gamma(0.5*q) cancels with term in H
    const3 = ( 3.0*sci.gamma(6-args.q)/(2.0*( (2.0*pi)**2.5 ) ) )
    def Fq(E, L):
    	assert E > 0
    	return const3  * (E**(3.5-args.q)) * H(0.0, 0.5*args.q, 4.5-args.q, 1.0, (L**2.0) / (2.0*E))

    sf = 1.1 	# increase fmax found on grid by sf
    steps = 100.0 	# step size in velocity space vmax/steps
    def DJ_maxfq(psi,r):
    	maxfq = 0.0
    	vmax = sqrt(2.0*psi)
    	incr = vmax/steps
    	for ev in np.arange(0,vmax,incr):
    		E=psi-0.5*ev**2
    		if E <= 0 and abs(E) < 1e-15:
    			continue
    		for jv in np.arange(0,vmax,incr):
    			if ( jv > ev ): #not realisitc
    				continue
    			L=jv*r
    			val = Fq(E,L)*jv
    			if (maxfq < val):
    				maxfq = val
    	return sf*maxfq

    psirange = np.linspace(0.000, 1.000, num=100, endpoint=True)
    maxfqarray = np.zeros(len(psirange))
    for i,psi in enumerate(psirange) :
    	if psi == 0.0:
    		maxfqarray[i] = 0.0
    		continue
    	ri = sqrt((1.0/psi)**2.0-1.0)
    	maxfqarray[i] = DJ_maxfq(psi,ri)

    # linear interpolation function to calculate bound
    psimax = interp1d(psirange, maxfqarray)

    # accept reject sampling: select a potential value
    r2 = np.power(r,2)
    psi = np.reciprocal(np.sqrt(r2 + 1.0))

    #
    fmaxv = psimax(psi)
    vmax = np.sqrt(2.0*psi)
    v = np.zeros(args.n)
    rvf =  np.random.rand(int(round(50*args.n)))
    rvc = 0
    for i in xrange(args.n):
    	fmax = fmaxv[i]
    	loopc = 0
    	while True:
    		#rv =  np.random.rand(3)
    		vr =  rvf[rvc+0]*vmax[i]
    		vt =  rvf[rvc+1]*vmax[i]
    		l = r[i]*vt
    		vsq = vr**2 + vt**2
    		E =  psi[i] - 0.5*vsq

    		if E < 0:
    			rvc+=2
    			continue

    		f1 =  rvf[rvc+2]*fmax
    		f = Fq(E,l)*vt

    		rvc += 3

    		if f >= f1:
    		#	print i,loopc
    			vrv,vtv = unitv(w[i,1:4],rbit,theta)
    			w[i,4:] = vrv*vr+vtv*vt
    			break

    		if rvc + 4 >= len(rvf):
    			rvf = np.random.rand(int(round(50*args.n)))
    			rvc = 0

    		loopc += 1
    		if loopc > 100000:
    			print(r[i], fmax, E, l)
    			raise NameError('Failed to sample')

	#print float(rvc)/float(len(rvf)),rvc,len(rvf),i
    return w
