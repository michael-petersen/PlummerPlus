

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


def OM_df(psi,r,vr,vt,ra):
	"""# pretabulate the DF"""
	E = psi - 0.5*(vr**2+vt**2)
	if E < 0:
		return 0.0
	q = -E + (r**2.0)*(vt**2.0)/(2.0*ra**2)
	if q >= 0:
		return 0.0

	# Plummer central velocity dispersion (scalar)
	sig0 = 1.0/6.0

    #                             is this ra supposed to be to the negative square?------|
    #                                                                                    v
	fi = (np.sqrt(2.0)/(378.0*(np.pi**3)*np.sqrt(sig0)))*((-q/sig0)**(7.0/2.0))*( 1.0-(ra**-2)+(63.0/4.0)*(ra**-2)*(-q/sig0)**(-2))

	assert fi >= 0, " DF negative! {0} r={1} vr,vt={2},{3} E,q={4},{5}".format(fi,r,vr,vt,E,q)

	return fi



def OM_maxfq(psi,r,ra,sf=1.1,steps=100.):
	"""find the maximum of the DF at each radius

	sf = 1.1 	# increase fmax found on grid by sf
	steps = 100.0 	# step size in velocity space vmax/steps

	"""
	maxfq = 0.0
	vmax = np.sqrt(2.0*psi)
	incr = vmax/steps
	for ev in np.arange(0,vmax,incr):
		E=psi-0.5*ev**2
		if E <= 0 and abs(E) < 1e-15:
			continue
		for jv in np.arange(0,vmax,incr):
			if ( jv > ev ): #not realisitc
				continue
			L=jv*r
			val = OM_df(psi,r,ev,jv,ra)*jv
			if (maxfq < val):
				maxfq = val
	return sf*maxfq


def generate_osipkov_merritt_radial_anisotropy(args,w,r,rbit,theta):
	assert args.ra >= +0.75, " ra value needs to be in range (+0.75,+inf) "

    # define the table for psi, from (0->1)
	psirange = np.linspace(0.000, 1.000, num=100, endpoint=True)
	maxfqarray = np.zeros(len(psirange))
	for i,psi in enumerate(psirange) :
		if psi == 0.0:
			maxfqarray[i] = 0.0
			continue
		ri = sqrt((1.0/psi)**2.0-1.0)
		maxfqarray[i] = OM_maxfq(psi,ri,args.ra)

	# linear interpolation function to calculate bound
	psimax = interp1d(psirange, maxfqarray)

	# accept reject sampling
	r2 = np.power(r,2)

    # draw the potential value from the cumulative potential function
	psi = np.reciprocal(np.sqrt(r2 + 1.0))
	fmaxv = psimax(psi)
	vmax = np.sqrt(2.0*psi)
	v = np.zeros(args.n)

    # select 40x the desired number of particles
	rvf =  np.random.rand(int(round(40*args.n)))
	rvc = 0
	for i in xrange(args.n):
		fmax = fmaxv[i]
		loopc = 0
		while True:
			rvc+=3
			vr = rvf[rvc]*vmax[i]
			vt = rvf[rvc+1]*vmax[i]
			l = r[i]*vt
			vsq = vr**2 + vt**2
			E =  psi[i] - 0.5*vsq

			if E < 0:
				continue

            # scale to the appropriate boundary
			f1 =  rvf[rvc+2]*fmax
			f = OM_df(psi[i],r[i],vr,vt,args.ra)*vt

			if f >= f1:
				# vrv,vtv create random vr vt unit vectors
				vrv,vtv = unitv(w[i,1:4],rbit,theta)
				w[i,4:] = vrv*vr+vtv*vt
				break
			loopc += 1
			if loopc > 10000:
				print(r[i], fmax, E, l)
				raise NameError('Failed to sample')
	#print float(rvc)/float(len(rvf))
	return w
