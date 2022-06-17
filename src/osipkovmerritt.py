

import numpy as np

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

	fi = (sqrt(2.0)/(378.0*(pi**3)*sqrt(sig0)))*((-q/sig0)**(7.0/2.0))*( 1.0-(ra**-2)+(63.0/4.0)*(ra**-2)*(-q/sig0)**(-2))

	assert fi >= 0, " DF negative! {0} r={1} vr,vt={2},{3} E,q={4},{5}".format(fi,r,vr,vt,E,q)

	return fi



def OM_maxfq(psi,r,ra,sf=1.1,steps=100.):
	"""find the maximum of the DF at each radius

	sf = 1.1 	# increase fmax found on grid by sf
	steps = 100.0 	# step size in velocity space vmax/steps

	"""
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
			val = OM_df(psi,r,ev,jv,ra)*jv
			if (maxfq < val):
				maxfq = val
	return sf*maxfq
