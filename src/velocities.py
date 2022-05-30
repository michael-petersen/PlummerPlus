import numpy as np


def isotropic_velocities(r,w,args):

    # for each realised radius, make the square
    r2 = np.power(r,2)

    # compute the maximum velocity for each radius
    vmax = np.sqrt(2.0*np.reciprocal(np.sqrt(r2 + 1.0)))

    # loop through all particles
    v = np.zeros(args.n)
    for i in range(args.n):
    	while True:

            # draw two random variates between zero and 1
    		xi = np.random.rand(2)

            # the the tolerance: 10% of the second random variate
    		f1 =  xi[1]*0.1

            # draw from the distribution function using the first random variate
    		f  =  xi[0]*xi[0]*(1.0 - xi[0]*xi[0])**3.5

            # reject the point if the random variate is larger than the guard
    		if f >=  f1:
    			break

        # realise the velocity
    	v[i] = vmax[i]*xi[0]

    # realise a random vector pointing on the surface of the sphere: is this the same set of points as before? are we exploiting that??
    ctheta = 2.0*np.random.rand(args.n)-1.0
    stheta = np.sin(np.arccos(ctheta))
    phi    = 2.0*np.pi*np.random.rand(args.n)

    # translate velocities onto the point
    w[:,4] = stheta*np.cos(phi)
    w[:,5] = stheta*np.sin(phi)
    w[:,6] = ctheta

    w[:,4:] = w[:,4:]*v[:,None]

    return w
