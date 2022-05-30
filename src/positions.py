import numpy as np


def generate_positions(args,w):

    # look for an outer cutoff radius
    if args.rcut == -1:
    	x = np.random.rand(args.n)
    else:
        # if there is an outer cutoff radius, then calculate the total mass of the cut model
        # in this formula, bc=1
    	mcut = args.rcut**3/(args.rcut**2 + 1.0)**1.5
    	x = np.random.uniform(0.0, mcut, args.n)

    # realise the radii of points by taking the cumulative mass curve for a Plummer sphere:
    # Menc/Mtot = (r^3)(1+r^2)^(-3/2)
    # x is a random variate between zero and 1: Menc/Mtot
    r = np.reciprocal(np.sqrt(np.power(x,-2.0/3.0)-1.0))

    # realise a random position on the surface of the sphere
    ctheta = np.random.uniform(-1.,1.0,args.n)
    stheta = np.sin(np.arccos(ctheta))
    #theta = np.arccos(theta)
    phi    = 2.0*np.pi*np.random.rand(args.n)

    # translate to Cartesian coordinates
    w[:,1] = stheta*np.cos(phi)
    w[:,2] = stheta*np.sin(phi)
    w[:,3] = ctheta

    w[:,1:4] = w[:,1:4]*r[:,None]

    return w,r
