import numpy as np



def flip_lcut(w,lcut,countflip):
    """lcut is zero to start

    so that's where the discontinuity is!
    """

    nparticles = w.shape[0]

    # create the angular momenta
    L = np.cross(w[:,1:4],w[:,4:])

    for i in xrange(nparticles):

        # normalise the angular momentum
    	crit = abs(L[i,2]/np.linalg.norm(L[i,:]))

        # lower branch on Lz cut (less than critical value)
    	if L[i,2] < 0.0 and crit < lcut[0]:

    		if lcut[1] >= np.random.rand():
    			w[i,4:] *= -1.0
    			countflip+=1

    	if L[i,2] < 0.0 and crit > lcut[0]:

    		if lcut[2] >= np.random.rand():
    			w[i,4:] *= -1.0
    			countflip+=1

    return w,countflip



def flip_icut(w,icut,countflip):

    nparticles = w.shape[0]

    L = np.cross(w[:,1:4],w[:,4:])

    for i in xrange(nparticles):

    	crit = acos(abs(L[i,2]/np.linalg.norm(L[i,:])))/(pi/2.)

    	if L[i,2] < 0.0 and crit < icut[0]:

    		if icut[1] >= np.random.rand():
    			w[i,4:] *= -1.0
    			countflip+=1

    	if L[i,2] < 0.0 and crit > icut[0]:

    		if icut[2] >= np.random.rand():
    			w[i,4:] *= -1.0
    			countflip+=1

    return w,countflip
