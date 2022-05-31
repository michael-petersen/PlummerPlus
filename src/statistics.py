import numpy as np
import sys


def report_statistics(w,args,countflip):
    # statistics
    if args.e:
    	name = " Einstein sphere"
    elif args.ra > 0:
    	name = " Osipkov-Merritt radially anisotropic (ra={})".format(args.ra)
    elif args.q != 0:
    	name = " Dejonghe (1987) anisotropic q={}".format(args.q)
    else:
     	name = " Isotropic"


    print("\n{} Plummer model with N = {}  (random seed {})".format(name, args.n, args.rs ))


    svr2 = 0.0
    svt2 = 0.0
    for i in range(args.n):
    	vi = w[i,4:]
    	xi = w[i,1:4]
    	xiu = xi/np.linalg.norm(xi)
    	vr = np.inner(xiu,vi)
    	vt = np.linalg.norm(vi - vr*xiu)
    	svr2 += w[i,0]*vr**2
    	svt2 += w[i,0]*vt**2

    if args.ms[0] == 0:
    	r2 = np.sum(np.power(w[:,1:4],2), axis=1)
    	rh = np.sqrt(np.median(r2))
    else:
    	r2 = np.sum(np.power(w[:,1:4],2), axis=1)
    	indx = sorted(range(args.n),key=lambda k: r2[k])
    	mc = 0.0
    	for i in range(args.n):
    		mc += w[indx[i],0]
    		if mc >= 0.5:
    			rh = np.sqrt(r2[indx[i]])
    			break

    print(" rh = {:.3e} K.E. = {:.3e} vt^2 = {:.3e} vr^2 = {:.3e} ".format(rh,0.5*(svt2+svr2) ,svt2,svr2))
    #print(" Ixx = {:.3e} Iyy = {:.3e} Izz = {:.3e}  ".format((w[:,1]**2).sum(),(w[:,2]**2).sum(),(w[:,3]**2).sum() ))
    #
    # crit val http://adsabs.harvard.edu/abs/1981SvA....25..533P
    # more recent
    if True:
    #if not args.e:
    	print(" 1-0.5*<vt^2>/<vr^2> = {:.3e} ".format(1.0 - 0.5*svt2/svr2))
    	print(" 2.0Tr/Tp = {:.3} (Polyachenko and Shukhman 1981, crit value 1.7 +/- 0.25) ".format(2.0*svr2/svt2))

    #if args.a > 0 or args.oa[0] > 0:
    	L = np.multiply(w[:,0,None],np.cross(w[:,1:4],w[:,4:]))
    	L = np.sum(L, axis=0)
    	print(" L = [{:.3e},{:.3e},{:.3e}]  |L|={:.3e} nf={}".format(L[0],L[1],L[2],np.linalg.norm(L),countflip ))

    #
    # bins size, number of rings nbin^2, for nbin=10 each ring contains 1% mass
    #
    def sign(x):
    	if x >= 0:
    		return +1.0
    	else:
    		return -1.0
    #hack use pre cal R2,|z| bins as rho fixed
    rindex = [0.0385891298799, 0.0866921555876, 0.148588985051, 0.230883647572, 0.345272808206,
    0.517942199325,0.807197126028, 1.38667629541, 3.10347289646]

    nbin = 10
    R2 = np.sum(np.power(w[:,1:3],2), axis=1)
    z = np.abs(w[:,3])
    indx = sorted(range(args.n),key=lambda k: R2[k])
    Rbins = np.array_split(indx,nbin)
    vi = w[:,4:]
    xi = np.zeros((args.n,3))
    xi[:,0] = w[:,2]
    xi[:,1] = -w[:,1]
    xi = np.divide(xi,np.sqrt(R2)[:,None])

    # note xi is unit vector in v_phi direction
    vphi = xi[:,0]*vi[:,0] + xi[:,1]*vi[:,1]
    vphike = 0.0
    for rl in Rbins:
    	zi =  z[rl]
    	indxi = sorted(range(len(zi)),key=lambda k: zi[k])
    	zbins = np.array_split(indxi,nbin)
    	#print R2[rl[0]],R2[rl[-1]]
    	for zl in zbins:

    		#print "{:.2e}".format(zi[zl[-1]]),
    		vphiavg = np.sum(vphi[rl[zl]])
    		vsign = sign(vphiavg)
    		#mass = np.sum(w[rl[zl],0])
    		vphike += (vphiavg**2)/len(zl)
    		# note vphiavg  accutally (ns*vphiavg)
    		# vphike = n*(vphiavg)**2
    		# divide by n before output
    	#print ""

    #http://adsabs.harvard.edu/abs/1973ApJ...186..467O
    print(" T_phi/|pot| = {:.3e}, (assuming pot = 0.5, see ostriker & peebles 1973, 0.14 +/- 0.03)".format(vphike/args.n))

    if args.oa[2] >= 0.0:
    	print(" Gamma = {:.3e} (|Ecut/Emin|) ".format(abs(1+ecut)))
    	print("De({:d},{:.2f},{:.3e}) {:.3e}".format(int(args.q),args.a,abs(1.+ecut),args.oa[2]) )

    #print command arguments
    commandstring = '\n ';
    for arg in sys.argv:
        if ' ' in arg:
            commandstring+= '"{}"  '.format(arg) ;
        else:
            commandstring+="{}  ".format(arg) ;
    print(commandstring+"\n");
