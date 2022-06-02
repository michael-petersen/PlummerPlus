#!/usr/bin/env python
#
# PROGRAM: PlummerPlus.py
#
# Description: Generates realisation of anisotropic and rotating models with plummer model density distribution
# 		Note on units G=M=a=1 for generation but values are rescales to Henon units for output
#		(see scale factors: lfact=3.0*pi/16.0, vfact=1/lfact. default output file fort.10)
#
# AUTHOR: Phil Breen, University of Edinburgh
#
# LICENSE:	MIT, This program is comes with ABSOLUTELY NO WARRANTY.
#
#
# CITATION: Please use Breen, Varri, Heggie 2017 (https://ui.adsabs.harvard.edu/abs/2017MNRAS.471.2778B/abstract)
#
# EXAMPLES:
#       use "chmod a+x PlummerPlus.py" to make script executable, alternativaly use "python PlummerPlus.py ...."
#
#
""" PlummerPlus.py: generates realisation of anisotropic and rotating models with plummer model density distribution """

import numpy as np
import scipy.special as sci
import sys
from scipy.interpolate import interp1d
from math import pi, sqrt, cos, sin, acos

# local imports
from src.quiet import *
from src.args import *
from src.positions import *
from src.velocities import *
from src.statistics import *

args = parse_all_args()

#-------- Python3 compatibility
try:
    xrange
except NameError:
    xrange = range
#----------------------

# set the random number seed
np.random.seed(args.rs)

if args.nk > 0:
	args.n = 1024*args.nk


# the workhorse here is the matrix w, which is
# [mass, xpos, ypos, zpos, xvel, yvel, zvel]

#--------- reduce n to n/qt
if args.qt > 1:
	#only for quiet starts
	m = 1.0/float(args.n)
	args.n = int(args.n/args.qt)
	w = np.zeros((args.n, 7))
	# w[:,0] = mass, w[:,1:4] = x,y,z, w[:,4:] = vx,vy,zx
	w[:,0] = m
else:
	# more general case
	w = np.zeros((args.n, 7))
	# w[:,0] = mass, w[:,1:4] = x,y,z, w[:,4:] = vx,vy,zx
	w[:,0] = 1.0/float(args.n)


#---------------------------------generate positions------------------------------------
w,r = generate_positions(args,w)



#---------------------------------generate velocities------------------------------------

#calucates orthogonal basis using r and returns random vr,vt units (needed for anisotropic models)
if args.q != 0 	or args.ra != 0 or args.e:
	sign = [-1.0,1.0]
	theta = 2.0*pi*np.random.rand(args.n)
	ctheta = np.cos(theta)
	stheta = np.sin(theta)
	rbit = np.random.randint(2, size=args.n)
	ui = 0
	def unitv(ri):
		ru = ri/np.linalg.norm(ri)
		e1 = [ri[1],-ri[0],0.0]
		e1 /= np.linalg.norm(e1)
		e2 = np.cross(ru,e1)
		vr = ru*sign[rbit[unitv.i]]
		vt = (e1*ctheta[unitv.i] + e2*stheta[unitv.i])
		unitv.i += 1
		return vr,vt
	unitv.i = 0



"""
In the anisotropic case, the magnitudes of the tangential and radial velocities (vt,vr) must be sampled from a joint probability distribution \propto f(r,vt,vr)vt. For acceptance-rejection sampling, a bound on the maximum value of the f(r,vt,vr)vt has to be known at each radius. This is calculated numerically in advance of sampling on a radial grid. During sampling, the bound a a particular r is calculated by interpolation. This technique is similar to Aarset, Hénon & Wielen (1974).

"""
#anisotropic plummer Dejonghe
if args.q != 0:
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
	def maxfq(psi,r):
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
		maxfqarray[i] = maxfq(psi,ri)

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
				vrv,vtv = unitv(w[i,1:4])
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


# anisotropoic plummer Osipkov-Merritt radial only Osipkov 1979; Merritt 1985
elif  args.ra != 0:
	assert args.ra >= +0.75, " ra value needs to be in range (+0.75,+inf) "

    # pretabulate the DF
	def df(psi,r,vr,vt):
		E = psi - 0.5*(vr**2+vt**2)
		if E < 0:
			return 0.0
		q = -E + (r**2.0)*(vt**2.0)/(2.0*args.ra**2)
		if q >= 0:
			return 0.0
		sig0 = 1.0/6.0
		fi = (sqrt(2.0)/(378.0*(pi**3)*sqrt(sig0)))*((-q/sig0)**(7.0/2.0))*( 1.0-(args.ra**-2)+(63.0/4.0)*(args.ra**-2)*(-q/sig0)**(-2))

		assert fi >= 0, " DF negative! {0} r={1} vr,vt={2},{3} E,q={4},{5}".format(fi,r,vr,vt,E,q)
		return fi

	sf = 1.1 	# increase fmax found on grid by sf
	steps = 100.0 	# step size in velocity space vmax/steps
	def maxfq(psi,r):
        """find the maximum of the DF at each radius"""
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
				val = df(psi,r,ev,jv)*jv
				if (maxfq < val):
					maxfq = val
		return sf*maxfq

    # define the table for psi, from (0->1)
	psirange = np.linspace(0.000, 1.000, num=100, endpoint=True)
	maxfqarray = np.zeros(len(psirange))
	for i,psi in enumerate(psirange) :
		if psi == 0.0:
			maxfqarray[i] = 0.0
			continue
		ri = sqrt((1.0/psi)**2.0-1.0)
		maxfqarray[i] = maxfq(psi,ri)

	# linear interpolation function to calculate bound
	psimax = interp1d(psirange, maxfqarray)

	# accept reject sampling
	r2 = np.power(r,2)
	psi = np.reciprocal(np.sqrt(r2 + 1.0))
	fmaxv = psimax(psi)
	vmax = np.sqrt(2.0*psi)
	v = np.zeros(args.n)
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

			f1 =  rvf[rvc+2]*fmax
			f = df(psi[i],r[i],vr,vt)*vt

			if f >= f1:
				# vrv,vtv random vr vt unit vectors
				vrv,vtv = unitv(w[i,1:4])
				w[i,4:] = vrv*vr+vtv*vt
				break
			loopc += 1
			if loopc > 10000:
				print(r[i], fmax, E, l)
				raise NameError('Failed to sample')
	#print float(rvc)/float(len(rvf))

#Einstein sphere
elif args.e:
	r2 = np.power(r,2)
	v = np.sqrt( np.divide(r2, np.power( np.sqrt( r2+1.0 ),3.0 ) ))

	for i in xrange(args.n):
		# create unit vectors e1,e2 in plane normal to r
		# assuming ri[0], ri[1] non-zero
		vr,vt = unitv(w[i,1:4])
		w[i,4:] = vt*v[i]


#isotropic plummer
else:
    isotropic_velocities(r,w,args)


#----------------- quiet start --------
if args.qt > 1:
    w = quiet_start(w,args)


#---------------------------------add rotation------------------------------------

countflip = 0

if args.lcut[0] >= 0.0:

	L = np.cross(w[:,1:4],w[:,4:])

	for i in xrange(args.n):

		crit = abs(L[i,2]/np.linalg.norm(L[i,:]))

		if L[i,2] < 0.0 and crit < args.lcut[0]:

			if args.lcut[1] >= np.random.rand():
				w[i,4:] *= -1.0
				countflip+=1

		if L[i,2] < 0.0 and crit > args.lcut[0]:

			if args.lcut[2] >= np.random.rand():
				w[i,4:] *= -1.0
				countflip+=1


elif args.icut[0] >= 0.0:

	L = np.cross(w[:,1:4],w[:,4:])

	for i in xrange(args.n):

		crit = acos(abs(L[i,2]/np.linalg.norm(L[i,:])))/(pi/2.)

		if L[i,2] < 0.0 and crit < args.icut[0]:

			if args.icut[1] >= np.random.rand():
				w[i,4:] *= -1.0
				countflip+=1

		if L[i,2] < 0.0 and crit > args.icut[0]:

			if args.icut[2] >= np.random.rand():
				w[i,4:] *= -1.0
				countflip+=1




elif args.oa[2] >= 0.0:
	# general trick energy cut and
	v2 = np.sum(np.power(w[:,4:],2.0),axis=1)
	E = 0.5*v2 - np.reciprocal(np.sqrt(np.power(r,2.0) + 1.0))
	#E = np.sort(E)
	#for i in range(args.n):
	#	print(" {} {} ".format(float(i)/float(args.n),E[i]))
	#exit()

	#poti = np.reciprocal(np.sqrt(np.power(r,2.0) + 1.0))
	ecut = np.percentile(E, 100.0*args.oa[2])
	#print(" {} {} {} {} ".format(max(r),min(r),max(E),min(E)))
	#exit()

#if args.oa[0] > 0:
	theta = pi*args.oa[0]/180.0
	ov = np.array([0.0,sin(theta),cos(theta)])
	L = np.cross(w[:,1:4],w[:,4:])

	for i in xrange(args.n):
		if L[i,2] < 0.0 and E[i] < ecut:
			if args.a >= np.random.rand():
				w[i,4:] *= -1.0
				countflip+=1
		if E[i] > ecut:
			if np.dot(L[i,:],ov) < 0.0:
				if args.oa[1] >= np.random.rand():
					w[i,4:] *= -1.0
					countflip+=1


elif args.hs != 0.0:
	v2 = np.sum(np.power(w[:,4:],2.0),axis=1)
	E = 0.5*v2 - np.reciprocal(np.sqrt(np.power(r,2.0) + 1.0))
	L = np.cross(w[:,1:4],w[:,4:])

	# calculation for Lz = 0 high shear models, no long used
	#Lz = abs(L[:,2])
	#Lztot = sum(Lz)
	#indx = sorted(range(args.n),key=lambda k: E[k])
	#Lzc = 0.0
	#for i in xrange(args.n):
	#	pid = indx[i]
	#	Lzc += Lz[pid]
	#	if Lzc >= 0.5*Lztot:
	#		ecut = E[pid]
	#		break
	#print(ecut)
	#exit()
	ecut = -1.0*args.hs #pot(0)=1.0

	for i in xrange(args.n):
		if L[i,2] < 0.0 and E[i] < ecut:
			if args.a > np.random.rand():
				w[i,4:] *= -1.0
				countflip+=1

		if L[i,2] > 0.0 and E[i] > ecut:
			if args.a > np.random.rand():
				w[i,4:] *= -1.0
				countflip+=1

elif args.a > 0:
	#basic LB trick
	countflip = 0
	L = np.cross(w[:,1:4],w[:,4:])
	for i in xrange(args.n):
		if L[i,2] < 0.0:
			if args.a > np.random.rand():
				w[i,4:] *= -1.0
				countflip += 1


#-----------------------------handle mass segregation---

if args.ms[0] > 0 and args.ms[1] > 0:

	if args.oa[2] == -1:
		v2 = np.sum(np.power(w[:,4:],2.0),axis=1)
		E = 0.5*v2 - np.reciprocal(np.sqrt(np.power(r,2.0) + 1.0))
	#print(max(E),min(E),len(x))
	#fexit()
	#E = np.sort(E)
	#for i in range(len(E)):
	#	if i % 1000:
	#		print(" {} {} ".format(float(i)/args.n,E[i]))
	#exit()
	ecut = np.percentile(E, 100.0*args.ms[0])

	nfold = args.ms[0]*args.n
	nkeep  = int(round(nfold/args.ms[1]))
	nc = 0

	rl = []
	for i in xrange(args.n):
		if E[i] < ecut:
			nc+=1
			if nc<=nkeep:
				w[i,0] *= args.ms[1]
			else:

				rl.append(i)
	print("\n Warning mass segregation particle number reduced n = {}!".format(args.n-len(rl)))

	w=np.delete(w,rl,axis=0)
	# reset n and renormlise mass to account for round error
	w[:,0] *= 1.0/sum(w[:,0])
	args.n = args.n-len(rl)
#



#--------------------------------scale to Henon units and save data--------------------------------

# scale to henon units and save data to output file "fort.10" (use -o to rename output)
if args.u == "HU":
	lfact=(3.0*pi)/16.0
	vfact = 1.0/sqrt(lfact)
	w[:,1:4] *= lfact
	w[:,4:] *= vfact

# apply zeroing?
if args.z:
    for i in range(1,7):
        w[:,i] -= np.nanmean(w[:,i])

# print the output file... we can almost certainly improve this.
if args.bods:
    with open(args.o, 'w') as f:
        print('{} 0 0'.format(args.n),file=f)
        for i in range(0,args.n):
            print(w[i,0],w[i,1],w[i,2],w[i,3],w[i,4],w[i,5],w[i,6],file=f)

else:
    np.savetxt(args.o, w)



# report statistics?
if args.v:
    report_statistics(w,args,countflip)
