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
from src.quiet              import *
from src.args               import *
from src.positions          import *
from src.velocities         import *
from src.statistics         import *

from src.osipkovmerritt     import *
from src.dejongheanisotropy import *
from src.einstein           import *

from src.lflip              import *
from src.io                 import *

# parse the input arguments
args = parse_all_args()

print(args.init)

#-------- Python3 compatibility: to be deprecated
try:
    xrange
except NameError:
    xrange = range
#----------------------

# set the random number seed
np.random.seed(args.rs)

# I don't think we need this feature anymore. Plan to deprecate.
if args.nk > 0:
    print("PlummerPlus: -nk will soon be deprecated, please default to using -n.")
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
# w is fully stored in memory now, so we can simply pass it around
w,r = generate_positions(args,w)



#---------------------------------generate velocities------------------------------------
# many different options here!

#calculates orthogonal basis using r and returns random vr,vt units (needed for anisotropic models)
if args.q != 0 	or args.ra != 0 or args.e:
	sign = [-1.0,1.0]

    # create random azimuthal angle
	theta = 2.0*pi*np.random.rand(args.n)

	rbit = np.random.randint(2, size=args.n)
	ui = 0


	unitv.i = 0



"""
In the anisotropic case, the magnitudes of the tangential and radial velocities (vt,vr) must be sampled from a joint probability distribution \propto f(r,vt,vr)vt. For acceptance-rejection sampling, a bound on the maximum value of the f(r,vt,vr)vt has to be known at each radius. This is calculated numerically in advance of sampling on a radial grid. During sampling, the bound a a particular r is calculated by interpolation. This technique is similar to Aarseth, Hénon & Wielen (1974).

"""
#anisotropic plummer Dejonghe
if args.q != 0:
    w = make_dejonghe_anistropy(args,w,r,rbit,theta)

# anisotropoic plummer Osipkov-Merritt radial only Osipkov 1979; Merritt 1985
elif  args.ra != 0:
    w = generate_osipkov_merritt_radial_anisotropy(args,w,r,rbit,theta)

#Einstein sphere
elif args.e:
    w = generate_einstein_sphere(args,w,r,rbit,theta)

#isotropic plummer
else:
    isotropic_velocities(r,w,args)


#----------------- quiet start --------
if args.qt > 1:
    w = quiet_start(w,args)


#---------------------------------add rotation------------------------------------

countflip = 0

if args.lcut[0] >= 0.0:

    w,countflip = flip_lcut(w,args.lcut,countflip)


elif args.icut[0] >= 0.0:

    w,countflip = flip_icut(w,args.icut,countflip)



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

elif args.a != 0:
    # basic LB trick: the primary workhorse
    countflip = 0

    # compute angular momentum for all particles
    L = np.cross(w[:,1:4],w[:,4:])

    # loop through particles and flip with some probability
    for i in xrange(args.n):

        # make positive spin
        if args.a > 0:
            if L[i,2] < 0.0:
                if args.a > np.random.rand():
                    w[i,4:] *= -1.0
                    countflip += 1

        # make negative spin
        else:
            if L[i,2] > 0.0:
                if np.abs(args.a) > np.random.rand():
                    w[i,4:] *= -1.0
                    countflip += 1

    print("Adding rotation. Countflip={}.".format(countflip))

#-----------------------------handle mass segregation-----------------------

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



#--------------------------------scale to Henon units--------------------------------

# scale to henon units?
if args.u == "HU":
	lfact=(3.0*pi)/16.0
	vfact = 1.0/sqrt(lfact)
	w[:,1:4] *= lfact
	w[:,4:] *= vfact


#--------------------------------mean-centre units--------------------------------

# apply zeroing?
if args.z:
    for i in range(1,7):
        w[:,i] -= np.nanmean(w[:,i])


#-------------------------------------save data--------------------------------

# print the output file... we can almost certainly improve this.
# save data to output file "output.txt" (use -o to rename output)

write_particles(w,args)

# report statistics?
if args.v:
    report_statistics(w,args,countflip)
