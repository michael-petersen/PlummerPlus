
import numpy as np


def rotation_matrix(axis, theta):
    # ref https://en.wikipedia.org/wiki/Euler%E2%80%93Rodrigues_formula
	axis = np.asarray(axis)
	axis = axis / sqrt(np.dot(axis, axis))
	a    = cos(theta / 2.0)
	b, c, d = -axis * sin(theta / 2.0)
	aa, bb, cc, dd = a * a, b * b, c * c, d * d
	bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
	return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                 [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                 [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


def quiet_start(w,args):

    qtl = [w]

    for i in range(1,args.qt):
    	wi = np.zeros_like(w)
    	wi[:,0] = w[:,0]
    	qtl.append(wi)

    for j in range(args.n):
    	ri = w[j,1:4]
    	vi = w[j,4:]
    	L = np.cross(ri,vi)
    	for i in range(1,args.qt):
    		angi = i*2.*pi/args.qt
    		rotmat = rotation_matrix(L, angi)
    		wi = qtl[i]
    		wi[j,1:4] = np.dot(rotmat, ri)
    		wi[j,4:] = np.dot(rotmat, vi)
    		#print(wi[i,1:4],wi[i,4:])
    		#print(w[i,1:4],w[i,4:])
    		#print(qtl[i][-1],wi[-1])
    		#exit()
    		#print(qtl[i][-1],wi[-1])

    args.n = args.n*args.qt
    w = np.concatenate(qtl)
    #print(w[-1])
    return w
