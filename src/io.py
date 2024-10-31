
import numpy as np

# consider adding this in the future
#import h5py


oldstyle = True

def write_particles(w,args):

    # multiply masses
    w[:,0] *= args.M

    # scale up velocities
    w[:,4:] *= np.sqrt(args.M)

    # apply the initial set shift (this needs to come last!)
    w[:,1:] += args.init

    # this version writes for EXP
    if args.bods:
        if oldstyle:
            with open(args.o, 'w') as f:
                print('{} 0 0'.format(args.n),file=f)
                for i in range(0,args.n):
                     #print(w[i,0],w[i,1]+args.init[0],w[i,2]+args.init[1],w[i,3]+args.init[2],w[i,4]+args.init[3],w[i,5]+args.init[4],w[i,6]+args.init[5],file=f)
                     print(w[i,0],w[i,1],w[i,2],w[i,3],w[i,4],w[i,5],w[i,6],file=f)
        else:
            with open(args.o, 'wb') as f:
                np.array([args.n]).tofile(f)
                np.array([0,0]).tofile(f)
                for i in range(0,args.n):
                    np.array([w[i,0],w[i,1],w[i,2],w[i,3],w[i,4],w[i,5],w[i,6]]).tofile(f)
    # this version writes for NBODY6++
    
    else:
        np.savetxt(args.o, w)
