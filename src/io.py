
import numpy as np

# consider adding this in the future
#import h5py


oldstyle = True

def write_particles(w,args):
    # this version writes for EXP
    if args.bods:
        if oldstyle:
            with open(args.o, 'w') as f:
                print('{} 0 0'.format(args.n),file=f)
                for i in range(0,args.n):
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
