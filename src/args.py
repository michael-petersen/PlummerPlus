
import argparse

def parse_all_args():
    parser = argparse.ArgumentParser(description="Generates anisotropic and rotating models with plummer model density distribution  ")

    npar = parser.add_mutually_exclusive_group()
    npar.add_argument("-n", help="number of particles (default: 1000)",
                        type=int,default=1000,metavar="")

    npar.add_argument("-nk", help="number of particles in units of 1024, e.g. 32K",
                        type=int,default=0,metavar="")

    parser.add_argument("-rs", help="random number generator seed (default: 101)",
                        type=int,default=101,metavar="rand_seed")

    parser.add_argument("-rcut", help="outer cutoff radius",
                        type=float,default=-1.0,metavar="")

    parser.add_argument("-o", help="name of output file (default: \"fort.10\")",
                        type=str,default="fort.10",metavar="")

    parser.add_argument("-u", help="units (default: \"Henon units\")",
                        type=str,default="HU",metavar="")

    # rotation parameters
    parser.add_argument("-a", help="Lynden-Bell trick flip fraction",
                        type=float,default=0,metavar="")

    parser.add_argument("-oa", help="Offset angle in degrees for particles below energy cut off and flip fraction for offset stars (e.g. -oa 90.0 0.5): OA offset angle, A2 flip fraction, MF mass fraction above which the offset is applied ",
                        type=float,default=[0,0,-1],metavar=("OA","A2","MF"),nargs=3)

    parser.add_argument("-ms", help="Mass segregated model, set fraction of total mass (MF) and mass ratio (MR). This will reduce number of particles see Readme file. ",type=float,default=[0,0],metavar=("MF","MR"),nargs=2)

    parser.add_argument("-hs", help=" defines energy cut hs*pot(0) for counter rotation ",
                       type=float,default=0,metavar="")

    parser.add_argument("-lcut", help=" cut on Lz/L  ",
                        type=float,default=[-1,0,0],metavar=("G","AU","AO"),nargs=3)

    parser.add_argument("-icut", help=" cut on inclination normalise to range 0-1 (0-90 degrees)",
                        type=float,default=[-1,0,0],metavar=("G","AU","AO"),nargs=3)

    # exclusice arguments for velocity space
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-q", help="q values of Dejonghe (1987) anisotropic plummer models, q<=+2 ",
                        type=float,default=0,metavar="")

    group.add_argument("-ra", help="anisotropic radius of Osipkov-Merritt radially anisotropic plummer model, ra >= 0.75",
                        type=float,default=0,metavar="")

    group.add_argument("-e", help=" Einstein sphere, plummer model with only circular orbits ",
                        action="store_true")

    parser.add_argument("-v", help=" print output statistics ",
                        action="store_true")

    parser.add_argument("-z", help=" zero the position and velocity of the model ",
                        action="store_true")

    parser.add_argument("-qt", help="Quiet start, place replicas of particles at 2*pi/qt intervals in plane of orbit, see Sellwood 1997",
                        type=int,default=0,metavar="")

    output = parser.add_mutually_exclusive_group()
    output.add_argument("-bods", help=" Output for EXP body file ",
                        action="store_true")

    args = parser.parse_args()

    return args
