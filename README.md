PlummerPlus
===========

Generates anisotropic and rotation Plummer models, see [Breen, Varri & Heggie (2019)](https://ui.adsabs.harvard.edu/abs/2017MNRAS.471.2778B/abstract) for details.

Modified by Michael Petersen for object-oriented maximalism.

Example:

```sh
$./PlummerPlus.py -n 10000 -q -6 -o outputfile -v

 Dejonghe (1987) anisotropic q=-6.0 Plummer model with N = 10000  (random seed 101)
 rh = 7.722e-01 K.E. = 2.482e-01 vt^2 = 4.138e-01 vr^2 = 8.267e-02
 1-0.5*<vt^2>/<vr^2> = -1.503e+00
 2.0Tr/Tp = 0.4 (Polyachenko and Shukhman 1981, crit value 1.7 +/- 0.25)
 L = [1.568e-03,-1.123e-03,1.566e-03]  |L|=2.484e-03 nf=0
 T_phi/|pot| = 2.642e-03, (assuming pot = 0.5, see ostriker & peebles 1973, 0.14 +/- 0.03)

```

A laundry list of examples:

```python
# 	isotropic plummmer with 8K particles
	./PlummerPlus.py -n 8192

#	8k anisotropic plummer with Dejonghe with q=-2 (see [Dejonghe 1987](http://adsabs.harvard.edu/full/1987MNRAS.224...13D))
	./PlummerPlus.py -n 8192 -q -2

# 	8k Osipkov-Merritt  radially anisotropic plummer with anisotropic radius ra=1.0 (see e.g. merritt, d. 1985. aj, 90, 1027)
	./PlummerPlus.py -n 8192 -ra 1.0

#	8k Einstien shpere i.e. plummer model with only circular orbits
	./PlummerPlus.py -n 8192 -e

#	8K isotropic plummmer with rotation via Lynden-Bell trick i.e. reverse velocities of 50% particles with L_z < 0
       ./PlummerPlus.py -n 8192 -a 0.5

#	8K isotropic plummmer with offset rotation, 50% of the most bound stars rotating about z-axis
# 	least bound 50% rotating is offset by 90 degrees (i.e. about y-axis).
#       Note most stars use -a for the fraction of stars reverse where as least bound stars use second value in -oa flag (i.e. -oa angle flipfraction)
       ./PlummerPlus.py -n 8192 -a 1.0  -oa 90.0 1.0 0.5

# 	for a high shear model, set offset to 180 degree (i.e. -z) and used mass fraction 0.67 (i.e. rotation model with 0 net L!)
	./PlummerPlus.py -n 10000  -a 1.0 -oa 180.0 1.0 0.67

#	create mass segregated model most bound 50% of the mass consisting of particles 10.0 times more massive then the least
	./PlummerPlus.py -n 10000 -ms 0.5 10.0
```

### Features
 - Ansiotropic models of Dejonghe (1987), tangential and radial velocity anisotropy
 - Osipkov-Merritt model radial velocity anisotropy
 - Rotation introduced via the Lyden-Bell trick (and generalisations)
 - use -h for full list of options

### Todos

 - include embedded Plummer models using Eddington Formula
 - Kroupa IMF (plus evolution with SSE)
