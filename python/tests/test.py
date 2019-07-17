# Import the modules
import sys
import MinVel as mv
import numpy as np

nMin = 1
nPT = 1
fl = 'tests/testInput.dat.DONOTCHANGE'
print('Reading input file {0:s}'.format(fl))
f = open(fl,"r")
if f.mode == "r":
    nPT = 0
    ln = 0
    for line in f:
        line = line.strip()
        columns = line.split(",")
        if ln < 2:
            nMin = len(columns)
        else:
            nPT = nPT + 1
        ln = ln + 1
f.close()
T = np.zeros((nPT),dtype=np.float64)
P = np.zeros((nPT),dtype=np.float64)
Cm = np.zeros((nMin),dtype=np.int8)
Cv = np.zeros((nMin),dtype=np.float64)
f = open(fl,"r")
if f.mode == "r":
    ln = 0
    jT = 0
    for line in f:
        line = line.strip()
        columns = line.split(",")
        if ln == 0:
            for j in range(0,len(columns)):
                Cm[j] = columns[j]
        elif ln == 1:
            for j in range(0,len(columns)):
                Cv[j] = columns[j]
        else:
            T[jT] = columns[0]
            P[jT] = columns[1]
            jT = jT + 1
        ln = ln + 1
f.close()
Kc = np.zeros((nPT),dtype=np.float64)
Gc = np.zeros((nPT),dtype=np.float64)
Ec = np.zeros((nPT),dtype=np.float64)
lc = np.zeros((nPT),dtype=np.float64)
vc = np.zeros((nPT),dtype=np.float64)
Vpc = np.zeros((nPT),dtype=np.float64)
Vsc = np.zeros((nPT),dtype=np.float64)
denc = np.zeros((nPT),dtype=np.float64)
fl = 'tests/results.npy'
print('Reading results file {0:s}'.format(fl))
Compare = np.load(fl)

# MAke sure volume fractions sum to 1
if sum(Cv) < 1:
    print('Composition does not sum to one. - Exiting')
    sys.exit()

Par, MinNames, nPar, nAllMin = mv.loadPar('../database/MineralPhysicsDatabase.nc')
MinIndex = Par[0,:];

print('{0:21s}{1:20s}'.format('Mineral','Volume fraction'))
for j in range(0,nMin):
    k = mv.find(MinIndex,Cm[j]);
    print(MinNames[:,k].tobytes().decode('utf-8'),'(',Cv[j],')')
if nPT > 1:
    print('There are',nPT,'temperature and pressure points')
else:
    print('Temperature',T)
    print('Pressure',P)
print('')

K, G, E, l, v, Vp, Vs, den = mv.CalcMV(Cm,Cv,T,P);
print('Difference between calculation and expectation for 8 parameters')
print(sum(K-Compare[0,:]))
print(sum(G-Compare[1,:]))
print(sum(E-Compare[2,:]))
print(sum(l-Compare[3,:]))
print(sum(v-Compare[4,:]))
print(sum(Vp-Compare[5,:]))
print(sum(Vs-Compare[6,:]))
print(sum(den-Compare[7,:]))
print('')

sm = sum(K-Compare[0,:]) + sum(G-Compare[1,:]) + sum(E-Compare[2,:]) + sum(l-Compare[3,:]) + \
     sum(v-Compare[4,:]) + sum(Vp-Compare[5,:]) + sum(Vs-Compare[6,:]) + sum(den-Compare[7,:])
if sm == 0:
    print('Passed')
else:
    print('Did not pass')


sys.exit()
