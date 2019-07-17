# Imp rt the modules
import sys
import numpy as np
from netCDF4 import Dataset 

def loadPar(fl):
    nc = Dataset(fl,'r')
    Par = nc.variables['Parameter Values'][:]
    MinNames = nc.variables['Mineral Name'][:]
    nPar = len(Par[:,0])
    nMin = len(Par[0,:])
    return Par,MinNames,nPar,nMin;

def find(v1,v2):
    for j in range(0,len(v1)):
        if v1[j] == v2:
            return j;
    print('Mineral not found {0:0}'.format(v2))
    sys.exit()

def CalcMV(Cm,Cv,T,P):
    nMin = len(Cm)
    nPT = len(T)

    p = np.zeros((nMin,nPT),dtype=np.float64);
    Vp = np.zeros((nMin,nPT),dtype=np.float64);
    Vs = np.zeros((nMin,nPT),dtype=np.float64);
    a = np.zeros((nMin,nPT),dtype=np.float64);
    K = np.zeros((nMin,nPT),dtype=np.float64);
    G = np.zeros((nMin,nPT),dtype=np.float64);
    V = np.zeros((nMin,nPT),dtype=np.float64);
    Cmp = np.zeros((nMin),dtype=np.float64);

    Par, MinNames, nPar, nAllMin = loadPar('../database/MineralPhysicsDatabase.nc')
    MinIndex = Par[0,:];
    for j in range(0,nMin):
        Cmp[j] = Cv[j];
        # Quartz
        if Cm[j] == 1:
            k = np.zeros((2),dtype=np.int8)
            k[0] = 0
            k[1] = 1
        else:
            k = find(MinIndex,Cm[j]);
        p[j,:], Vp[j,:], Vs[j,:], a[j,:], K[j,:], G[j,:], V[j,:] = GetV(P,T, \
                Par[22,k], Par[24,k], Par[18,k], \
                Par[5,k], Par[8,k], \
                Par[15,k], Par[20,k], \
                Par[10,k], Par[13,k], Par[3,k]);
    K, G, E, l, v, Vp, Vs, den, Vpv, Vpr, Vsv, Vsr, a = GetMV(p, K, G, a, Cmp);
    return K, G, E, l, v, Vp, Vs, den, Vpv, Vpr, Vsv, Vsr, a;

def GetV(P, T, th, dt, gamma, p, a, G, dGdPr, K, dKdPr, Vo):

    # GetV returns the density, p, P and S velocity, Vp and Vs, and thermal
    # expansivity, a, at a given and temperature and pressure, T and P. The postscript
    # TP stands for temperature and pressure.
    # It requires the following parameters:
    # 
    #        P - Pressure
    #        T - Temperature
    #        a0 - Thermal expansivity
    #        p - Density
    #        G - Shear Modulus
    #        K - Bulk Modulus
    #        d.dP - Derivative with respect to pressure
    #        dt = isothermal second gruneisen parameter
    #        th - first gruneisen parameter
    #        gamma - dlnGdp at constant pressure
    #        Vo - Volume at STP
    #
    # Equations based on Hacker and Abers (2004) with modifications to equations in GetF by Boyd et al (2004).
    # Further modifcation for quartz adapted from Abers and Hacker (2016)
    #
    nT = len(T)
    nMin = K.size

    phi = np.zeros((nT),dtype=np.float64)
    aT = np.zeros((nT),dtype=np.float64)
    pT = np.zeros((nT),dtype=np.float64)
    VT = np.zeros((nT),dtype=np.float64)
    KT = np.zeros((nT),dtype=np.float64)
    GT = np.zeros((nT),dtype=np.float64)
    dKdP = np.zeros((nT),dtype=np.float64)
    dGdP = np.zeros((nT),dtype=np.float64)
    thg = np.zeros((nT),dtype=np.float64)
    f = np.zeros((nT),dtype=np.float64)
    p0 = np.zeros((nT),dtype=np.float64)
    pP = np.zeros((nT),dtype=np.float64)
    GTP = np.zeros((nT),dtype=np.float64)
    if nMin == 2: # Dealing with quartz
        Tab = 273 + 574.3 + 0.2559*P/1e6 - 6.406e-6*(P/1e6)**2;

        dt = dt[0]*np.ones((T.shape));
        gamma = gamma[0]*np.ones((T.shape));

        for j in range(0,nT):
            # Alpha quarta
            if T[j] <= Tab[j]:
                #### Expansion fit: form of Carpenter et al. 1998
                Ttr = Tab[j];
                Tc = Ttr - 3.1;    # values from Carpenter et al.
                Vso = -5.4E-3;     # Slightly adjusted from -5.1E-3, to better fit Ohno density data
                Vs = (2/3)*Vso*(1 + np.sqrt(1 - 0.75*(T[j] - Tc)/(Ttr - Tc)));  # Carpenter et al eq 25
                Vsstp = (2/3)*Vso*(1 + np.sqrt(1 - 0.75*(298 - Tc)/(Ttr - Tc)));
                phi[j] = np.log(1 + Vs) - np.log(1 + Vsstp);  # integrated therm. expans.;
                # Derivative gives alpha(T) thermal expansion
                aT[j] = -0.25*Vso/((1 + Vs)*(Ttr - Tc)*np.sqrt(1 - 0.75*(T[j] - Tc)/(Ttr - Tc)));
        
                pT[j] = p[0]*np.exp(-phi[j]);
                VT[j] = Vo[0]*np.exp(phi[j]);

                 # regressions by Abers and Hacker (2016) on Ohno et al. 2006 data
                aktp0 = (11.384 + 0.6478*np.log(Tab[j] - 298) + 0.5878*(np.log(Tab[j] - 298)**2)) ; #Ks(P,25C)
                amup0 = (41.721 - 0.034879*np.log(Tab[j] - 298) + 0.08424*(np.log(Tab[j] - 298)**2)); #G(P,25C)
                KT[j] = (K[0]/aktp0)*(11.384 + 0.6478*np.log(Tab[j] - T[j]) + 0.5878*(np.log(Tab[j] - T[j])**2)) ; #Ks(0,T)
                GT[j] = (G[0]/amup0)*(41.721 - 0.034879*np.log(Tab[j] - T[j]) + 0.08424*(np.log(Tab[j] - T[j])**2)); # G(0,T)

                dKdP[j] = dKdPr[0];
                dGdP[j] = dGdPr[0];
                f[j] = GetF(P[j],KT[j],dKdPr[0]);
                thg[j] = th[0];
                p0[j] = p[0];
                pP[j] = p[0]*(1 + 2*f[j])**(3/2);
            else:

            #Beta quartz
                # I was uncomfortable assuming phi = 0 and aT = 0. I therefor maintained the original formulation
                # but normalized dt and gamma by 1/phi to get the same fit.

                # This is the orginal method in Abers and Hacker (2016)
                #dt[j] = -0.0963*log(T[j] - Tab[j]) + 0.0072; #fit to Ohno et al. (2006) & Ackerman et al. (1974)
                #gamma[j] = 0.0002 - 0.0002*(T[j] - Tab[j]); #fit to Ohno et al. (2006)   dGdT in Excel macro

                phi[j] = a[1]*((T[j]-298) - 20*(np.sqrt(T[j]) - np.sqrt(298)));
                aT[j] = a[1]*(1 - 10/np.sqrt(T[j]));
                dt[j] = (1/phi[j])*(-0.0963*np.log(T[j] - Tab[j]) + 0.0072); #fit to Ohno et al. (2006) & Ackerman et al. (1974)  
                gamma[j] = (1/phi[j])*(0.0002 - 0.0002*(T[j] - Tab[j])); #fit to Ohno et al. (2006) 
                pT[j] = p[1]*np.exp(-phi[j]);
                VT[j] = Vo[1]*np.exp(phi[j]);
                KT[j] = K[1]*np.exp(-dt[j]*phi[j]);
                GT[j] = G[1]*np.exp(-gamma[j]*phi[j]);
                dKdP[j] = dKdPr[1];
                dGdP[j] = dGdPr[1];
                f[j] = GetF(P[j],KT[j],dKdPr[1]);
                thg[j] = th[1];
                p0[j] = p[1];
                pP[j] = p[1]*(1 + 2*f[j])**(3/2);

    else:
        phi = a*((T-298) - 20*(np.sqrt(T) - np.sqrt(298)));
        aT = a*(1 - 10./np.sqrt(T));
        pT = p*np.exp(-phi);
        VT = Vo*np.exp(phi);
        KT = K*np.exp(-dt*phi);
        GT = G*np.exp(-gamma*phi);
        dKdP = dKdPr*np.ones((nT),dtype=np.float64);
        dGdP = dGdPr*np.ones((nT),dtype=np.float64);
        f = GetF(P,KT,dKdP);
        thg = th;
        p0 = p;
        pP = p*(1 + 2*f)**(3/2);

    VTP = VT*(1 + 2*f)**(-3/2);

    #Assume ddKddP is 1/KT and ddGddP is 1/GT as suggested by Bina and Helffrich, 1992
    KTP = KT*((1+2*f)**(5/2))*(1 - f*(5 - 3*dKdP) + (f**2)*(9 + (3*dKdP - 7)*(3*dKdP - 5))/2);
    for j in range(0,nT):
        if GT[j] > 0:
            GTP[j] = GT[j]*((1+2*f[j])**(5/2))*(1 - f[j]*(5 - 3*dGdP[j]*KT[j]/GT[j]) + \
                    (f[j]**2)*(9*KT[j]**2/GT[j]**2 + 9*(dKdP[j] - 4)*dGdP[j]*KT[j]/GT[j] + 35)/2);
        else:
            GTP[j] = 0;

    aTP = aT*((pP/p0)**-dt);
    pTP = pP*pT/p0;
    Ks = KTP*(1 + T*thg*aTP);

    Vp = np.sqrt((Ks + 4*GTP/3)/pTP);
    Vs = np.sqrt(GTP/pTP);

    return pTP, Vp, Vs, aTP, Ks, GTP, VTP;

def GetF(P,KT,dKdP):
    #function f = GetF(P,KT,dKdP)

    # From Bina and Helffrich, 1992. Assume second derivative of K with respect to P,
    # ddKddP, is 1/dKdP as suggested by Bina and Helffrich, else
    #
    # y = (9*dKdP*ddKddP + 4*e*(4 - 3*dKdP) + 5*(3*dKdP - 5))/6;
    #
    # Took Eq from Bina and Helfrich, Turned into polynomial, and used matlab roots
    # command to get the roots of the polynomial. I use solution closest to First
    # order solution.
    #

    e = 0.75*(4-dKdP);
    y = (9 + 4*e*(4 - 3*dKdP) + 5*(3*dKdP - 5))/6;

    y2 = y**2;
    e2 = e**2;
    ey = e*y;

    if P.size == 1:
        P = P*np.ones((1),dtype=np.float64)
        KT = KT*np.ones((1),dtype=np.float64)

    nPT = len(P)
    C = np.zeros((12,nPT),dtype=np.float64)
    f = np.zeros(((nPT)),dtype=np.float64)

    C[0,:] = 32*y2;
    C[1,:] = 80*y2 - 128*ey;
    C[2,:] = 80*y2 - 320*ey + 64*y + 128*e2;
    C[3,:] = 40*y2 - 320*ey + 160*y + 320*e2 - 128*e;
    C[4,:] = 10*y2 - 160*ey + 160*y + 320*e2 - 320*e + 32;
    C[5,:] = y2 - 40*ey + 80*y + 160*e2 - 320*e + 80;
    C[6,:] = -4*ey + 20*y + 40*e2 - 160*e + 80;
    C[7,:] = 2*y + 4*e2 - 40*e + 40;
    C[8,:] = 10 - 4*e;
    C[9,:] = 1;
    C[10,:] = 0;
    C[11,:] = -((P/KT)**2)/9;
    for j in range(0,nPT):
        F1 = np.sqrt(-C[11,j]);
        F = np.roots(C[:,j]);
        abF = np.abs(F-F1)
        mndF = np.min(abF)
        imF = np.imag(F)
        for l in range(0,12):
            if (imF[l] == 0 and abF[l] == mndF):
                f[j] = F[l].real;
                break

    return f;

def GetMV(p, k, g, a, Cmp):

    # GetMV calculates the Voigt-Reuss averages (of velocities but should be moduli
    # though VRH will be about the same) of the s and p-wave velocity (Anderson 1997).

    nMin,nPT = k.shape

    Kv = np.zeros((nPT),dtype=np.float64)
    Kr = np.zeros((nPT),dtype=np.float64)
    Gv = np.zeros((nPT),dtype=np.float64)
    Gr = np.zeros((nPT),dtype=np.float64)
    iKr = np.zeros((nPT),dtype=np.float64)
    iGr = np.zeros((nPT),dtype=np.float64)
    K = np.zeros((nPT),dtype=np.float64)
    G = np.zeros((nPT),dtype=np.float64)
    l = np.zeros((nPT),dtype=np.float64)
    v = np.zeros((nPT),dtype=np.float64)
    E = np.zeros((nPT),dtype=np.float64)
    Vpv = np.zeros((nPT),dtype=np.float64)
    Vpr = np.zeros((nPT),dtype=np.float64)
    Vsv = np.zeros((nPT),dtype=np.float64)
    Vsr = np.zeros((nPT),dtype=np.float64)
    Vp = np.zeros((nPT),dtype=np.float64)
    Vs = np.zeros((nPT),dtype=np.float64)
    den = np.zeros((nPT),dtype=np.float64)
    al = np.zeros((nPT),dtype=np.float64)

    Kv = Cmp[0]*k[0,:]
    Gv = Cmp[0]*g[0,:]
    iKr = Cmp[0]/k[0,:]
    iGr = Cmp[0]/g[0,:]
    den = Cmp[0]*p[0,:]
    al = Cmp[0]*a[0,:]
    for j in range(1,nMin):
        den = den + Cmp[j]*p[j,:]
        al = al + Cmp[j]*a[j,:]
        # Voigt Average
        Kv = Kv + Cmp[j]*k[j,:]
        Gv = Gv + Cmp[j]*g[j,:]
        # Reuss Average
        iKr = iKr + Cmp[j]/k[j,:]
        iGr = iGr + Cmp[j]/g[j,:]

    Vpv = np.sqrt((Kv + 4*Gv/3)/den)
    Vsv = np.sqrt(Gv/den)
    Kr = 1/iKr
    Gr = 1/iGr
    Vpr = np.sqrt((Kr + 4*Gr/3)/den)
    Vsr = np.sqrt(Gr/den)

    # Voigt Reuss Hill Average
    Vp = (Vpv + Vpr)/2
    Vs = (Vsv + Vsr)/2
    K = (Kv + Kr)/2
    G = (Gv + Gr)/2

    l = K-2*G/3
    v = l/(2*(l+G))
    E = 2*G*(1+v)

    return K, G, E, l, v, Vp, Vs, den, Vpv, Vpr, Vsv, Vsr, al;
