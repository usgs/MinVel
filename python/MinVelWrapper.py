# Import the modules
import sys
import MinVel as mv
import numpy as np

# NOTES: May want to update temperature dependence of thermal expansivity using Holland and Powell's (2011) 
#        new revised equations (see figure 1 in that article). This will necessitate recalculating the first
#         Gruneisen parameters. This could provide more realistic temperature dependence of material 
#         properties within the mantle.

if len(sys.argv) > 1:
    if sys.argv[1] == "-h":
        print('MinVel -- Program to calculate mineral aggregate moduli and density')
        print('')
        print('        Written by Oliver Boyd')
        print('')
        print('        This program calculates the velocity and density of a mineral assemblage ')
        print('        at a given pressure and temperature (which may be vectors).')
        print('        The velocities are expressed as Voigt, Reuss, and Voigt-Reuss-Hill averages.')
        print('')
        print('        The data required for this analysis is taken from Hacker and Abers (2003),')
        print('        updated by Abers and Hacker in 2016, and expanded by Boyd in 2018.')
        print('        The moduli at pressure and temperature are calculated based on the')
        print('        procedures of Hacker and Abers (2004), Bina and Helffrich (1992) and')
        print('        Holland and Powell (1998) as outlined in the supplementary section of ')
        print('        Boyd et al. (2004) with updates by Abers and Hacker (2016) for quartz.')
        print('')
        print('        OUTPUT (SI Units)')
        print('                results.npy - numpy binary file containing the following vectors:')
        print('                    Voigt-Reuss-Hill averages')
        print('                        K   - Bulk modulus')
        print('                        G   - Shear modulus')
        print('                        E   - Youngs modulus')
        print('                        l   - Lambda')
        print('                        v   - Poissons ratio')
        print('                        Vp  - P-wave velocity')
        print('                        Vs  - S-wave velocity')
        print('                        p   - Density')
        print('                        a   - Thermal Expansivity')
        print('                    Voigt(v) and Reuss(r) bounds on velocity')
        print('                        Vpv - P-wave velocity')
        print('                        Vpr - P-wave velocity')
        print('                        Vsv - S-wave velocity')
        print('                        Vsr - S-wave velocity')
        print('')
        print('        INPUTS')
        print('        Command line options')
        print('        -h Help about this program.')
        print('')
        print('        -f InputFile - File containing composition, temperature, and pressure ')
        print('              information with the following format')
        print('           MinIndx 1, MinIndx 2, ..., MinIndx N')
        print('           VolFrac 1, VolFrac 2, ..., VolFrac N')
        print('           T1, P1')
        print('           T2, P2')
        print('           ...')
        print('           TN, PN')
        print('')
        print('        -p Pressure - desired pressure or comma separated vector of pressures (Pa)')
        print('        -t Temperature - desired temperature or comma separated vector of temperatures (K)')
        print('')
        print('        Composition parmeters - a composition structure with the following fields: ')
        print('        -cm Min - The mineral index comma separated vector.')
        print('        -cv Fr - Volume fraction for each mineral in Min (0 to 1), comma separated.')
        print('')
        print('          Mineral Indexes')
        print('            Quartz')
        print('                 1. Alpha Quartz                  ')
        print('                 2. Beta Quartz                   ')
        print('                 3. Coesite                       ')
        print('               Feldspar group')
        print('             Plagioclase')
        print('                 4. High Albite                   ')
        print('                 5. Low Albite                    ')
        print('                 6. Anorthite                     ')
        print('')
        print('                 7. Orthoclase                    ')
        print('                 8. Sanidine                      ')
        print('            Garnet structural group')
        print('                 9. Almandine                     ')
        print('                10. Grossular                     ')
        print('                11. Pyrope                        ')
        print('            Olivine group')
        print('                12. Forsterite                    ')
        print('                13. Fayalite                      ')
        print('            Pyroxene group')
        print('                14. Diopside                      ')
        print('                15. Enstatite                     ')
        print('                16. Ferrosilite                   ')
        print('                79. Mg-Tschermak                  ')
        print('                17. Jadeite                       ')
        print('                18. Hedenbergite                  ')
        print('                80. Acmite                        ')
        print('               81. Ca-Tschermak                  ')
        print('            Amphibole supergroup')
        print('                19. Glaucophane                   ')
        print('                20. Ferroglaucophane              ')
        print('                21. Tremolite                     ')
        print('                22. Ferroactinolite               ')
        print('                23. Tshermakite                   ')
        print('                24. Pargasite                     ')
        print('                25. Hornblende                    ')
        print('                26. Anthophyllite                 ')
        print('            Mica group')
        print('                27. Phlogopite                    ')
        print('                28. Annite                        ')
        print('                29. Muscovite                     ')
        print('                30. Celadonite                    ')
        print('            Other')
        print('                31. Talc                          ')
        print('                32. Clinochlore                   ')
        print('                33. Daphnite                      ')
        print('                34. Antigorite                    ')
        print('                35. Zoisite                       ')
        print('                36. Clinozoisite                  ')
        print('                37. Epidote                       ')
        print('                38. Lawsonite                     ')
        print('                39. Prehnite                      ')
        print('                40. Pumpellyite                   ')
        print('                41. Laumontite                    ')
        print('                42. Wairakite                     ')
        print('                43. Brucite                       ')
        print('                44. Clinohumite                   ')
        print('                45. Phase A                       ')
        print('                46. Sillimanite                   ')
        print('                47. Kyanite                       ')
        print('                48. Spinel                        ')
        print('                49. Hercynite                     ')
        print('                50. Magnetite                     ')
        print('                51. Calcite                       ')
        print('                52. Aragonite                     ')
        print('                82. Magnesite                     ')
        print('                83. En79Fs09Ts12                  ')
        print('                84. Di75He9Jd3Ts12                ')
        print('                85. ilmenite                      ')
        print('                86. cordierite                    ')
        print('                87. scapolite (meionite)          ')
        print('                88. rutile                        ')
        print('                89. sphene                        ')
        print('                53. Corundum                      ')
        print('                54. Dolomite                      ')
        print('                74. Halite                        ')
        print('                77. Pyrite                        ')
        print('                78. Gypsum   ')
        print('                90. Anhydrite   ')
        print('                 0. Water   ')
        print('                -1. Ice   ')
        print('            Clays')
        print('                55. Montmorillonite (Saz-1)')
        print('                56. Montmorillonite (S Wy-2)')
        print('                57. Montmorillonite (STX-1)')
        print('                58. Montmorillonite (S Wy-1)')
        print('                59. Montmorillonite (Shca-1)')
        print('                60. Kaolinite (Kga-2)')
        print('                61. Kaolinite (Kga-1b)')
        print('                62. Illite (IMT-2)')
        print('                63. Illite (ISMT-2)')
        print('                66. Smectite (S Wa-1)')
        print('                70. Montmorillonite (S YN-1)')
        print('                71. Chrysotile                    ')
        print('                72. Lizardite                     ')
        print('                76. Dickite                       ')
        print('')
        print('     Example:');
        print('     Geophysical parameters for 20% Quartz, 20% low Albite, 30% Forsterite, and 30% Fayalite at')
        print('     300, 400, and 500K and 0.1, 0.3, and 0.5 MPa')
        print('     > python MinVelWrapper.py -t 300,400,500 -p 0.1e6,0.3e6,0.5e6 -cm 1,5,12,13 -cv 0.2,0.2,0.3,0.3')
        print('')
        sys.exit()

nMin = 1
nPT = 1
nT = 0
nP = 0
if len(sys.argv) > 1:
    for j in range(1,len(sys.argv),2):
        if sys.argv[j] == "-t":
            entries = sys.argv[j+1].split(",")
            nT = len(entries)
            T = np.zeros((nT),dtype=np.float64)
            for k in range(0,nT):
                T[k] = entries[k]
        if sys.argv[j] == "-p":
            entries = sys.argv[j+1].split(",")
            nP = len(entries)
            P = np.zeros((nP),dtype=np.float64)
            for k in range(0,nP):
                P[k] = entries[k]
        if sys.argv[j] == "-cm":
            entries = sys.argv[j+1].split(",")
            nMin = len(entries)
            Cm = np.zeros((nMin),dtype=np.int8)
            for k in range(0,nMin):
                Cm[k] = entries[k]
        if sys.argv[j] == "-cv":
            entries = sys.argv[j+1].split(",")
            nFr = len(entries)
            Cv = np.zeros((nFr),dtype=np.float64)
            for k in range(0,nFr):
                Cv[k] = entries[k]
        if sys.argv[j] == "-f":
            fl = sys.argv[j+1]
            print('Reading {0:s}'.format(fl))
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
            nT = nPT
            nP = nPT
            nFr = nMin
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

# MAke sure volume fractions sum to 1
if sum(Cv) < 1:
    print('Composition does not sum to one. - Exiting')
    sys.exit()
if nT != nP:
    print('Number of temperature inputs must be equal to the number of pressure inputs')
    sys.exit()
else:
    nPT = nT
if nMin != nFr:
    print('Number of minerals types must be equal to the number of mineral fractional volumes')
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

K, G, E, l, v, Vp, Vs, den, Vpv, Vpr, Vsv, Vsr, a = mv.CalcMV(Cm,Cv,T,P);
print('K  ',K)
print('G  ',G)
print('E  ',E)
print('l  ',l)
print('v  ',v)
print('Vp ',Vp)
print('Vs ',Vs)
print('den',den)
print('a  ',a)
print('')
print('Voigt(v) and Reuss(r) bounds on velocity')
print('Vpv',Vpv)
print('Vpr',Vpr)
print('Vsv',Vsv)
print('Vsr',Vsr)
print('')

res = np.zeros((13,nPT),dtype=np.float64)
res[0,:] = K
res[1,:] = G
res[2,:] = E
res[3,:] = l
res[4,:] = v
res[5,:] = Vp
res[6,:] = Vs
res[7,:] = den
res[8,:] = a
res[9,:] = Vpv
res[10,:] = Vpr
res[11,:] = Vsv
res[12,:] = Vsr

f = 'results.npy'
np.save(f,res)

sys.exit()
