import math
import random

# -------------------------------------------------
# Helper function to mimic Fortran rndm(x)
def rndm():
    return random.random()

# -------------------------------------------------
# Constants and parameters

drlong = 205.0       # PM drift length in mm
drshort = 45.0
rainn = 15.0         # inner anode radius
raout = 30.0         # outer anode radius

resINNL = 0.08       # MeV
resOUTL = 0.08
resINNS = 0.08
resOUTS = 0.08

dtave = 2 * 1000.0 / (18.0 + 1.4)  # average time between events (us)
srcfrac = 1.4 / (18.0 + 1.4)
drfiv = 1.0 / 1.5                  # us/mm
attdist = 3000.0 / drfiv           # attenuation distance (mm)
thr = 0.26

# A source probabilities
Aprob1 = 0.0015
Aprob2 = 0.0044 + Aprob1
Aprob3 = 0.0152 + Aprob2

EleEne1 = 0.482
EleEne2 = 0.556
EleEne3 = 0.566
GamEne1 = 0.570
GamInt1 = 80.0

# B source probabilities
Bprob1 = 0.0054             # 1060 keV
Bprob2 = 0.0184 + Bprob1    # 1049 keV
Bprob3 = 0.0703 + Bprob2    #  976 keV
Bprob4 = 0.7458 + Bprob3    # 1064 keV gam
Bprob5 = 0.0687 + Bprob4    # 1770 keV gam
Bprob6 = 0.0002 + Bprob5    # 1682 keV
Bprob7 = 0.0013 + Bprob6    # 1440 keV gam

EleEne4 = 1.060
EleEne5 = 1.049
EleEne6 = 0.976
GamEne2 = 1.064
GamInt2 = 120.0
GamEne3 = 1.77
GamInt3 = 160.0
EleEne7 = 1.682
GamEne4 = 1.440
GamInt4 = 140.0

print(Aprob3, Bprob7)

pi2 = 2.0 * math.pi

# -------------------------------------------------
# Output file
fout = open("bi207stream.txt", "w")

ilosh = 0
timtot = 0.0
timprv = 0.0

# -------------------------------------------------
# Main loop
for i in range(1, 1_000_001):

    # Long / short drift
    if rndm() > srcfrac:
        drmax = drlong
        ilosh = 1
    else:
        drmax = drshort
        ilosh = 0

    # Time increment
    dtim = -dtave * math.log(rndm())
    timtot += dtim

    # -------------------------------
    # A decay
    Aprob = rndm()
    if Aprob <= Aprob1:
        DinterA = 0.0
        EnergyA = EleEne1
        ieleA = 1
    elif Aprob <= Aprob2:
        DinterA = 0.0
        EnergyA = EleEne2
        ieleA = 1
    elif Aprob <= Aprob3:
        DinterA = 0.0
        EnergyA = EleEne3
        ieleA = 1
    else:
        DinterA = GamInt1
        EnergyA = GamEne1
        ieleA = 0

    # Random starting point
    while True:
        x0 = 5.0 * rndm() - 2.5
        y0 = 5.0 * rndm() - 2.5
        if x0 * x0 + y0 * y0 < 6.25:
            break
    z0 = 0.0

    d = -DinterA * math.log(rndm())
    ctheta = rndm()
    phi = rndm() * pi2
    stheta = math.sqrt(1.0 - ctheta * ctheta)

    x1 = x0 + d * stheta * math.sin(phi)
    y1 = y0 + d * stheta * math.cos(phi)
    z1A = z0 + d * ctheta
    r1A = math.sqrt(x1 * x1 + y1 * y1)

    IokA = 0
    rrA = 0.0
    ElectEA = 0.0

    if EnergyA > 0.0:
        if 0.0 <= z1A < drmax and r1A < raout:

            if ieleA == 0:
                gg = EnergyA / 0.511
                while True:
                    ct = 2.0 * rndm() - 1.0
                    eps = 1.0 / (1.0 + gg * (1.0 - ct))
                    sctprob = 0.5 * eps**2 * (eps + 1/eps - (1.0 - ct**2))
                    if rndm() <= sctprob:
                        GammaE = EnergyA * eps
                        ElectEA = EnergyA - GammaE
                        break
            else:
                ElectEA = EnergyA

            ElectEA *= math.exp((z1A - drmax) / attdist)

            if r1A <= rainn:
                IokA = 2
                resol = resINNL if ilosh == 1 else resINNS
            else:
                IokA = 1
                resol = resOUTL if ilosh == 1 else resOUTS

            res = sum(rndm() for _ in range(12))
            res = resol * (res - 6.0)
            rrA = ElectEA + res

    # -------------------------------
    # B decay (same structure)
    Bprob = rndm()
    if Bprob <= Bprob1:
        DinterB = 0.0
        EnergyB = EleEne4
        ieleB = 1
    elif Bprob <= Bprob2:
        DinterB = 0.0
        EnergyB = EleEne5
        ieleB = 1
    elif Bprob <= Bprob3:
        DinterB = 0.0
        EnergyB = EleEne6
        ieleB = 1
    elif Bprob <= Bprob4:
        DinterB = GamInt2
        EnergyB = GamEne2
        ieleB = 0
    elif Bprob <= Bprob5:
        DinterB = GamInt3
        EnergyB = GamEne3
        ieleB = 0
    elif Bprob <= Bprob6:
        DinterB = 0.0
        EnergyB = EleEne7
        ieleB = 1
    elif Bprob <= Bprob7:
        DinterB = GamInt4
        EnergyB = GamEne4
        ieleB = 0
    else:
        DinterB = 0.0
        EnergyB = 0.0
        ieleB = 0

    while True:
        x0 = 5.0 * rndm() - 2.5
        y0 = 5.0 * rndm() - 2.5
        if x0 * x0 + y0 * y0 < 6.25:
            break
    z0 = 0.0

    d = -DinterB * math.log(rndm())
    ctheta = rndm()
    phi = rndm() * pi2
    stheta = math.sqrt(1.0 - ctheta * ctheta)

    x1 = x0 + d * stheta * math.sin(phi)
    y1 = y0 + d * stheta * math.cos(phi)
    z1B = z0 + d * ctheta
    r1B = math.sqrt(x1 * x1 + y1 * y1)

    IokB = 0
    rrB = 0.0
    ElectEB = 0.0

    if EnergyB > 0.0:
        if 0.0 <= z1B < drmax and r1B < raout:

            if ieleB == 0:
                gg = EnergyB / 0.511
                while True:
                    ct = 2.0 * rndm() - 1.0
                    eps = 1.0 / (1.0 + gg * (1.0 - ct))
                    sctprob = 0.5 * eps**2 * (eps + 1/eps - (1.0 - ct**2))
                    if rndm() <= sctprob:
                        GammaE = EnergyB * eps
                        ElectEB = EnergyB - GammaE
                        break
            else:
                ElectEB = EnergyB

            ElectEB *= math.exp((z1B - drmax) / attdist)

            if r1B <= rainn:
                IokB = 2
                resol = resINNL if ilosh == 1 else resINNS
            else:
                IokB = 1
                resol = resOUTL if ilosh == 1 else resOUTS

            res = sum(rndm() for _ in range(12))
            res = resol * (res - 6.0)
            rrB = ElectEB + res

    # -------------------------------
    # Threshold and output
    if rrA < thr:
        IokA = 0
    if rrB < thr:
        IokB = 0

    dtimA = timtot + (drmax - z1A) * drfiv
    dtimB = timtot + (drmax - z1B) * drfiv

    if IokA != 0 and IokB != 0:
        if z1B > z1A:
            fout.write(f"{i} {ilosh} {IokB} {ieleB} {dtimB} {r1B} {z1B} {EnergyB} {ElectEB} {rrB} 1\n")
            fout.write(f"{i} {ilosh} {IokA} {ieleA} {dtimA} {r1A} {z1A} {EnergyA} {ElectEA} {rrA} 0\n")
        else:
            fout.write(f"{i} {ilosh} {IokA} {ieleA} {dtimA} {r1A} {z1A} {EnergyA} {ElectEA} {rrA} 0\n")
            fout.write(f"{i} {ilosh} {IokB} {ieleB} {dtimB} {r1B} {z1B} {EnergyB} {ElectEB} {rrB} 1\n")
    elif IokA != 0:
        fout.write(f"{i} {ilosh} {IokA} {ieleA} {dtimA} {r1A} {z1A} {EnergyA} {ElectEA} {rrA} 0\n")
    elif IokB != 0:
        fout.write(f"{i} {ilosh} {IokB} {ieleB} {dtimB} {r1B} {z1B} {EnergyB} {ElectEB} {rrB} 1\n")

# -------------------------------------------------
fout.close()
print(i, timtot)
