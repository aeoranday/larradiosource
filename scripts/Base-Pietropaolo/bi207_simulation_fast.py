import math
import random
import numpy as np


TOTAL_EMISSIONS: int = 0
TOTAL_ELECTRON_EMISSIONS: int = 0

TOTAL_OUTSIDE_EMISSIONS: int = 0

# -------------------------------------------------
# Local bindings (major speedup)
np.random.seed(42)
rand = np.random.rand
log = np.log
sqrt = np.sqrt
sin = np.sin
cos = np.cos
exp = np.exp
pi2 = 2.0 * np.pi

# -------------------------------------------------
# Constants

drlong = 205.0
drshort = 45.0
rainn = 15.0
raout = 30.0

resINNL = resOUTL = resINNS = resOUTS = 0.08

dtave = 2 * 1000.0 / (18.0 + 1.4)
srcfrac = 1.4 / (18.0 + 1.4)
drfiv = 1.0 / 1.5
attdist = 3000.0 / drfiv
thr = 0.26

# -------------------------------------------------
# Sources

Abranches = (
    (0.0015, 0.0, 0.482, 1),
    (0.0059, 0.0, 0.556, 1),
    (0.0211, 0.0, 0.566, 1),
    (1.0,   80.0, 0.570, 0),
)

Bbranches = (
    (0.0054, 0.0,   1.060, 1),
    (0.0238, 0.0,   1.049, 1),
    (0.0941, 0.0,   0.976, 1),
    (0.8399, 120.0, 1.064, 0),
    (0.9086, 160.0, 1.77,  0),
    (0.9088, 0.0,   1.682, 1),
    (1.0,    140.0, 1.440, 0),
)

# -------------------------------------------------
def pick_branch(branches):
    r = rand()
    for p, d, e, ie in branches:
        if r <= p:
            return d, e, ie
    return 0.0, 0.0, 0


def random_disk():
    while True:
        x = 5.0 * rand() - 2.5
        y = 5.0 * rand() - 2.5
        if x * x + y * y < 6.25:
            return x, y


def gaussian_noise(sigma):
    s = 0.0
    for _ in range(12):
        s += rand()
    return sigma * (s - 6.0)


def compton(E):
    gg = E / 0.511
    while True:
        ct = 2.0 * rand() - 1.0
        eps = 1.0 / (1.0 + gg * (1.0 - ct))
        if rand() <= 0.5 * eps * eps * (eps + 1/eps - (1.0 - ct * ct)):
            return E * (1.0 - eps)


def simulate_hit(Dinter, Energy, isele, drmax, ilosh):
    if Energy <= 0.0:
        return 0, 0, 0, 0, 0

    x0, y0 = random_disk()

    d = -Dinter * log(rand()) if Dinter > 0 else 0.0
    cth = rand()
    phi = pi2 * rand()
    sth = sqrt(1.0 - cth * cth)

    x1 = x0 + d * sth * sin(phi)
    y1 = y0 + d * sth * cos(phi)
    z1 = d * cth
    r1 = sqrt(x1 * x1 + y1 * y1)

    if not (0.0 <= z1 < drmax and r1 < raout):
#        if r1 > 1000:
#            print(z1, r1)
#            input()
        global TOTAL_OUTSIDE_EMISSIONS
        TOTAL_OUTSIDE_EMISSIONS += 1
        return 0, 0, 0, 0, 0

    Ee = Energy if isele else compton(Energy)
    Ee *= exp((z1 - drmax) / attdist)

    if r1 <= rainn:
        Iok = 2
        resol = resINNL if ilosh else resINNS
    else:
        Iok = 1
        resol = resOUTL if ilosh else resOUTS

    rr = Ee + gaussian_noise(resol)
    return Iok, z1, r1, Ee, rr


# -------------------------------------------------
# Main loop

with open("bi207stream.txt", "w") as fout:

    timtot = 0.0

    for i in range(1, 100_001):

        if rand() > 0:#srcfrac:
            drmax = drlong
            ilosh = 1
        else:
            drmax = drshort
            ilosh = 0

        timtot += -dtave * log(rand())

        DA, EA, ieA = pick_branch(Abranches)
        DB, EB, ieB = pick_branch(Bbranches)
        TOTAL_EMISSIONS += 2
        if ieB:
            TOTAL_ELECTRON_EMISSIONS += 1
        if ieA:
            TOTAL_ELECTRON_EMISSIONS += 1

        IA, zA, rA, EeA, rrA = simulate_hit(DA, EA, ieA, drmax, ilosh)
        IB, zB, rB, EeB, rrB = simulate_hit(DB, EB, ieB, drmax, ilosh)

        if rrA < thr:
            IA = 0
        if rrB < thr:
            IB = 0

        if IA:
            tA = timtot + (drmax - zA) * drfiv
        if IB:
            tB = timtot + (drmax - zB) * drfiv

        if IA and IB:
            if zB > zA:
                fout.write(f"{i} {ilosh} {IB} {ieB} {tB} {rB} {zB} {EB} {EeB} {rrB} 1\n")
                fout.write(f"{i} {ilosh} {IA} {ieA} {tA} {rA} {zA} {EA} {EeA} {rrA} 0\n")
            else:
                fout.write(f"{i} {ilosh} {IA} {ieA} {tA} {rA} {zA} {EA} {EeA} {rrA} 0\n")
                fout.write(f"{i} {ilosh} {IB} {ieB} {tB} {rB} {zB} {EB} {EeB} {rrB} 1\n")
        elif IA:
            fout.write(f"{i} {ilosh} {IA} {ieA} {tA} {rA} {zA} {EA} {EeA} {rrA} 0\n")
        elif IB:
            fout.write(f"{i} {ilosh} {IB} {ieB} {tB} {rB} {zB} {EB} {EeB} {rrB} 1\n")


print(f"Total Emissions: {TOTAL_EMISSIONS}. Total Electron Emissions: {TOTAL_ELECTRON_EMISSIONS}.")
print(f"Outside Emissions: {TOTAL_OUTSIDE_EMISSIONS}.")
print("Done")
