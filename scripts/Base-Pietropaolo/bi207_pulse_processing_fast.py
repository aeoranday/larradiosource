import numpy as np
import math

# =================================================
# Constants
# =================================================

NTRACE = 4000
NSIGNAL = 500
NSPC = 400

thr = 0.26
ioffset0 = 1500
WINDOW = 200.0
DT_SCALE = 10.0

# =================================================
# Precompute shaping signal (vectorized)
# =================================================

x = (np.arange(1, NSIGNAL + 1) * 0.1)
signal = np.exp(-(x - 10.0) / 5.0) / (0.56988 * (1.0 + np.exp(-(x - 10.0) / 1.25)))

# =================================================
# Traces & spectra
# =================================================

trace1 = np.zeros(NTRACE)  # outer long
trace2 = np.zeros(NTRACE)  # inner long
trace3 = np.zeros(NTRACE)  # outer short
trace4 = np.zeros(NTRACE)  # inner short

spcOUTL = np.zeros(NSPC, dtype=np.int32)
spcEINL = np.zeros(NSPC, dtype=np.int32)
spcGINL = np.zeros(NSPC, dtype=np.int32)
spcOUTS = np.zeros(NSPC, dtype=np.int32)
spcEINS = np.zeros(NSPC, dtype=np.int32)
spcGINS = np.zeros(NSPC, dtype=np.int32)

# =================================================
# File handles
# =================================================

fin = open("bi207stream.txt")
fOUTL = open("outer_anode_long.txt", "w")
fINNL = open("inner_anode_long.txt", "w")
fOUTS = open("outer_anode_short.txt", "w")
fINNS = open("inner_anode_short.txt", "w")
fspec = open("bi207spectra.txt", "w")

# =================================================
# State variables
# =================================================

timtrig = 0.0
ioffset = ioffset0
imult = 0

phmaxOUTL = phmaxINNL = 0.0
phmaxOUTS = phmaxINNS = 0.0

ieleOUTL = ieleINNL = -1
ieleOUTS = ieleINNS = -1

# =================================================
# Helper: add pulse (vectorized)
# =================================================

def add_pulse(trace, kstart, amplitude):
    k0 = max(kstart, 0)
    k1 = min(kstart + NSIGNAL, NTRACE)
    s0 = k0 - kstart
    s1 = s0 + (k1 - k0)
    trace[k0:k1] += signal[s0:s1] * amplitude

# =================================================
# Main loop
# =================================================

for line in fin:
    i, ilosh, iok, iele, dtim, r1, z1, energy, electe, rr = line.split()
    i = int(i)
    ilosh = int(ilosh)
    iok = int(iok)
    iele = int(iele)
    dtim = float(dtim)
    rr = float(rr)

    if i > 999999:
        break

    kstart = int(DT_SCALE * (dtim - timtrig) + ioffset)

    # ---------------------------------------------
    # Same trigger window
    # ---------------------------------------------
    if dtim - timtrig <= WINDOW:
        imult += 1

        if ilosh == 1:
            if iok == 1:
                if rr > phmaxOUTL:
                    phmaxOUTL, ieleOUTL = rr, iele
                add_pulse(trace1, kstart, rr)
            elif iok == 2:
                if rr > phmaxINNL:
                    phmaxINNL, ieleINNL = rr, iele
                add_pulse(trace2, kstart, rr)

        else:
            if iok == 1:
                if rr > phmaxOUTS:
                    phmaxOUTS, ieleOUTS = rr, iele
                add_pulse(trace3, kstart, rr)
            elif iok == 2:
                if rr > phmaxINNS:
                    phmaxINNS, ieleINNS = rr, iele
                add_pulse(trace4, kstart, rr)

    # ---------------------------------------------
    # New trigger
    # ---------------------------------------------
    if dtim - timtrig <= 0.0:
        timtrig = dtim
        ioffset = kstart

    # ---------------------------------------------
    # Event ends → analyze
    # ---------------------------------------------
    if dtim - timtrig > WINDOW:

        sl = slice(ioffset, min(ioffset + 2200, NTRACE))
        phOUTL = trace1[sl].max()
        phINNL = trace2[sl].max()
        phOUTS = trace3[sl].max()
        phINNS = trace4[sl].max()

        imult = 1

        if phmaxOUTL > thr:
            fOUTL.write(f"{i} {ieleOUTL} {timtrig} {phmaxOUTL} {phOUTL}\n")
            k = int(phOUTL * 200)
            if k < NSPC:
                spcOUTL[k] += 1

        if phmaxINNL > thr:
            fINNL.write(f"{i} {ieleINNL} {timtrig} {phmaxINNL} {phINNL}\n")
            k = int(phINNL * 200)
            if k < NSPC:
                (spcEINL if ieleINNL == 1 else spcGINL)[k] += 1

        if phmaxOUTS > thr:
            fOUTS.write(f"{i} {ieleOUTS} {timtrig} {phmaxOUTS} {phOUTS}\n")
            k = int(phOUTS * 200)
            if k < NSPC:
                spcOUTS[k] += 1

        if phmaxINNS > thr:
            fINNS.write(f"{i} {ieleINNS} {timtrig} {phmaxINNS} {phINNS}\n")
            k = int(phINNS * 200)
            if k < NSPC:
                (spcEINS if ieleINNS == 1 else spcGINS)[k] += 1

        # Reset
        timtrig = dtim
        ioffset = ioffset0
        trace1.fill(0)
        trace2.fill(0)
        trace3.fill(0)
        trace4.fill(0)
        phmaxOUTL = phmaxINNL = phmaxOUTS = phmaxINNS = 0.0
        ieleOUTL = ieleINNL = ieleOUTS = ieleINNS = -1

# =================================================
# Write spectra
# =================================================

for k in range(NSPC):
    E = k * 0.005 - 0.0025
    fspec.write(
        f"{E} {spcOUTL[k]} {spcEINL[k]} {spcGINL[k]} {spcEINL[k]+spcGINL[k]} "
        f"{spcOUTS[k]} {spcEINS[k]} {spcGINS[k]} {spcEINS[k]+spcGINS[k]}\n"
    )

fin.close()
fOUTL.close()
fINNL.close()
fOUTS.close()
fINNS.close()
fspec.close()
