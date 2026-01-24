import math
import numpy as np

# =================================================
# Configuration
# =================================================

NTRACE = 4000
NSIGNAL = 500
NSPC = 400

thr = 0.26
ioffset0 = 1500

# =================================================
# Signal shaping function (precomputed)
# =================================================

signal = np.zeros(NSIGNAL)
for k in range(NSIGNAL):
    x = (k + 1) * 0.1
    signal[k] = math.exp(-(x - 10.0) / 5.0) / (
        0.56988 * (1.0 + math.exp(-(x - 10.0) / 1.25))
    )

# =================================================
# Traces and spectra
# =================================================

trace1 = np.zeros(NTRACE)  # outer long
trace2 = np.zeros(NTRACE)  # inner long
trace3 = np.zeros(NTRACE)  # outer short
trace4 = np.zeros(NTRACE)  # inner short

spcOUTL = np.zeros(NSPC, dtype=int)
spcEINL = np.zeros(NSPC, dtype=int)
spcGINL = np.zeros(NSPC, dtype=int)

spcOUTS = np.zeros(NSPC, dtype=int)
spcEINS = np.zeros(NSPC, dtype=int)
spcGINS = np.zeros(NSPC, dtype=int)

# =================================================
# Open files
# =================================================

fin = open("bi207stream.txt", "r")

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

ifil = 80  # diagnostic waveform output index

# =================================================
# Main loop over hits
# =================================================

for line in fin:

    vals = line.split()
    if len(vals) < 10:
        continue

    i = int(vals[0])
    ilosh = int(vals[1])
    iok = int(vals[2])
    iele = int(vals[3])
    dtim = float(vals[4])
    r1 = float(vals[5])
    z1 = float(vals[6])
    energy = float(vals[7])
    electe = float(vals[8])
    rr = float(vals[9])

    if i > 999999:
        break

    kstart = int(10 * (dtim - timtrig) + ioffset)

    # ---------------------------------------------
    # Same event (pile-up window)
    # ---------------------------------------------
    if dtim - timtrig <= 200.0:

        imult += 1

        if ilosh == 1:  # long drift
            if iok == 1:  # outer
                if rr > phmaxOUTL:
                    phmaxOUTL = rr
                    ieleOUTL = iele
                kk = kstart
                for k in range(NSIGNAL):
                    kk += 1
                    if kk < NTRACE:
                        trace1[kk] += signal[k] * rr

            elif iok == 2:  # inner
                if rr > phmaxINNL:
                    phmaxINNL = rr
                    ieleINNL = iele
                kk = kstart
                for k in range(NSIGNAL):
                    kk += 1
                    if kk < NTRACE:
                        trace2[kk] += signal[k] * rr

        else:  # short drift
            if iok == 1:
                if rr > phmaxOUTS:
                    phmaxOUTS = rr
                    ieleOUTS = iele
                kk = kstart
                for k in range(NSIGNAL):
                    kk += 1
                    if kk < NTRACE:
                        trace3[kk] += signal[k] * rr

            elif iok == 2:
                if rr > phmaxINNS:
                    phmaxINNS = rr
                    ieleINNS = iele
                kk = kstart
                for k in range(NSIGNAL):
                    kk += 1
                    if kk < NTRACE:
                        trace4[kk] += signal[k] * rr

    # ---------------------------------------------
    # Trigger reset
    # ---------------------------------------------
    if dtim - timtrig <= 0.0:
        timtrig = dtim
        ioffset = kstart

    # ---------------------------------------------
    # Event ends → analyze traces
    # ---------------------------------------------
    if dtim - timtrig > 200.0:

        phOUTL = np.max(trace1[ioffset:ioffset + 2200])
        phINNL = np.max(trace2[ioffset:ioffset + 2200])
        phOUTS = np.max(trace3[ioffset:ioffset + 2200])
        phINNS = np.max(trace4[ioffset:ioffset + 2200])

        # Optional waveform dump
        if imult > 4 and ifil < 90:
            for k in range(NTRACE):
                pass  # waveform dump omitted for clarity
            ifil += 1

        imult = 1

        # -----------------------------------------
        # Spectra and output
        # -----------------------------------------
        if phmaxOUTL > thr:
            fOUTL.write(f"{i} {ieleOUTL} {timtrig} {phmaxOUTL} {phOUTL}\n")
            k = int(phOUTL * 200)
            if k < NSPC:
                spcOUTL[k] += 1

        if phmaxINNL > thr:
            fINNL.write(f"{i} {ieleINNL} {timtrig} {phmaxINNL} {phINNL}\n")
            k = int(phINNL * 200)
            if k < NSPC:
                if ieleINNL == 1:
                    spcEINL[k] += 1
                else:
                    spcGINL[k] += 1

        if phmaxOUTS > thr:
            fOUTS.write(f"{i} {ieleOUTS} {timtrig} {phmaxOUTS} {phOUTS}\n")
            k = int(phOUTS * 200)
            if k < NSPC:
                spcOUTS[k] += 1

        if phmaxINNS > thr:
            fINNS.write(f"{i} {ieleINNS} {timtrig} {phmaxINNS} {phINNS}\n")
            k = int(phINNS * 200)
            if k < NSPC:
                if ieleINNS == 1:
                    spcEINS[k] += 1
                else:
                    spcGINS[k] += 1

        # -----------------------------------------
        # Reset for next event
        # -----------------------------------------
        timtrig = dtim
        ioffset = ioffset0

        trace1[:] = 0.0
        trace2[:] = 0.0
        trace3[:] = 0.0
        trace4[:] = 0.0

        phmaxOUTL = phmaxINNL = 0.0
        phmaxOUTS = phmaxINNS = 0.0

        ieleOUTL = ieleINNL = -1
        ieleOUTS = ieleINNS = -1

# =================================================
# Write final spectra
# =================================================

for k in range(NSPC):
    energy = k * 0.005 - 0.0025
    fspec.write(
        f"{energy} "
        f"{spcOUTL[k]} {spcEINL[k]} {spcGINL[k]} {spcEINL[k] + spcGINL[k]} "
        f"{spcOUTS[k]} {spcEINS[k]} {spcGINS[k]} {spcEINS[k] + spcGINS[k]}\n"
    )

# =================================================
# Close files
# =================================================

fin.close()
fOUTL.close()
fINNL.close()
fOUTS.close()
fINNS.close()
fspec.close()
