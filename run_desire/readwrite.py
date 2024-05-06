# Functions to interact with SIR/DeSIRe input and output files

# Author: Carlos Quintero Noda
# Revised by Basilio Ruiz Cobo
# Version: 1.0
# Date: 23 - June - 2022

# We need to load some packages
import numpy as np

# -------------------------------------------------------------------------------------------------
# -------------------------------- Reading an atmosphere ------------------------------------------
# -------------------------------------------------------------------------------------------------

# This part will read an atmosphere in SIR/DeSIRe format and will return the header
# (atmospheric parameters that do not change with height), and the height dependent atmospheric
# parameters.

# Running example:
# >>> from readwrite import readmod
# >>> atm0,atm1 = readmod('FALC.mod')

def readmod(name):

    header = np.loadtxt(name,skiprows=0,max_rows=1)
    atm = np.loadtxt(name,skiprows=1)

    return header,atm

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------
# -------------------------------- Writing an atmosphere ------------------------------------------
# -------------------------------------------------------------------------------------------------
# This function will write an atmosphere in SIR/DeSIRe format using as input the
# header and the height dependent atmospheric parameters

# Running example:
# >>> from readwrite import writemod
# >>> writemod(atm0,atm1,'FALCnb.mod')

def writemod(headeratm,atm1,nameatm):
    
    format = ['%.3f','%.4f','%.6e','%.6e','%.4f','%.4e','%.4f','%.4f','%.6f','%.6e','%.9e']

    with open(nameatm, 'wb') as f:
        np.savetxt(f,headeratm[None,:], delimiter='    ', fmt='%.5f')
        np.savetxt(f, atm1, delimiter='    ',fmt=format)

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------        


# -------------------------------------------------------------------------------------------------
# ------------------------------------ Reading a profile ------------------------------------------
# -------------------------------------------------------------------------------------------------

# The function reads a profile in SIR/DeSIRe format and will return an array containing
# the line index used in the *.grid file, the wavelength vector and the four Stokes parameters. 

# Running example:
# >>> from readwrite import readprof
# >>> profiles = readprof('profiles.per')

def readprof(name):

    prof = np.loadtxt(name,skiprows=0)

    return prof

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------
# ------------------------------------ Writing a profile ------------------------------------------
# -------------------------------------------------------------------------------------------------
# This part will write a profile in SIR/DeSIRe format using as input the
# 6 x wavelength array containing the line index used in the *.grid file, the wavelength
# vector and the four Stokes parameters. 

# Running example:
# >>> from readwrite import writeprof
# >>> writeprof(profiles,'profilesb.per')

def writeprof(prof,nameatm):
    
    format = ['%.1f','%.4f','%.8e','%.8e','%.8e','%.8e']

    with open(nameatm, 'wb') as f:
        np.savetxt(f, prof, delimiter='    ',fmt=format)          

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------
# ------------------------------ Reading Response Functions ---------------------------------------
# -------------------------------------------------------------------------------------------------
# This function will read the Response Functions computed with SIR/DeSIRe. The output will be four
# 2D arrays containing the Stokes (I,Q,U,V) RFs that depend on height and wavelength.
# 

# Running example:
# >>> from readwrite import readrf
# >>> SiRF, SqRF, SuRF, SvRF = readrf(name)

def readrf(name):

    # The information about n_tau and n_wave are in the first row

    header = np.loadtxt(name,skiprows=0,max_rows=1)
    ntau = header[0]
    nwave = header[1] 

    # RFs are written in specific way where we have that rows correspond to the Stokes I, Q, U and
    # V response functions, one after the other for all optical depths and wavelength points. Thus,
    # reading them is not as straightforwards as with the atmospheres and the profiles. Still it
    # can be done rather easily.

    RF = np.loadtxt(name,skiprows=1)
    RF2 = RF.reshape(int(ntau),int(nwave),order='C')

    nwaveS = int(nwave/4.)

    StokesIRF = RF2[:,int(nwaveS*0.):int(nwaveS*1.)]
    StokesQRF = RF2[:,int(nwaveS*1.):int(nwaveS*2.)]
    StokesURF = RF2[:,int(nwaveS*2.):int(nwaveS*3.)]
    StokesVRF = RF2[:,int(nwaveS*3.):int(nwaveS*4.)]

    return StokesIRF, StokesQRF, StokesURF, StokesVRF

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------
# ------------------------------ Reading Response Functions when B = 0 ----------------------------
# -------------------------------------------------------------------------------------------------
# This function will read the Response Functions computed with SIR/DeSIRe. The output will be only 
# one 2D array containing the Stokes I RFs that depend on height and wavelength.
# 

# Running example:
# >>> from readwrite import readrfi
# >>> SiRF = readrfi(name)

def readrfi(name):

    # The information about n_tau and n_wave are in the first row

    header = np.loadtxt(name,skiprows=0,max_rows=1)
    ntau = header[0]
    nwave = header[1] 

    RF = np.loadtxt(name,skiprows=1)
    RF2 = RF.reshape(int(ntau),int(nwave),order='C')

    nwaveS = int(nwave) # We read all of them instead of every 4 points (see readrf)

    StokesIRF = RF2[:,int(nwaveS*0.):int(nwaveS*1.)]

    return StokesIRF

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
