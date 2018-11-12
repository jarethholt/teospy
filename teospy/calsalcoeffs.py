"""Calculate coefficients for the salt Gibbs function.

Calculate modified coefficients for the Gibbs energy function of salt in
seawater. These coefficients need to be modified from those given in the
reference documents in a few ways. First, the implementation here uses
the logarithm of salinity rather than the logarithm of the square root
of salinity, so the corresponding coefficients have to be divided by
two. Second, the implementation here allows the primary routines to use
salinity directly by scaling the coefficients here by the standard
salinity value. Third, there are two coefficients that can be freely
chosen, and are only fixed by a choice of reference state. The reference
state chosen by IAPWS is that the enthalpy and entropy of seawater under
standard conditions is zero, where the liquid water component is
calculated using the Feistel (2003) polynomial formulation.

This script is intended to be called as::

    python calsalcoeffs.py [shift]

which prints the modified coefficients in a format that can be directly
copied and pasted into the module :mod:`sal1`. This is not recommended
unless the module behavior does not seem to match the triple point and
standard parameter conditions of seawater. If the optional keyword
``shift`` is given, the coefficients ``g200`` and ``g210`` are adjusted
so that the enthalpy and entropy of seawater under standard conditions
is zero when the Feistel (2003) Gibbs energy of liquid water is used.
This currently produces more inconsistencies with the Fortran library
than it resolves.
"""

import copy
import math
import constants0
SAL1 = constants0.SAL1
TCELS = constants0.TCELS
PATM = constants0.PATM
SU = SAL1 * 40./35.
LSU = math.log(SU)
RSU = SU**.5

# Standard coefficients
GSC_ORIG = {
    (1,0,0):  0.581281456626732e+4, (1,1,0):  0.851226734946706e+3,
    (2,0,0):  0.141627648484197e+4, (2,1,0):  0.168072408311545e+3,
    (3,0,0): -0.243214662381794e+4, (4,0,0):  0.202580115603697e+4,
    (5,0,0): -0.109166841042967e+4, (6,0,0):  0.37460123787784e+3 ,
    (7,0,0): -0.485891069025409e+2, (3,1,0): -0.493407510141682e+3,
    (4,1,0):  0.543835333000098e+3, (5,1,0): -0.196028306689776e+3,
    (6,1,0):  0.367571622995805e+2, (2,2,0):  0.880031352997204e+3,
    (3,2,0): -0.430664675978042e+2, (4,2,0): -0.685572509204491e+2,
    (2,3,0): -0.225267649263401e+3, (3,3,0): -0.100227370861875e+2,
    (4,3,0):  0.493667694856254e+2, (2,4,0):  0.914260447751259e+2,
    (3,4,0):  0.875600661808945e0 , (4,4,0): -0.171397577419788e+2,
    (2,5,0): -0.216603240875311e+2, (4,5,0):  0.249697009569508e+1,
    (2,6,0):  0.213016970847183e+1, (2,0,1): -0.331049154044839e+4,
    (3,0,1):  0.199459603073901e+3, (4,0,1): -0.547919133532887e+2,
    (5,0,1):  0.360284195611086e+2, (2,1,1):  0.729116529735046e+3,
    (3,1,1): -0.175292041186547e+3, (4,1,1): -0.226683558512829e+2,
    (2,2,1): -0.860764303783977e+3, (3,2,1):  0.383058066002476e+3,
    (2,3,1):  0.694244814133268e+3, (3,3,1): -0.460319931801257e+3,
    (2,4,1): -0.297728741987187e+3, (3,4,1):  0.234565187611355e+3,
    (2,0,2):  0.384794152978599e+3, (3,0,2): -0.522940909281335e+2,
    (4,0,2): -0.408193978912261e+1, (2,1,2): -0.343956902961561e+3,
    (3,1,2):  0.831923927801819e+2, (2,2,2):  0.337409530269367e+3,
    (3,2,2): -0.541917262517112e+2, (2,3,2): -0.204889641964903e+3,
    (2,4,2):  0.74726141138756e+2 , (2,0,3): -0.965324320107458e+2,
    (3,0,3):  0.680444942726459e+2, (4,0,3): -0.301755111971161e+2,
    (2,1,3):  0.124687671116248e+3, (3,1,3): -0.29483064349429e+2 ,
    (2,2,3): -0.178314556207638e+3, (3,2,3):  0.256398487389914e+2,
    (2,3,3):  0.113561697840594e+3, (2,4,3): -0.364872919001588e+2,
    (2,0,4):  0.158408172766824e+2, (3,0,4): -0.341251932441282e+1,
    (2,1,4): -0.31656964386073e+2 , (2,2,4):  0.442040358308e+2   ,
    (2,3,4): -0.111282734326413e+2, (2,0,5): -0.262480156590992e+1,
    (2,1,5):  0.704658803315449e+1, (2,2,5): -0.792001547211682e+1
}
GSCEXT_ORIG = {
    (2,1,1): -0.347600838235511e2, (3,1,1):  0.935178208024272e2,
    (4,1,1): -0.603483495593212e2, (2,2,1):  0.228081199116236e2,
    (4,2,1): -0.242869748847311e2, (2,3,1): -0.560725782144008e2,
    (3,3,1): -0.143523729429211e2, (4,3,1):  0.926367388049097e2,
    (4,4,1): -0.416658900599273e2, (2,5,1):  0.645288813326254e2,
    (3,5,1): -0.403505133068118e2, (2,6,1): -4.32746069361075e0 ,
    (4,6,1):  2.05954716712622e0 
}

GSCTOT_ORIG = copy.copy(GSCEXT_ORIG)
for ((i,j,k),gc) in GSC_ORIG.items():
    ge = GSCEXT_ORIG.get((i,j,k),0.)
    GSCTOT_ORIG[(i,j,k)] = gc + ge


## Calculate the modified coefficients
def _prntformatted(gscdict):
    """Print a dictionary of coefficients.
    
    Print a dictionary of salinity coefficients in a format that can be
    directly copied into :mod:`sal1`.
    """
    COEFFMT = '({0:1d},{1:1d},{2:1d},{3: 21.14e})'
    line = '\t'
    ncoeff = len(gscdict)
    inds = list(gscdict.keys())
    inds.sort()
    for (icoeff,ind) in enumerate(inds):
        i, j, k = ind
        gi = gscdict[ind]
        line += COEFFMT.format(i,j,k,gi)
        if icoeff == ncoeff-1:
            line += '\n)'
        elif icoeff % 2 == 0:
            line += ', '
        else:
            line += ',\n\t'
    return line

def _calcoeffs():
    """Calculate the modified salinity coefficients.
    
    Calculate the modified salinity coefficients from the original ones.
    The modified coefficients are scaled such that the Gibbs free energy
    of salt in seawater is calculated directly from the absolute
    salinity without scaling or calculating the square root.
    
    :returns: Dictionaries of the modified coefficients with and without
        the high-temperature, high-salinity extension.
    """
    gsc_mod = dict()
    for ((i,j,k),go) in GSC_ORIG.items():
        if i==1:
            gm = go * .5/SU
        elif i==2:
            g1jk = GSC_ORIG.get((1,j,k),0.)
            gm = (go - .5*g1jk*LSU)/SU
        else:
            gm = go / RSU**i
        gsc_mod[(i,j,k)] = gm
    
    gsctot_mod = dict()
    for ((i,j,k),go) in GSCTOT_ORIG.items():
        if i==1:
            gm = go * .5/SU
        elif i==2:
            g1jk = GSCTOT_ORIG.get((1,j,k),0.)
            gm = (go - .5*g1jk*LSU)/SU
        else:
            gm = go / RSU**i
        gsctot_mod[(i,j,k)] = gm
    return (gsc_mod, gsctot_mod)

def _calcoeffs_shifted():
    """Calculate the modified salinity coefficients.
    
    Calculate the modified salinity coefficients from the original ones.
    The modified coefficients are scaled such that the Gibbs free energy
    of salt in seawater is calculated directly from the absolute
    salinity without scaling or calculating the square root.
    
    In addition, the coefficients g200 and g210 are adjusted so that the
    enthalpy and entropy of seawater at the standard ocean conditions
    are exactly zero when the Feistel (2003) formulation for the Gibbs
    energy of seawater is used. This choice currently produces more
    inconsistencies with the Fortran library than it resolves.
    
    :returns: Dictionaries of the modified coefficients with and without
        the high-temperature, high-salinity extension.
    """
    import sal1
    import liq5_f03
    TTP = constants0.TTP
    PTPI = constants0.PTPI
    UWTP = liq5_f03.internalenergy(TTP,PTPI)
    HW0 = liq5_f03.enthalpy(TCELS,PATM)
    SWTP = liq5_f03.entropy(TTP,PTPI)
    SW0 = liq5_f03.entropy(TCELS,PATM)
    FWTP = UWTP - TTP*SWTP
    GW0 = HW0 - TCELS*SW0
    TRED = sal1._TRED
    LS1 = math.log(SAL1)
    RS1 = SAL1**.5
    
    gsc_mod = dict()
    for ((i,j,k),go) in GSC_ORIG.items():
        if i==1:
            gm = go * .5/SU
        elif i==2:
            g1jk = GSC_ORIG.get((1,j,k),0.)
            gm = (go - .5*g1jk*LSU)/SU
        else:
            gm = go / RSU**i
        gsc_mod[(i,j,k)] = gm
    
    g200 = (FWTP - GW0)/SAL1 - gsc_mod[(1,0,0)]*LS1
    g210 = TRED/SAL1*(SW0 - SWTP) - gsc_mod[(1,1,0)]*LS1
    for ((i,j,k),gm) in gsc_mod.items():
        if i>2 and k==0:
            if j==0:
                g200 -= gm * RS1**(i-2)
            if j==1:
                g210 -= gm * RS1**(i-2)
    gsc_mod[(2,0,0)] = g200
    gsc_mod[(2,1,0)] = g210
    
    gsctot_mod = dict()
    for ((i,j,k),go) in GSCTOT_ORIG.items():
        if i==1:
            gm = go * .5/SU
        elif i==2:
            g1jk = GSCTOT_ORIG.get((1,j,k),0.)
            gm = (go - .5*g1jk*LSU)/SU
        else:
            gm = go / RSU**i
        gsctot_mod[(i,j,k)] = gm
    
    gt200 = (FWTP - GW0)/SAL1 - gsctot_mod[(1,0,0)]*LS1
    gt210 = TRED/SAL1*(SW0 - SWTP) - gsctot_mod[(1,1,0)]*LS1
    for ((i,j,k),gm) in gsctot_mod.items():
        if i>2 and k==0:
            if j==0:
                gt200 -= gm * RS1**(i-2)
            if j==1:
                gt210 -= gm * RS1**(i-2)
    gsctot_mod[(2,0,0)] = gt200
    gsctot_mod[(2,1,0)] = gt210
    return (gsc_mod, gsctot_mod)


## Command line interface
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == 'shift':
        gsc_mod, gsctot_mod = _calcoeffs_shifted()
    else:
        gsc_mod, gsctot_mod = _calcoeffs()
    
    print('_GSCOEFFS = (')
    line = _prntformatted(gsc_mod)
    print(line)
    print('_GSCOEFFS_EXT = (')
    line = _prntformatted(gsctot_mod)
    print(line)

