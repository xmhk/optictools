import numpy as np


def gauss_peak_power( nurep, pmean, taufwhm):
    t0 = taufwhm / np.sqrt(2 * np.log(2))
    return pmean / ( nurep * t0 * np.sqrt( np.pi/2)) 


def optical_density_from_nu_to_lam( nuvec, Snu):
    c = 2.99792458e8
    lamvec = c / nuvec
    sortindx = sorted( range(len(nuvec)), key=lambda k:lamvec[k])
    Slam = Snu * c / lamvec**2
    return lamvec[sortindx], Slam[sortindx]

def optical_density_from_lam_to_nu( lamvec, Slam):
    c = 2.99792458e8
    nuvec = c / lamvec
    Snu = Slam * c / nuvec**2
    sortindx = sorted( range(len(nuvec)), key=lambda k:nuvec[k])
    return nuvec[sortindx], Snu[sortindx]
    
    



