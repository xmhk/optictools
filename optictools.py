
import numpy as np
from matplotlib import pyplot as plt
from scipy.misc import factorial
from scipy.interpolate import interp1d
import ctypes

def beta2_curve(omega , omega0, betas):
    """
    calculate the group velocity curve from a beta series

    INPUT:
    - omegas: a vector of angular frequencies
    - omega0: the center angular frequency
    - betas: a list with the beta coefficients: [beta0, beta1, beta2, ...]

    OUTPUT: 
    - the beta2(omega) curve
    """
    domega  = np.array(omega)-omega0
    b2k = np.zeros(len(omega))
    for i in range(2,len(betas)):
        b2k = b2k + betas[i]/factorial(i-2) * domega**(i-2)
    return b2k

def convert_b2curve_to_dcurve(b2vals, omegas):
    """
    convert a GVD curve in the frequency representation (beta2) to a GVD curve (D) in the wavelength representation

    INPUT:
    -b2vals: (array of beta2 values)
    -omegas: angular frequency vector
    
    OUTPUT:
    -D:      GVD in wavelength representation
    """    
    c = 2.99792458e8
    lams = 2 * np.pi * c / omegas
    d = -1.0 * omegas/lams * b2vals
    return d

def beta0_curve(omegas, omega0, betas):
    """
    calculate the dispersion curve from a beta series

    INPUT:
    - omegas: a vector of angular frequencies
    - omega0: the center angular frequency
    - betas: a list with the beta coefficients: [beta0, beta1, beta2, ...]

    OUTPUT: 
    - the beta(omega) curve
    """
    bc = np.zeros(len(omegas))
    for i in range(len(betas)):
        bc = bc + betas[i]/factorial(i) * (omegas-omega0)**i
    return bc

def beta0_curve_from_b2data( omegas, b2data, omega0):
    """
    get the beta curve from (numerical beta data

    this is done via integration
    
    INPUT:
    - omegas: aequidistant(!) vector of angular frequencies
    - b2data: beta2 data for omegas
    - omega0: center frequency (result will be zero here)
    
    OUTPUT:
    - betacurve(omegas)
    
    """
    dom = omegas[2]-omegas[1]
    b = np.cumsum(np.cumsum(b2data) ) * dom**2
    bk = interp1d( omegas, b)
    bb = b - bk(omega0)
    return bb

def get_even_part( omvec, om0, k_curve):
    """
    extract the even part of a (dispersion) curve

    INPUT:
    -omvec - angular frequency vector
    -om0   - center frequency (reference)
    -k_curve - dispersion curve

    OUTPUT:
    -even part of the dispersion curve
    """
    f1 = interp1d( omvec, k_curve,'linear')
    k_even = np.zeros( np.shape(omvec))
    for i in range(len(omvec)):
        deltao = om0-omvec[i]
        if om0-deltao<omvec[0]:
            k_even[i] = (k_curve[0]+k_curve[-1])/2.0
        elif om0 + deltao > omvec[-1]:
            k_even[i] = (k_curve[-1]+k_curve[0])/2.0
        else:
            k_even[i] =  ( f1( om0-deltao) + f1(om0+deltao) )/2.0
    return k_even

def poly2beta(p,omega0):
    """
    convert the polynomial coefficients of a beta2(omega) fit to a beta series
    
    INPUT:
    - p: polynomcoefficients from polyfit beta2(omega)
    - omega0: center angular frequency for extension

    OUTPUT: 
    -betas: a list containing [beta0, beta1, beta2, beta3, ...]    
    """
    betas = [0,0]
    reducedpoly=p
    for i in range(len(p)):
        bvalue = np.polyval(reducedpoly,omega0)
        betas.append(bvalue)
        reducedpoly = np.polyder(reducedpoly)
    return betas

def poly2beta_B(p,omega0):
    """
    convert the polynomial coefficients of a beta(omega) fit to a beta series
    
    INPUT:
    - p: polynomcoefficients from polyfit of beta(omega)
    - omega0: center angular frequency for extension

    OUTPUT: 
    -betas: a list containing [beta0, beta1, beta2, beta3, ...]
    """
    betas = []
    reducedpoly=p
    for i in range(len(p)):
        bvalue = np.polyval(reducedpoly,omega0)
        betas.append(bvalue)
        reducedpoly = np.polyder(reducedpoly)
    return betas


def fwhm3(valuelist, peakpos=-1):
    """calculates the full width at half maximum (fwhm) of some curve.

    the function will return the fwhm with sub-pixel interpolation. It will start at the maximum position and 'walk' left and right until it approaches the half values.

    INPUT: 
    - valuelist: e.g. the list containing the temporal shape of a pulse 

    OPTIONAL INPUT: 
    -peakpos: position of the peak to examine (list index)
    the global maximum will be used if omitted.

    OUTPUT:
    -fwhm (value)
    """
    if peakpos== -1: #no peakpos given -> take maximum
        peak = np.max(valuelist)
        peakpos = np.min( np.nonzero( valuelist==peak  )  )

    peakvalue = valuelist[peakpos]
    phalf = peakvalue / 2.0

    # go left and right, starting from peakpos
    ind1 = peakpos
    ind2 = peakpos   

    while ind1>2 and valuelist[ind1]>phalf:
        ind1=ind1-1
    while ind2<len(valuelist)-1 and valuelist[ind2]>phalf:
        ind2=ind2+1  
    #ind1 and 2 are now just below phalf
    grad1 = valuelist[ind1+1]-valuelist[ind1]
    grad2 = valuelist[ind2]-valuelist[ind2-1]
    #calculate the linear interpolations
    p1interp= ind1 + (phalf -valuelist[ind1])/grad1
    p2interp= ind2 + (phalf -valuelist[ind2])/grad2
    #calculate the width
    width = p2interp-p1interp
    return width



def pyfindpeaks( environment, valuelist , thresh):
    """
    find peak positions in a list of values

    INPUT:
    - environment: (INT) a maxima has to be the local maximum in this environment of points
    - valuelist: list or array of points to find the maxima in
    - thresh: a maximum has to be larger than this value

    OUTPUT:
    - listindices: positions of the peaks found
    """
    def limitss(diff,length,pos):
    #this prevents hitting the borders of the array
        mi = np.max( [0, pos-diff])
        ma = np.min( [length, pos+diff])
        return mi,ma
    #range left/right
    half = int( np.floor( environment/2))
    valuelistlength = len(valuelist)
    #pre-filter the peaks above threshold
    abovethresh = np.nonzero( valuelist>thresh )[0]
    i = 0
    peakpos =np.array([],int)  
    # circle through the candidates
    while (i < len(abovethresh)):
        mi,ma = limitss(half, valuelistlength, abovethresh[i])
        partialvaluelist = valuelist[mi:ma]
    # is the valuelist value of the actual position the maximum of the environment?
        if valuelist[abovethresh[i]] == max(partialvaluelist):
            peakpos=np.append(peakpos,abovethresh[i])
            i = i+half-1 #skip the non-maxima
        else:
            i = i+1
    return peakpos

def cfindpeaks(env, valuelist, threshval):
    """
    find peaks in a list or an array of value
    
    this is a python wrapper for the C-lib libfindpeaks (github/xhmk)

    INPUT:
    - environment: (INT) a maxima has to be the local maximum in this environment of points
    - valuelist: list or array of points to find the maxima in
    - thresh: a maximum has to be larger than this value

    OUTPUT:
    - listindices: positions of the peaks found
    """
    libfp = ctypes.cdll.LoadLibrary("libfindpeaks.so")
    c_valuelist = (ctypes.c_double * len(valuelist))()
    c_peakposes  = (ctypes.c_long * len(valuelist))()
    c_listlength = ctypes.c_long()
    c_peaknumber = ctypes.c_long()
    c_env   = ctypes.c_long()
    c_valuelist[:] = valuelist[:]
    c_listlength = len(valuelist)
    c_env = env
    c_threshval = (ctypes.c_double  *1)()
    c_threshval[0] = threshval
    libfp.findpeaks(c_env, 
                     c_listlength,
                     ctypes.pointer(c_valuelist), 
                     ctypes.pointer(c_peakposes),
                     ctypes.pointer(c_peaknumber),
                     c_threshval)
    pindx = np.array(c_peakposes[0:c_peaknumber.value])
    return pindx


def sechfield( p0, width, tvec,mode):
    """ returns the field of a temporal sech pulse
    
    INPUT:
    - p0: optical power in W
    - width: temporal width in s
    - tvec: time vector
    - mode: can be either 
            'fwhm'  (full width at half maximum of intensity)
            or 't0' (argument of sech)  

    OUTPUT:
    - temporal sech field
    """
    if mode.lower() == 'fwhm':
        t0 = width/ 2 / np.arcsinh(1)
    elif mode.lower() == 't0':
        t0 = width
    else:
        print("sechfield error! no valid width given!!")
        t0 = 0
    return np.sqrt(p0) * 1/np.cosh( tvec / t0)

def gaussfieldA( p0,width, tvec,mode):
    """ 
    returns the field of a gaussian pulse (A)
    type A: Intensity = 1/e**2 at T0
                                                                                                                   
    INPUT:
    - p0: optical power in W                                                                                        
    - width: temporal width in s                                                                                    
    - tvec: time vector                                                                                             
    - mode: can be either                                                                                           
            'fwhm'  (full width at half maximum of intensity)                                                       
            or 't0' (argument of exp)          
    OUTPUT:
    - temporal gaussian field (A)                                                                     
    """
    if mode.lower() == 'fwhm':
        t0 = width / 2 / np.sqrt(np.log(2))
    elif mode.lower() == 't0':
        t0 = width
    else:
        print( "gaussfieldA error! no valid width given!!")
        t0 = 0
    return np.sqrt(p0) * np.exp( -0.5* (tvec/t0)**2)

def gaussfieldB( p0,width, tvec,mode):
    """ 
    returns the field of a gaussian pulse (B)
    type B: Intensity = 1/e at T0
                                                                                                                   
    INPUT:
    - p0: optical power in W                                                                                        
    - width: temporal width in s                                                                                    
    - tvec: time vector                                                                                             
    - mode: can be either                                                                                           
            'fwhm'  (full width at half maximum of intensity)                                                       
            or 't0' (argument of exp)                                                                             
    OUTPUT:
    - temporal gaussian field (B)  
    """
    if mode.lower() == 'fwhm':
        t0 = width / 2 / np.sqrt(np.log(2)/2)
    elif mode.lower() == 't0':
        t0 = width
    else:
        print( "gaussfieldB error! no valid width given!!")
        t0 = 0
    return np.sqrt(p0) * np.exp( - (tvec/t0)**2)

def gauss_peak_power( nurep, pmean, taufwhm):
    """ 
    calculate the peak power of gaussian pulses (A) from 
    the repetition frequency, the mean power and the fwhm

    INPUT:
    - nurep: repetition frequency
    - pmean: mean power
    - taufwhm: temporal fwhm (intensity) of the pulses
    
    OUTPUT:
    - peak power of the pulses
    """
    t0 = taufwhm / np.sqrt(2 * np.log(2))
    return pmean / ( nurep * t0 * np.sqrt( np.pi/2)) 

def sech_peak_power( nurep, pmean, taufwhm):
    """ 
    calculate the peak power of sech from 
    the repetition frequency, the mean power and the fwhm

    INPUT:
    - nurep: repetition frequency
    - pmean: mean power
    - taufwhm: temporal fwhm (intensity) of the pulses
    
    OUTPUT:
    - peak power of the pulses
    """
    t0 = taufwhm / 2 / np.arcsinh(1)
    return pmean / ( nurep * t0 * np.sqrt( np.pi/2)) 


def optical_density_from_nu_to_lam( nuvec, Snu):
    """
    convert a spectral density from frequency to wavelength representation

    INPUT:
    - nuvec: vector of optical frequency
    - Snu: spectral density Snu(nuvec)

    OUTPUT:
    - lamvec: vector of wavelengths
    - Slam: spectral density Slam(lamvec)
    """
    c = 2.99792458e8
    lamvec = c / nuvec
    sortindx = sorted( range(len(nuvec)), key=lambda k:lamvec[k])
    Slam = Snu * c / lamvec**2
    return lamvec[sortindx], Slam[sortindx]

def optical_density_from_lam_to_nu( lamvec, Slam):
    """
    convert a spectral density from wavlength to frequency representation

    INPUT:
    - lamvec: vector of wavelengths
    - Slam: spectral density Slam(lamvec)

    OUTPUT:
    - nuvec: vector of optical frequency
    - Snu: spectral density Snu(nuvec)
    """
    c = 2.99792458e8
    nuvec = c / lamvec
    Snu = Slam * c / nuvec**2
    sortindx = sorted( range(len(nuvec)), key=lambda k:nuvec[k])
    return nuvec[sortindx], Snu[sortindx]

def passnotch(vec,n1,n2,mode="pass"):
    """
    a very simple bandpass or notch binary filter
    
    INPUT:
    - vec: a monotonously growing vector e.g of frequency, time
    - n1,n2: border frequency of notch/pass
    - mode: either
          'notch': create a filter with zeros between n1 and n2, ones otherwise
       or 'pass': create a filter with ones between n1 and n2, zeros otherwise
     
    OUTPUT:
    - a vector with the same length as vec. Filled with zeros and ones    
    """
    if mode == "pass":
        return np.multiply( vec > np.min([n1,n2]), vec < np.max([n1,n2]))
    elif mode == "notch":
        return 1-np.multiply( vec > np.min([n1,n2]), vec < np.max([n1,n2]))
    else:
        print( "error: mode should be 'notch' or 'pass'")
        return None


#
# tools for dispersion measurement
#
def heneint(henespur):
    """
    helper for fourier transform white light interferometry
    
    INPUT: 
    - henespur - voltage from reference interferometer
    OUTPUT:
    - xn - integer positions of zeros
    - xnint - (float) interpolated positions of zeros
    """
    lh = len(henespur)
    henespur = henespur - np.mean(henespur)
    # wo unterscheiden sich die Vorzeichen von zwei aufeinanderfolgenden Punkten?
    diffprod = np.multiply( henespur[0:lh-1],henespur[1:lh])    
    nsfilter = diffprod<0
    #interpoliere die Nullstelle
    xn = np.array(np.nonzero(nsfilter))
    yn = np.array([henespur[n] for n in xn])
    yn1 = np.array([henespur[n+1] for n in xn])
    xnint = np.array(-1.0*yn/(yn1-yn)+xn) #interpolierte Nullstellen
    return xn[0],xnint[0]


def reduziertes_interferogramm( xn, xnint, interferogrammspur ):
    """
    reduce a interferogramm with interpolated reference points
    
    INPUT:
    -xn, xnint - interpolated reference points from henespur
    -interferogrammspur: recorded white light interferogramm
    OUTPUT:
    -interpolated interferogramm
    """
    yn = interferogrammspur[xn]
    yn1 = interferogrammspur[xn+1]
    yint = np.multiply(  (yn1-yn), xnint-xn)+yn   
    return np.array(yint)

def unwrap2(phase):
    """
    unwrap a phase

    compared to standard unwrap, this also ensures that the first derivative
    of the phase is steady

    INPUT:
    -phase curve

    OUTPUT:
    -unwrapped phase cure
    """
    p1 = np.unwrap(phase)
    deltal=[]
    for i in range(1,len(phase)-1):
        delta = p1[i+1]-2*p1[i]+p1[i-1]
        kfak = delta/( np.pi)             
        if np.abs(kfak)>1.0:            
            #kk = -1* np.sign(kfak)*np.floor(np.abs(kfak)+1)
            kk=-1*np.sign(kfak)*np.ceil(np.abs(kfak))
            p1[i+1::]=p1[i+1::]+(kk)*np.pi
            #print t[i], delta,kfak,kk,delta+kk*np.pi
        else:
            kk = 0.0      
    return p1

def ge_index(valuelist, val):
    """ greater-equal index 

    returns the index of the first value in valuelist that is greater or equal val
    INPUT:
    - valuelist 
    - val
    
    OUTPUT:
    - index
    """
    arra = np.array(valuelist)        
    return np.min( np.nonzero( arra>=val))

def le_index(valuelist, val):
    """ lower-equal index

    returns the index of the first value in valuelist that is smaller or equal val
    INPUT:
    - valuelist 
    - val
    
    OUTPUT:
    - index
    """
    arra = np.array(valuelist)        
    return np.max( np.nonzero( arra<=val))

def db_abs2(y):
    """
    return the decadic logarithm of the absolut square of a value (Decibel)
    
    INPUT:
    -y: value
    
    OUTPUT:
    -dby = 10 * np.log10( np.abs(y)**2)
    """
    return 10 * np.log10( np.abs(y)**2)

def db_abs(y):
    """
    return the decadic logarithm of the abs of a value (Decibel)
    
    INPUT:
    -y: value
    
    OUTPUT:
    -dby = 10 * np.log10( np.abs(y))
    """
    return 10 * np.log10( np.abs(y))
