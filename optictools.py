
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



def convert_dcurve_to_b2curve(dvals,omegas):
    """
    convert a GVD curve in the  wavelength representation (D) to a GVD curve (beta2) in the frequency representation

    INPUT:
    -dvals: (array of D values)
    -omegas: angular frequency vector
    
    OUTPUT:
    -beta2 :      GVD in frequency representation
    """    


    c = 2.99792458e8
    lams = 2 * np.pi * c / omegas
    return -lams/omegas*dvals



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

#def poly2beta(p,omega0):
#    """
#    convert the polynomial coefficients of a beta2(omega) fit to a beta series
#    
#    INPUT:
#    - p: polynomcoefficients from polyfit beta2(omega)
#    - omega0: center angular frequency for extension
#
#    OUTPUT: 
#    -betas: a list containing [beta0, beta1, beta2, beta3, ...]    
#    """
#    betas = [0,0]
#    reducedpoly=p
#    for i in range(len(p)):
#        bvalue = np.polyval(reducedpoly,omega0)
#        betas.append(bvalue)
#        reducedpoly = np.polyder(reducedpoly)
#    return betas

def poly2beta(p,omega0):
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

def beta2poly(betas):
    """
    convert the beta  coefficients to a  polynom
    
    INPUT:
    - betas:  beta coefficients

    OUTPUT: 
    -p: a polynom
    """
    p = []
    print("-----")
    for i in range(len(betas),0,-1):
        ii = i-1
        #print i, ii ,    betas[ii], factorial(ii), betas[ii]/max(1.0,factorial(ii))
        p.append(betas[ii]/max(1.0,factorial(ii)))
        
    return np.array( p )



def beta_change_base(betas, oldom0, newom0):
    """
    convert the beta series from one center frequency to the other
    
    INPUT:
    - betas:  beta coefficients
    - oldom0: old center frequency
    - newom0: new center frequency

    OUTPUT: 
    -betas : beta series for newom0
    """
    dom = newom0-oldom0
    ppoly = beta2poly(betas)
    newbetas = []
    for i in range(len(ppoly)+1 ):
        #print i, dom, np.polyval( ppoly, dom), ppoly        
        newbetas.append( np.polyval( ppoly, dom))
        ppoly = np.polyder(ppoly)
    return newbetas
                   
    



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


class pyfindpeaks2():   
    """objectified version of pyfindpeaks. use pfo = pyfindpeaks( environment, valuelist, thresh;
    peaks are in pfo.peakpos;
    
    plot with:
    plt.figure();pfo.show();plt.show()
    """
    
    
    def limitss(self, diff,length,pos):
    #this prevents hitting the borders of the array
        mi = np.max( [0, pos-diff])
        ma = np.min( [length, pos+diff])        
        return mi,ma
    
    
    def __init__(self, environment, valuelist, thresh):
        half = int( np.floor( environment/2))
        self.valuelist = valuelist
        self.thresh = thresh
        self.valuelistlength = len(valuelist)
        
        #pre-filter the peaks above threshold
        self.abovethresh = np.nonzero( self.valuelist>self.thresh )[0]
        self.message = "# above thresh : %d"%len(self.abovethresh)
        i = 0
        self.peakpos =np.array([],int)  
        self.milist = []
        self.malist = []
        while (i < len(self.abovethresh)):
            mi,ma = self.limitss(half, self.valuelistlength, self.abovethresh[i])  
            self.milist.append(mi)
            self.malist.append(ma)
            self.message += "\nabove thresh : %d     start = %d   stop=%d"%(i, mi, ma)
            
            partialvaluelist = valuelist[mi:ma]
            # is the valuelist value of the actual position the maximum of the environment?
            self.message+="\n--\nindex = %d,  value = %f    (max val = %f)"%(
                        self.abovethresh[i],
                        self.valuelist[self.abovethresh[i]],
                        max(partialvaluelist ))
            if self.valuelist[self.abovethresh[i]] == max(partialvaluelist):
                self.peakpos=np.append(self.peakpos,self.abovethresh[i])
                #skip forward to borders
                #print(self.abovethresh[i], ma)
                while ( self.abovethresh[i]<ma-2) and (i<len(self.abovethresh)-1):
                    i+=1
            #else :               
            i = i+1
                
            #else:
            #    i = i+1            
    def show(self):        
        plt.plot( self.valuelist)
        plt.axhline(y=self.thresh, dashes = (3,6), color="0.7")
        for i in range( len(self.milist)):
            plt.axvspan(self.milist[i],self.malist[i], alpha=0.1, color='blue')
        plt.plot( self.abovethresh, self.valuelist[self.abovethresh],'s', 
                 mfc='1.0', label="candidates", ms=8)
        plt.plot( self.peakpos, self.valuelist[self.peakpos],'o',
                 mfc='red',label='verified',
                ms=4)
        plt.legend(loc=0)    
    

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
        i = i+1        
    return peakpos

#def cfindpeaks(env, valuelist, threshval):
#    """
#    find peaks in a list or an array of value
#    
#    this is a python wrapper for the C-lib libfindpeaks (github/xhmk)
#
#   INPUT:
#    - environment: (INT) a maxima has to be the local maximum in this environment of points
#    - valuelist: list or array of points to find the maxima in
#    - thresh: a maximum has to be larger than this value
#
#    OUTPUT:
#    - listindices: positions of the peaks found
#    """
#    libfp = ctypes.cdll.LoadLibrary("libfindpeaks.so")
#    c_valuelist = (ctypes.c_double * len(valuelist))()
#    c_peakposes  = (ctypes.c_long * len(valuelist))()
#    c_listlength = ctypes.c_long()
#    c_peaknumber = ctypes.c_long()
#    c_env   = ctypes.c_long()
#    c_valuelist[:] = valuelist[:]
#    c_listlength = len(valuelist)
#    c_env = env
#    c_threshval = (ctypes.c_double  *1)()
#    c_threshval[0] = threshval
#    libfp.findpeaks(c_env, 
#                     c_listlength,
#                     ctypes.pointer(c_valuelist), 
#                     ctypes.pointer(c_peakposes),
#                     ctypes.pointer(c_peaknumber),
#                     c_threshval)
#    pindx = np.array(c_peakposes[0:c_peaknumber.value])
#    return pindx


def sechfield( p0, width, tvec,mode):
    """ returns the field of a temporal sech pulse
    
    INPUT:
    - p0: optical power in W
    - width: temporal width in s
    - tvec: time vector
    - mode: can be either 
            'fwhm'  (full width at half maximum of power (field squared))
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

def N2sol(t,z,b2,gamma, t0):
    """
    complete field for the N=2 soliton
    """
    p01 = np.abs(b2) / gamma / t0**2
    ld = t0**2 / np.abs(b2)  
    return np.sqrt(p01) * 4 * (np.cosh(3*t/t0) + 3 * np.exp(1.0j * 4 * z/ld ) * np.cosh(t/t0)) / (np.cosh(4*t/t0)+4*np.cosh(2*t/t0)+ 3*np.cos(4*z/ld)) * np.exp( 1.0j * z/ld /2.)  
    
def N2solnorm( t, z ): 
    """
    analytical field for an N=2-Soliton. z0 in this units is pi/2 
    """
    return 4 * (np.cosh(3*t) + 3 * np.exp(1.0j * 4 * z ) * np.cosh(t)) / (np.cosh(4*t)+4*np.cosh(2*t)+ 3*np.cos(4*z)) * np.exp( 1.0j * z /2.)    
    
def gaussfieldA( p0,width, tvec,mode):
    """ 
    returns the field of a gaussian pulse (A)
    type A: Intensity = 1/e**2 at T0
                                                                                                                   
    INPUT:
    - p0: optical power in W                                                                                        
    - width: temporal width in s                                                                                    
    - tvec: time vector                                                                                             
    - mode: can be either                                                                                           
            'fwhm'  (full width at half maximum of power (field squared))                                                       
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



def nuvec_from_tvec(tvec):
    """ 
    calculate the relative frequency vector 
    for a given time vector, e.g. for fft.

    INPUT:
    - tvec time vector (aequidistant)
    OUTPUT:
    - nuvec relative freq vector
    """
    N = len(tvec)
    dt = tvec[2]-tvec[1]
    dnu =  1/((N)*dt)
    # range(1,N) means [1,2,...,(N-1)] !!!
    nuvec = np.array( range(-N/2,N/2))*dnu
    return nuvec


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


def specnorm( nuv, snuv, pnorm = 1.0):
    """
    normalize spectra to have an integral power of pnorm,
    which is 1.0 as standard
    
    INPUT: 
        -nuv nu-vector (Hz)
        -snu S(nu) unit/Hz
        -OPTIONAL pnorm = 

    OUTPUT:
        -nunorm (interpolated to be aequidistant)
        -Snunorm
    
    """
    from scipy.interpolate import interp1d
    nu2 = np.linspace( np.min(nuv), np.max(nuv),len(nuv))
    sf = interp1d( nuv, snuv)
    snu2 = sf(nu2)
    dnu = nu2[2]-nu2[1]
    tscale = dnu *  np.trapz( np.abs(snu2))
    snu2 = snu2/tscale * pnorm
    return nu2, snu2




def calc_g12( fieldlist, verbose=False ):
    """
    calculate the absolute value of the complex spectral coherence function g12

    as defined in J. M. Dudley,  S. Coen: Opt. Lett 27, 1180 (2002).

    INPUT:
    - a list of fields [field1, field2, ...]
    - [optional] verbose = True/False 

    OUTPUT:
    - a dict containing the fields:
       - fsqmean : mean value of the intensity (|A|**2)
       - fsqstd  : standard deviation of the intensity (|A|**2)
       - g12 : absolute value of g12

    """
    rv = {}
    rv['fsqmean'] = np.mean( np.abs( np.array(fieldlist))**2, axis=0)
    rv['fsqstd'] = np.std( np.abs( np.array(fieldlist))**2, axis=0)
    numerator = np.zeros(len(fieldlist[0]))+0.0j
    denominator1 = np.zeros(len(fieldlist[0]))+0.0j
    denominator2 = np.zeros(len(fieldlist[0]))+0.0j
    pairs = 0
    for i in range(len(fieldlist)):
        for k in range(i+1,len(fieldlist)):
            pairs +=1
            if verbose==True:
                print("i = %d   k = %d"%(i, k))
            numerator +=  np.multiply( np.conj(fieldlist[i]), fieldlist[k])
            denominator1 +=  np.abs(fieldlist[i])**2
            denominator2 +=  np.abs(fieldlist[k])**2
    numerator = numerator / pairs
    denominator1 = denominator1 / pairs
    denominator2 = denominator2 / pairs
    rv['g12'] = np.abs( numerator / np.sqrt( np.multiply( denominator1, denominator2)))
    return rv


def frog_corr(field1, field2):
    """
    calculate the spectrogramm of two fields.
    
    calculated results spectrogramm that would be yielded by
    cross-correlating the fields and than measure the fourier
    transform.

    INPUT: field1, field2
    (SHG-FROG-trace when field1 == field2)

    OUTPUT: MF
    with MF = | Fourier Transform ( int dt field(t) * field2(t-tau) | **2
    (first axis-> tau, second axis->nu)
    """
    def FROG_row_rotate(M):
        for r in range(np.shape(M)[0]):
            M[r,:]=np.roll( M[r,:],-r)
        return M    
    TDfrog = FROG_row_rotate( np.outer(field1, field2))
    FDfrog =  np.fft.fftshift( np.fft.fft( TDfrog,axis=0 ))
    FI = np.abs(FDfrog)**2
    return FI / np.max(FI)



def sg( xvec, mu, width, order):
    """
    brief version of a supergauss - function
    INPUT:
    -xvec time( or x-)-vector
    -mu time offset
    -width
    -order
    """
    return np.exp( - ((xvec-mu)/width)**order)

def supergauss( xvec, mu, width, order):
    """
    supergauss - function
    INPUT:
    -xvec time( or x-)-vector
    -mu time offset
    -width
    -order
    """
    return sg( xvec, mu, width, order) 

def specfilter_on_timefield( field, filt ):
    """
    apply a spectral filter on a timefield,
    return the filtered timefield
    (IFFT type)
    """
    return np.fft.fft( np.fft.ifft(field )* np.fft.fftshift(filt))


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
# little helpers
#

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

def get_slice(signal, xvec, lower,upper):
    """
    cuts a slice from an array
    
    INPUT:
    - signale
    - xvec
    - lower x
    - upper x

    OUTPUT:
    - slice of input array
    """
    llind = ge_index(xvec, lower)
    uuind = ge_index(xvec, upper)
    return signal[llind:uuind]

def moving_average(somearray, environment):
    """
    moving average
    
    dumb implementation (slow!)

    INPUT:
    somearray - some array
    environment - number of points to include in the average

    OUTPUT:
    movav : array with moving average
    movavstd : standard deviation

    length of movav, movavstd is the same as the input array. 
    the borders are filled up with the first/last calculated value
    """
    envhalf = int(np.floor(environment/2))    
    movav =  np.zeros(np.shape(somearray))
    movavstd = np.zeros(np.shape(somearray))
    #move throught the array and calculate mean, std
    for i in range(envhalf+1,len(somearray)-envhalf):
        movav[i]     = np.mean( somearray[i-envhalf:i+envhalf] )
        movavstd[i] = np.std( somearray[i-envhalf:i+envhalf])
    #fill up the borders
    movavstd[len(somearray)-envhalf-1:len(somearray)]=movavstd[len(somearray)-envhalf-3]
    movavstd[0:envhalf+1]=movavstd[envhalf+2]

    movav[len(somearray)-envhalf-1:len(somearray)]=movav[len(somearray)-envhalf-3]
    movav[0:envhalf+1]=movav[envhalf+2]

    return movav,movavstd
