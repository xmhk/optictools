
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

def beta0_curve(omvec, om0, betas):
    bc = np.zeros(len(omvec))
    for i in range(len(betas)):
        bc = bc + betas[i]/factorial(i) * (omvec-om0)**i
    return bc

def beta0_curve_from_b2data( omvec, b2data, om0):
    dom = omvec[2]-omvec[1]
    b = np.cumsum(np.cumsum(b2data) ) * dom**2
    bk = interp1d( omvec, b)
    bb = b - bk(om0)
    return bb

def remove_b1_slope_from_betacurve( betacurve, omvec, om0, deltaom, fignr=0):
    pointstoconsider = np.multiply( omvec > om0-deltaom, omvec<om0+deltaom)
    p = np.polyfit( omvec[pointstoconsider], betacurve[pointstoconsider], 1)
    betareturn = betacurve - np.polyval(p, omvec)
    if fignr>0:
        plt.figure(fignr)
        plt.plot( omvec, betacurve)
        plt.plot( omvec, betareturn)
        plt.xlim( [om0-1.1*deltaom, om0+1.1*deltaom])
        plt.ylim( [ min( betacurve[pointstoconsider]),max(betacurve[pointstoconsider])])
        plt.axvline(om0,c="#777777")
        plt.legend(["before","after","om0"])
    return betareturn


def get_even_part( omvec, om0, k_curve):
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

def poly2beta(p,xo):
    betas = [0,0]
    reducedpoly=p
    for i in range(len(p)):
        bvalue = np.polyval(reducedpoly,xo)
        betas.append(bvalue)
        reducedpoly = np.polyder(reducedpoly)
    return betas

def poly2beta_B(p,xo):
#    betas = [0,0]
    betas = []
    reducedpoly=p
    for i in range(len(p)):
        bvalue = np.polyval(reducedpoly,xo)
        betas.append(bvalue)
        reducedpoly = np.polyder(reducedpoly)
    return betas


#full width at half maximum with linear interpolation
def fwhm3(list, peakpos=-1):
  if peakpos== -1: #no peakpos given -> take maximum
    peak = np.max(list)
    peakpos = np.min( np.nonzero( list==peak  )  )

  peakvalue = list[peakpos]
  phalf = peakvalue / 2.0

  # go left and right, starting from peakpos
  ind1 = peakpos
  ind2 = peakpos   

  while ind1>2 and list[ind1]>phalf:
    ind1=ind1-1
  while ind2<len(list)-1 and list[ind2]>phalf:
    ind2=ind2+1
  
  #ind1 and 2 are now just below phalf
  grad1 = list[ind1+1]-list[ind1]
  grad2 = list[ind2]-list[ind2-1]
  #calculate the linear interpolations
  p1interp= ind1 + (phalf -list[ind1])/grad1
  p2interp= ind2 + (phalf -list[ind2])/grad2
  #calculate the width
  width = p2interp-p1interp
  return width


#-------------------------------------------------------------------------------------------
# find maxima in a list
def pyfindpeaks( environment, list , thresh):
  def limitss(diff,length,pos):
    #this prevents hitting the borders of the array
    mi = np.max( [0, pos-diff])
    ma = np.min( [length, pos+diff])
    return mi,ma
  #range left/right
  half = int( np.floor( environment/2))
  listlength = len(list)
  #pre-filter the peaks above threshold
  abovethresh = np.nonzero( list>thresh )[0]
  i = 0
  peakpos =np.array([],int)  
  # circle through the candidates
  while (i < len(abovethresh)):
    mi,ma = limitss(half, listlength, abovethresh[i])
    partiallist = list[mi:ma]
    # is the list value of the actual position the maximum of the environment?
    if list[abovethresh[i]] == max(partiallist):
      peakpos=np.append(peakpos,abovethresh[i])
      i = i+half-1 #skip the non-maxima
    else:
      i = i+1
  return peakpos





#-------------------------------------------------------------------------------------------
# find maxima in a list (c version)

def cfindpeaks(env, liste, threshval):
    libfp = ctypes.cdll.LoadLibrary("libfindpeaks.so")
    c_liste = (ctypes.c_double * len(liste))()
    c_peakposes  = (ctypes.c_long * len(liste))()
    c_listlength = ctypes.c_long()
    c_peaknumber = ctypes.c_long()
    c_env   = ctypes.c_long()
    c_liste[:] = liste[:]
    c_listlength = len(liste)
    c_env = env
    c_threshval = (ctypes.c_double  *1)()
    c_threshval[0] = threshval
    libfp.findpeaks(c_env, 
                     c_listlength,
                     ctypes.pointer(c_liste), 
                     ctypes.pointer(c_peakposes),
                     ctypes.pointer(c_peaknumber),
                     c_threshval)
    pindx = np.array(c_peakposes[0:c_peaknumber.value])
    return pindx




# -----------------------------------------------------------------------------


def sechfield( p0, width, tvec,mode):
    if mode.lower() == 'fwhm':
        t0 = width/ 2 / np.arcsinh(1)
    elif mode.lower() == 't0':
        t0 = width
    else:
        print("sechfield error! no valid width given!!")
        t0 = 0
    return np.sqrt(p0) * 1/np.cosh( tvec / t0)

def gaussfieldA( p0,width, tvec,mode):
    if mode.lower() == 'fwhm':
        t0 = width / 2 / np.sqrt(np.log(2))
    elif mode.lower() == 't0':
        t0 = width
    else:
        print( "gaussfieldA error! no valid width given!!")
        t0 = 0
    return np.sqrt(p0) * np.exp( -0.5* (tvec/t0)**2)

def gaussfieldB( p0,width, tvec,mode):
    if mode.lower() == 'fwhm':
        t0 = width / 2 / np.sqrt(np.log(2)/2)
    elif mode.lower() == 't0':
        t0 = width
    else:
        print( "gaussfieldB error! no valid width given!!")
        t0 = 0

    return np.sqrt(p0) * np.exp( - (tvec/t0)**2)



def gauss_peak_power( nurep, pmean, taufwhm):
    t0 = taufwhm / np.sqrt(2 * np.log(2))
    return pmean / ( nurep * t0 * np.sqrt( np.pi/2)) 

def sech_peak_power( nurep, pmean, taufwhm):
    t0 = taufwhm / 2 / np.arcsinh(1)
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
    
    



#
# passnotch very simple bandpass or notch binary filters
#

def passnotch(vec,n1,n2,mode="pass"):
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
    lh = len(henespur)
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
    yn = interferogrammspur[xn]
    yn1 = interferogrammspur[xn+1]
    yint = np.multiply(  (yn1-yn), xnint-xn)+yn   
    return np.array(yint)


#
# list / array handling
#


def ge_index(liste, val):
    """ greater-equal index """
    arra = np.array(liste)        
    return np.min( np.nonzero( arra>=val))

def le_index(liste, val):
    """ lower-equal index """
    arra = np.array(liste)        
    return np.max( np.nonzero( arra<=val))
    

#
# little helper
#

def db_abs2(y):
    """
    return logarithmic absolut square of a value in Decibel 
    """
    return 10 * np.log10( np.abs(y)**2)

def db_abs(y):
    return 10 * np.log10( np.abs(y))
