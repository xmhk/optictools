import numpy as np
from scipy.misc import factorial
from scipy.interpolate import interp1d

def beta2_curve(om, om0, betas):
    dom = np.array( om) -om0
    b2k = np.zeros(len(om))
    for i in range(2,len(betas)):
#        print i
        b2k = b2k + betas[i]/factorial(i-2) * dom**(i-2)
    return b2k

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

# -----------------------------------------------------------------------------


def sechfield( p0, width, tvec,mode):
    if mode.lower() == 'fwhm':
        t0 = width/ 2 / np.arcsinh(1)
    elif mode.lower() == 't0':
        t0 = width
    else:
        print "sechfield error! no valid width given!!"
        t0 = 0
    return np.sqrt(p0) * 1/np.cosh( tvec / t0)

def gaussfieldA( p0,width, tvec,mode):
    if mode.lower() == 'fwhm':
        t0 = width / 2 / np.sqrt(np.log(2))
    elif mode.lower() == 't0':
        t0 = width
    else:
        print "gaussfieldA error! no valid width given!!"
        t0 = 0
    return np.sqrt(p0) * np.exp( -0.5* (tvec/t0)**2)

def gaussfieldB( p0,width, tvec,mode):
    if mode.lower() == 'fwhm':
        t0 = width / 2 / np.sqrt(np.log(2)/2)
    elif mode.lower() == 't0':
        t0 = width
    else:
        print "gaussfieldB error! no valid width given!!"
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
    
    



