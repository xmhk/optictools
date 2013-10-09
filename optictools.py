import numpy as np
from scipy.misc import factorial

def beta2_curve(om, om0, betas):
    dom = np.array( om) -om0
    b2k = np.zeros(len(om))
    for i in range(2,len(betas)):
#        print i
        b2k = b2k + betas[i]/factorial(i-2) * dom**(i-2)
    return b2k



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
    
    



