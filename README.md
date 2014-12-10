#optictools.py
Rev 35, Dec 10, 2014

## Dispersion tools

###beta2\_curve(om, om0, betas)
   calculate the group velocity curve from a beta series
        
        INPUT:
        - omegas: a vector of angular frequencies
        - omega0: the center angular frequency
        - betas: a list with the beta coefficients: [beta0, beta1, beta2, ...]
        
        OUTPUT: 
        - the beta2(omega) curve

###convert\_b2curve\_to\_dcurve(b2k, omvec)
  convert a GVD curve in the frequency representation (beta2) to a GVD curve (D) in the wavelength representation
        
        INPUT:
        -b2vals: (array of beta2 values)
        -omegas: angular frequency vector
        
        OUTPUT:
        -D:      GVD in wavelength representation


### convert\_dcurve\_to\_b2curve(dvals,omegas)
   
    convert a GVD curve in the  wavelength representation (D) to a GVD curve (beta2) in the frequency representation

    INPUT:
    -dvals: (array of D values)
    -omegas: angular frequency vector
    
    OUTPUT:
    -beta2 :      GVD in frequency representation

### beta0\_curve(om, om0, betas)
calculate the dispersion curve from a beta series

        INPUT:
        - omegas: a vector of angular frequencies
        - omega0: the center angular frequency
        - betas: a list with the beta coefficients: [beta0, beta1, beta2, ...]
        
        OUTPUT: 
        - the beta(omega) curve

### beta0\_curve\_from\_b2data( omvec, b2data, om0)


get the beta curve from (numerical beta data
        
        this is done via integration
        
        INPUT:
        - omegas: aequidistant(!) vector of angular frequencies
        - b2data: beta2 data for omegas
        - omega0: center frequency (result will be zero here)
        
        OUTPUT:
        - betacurve(omegas)
    

### get\_even\_part( omvec, om0, k_curve)

extract the even part of a (dispersion) curve
        
        INPUT:
        -omvec - angular frequency vector
        -om0   - center frequency (reference)
        -k_curve - dispersion curve
	OUTPUT:
	-even part of the dispersion curve

### poly2beta(p,xo)
convert the polynomial coefficients of a beta2(omega) fit to a beta series
        
        INPUT:
        - p: polynomcoefficients from polyfit beta2(omega)
        - omega0: center angular frequency for extension
        
        OUTPUT: 
        -betas: a list containing [beta0, beta1, beta2, beta3, ...]



### poly2beta\_B(p,xo)

convert the polynomial coefficients of a beta(omega) fit to a beta series
        
        INPUT:
        - p: polynomcoefficients from polyfit of beta(omega)
        - omega0: center angular frequency for extension
        
        OUTPUT: 
        -betas: a list containing [beta0, beta1, beta2, beta3, ...]

## functions for generate / measure pulses

### fwhm3(list, peakpos=-1):
calculates the full width at half maximum (fwhm) of some curve.
        
the function will return the fwhm with sub-pixel interpolation. 
it will start at the maximum position and walk left and right until it approaches the half values.
        
        INPUT: 
        - valuelist: e.g. the list containing the temporal shape of a pulse 
        
        OPTIONAL INPUT: 
        -peakpos: position of the peak to examine (list index)
        the global maximum will be used if omitted.
        
        OUTPUT:
        -fwhm (value)
    


### pyfindpeaks( environment, list , thresh)
find peak positions in a list of values
        
        INPUT:
        - environment: (INT) a maxima has to be the local maximum in this environment of points
        - valuelist: list or array of points to find the maxima in
        - thresh: a maximum has to be larger than this value
        
        OUTPUT:
        - listindices: positions of the peaks found


### cfindpeaks(env, liste, threshval)
find peaks in a list or an array of value
        
        this is a python wrapper for the C-lib libfindpeaks (github/xhmk)
        
        INPUT:
        - environment: (INT) a maxima has to be the local maximum in this environment of points
        - valuelist: list or array of points to find the maxima in
        - thresh: a maximum has to be larger than this value
        
        OUTPUT:
        - listindices: positions of the peaks found

### sechfield( p0, width, tvec,mode)
 returns the field of a temporal sech pulse
        
        INPUT:
        - p0: optical power in W
        - width: temporal width in s
        - tvec: time vector
        - mode: can be either 
                'fwhm'  (full width at half maximum of intensity)
                or 't0' (argument of sech)
        OUTPUT:
        - temporal sech field


### gaussfieldA( p0,width, tvec,mode):

 returns the field of a gaussian pulse (A)
        type A: Intensity = 1/e**2 at T0
                                                                                                                    
   
        INPUT:
        - p0: optical power in W                                                                                        
        - width: temporal width in s                                                                                    
        - tvec: time vector                                                                                             
        - mode: can be either                                                                                     
    
                -'fwhm'  (full width at half maximum of intensity)                                                       
                -or 't0' (argument of exp)          
        OUTPUT:
        - temporal gaussian field (A)

### gaussfieldB( p0,width, tvec,mode):
returns the field of a gaussian pulse (B)
        type B: Intensity = 1/e at T0
                                                                                                                       
        INPUT:
        - p0: optical power in W                                                                                        
        - width: temporal width in s                                                                                    
        - tvec: time vector                                                                                             
        - mode: can be either                                                                                           
                -'fwhm'  (full width at half maximum of intensity)                                                       
                -or 't0' (argument of exp)                                                                             
        OUTPUT:
        - temporal gaussian field (B)

### sech\_peak\_power( nurep, pmean, taufwhm):

calculate the peak power of sech from the repetition frequency, the mean power and the fwhm
        
        INPUT:
        - nurep: repetition frequency
        - pmean: mean power
        - taufwhm: temporal fwhm (intensity) of the pulses
        
        OUTPUT:
        - peak power of the pulses


###gauss\_peak\_power( nurep, pmean, taufwhm)
calculate the peak power of gaussian pulses (A) from the repetition frequency, the mean power and the fwhm
        
        INPUT:
        - nurep: repetition frequency
        - pmean: mean power
        - taufwhm: temporal fwhm (intensity) of the pulses
        
        OUTPUT:
        - peak power of the pulses

## spectral conversion

### nuvec\_from\_tvec(tvec)
 calculate the relative frequency vector 
 for a given time vector, e.g. for fft.

    INPUT:
    - tvec time vector (aequidistant)
    OUTPUT:
    - nuvec relative freq vector
 



### lamvec,Slam = optical\_density\_from\_nu\_to\_lam( nuvec, Snu)
 convert a spectral density from frequency to wavelength representation
        
        INPUT:
        - nuvec: vector of optical frequency
        - Snu: spectral density Snu(nuvec)
        
        OUTPUT:
        - lamvec: vector of wavelengths
        - Slam: spectral density Slam(lamvec)


### nuvec, Snu = optical\_density\_from\_lam\_to\_nu( lamvec, Slam)
convert a spectral density from wavlength to frequency representation
        
        INPUT:
        - lamvec: vector of wavelengths
        - Slam: spectral density Slam(lamvec)
        
        OUTPUT:
        - nuvec: vector of optical frequency
        - Snu: spectral density Snu(nuvec)


### specnorm( nuv, snuv, pnorm = 1.0)

 normalize spectra to have an integral power of pnorm,
    which is 1.0 as standard
    
    INPUT: 
        -nuv nu-vector (Hz)
        -snu S(nu) unit/Hz
        -OPTIONAL pnorm = 

    OUTPUT:
        -nunorm (interpolated to be aequidistant)
        -Snunorm



### frog\_corr(field1, field2)
    
    calculate the spectrogramm of two fields.
    
    calculated results spectrogramm that would be yielded by
    cross-correlating the fields and than measure the fourier
    transform.

    INPUT: field1, field2
    (SHG-FROG-trace when field1 == field2)

    OUTPUT: MF,
    with MF = | Fourier Transform ( int dt field(t) * field2(t-tau) | **2
    (first axis-> nuvec, second axis->tau)
    

## coherence

### calc_g12( fieldlist, verbose=False ) 
calculate the absolute value of the complex spectral coherence function g12

as defined in J. M. Dudley,  S. Coen: Opt. Lett. 27, 1180–1182 (2002).

    INPUT:
    - a list of fields [field1, field2, ...]
    - [optional] verbose = True/False 

    OUTPUT:
    - a dict containing the fields:
       - fsqmean : mean value of the intensity (|A|**2)
       - fsqstd  : standard deviation of the intensity (|A|**2)
       - g12 : absolute value of g12


## little helpers



### specfilter\_on\_timefield( field, filt )
apply a spectral filter on a timefield, return the filtered timefield (IFFT type)



### supergauss( xvec, mu, width, order)
    supergauss - function
    INPUT:
        - xvec time( or x-)-vector
    	- mu time offset
    	- width
	- order


### sg( xvec, mu, width, order)
brief version of a supergauss - function
   

 

### passnotch(vec,n1,n2,mode="pass")
a very simple bandpass or notch binary filter
        
        INPUT:
        - vec: a monotonously growing vector e.g of frequency, time
        - n1,n2: border frequency of notch/pass
        - mode: either
              'notch': create a filter with zeros between n1 and n2, ones otherwise
           or 'pass': create a filter with ones between n1 and n2, zeros otherwise
         
        OUTPUT:
        - a vector with the same length as vec. Filled with zeros and ones


###ge\_index(liste, val)
greater-equal index 
        
        returns the index of the first value in valuelist that is greater or equal val
        INPUT:
        - valuelist 
        - val
        
        OUTPUT:
        - index

###le\_index(liste, val)
 lower-equal index
        
        returns the index of the first value in valuelist that is smaller or equal val
        INPUT:
        - valuelist 
        - val
        
        OUTPUT:
        - index


###db\_abs2(y)
return the decadic logarithm of the absolut square of a value (Decibel)
        
        INPUT:
        -y: value
        
        OUTPUT:
        -dby = 10 * np.log10( np.abs(y)**2)

###db\_abs(y)
return the decadic logarithm of the absolut of a value (Decibel)
        
        INPUT:
        -y: value
        
        OUTPUT:
        -dby = 10 * np.log10( np.abs(y) )



### get\_slice(signal, xvec, lower,upper)
    cuts a slice from an array
    
    INPUT:
    - signal
    - xvec
    - lower x
    - upper x

    OUTPUT:
    - slice of input array

### moving\_average(somearray, environment)
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
