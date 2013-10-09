#optictools.py
Rev 3 , oct 09 2013

###beta2\_curve(om, om0, betas)
* returns the group-velocity dispersion curve for a given omega vector **om**
* **om0** is the reference frequency
* **betas** has to be in the form [beta0, beta1, beta2, beta3, ...]

###fwhm3(list, peakpos=-1):
* get the full width at half maximum of a strutured list
* returns the fwhm of the global maximum when no peakpos is given


### sechfield( p0, width, tvec,mode)
* returns the **field** of a sech shaped pulse np.sqrt(p0) * 1/np.cosh( tvec / t0)
* **width** can be given as t0 width **mode='t0'** or as fwhm (power) width **mode="fwhm"**

### gaussfieldA( p0,width, tvec,mode):
* returns the **field** of a gaussian type A (P(t0)=1/e
* functional form np.sqrt(p0) * np.exp( -0.5* (tvec/t0)**2)
* **width** can be given as t0 width **mode='t0'** or as fwhm (power) width **mode="fwhm"**

### gaussfieldB( p0,width, tvec,mode):
* returns the **field** of a gaussian type B (P(t0)=1/e**2
* functional form np.sqrt(p0) * np.exp( - (tvec/t0)**2)
* **width** can be given as t0 width **mode='t0'** or as fwhm (power) width **mode="fwhm"**

### sech_peak_power( nurep, pmean, taufwhm):
* returns the peak power of a sech**2 pulse from the repetition rate *nurep*, the mean power *pmean* and the full width at half maximum *taufwhm*


###gauss\_peak\_power( nurep, pmean, taufwhm)
* !!! unchecked ! 
* returns the peak power of a gaussian pulse from the repetition rate *nurep*, the mean power *pmean* and the full width at half maximum *taufwhm*



### lamvec,Slam = optical\_density\_from\_nu\_to\_lam( nuvec, Snu)
* converts the spectral (e.g. power) density from frequency to wavelength representation.
* return the wavelength vector (**not** aequidistant!) and the spectral density (normalized to length)
* lamvec-Slam pairs are already sorted by wavelength (ascending)

### nuvec, Snu = optical\_density\_from\_lam\_to\_nu( lamvec, Slam)
* converts the spectral (e.g. power) density from wavelength to frequency representation.
* return the frequency vector (**not** aequidistant!) and the spectral density (normalized to frequency)
* nuvec-Snu pairs are already sorted by frequency (ascending)