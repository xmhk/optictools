#optictools.py
Rev 1 , oct 02 2013

###gauss\_peak\_power( nurep, pmean, taufwhm)
* returns the peak power of a gaussian pulse from the repetition rate *nurep*, the mean power *pmean* and the full width at half maximum *taufwhm*

### lamvec,Slam = optical\_density\_from\_nu\_to\_lam( nuvec, Snu)
* converts the spectral (e.g. power) density from frequency to wavelength representation.
* return the wavelength vector (**not** aequidistant!) and the spectral density (normalized to length)

### nuvec, Snu = optical\_density\_from\_lam\_to\_nu( lamvec, Slam)
* converts the spectral (e.g. power) density from wavelength to frequency representation.
* return the frequency vector (**not** aequidistant!) and the spectral density (normalized to frequency)