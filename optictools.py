import numpy as np


def gauss_peak_power( nurep, pmean, taufwhm):
    t0 = taufwhm / np.sqrt(2 * np.log(2))
    return pmean / ( nurep * t0 * np.sqrt( np.pi/2)) 




