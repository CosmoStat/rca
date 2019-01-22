""" Defines proximal operators to be fed to ModOpt algorithm that are specific to RCA
(or rather, not currently in `modopt.opt.proximity`)."""

from utils import lineskthresholding, reg_format, rca_format, SoftThresholding
from modopt.signal.noise import thresh
from modopt.signal.wavelet import filter_convolve
import numpy as np


class KThreshold(object):
    """ KThreshold proximity operator

    This class defines linewise hard threshold operator with variable threshold

    Parameters
    ----------
    iter_func : function
        Input function that calcultates the number of non-zero values to keep in each line at each iteration
        
    Calls:
    
    * :func:`utils.lineskthresholding`

    """
    def __init__(self, iter_func):

        self.iter_func = iter_func
        self.iter = 0

    def reset_iter(self):
        """Reset iter

        This method sets the iterations counter to zero

        """
        self.iter = 0


    def op(self, data, extra_factor=1.0):
        """Operator

        This method returns the input data thresholded

        Parameters
        ----------
        data : np.ndarray
            Input data array
        extra_factor : float
            Additional multiplication factor

        Returns
        -------
        np.ndarray thresholded data

        """


        self.iter += 1

        return lineskthresholding(data,self.iter_func(self.iter))
        
        
class StarletThreshold(object):
    """Apply soft thresholding in Starlet domain.
    """

    def __init__(self, filters, threshold, thresh_type='soft'):
        self.threshold = threshold
        # get Starlet filters
        self._filters = filters
        self._thresh_type = thresh_type

    def op(self, data, **kwargs):
        """Operator

        Parameters
        ----------
        data : np.ndarray
            Input data array

        Returns
        -------
        np.ndarray thresholded and recomposed data

        """
        data = reg_format(data)
        transf_data = np.array([filter_convolve(im, self._filters) for im in data])
        # Threshold all scales but the coarse
        transf_data[:,:-1] = SoftThresholding(transf_data[:,:-1], self.threshold[:,:-1])
        thresh_data = np.sum(transf_data, axis=1)
        data, thresh_data = rca_format(data), rca_format(thresh_data)
        return thresh_data

    def cost(self, x, y):
        return 0 #TODO
