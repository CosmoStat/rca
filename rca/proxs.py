""" Defines proximal operators to be fed to ModOpt algorithm that are specific to RCA
(or rather, not currently in ``modopt.opt.proximity``)."""

from __future__ import absolute_import, print_function
import numpy as np
from modopt.signal.noise import thresh
from modopt.opt.linear import LinearParent
from modopt.signal.wavelet import filter_convolve
from rca.utils import lineskthresholding, reg_format, rca_format, SoftThresholding, apply_transform
     
class LinRecombine(object):
    """ Multiply eigenvectors and (factorized) weights."""
    def __init__(self, A, filters, compute_norm=False):
        self.A = A
        self.op = self.recombine
        self.adj_op = self.adj_rec
        self.filters = filters
        if compute_norm:
            U, s, Vt = np.linalg.svd(self.A.dot(self.A.T),full_matrices=False)
            self.norm = np.sqrt(s[0])
        
    def recombine(self, transf_S):
        S = np.array([filter_convolve(transf_Sj, self.filters, filter_rot=True)
                      for transf_Sj in transf_S])
        return rca_format(S).dot(self.A)
        
    def adj_rec(self, Y):
        return apply_transform(Y.dot(self.A.T), self.filters)
        
    def update_A(self, new_A, update_norm=True):
        self.A = new_A
        if update_norm:
            U, s, Vt = np.linalg.svd(self.A.dot(self.A.T),full_matrices=False)
            self.norm = np.sqrt(s[0])

class KThreshold(object):
    """This class defines linewise hard-thresholding operator with variable thresholds.

    Parameters
    ----------
    iter_func: function
        Input function that calcultates the number of non-zero values to keep in each line at each iteration.
    """
    def __init__(self, iter_func):

        self.iter_func = iter_func
        self.iter = 0

    def reset_iter(self):
        """Set iteration counter to zero.
        """
        self.iter = 0


    def op(self, data, extra_factor=1.0):
        """Return input data after thresholding.
        """
        self.iter += 1
        return lineskthresholding(data,self.iter_func(self.iter))
        
    def cost(self, x):
        """Returns 0 (Indicator of :math:`\Omega` is either 0 or infinity).
        """
        return 0 
        
class StarletThreshold(object):
    """Apply soft thresholding in Starlet domain.

    Parameters
    ----------
    threshold: np.ndarray
        Threshold levels.
    thresh_type: str
        Whether soft- or hard-thresholding should be used. Default is ``'soft'``.
    """
    def __init__(self, threshold, thresh_type='soft'):
        self.threshold = threshold
        self._thresh_type = thresh_type

    def update_threshold(self, new_threshold, new_thresh_type=None):
        self.threshold = new_threshold
        if new_thresh_type in ['soft', 'hard']:
            self._thresh_type = new_thresh_type

    def op(self, transf_data, **kwargs):
        """Applies Starlet transform and perform thresholding.
        """
        # Threshold all scales but the coarse
        transf_data[:,:-1] = SoftThresholding(transf_data[:,:-1], self.threshold[:,:-1])
        return transf_data

    def cost(self, x, y):
        return 0 #TODO
