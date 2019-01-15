import numpy as np
from modopt.opt.gradient import GradParent, GradBasic
from modopt.math.matrix import PowerMethod
import utils
from scipy.signal import fftconvolve
        
        
class SourceGrad(GradParent, PowerMethod):
    """Parameters
    ----------
    data : np.ndarray
        Input data array, a array of 2D observed images (i.e. with noise)
    S : np.ndarray
        Current estimation of eigenPSFs
    A : np.ndarray
        Current estimation of corresponding coefficients
    flux : np.ndarray
        Per-object flux value
    sig : np.ndarray
        Noise levels
    ker : np.ndarray
        Shifting kernels
    ker_rot : np.ndarray
        Inverted shifting kernels
    D : float
        Upsampling factor
    data_type : type
        
    Notes
    -----
    The properties of `GradBasic` and `PowerMethod` are inherited in this class
    """

    def __init__(self, data, A, flux, sig, ker, ker_rot, D, data_type=float):
        self._grad_data_type = data_type
        self.obs_data = data
        self.op = self.MX 
        self.trans_op = self.MtX 
        shap = data.shape
        self.shape = (shap[0]*D,shap[1]*D)
        self.D = D
        self.A = np.copy(A)
        self.flux = flux
        self.sig = sig
        self.ker = ker
        self.ker_rot  = ker_rot
        #PowerMethod.__init__(self, self.trans_op_op, (np.prod(self.shape),np.prod(self.shape),A.shape[0]))
        #print " > SPECTRAL RADIUS:\t{}".format(self.spec_rad)
        #TODO: uncomment
        
        self._current_rec = None # stores latest application of self.MX

    def set_A(self,A_new,pwr_en=True):
        self.A = np.copy(A_new)
        if pwr_en:
            PowerMethod.__init__(self, self.trans_op_op, (np.prod(self.shape),np.prod(self.shape),self.A.shape[0]))

    def set_flux(self,flux_new,pwr_en=False):
        self.flux = np.copy(flux_new)
        if pwr_en:
            PowerMethod.__init__(self, self.trans_op_op, (np.prod(self.shape),np.prod(self.shape),self.A.shape[0]))

    def get_flux(self):
        return self.flux

    def degradation_op(self, S, A_i, shift_ker):
        """ Shift and decimate reconstructed PSF."""
        return utils.decim(fftconvolve(S.dot(A_i),
                                       shift_ker,mode='same'),self.D,av_en=0)

    def adjoint_degradation_op(self, S, A_i, shift_ker):
        """ Apply adjoint of the degradation operator."""
        return fftconvolve(utils.transpose_decim(S,self.D),
                           shift_ker,mode='same').dot(A_i)

    def MX(self, S):
        """Apply degradation operator and renormalize.

        Parameters
        ----------
        S : np.ndarray
            Input data array, a set of eigenPSFs

        Returns
        -------
        np.ndarray result

        """
        normfacs = self.flux / (np.median(self.flux)*self.sig)
        dec_rec = np.array([nf * self.degradation_op(S,A_i,shift_ker) for nf,A_i,shift_ker 
                       in zip(normfacs, self.A.T, utils.reg_format(self.ker))])
        self._current_rec = utils.rca_format(dec_rec)
        return self._current_rec

    def MtX(self, x):
        """MtX

        This method calculates the action of the transpose of the matrix Mt on
        the data X

        Parameters
        ----------
        x : np.ndarray
            Input data array, a cube of 2D images

        Returns
        -------
        np.ndarray result

        """
        normfacs = self.flux / (np.median(self.flux)*self.sig)
        upsamp_x = np.array([nf * self.adjoint_degradation_op(x,A_i,shift_ker) for nf,A_i,shift_ker 
                       in zip(normfacs, self.A.T, utils.reg_format(self.ker_rot))])
        return utils.rca_format(upsamp_x)
                
    def cost(self, x, y=None, verbose=False):
        """ Compute data fidelity term. ``y`` is unused (it's just so ``modopt.opt.algorithms.Condat`` can feed
        the dual variable.)
        """
        if isinstance(self._current_rec, type(None)):
            self._current_rec = self.MX(x)
        cost_val = 0.5 * np.linalg.norm(self._current_rec - self.obs_data) ** 2
        if verbose:
            print " > MIN(X):\t{}".format(np.min(x))
            print " > Current cost: {}".format(cost_val)
        return cost_val
                
    def get_grad(self, x):
        """Get the gradient step
        This method calculates the gradient step from the input data
        Parameters
        ----------
        x : np.ndarray
            Input data array
        Returns
        -------
        np.ndarray gradient value
        Notes
        -----
        Calculates M^T (MX - Y)
        """

        self.grad = self.MtX(self.MX(x) - self.obs_data)



