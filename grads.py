import numpy as np
from modopt.opt.gradient import GradParent, GradBasic
from modopt.math.matrix import PowerMethod
import utils
from scipy.signal import fftconvolve


def degradation_op(S, A_i, shift_ker, D):
    """ Shift and decimate reconstructed PSF."""
    return utils.decim(fftconvolve(S.dot(A_i),shift_ker,mode='same'),
                       D,av_en=0)

def adjoint_degradation_op(x_i, shift_ker, D):
    """ Apply adjoint of the degradation operator."""
    return fftconvolve(utils.transpose_decim(x_i,D),
                       shift_ker,mode='same')
        
class CoeffGrad(GradParent, PowerMethod):
    def __init__(self, data, S, VT, flux, sig, ker, ker_rot, D, data_type='float'):
        self._grad_data_type = data_type
        self.obs_data = data
        self.op = self.MX 
        self.trans_op = self.MtX 
        self.S = np.copy(S)
        self.VT = VT
        self.flux = flux
        self.sig = sig
        self.ker = ker
        self.ker_rot  = ker_rot
        self.D = D
        # save number of (superresolved) pixels to flatten stacks of images
        self.hr_npix = np.prod(np.array(data.shape[:2])*D)
        # initialize Power Method to compute spectral radius
        PowerMethod.__init__(self, self.trans_op_op, 
                        (S.shape[-1],VT.shape[0]), auto_run=False)
        
        self._current_rec = None # stores latest application of self.MX

    def update_S(self, new_S, update_spectral_radius=True):
        self.S = new_S
        if update_spectral_radius:
            PowerMethod.get_spec_rad(self)

    def MX(self, alph):
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
        A = alph.dot(self.VT) 
        dec_rec = np.array([nf * degradation_op(self.S,A_i,shift_ker,self.D) for nf,A_i,shift_ker 
                       in zip(normfacs, A.T, utils.reg_format(self.ker))])
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
        x = utils.reg_format(x)
        upsamp_x = np.array([nf * adjoint_degradation_op(x_i,shift_ker,self.D) for nf,x_i,shift_ker 
                       in zip(normfacs, x, utils.reg_format(self.ker_rot))])
        x, upsamp_x = utils.rca_format(x), utils.rca_format(upsamp_x)
        STx = self.S.T.reshape(-1,self.hr_npix).dot(upsamp_x.reshape(self.hr_npix,-1))
        return STx.dot(self.VT.T) #aka... "V"
                
    def cost(self, x, y=None, verbose=False):
        """ Compute data fidelity term. ``y`` is unused (it's just so ``modopt.opt.algorithms.Condat`` can feed
        the dual variable.)
        """
        if isinstance(self._current_rec, type(None)):
            self._current_rec = self.MX(x)
        cost_val = 0.5 * np.linalg.norm(self._current_rec - self.obs_data) ** 2
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
      

class SourceGrad(GradParent, PowerMethod):
    """Parameters
    ----------
    data : np.ndarray
        Input data array, a array of 2D observed images (i.e. with noise)
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
        
    Notes
    -----
    The properties of `GradParent` and `PowerMethod` are inherited in this class
    """

    def __init__(self, data, A, flux, sig, ker, ker_rot, D, data_type='float'):
        self._grad_data_type = data_type
        self.obs_data = data
        self.op = self.MX 
        self.trans_op = self.MtX 
        self.A = np.copy(A)
        self.flux = flux
        self.sig = sig
        self.ker = ker
        self.ker_rot  = ker_rot
        self.D = D
        # initialize Power Method to compute spectral radius
        hr_shape = np.array(data.shape[:2])*D
        PowerMethod.__init__(self, self.trans_op_op, 
                        tuple(hr_shape)+(A.shape[0],), auto_run=False)
        
        self._current_rec = None # stores latest application of self.MX

    def update_A(self, new_A, update_spectral_radius=True):
        self.A = new_A
        if update_spectral_radius:
            PowerMethod.get_spec_rad(self)

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
        dec_rec = np.array([nf * degradation_op(S,A_i,shift_ker,self.D) for nf,A_i,shift_ker 
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
        x = utils.reg_format(x)
        upsamp_x = np.array([nf * adjoint_degradation_op(x_i,shift_ker,self.D) for nf,x_i,shift_ker 
                       in zip(normfacs, x, utils.reg_format(self.ker_rot))])
        x, upsamp_x = utils.rca_format(x), utils.rca_format(upsamp_x)
        return upsamp_x.dot(self.A.T)
                
    def cost(self, x, y=None, verbose=False):
        """ Compute data fidelity term. ``y`` is unused (it's just so ``modopt.opt.algorithms.Condat`` can feed
        the dual variable.)
        """
        if isinstance(self._current_rec, type(None)):
            self._current_rec = self.MX(x)
        cost_val = 0.5 * np.linalg.norm(self._current_rec - self.obs_data) ** 2
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



