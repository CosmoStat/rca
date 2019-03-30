import numpy as np
from modopt.opt.gradient import GradParent, GradBasic
from modopt.math.matrix import PowerMethod
import utils
from scipy.signal import fftconvolve
import psf_toolkit as tk

def degradation_op(X, shift_ker, D):
    """ Shift and decimate fine-grid image."""
    return utils.decim(fftconvolve(X,shift_ker,mode='same'),
                       D,av_en=0)

def adjoint_degradation_op(x_i, shift_ker, D):
    """ Apply adjoint of the degradation operator."""
    return fftconvolve(utils.transpose_decim(x_i,D),
                       shift_ker,mode='same')
        
class CoeffGrad(GradParent, PowerMethod):
    """ Gradient class for the coefficient update.
    
    Parameters
    ----------
    data: np.ndarray
        Observed data.
    S: np.ndarray
        Current eigenPSFs :math:`S`.
    VT: np.ndarray
        Matrix of concatenated graph Laplacians.
    flux: np.ndarray
        Per-object flux value.
    sig: np.ndarray
        Noise levels.
    ker: np.ndarray
        Shifting kernels.
    ker_rot: np.ndarray
        Inverted shifting kernels.
    D: float
        Upsampling factor.
    """
    def __init__(self, data, S, VT, flux, sig, ker, ker_rot, D, data_type='float'):
        self._grad_data_type = data_type
        self.obs_data = data
        self.op = self.MX 
        self.trans_op = self.MtX 
        self.VT = VT
        self.flux = flux
        self.sig = sig
        self.ker = ker
        self.ker_rot  = ker_rot
        self.D = D
        # initialize Power Method to compute spectral radius
        PowerMethod.__init__(self, self.trans_op_op, 
                        (S.shape[-1],VT.shape[0]), auto_run=False)
        self.update_S(np.copy(S), update_spectral_radius=False)
        
        self._current_rec = None # stores latest application of self.MX

    def update_S(self, new_S, update_spectral_radius=True):
        """ Update current eigenPSFs."""
        self.S = new_S
        # Apply degradation operator to components
        normfacs = self.flux / (np.median(self.flux)*self.sig)
        self.FdS = np.array([[nf * degradation_op(S_j,shift_ker,self.D) 
                              for nf,shift_ker in zip(normfacs, utils.reg_format(self.ker))] 
                              for S_j in utils.reg_format(self.S)])
        if update_spectral_radius:
            PowerMethod.get_spec_rad(self)

    def MX(self, alpha):
        """Apply degradation operator and renormalize.

        Parameters
        ----------
        alpha: np.ndarray
            Current weights (after factorization by :math:`V^\\top`).
        """
        A = alpha.dot(self.VT) 
        dec_rec = np.empty(self.obs_data.shape)
        for j in range(dec_rec.shape[-1]):
            dec_rec[:,:,j] = np.sum(A[:,j].reshape(-1,1,1)*self.FdS[:,j],axis=0)
        self._current_rec = dec_rec
        return self._current_rec

    def MtX(self, x):
        """Adjoint to degradation operator :func:`MX`.

        Parameters
        ----------
        x : np.ndarray
            Set of finer-grid images.
        """ 
        x = utils.reg_format(x)
        STx = np.array([np.sum(FdS_i*x, axis=(1,2)) for FdS_i in self.FdS])
        return STx.dot(self.VT.T) #aka... "V"
                
    def cost(self, x, y=None, verbose=False):
        """ Compute data fidelity term. ``y`` is unused (it's just so ``modopt.opt.algorithms.Condat`` 
        can feed the dual variable.)
        """
        if isinstance(self._current_rec, type(None)):
            self._current_rec = self.MX(x)
        cost_val = 0.5 * np.linalg.norm(self._current_rec - self.obs_data) ** 2
        return cost_val
                
    def get_grad(self, x):
        """Compute current iteration's gradient.
        """

        self.grad = self.MtX(self.MX(x) - self.obs_data)
      

class SourceGrad(GradParent, PowerMethod):
    """Gradient class for the eigenPSF update.
    
    Parameters
    ----------
    data: np.ndarray
        Input data array, a array of 2D observed images (i.e. with noise).
    A: np.ndarray
        Current estimation of corresponding coefficients.
    flux: np.ndarray
        Per-object flux value.
    sig: np.ndarray
        Noise levels.
    ker: np.ndarray
        Shifting kernels.
    ker_rot: np.ndarray
        Inverted shifting kernels.
    D: float
        Upsampling factor.
    filters: np.ndarray
        Set of filters.
    """

    def __init__(self, data, A, flux, sig, ker, ker_rot, D, filters, data_type='float'):
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
        self.filters = filters
        # initialize Power Method to compute spectral radius
        hr_shape = np.array(data.shape[:2])*D
        PowerMethod.__init__(self, self.trans_op_op, 
                        (A.shape[0],filters.shape[0])+tuple(hr_shape), auto_run=False)
        
        self._current_rec = None # stores latest application of self.MX

    def update_A(self, new_A, update_spectral_radius=True):
        """Update current weights.
        """
        self.A = new_A
        if update_spectral_radius:
            PowerMethod.get_spec_rad(self)

    def MX(self, transf_S):
        """Apply degradation operator and renormalize.

        Parameters
        ----------
        transf_S : np.ndarray
            Current eigenPSFs in Starlet space.

        Returns
        -------
        np.ndarray result

        """
        normfacs = self.flux / (np.median(self.flux)*self.sig)
        S = utils.rca_format(np.sum(transf_S, axis=1)) #[SCALESUMTAG]
        dec_rec = np.array([nf * degradation_op(S.dot(A_i),shift_ker,self.D) for nf,A_i,shift_ker 
                       in zip(normfacs, self.A.T, utils.reg_format(self.ker))])
        self._current_rec = utils.rca_format(dec_rec)
        return self._current_rec

    def MtX(self, x):
        """Adjoint to degradation operator :func:`MX`.

        """
        normfacs = self.flux / (np.median(self.flux)*self.sig)
        x = utils.reg_format(x)
        upsamp_x = np.array([nf * adjoint_degradation_op(x_i,shift_ker,self.D) for nf,x_i,shift_ker 
                       in zip(normfacs, x, utils.reg_format(self.ker_rot))])
        x, upsamp_x = utils.rca_format(x), utils.rca_format(upsamp_x)
        return utils.apply_transform(upsamp_x.dot(self.A.T), self.filters) #[SCALESUMTAG]
                                                                           # ok not really but same difference
                
    def cost(self, x, y=None, verbose=False):
        """ Compute data fidelity term. ``y`` is unused (it's just so 
        ``modopt.opt.algorithms.Condat`` can feed the dual variable.)
        """
        if isinstance(self._current_rec, type(None)):
            self._current_rec = self.MX(x)
        cost_val = 0.5 * np.linalg.norm(self._current_rec - self.obs_data) ** 2
        return cost_val
                
    def get_grad(self, x):
        """Compute current iteration's gradient.
        """

        self.grad = self.MtX(self.MX(x) - self.obs_data)



