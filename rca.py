import numpy as np
import utils
from scipy.signal import fftconvolve
from modopt.signal.wavelet import get_mr_filters, filter_convolve
from modopt.opt.cost import costObj
from modopt.opt.proximity import Positivity
import modopt.opt.algorithms as optimalg
import proxs as rca_prox
import grads
from modopt.opt.reweight import cwbReweight
from scipy.interpolate import Rbf

class RCA(object):
    """ Resolved Components Analysis.
    
    Parameters
    ----------
    n_comp: int
        Number of components to learn.
    upfact: int
        Upsampling factor. Default is 1 (no superresolution).
    ksig: float
        Value of :math:`k` for the thresholding in Starlet domain (taken to be 
        :math:`k\sigma`, where :math:`\sigma` is the estimated noise standard deviation.
    n_scales: int
        Number of Starlet scales to use for the sparsity constraint. Default is 3.
    ksig_init: float
        Similar to ``ksig``, for use when estimating shifts and noise levels, as it might 
        be desirable to have it set higher than ``ksig``. Unused if ``shifts`` are provided 
        when running :func:`RCA.fit`. Default is 5.
    n_scales_init: int
        Similar to ``n_scales``, for use when estimating shifts and noise levels, as it might 
        be sufficient to use fewer scales when initializing. Unused if ``sigs`` are provided
        when running :func:`RCA.fit`. Default is 2.
    verbose: bool or int
        If True, will only output RCA-specific lines to stdout. #TODO: If verbose is set to 2,
        will run ModOpt's optimization algorithms in verbose mode. 
        
    """
    def __init__(self, n_comp, upfact=1, ksig=4, n_scales=3,
                 ksig_init=5, n_scales_init=2, verbose=True):
        self.n_comp = n_comp
        self.upfact = upfact
        self.ksig = ksig
        self.ksig_init = ksig_init
        
        # option strings for mr_transform
        self.opt_sig_init = ['-t2', '-n{}'.format(n_scales_init)]
        self.opt = ['-t2', '-n{}'.format(n_scales)]
        self.verbose = verbose
        self.is_fitted = False
        
    def fit(self, obs_data, obs_pos, S=None, VT=None, alpha=None,
            shifts=None, centroids=None, sigs=None, flux=None,
            nb_iter=2, nb_subiter_S=300, nb_reweight=0, 
            nb_subiter_weights=None, n_eigenvects=None, graph_kwargs={}):
        """ Fits RCA to observed star field.
        
        Parameters
        ----------
        obs_data: np.ndarray
            Observed data.
        obs_pos: np.ndarray
            Corresponding positions.
        S: np.ndarray
            First guess (or warm start) eigenPSFs :math:`S`. Default is ``None``.
        VT: np.ndarray
            Matrix of concatenated graph Laplacians. Default is ``None``.
        alpha: np.ndarray
            First guess (or warm start) weights :math:`\\alpha`, after factorization by ``VT``. Default is ``None``.
        shifts: np.ndarray
            Corresponding sub-pixel shifts. Default is ``None``; will be estimated from
            observed data if not provided.
        centroids: np.ndarray
            Corresponding image centroids. Default is ``None``; will be estimated from
            observed data if not provided.
        sigs: np.ndarray
            Estimated noise levels. Default is ``None``; will be estimated from data
            if not provided.
        flux: np.ndarray
            Flux levels. Default is ``None``; will be estimated from data if not provided.
        nb_iter: int
            Number of overall iterations (i.e. of alternations). Note the weights do not
            get updated the last time around, so they actually get ``nb_iter-1`` updates.
            Default is 300.
        nb_subiter_S: int
            Maximum number of iterations for :math:`S` updates. If ModOpt's optimizers achieve 
            internal convergence, that number may (and often is) not reached. Default is
            300.
        nb_reweight: int 
            Number of reweightings to apply during :math:`S` updates. See equation (33) in RCA paper. 
            Default is 0.
        nb_subiter_weights: int
            Maximum number of iterations for :math:`\\alpha` updates. If ModOpt's optimizers achieve 
            internal convergence, that number may (and often is) not reached. Default is None;
            if not provided, will be set to ``2*nb_subiter_S`` (as it was in RCA v1). 
        n_eigenvects: int
            Maximum number of eigenvectors to consider per :math:`(e,a)` couple. Default is ``None``;
            if not provided, *all* eigenvectors will be considered, which can lead to a poor
            selection of graphs, especially when data is undersampled. Ignored if ``VT`` and
            ``alpha`` are provided.
        graph_kwargs: dictionary
            List of optional kwargs to be passed on to the :func:`utils.GraphBuilder`.
        """
        
        self.obs_data = np.copy(obs_data)
        self.shap = self.obs_data.shape
        self.im_hr_shape = (self.upfact*self.shap[0],self.upfact*self.shap[1],self.shap[2])
        self.obs_pos = obs_pos
        if S is None:
            self.S = np.zeros(self.im_hr_shape[:2] + (self.n_comp,))
        else:
            self.S = S
        self.VT = VT
        self.alpha = alpha
        self.shifts = shifts
        self.centroids = centroids #TODO: compute them from shifts
        self.sigs = sigs
        self.flux = flux
        self.nb_iter = nb_iter
        self.nb_subiter_S = nb_subiter_S
        if nb_subiter_weights is None:
            nb_subiter_weights = 2*nb_subiter_S
        self.nb_subiter_weights = nb_subiter_weights
        self.nb_reweight = nb_reweight
        self.n_eigenvects = n_eigenvects
        self.graph_kwargs = graph_kwargs
            
        if self.verbose:
            print 'Running basic initialization tasks...'
        self._initialize()
        if self.verbose:
            print '... Done.'
        if self.VT is None or self.alpha is None:
            if self.verbose:
                print 'Constructing graph constraint...'
            self._initialize_graph_constraint()
            if self.verbose:
                print '... Done.'
        else:
            self.weights = self.alpha.dot(self.VT)
        self._fit()
        self.is_fitted = True
        return self.S, self.weights
        
    def estimate_psf(self, test_pos, n_neighbors=15, rbf_function='thin_plate'):
        """ Estimate and return PSF at desired positions.
        
        Parameters
        ----------
        test_pos: np.ndarray
            Positions where the PSF should be estimated. Should be in the same format (units,
            etc.) as the ``obs_pos`` fed to :func:`RCA.fit`.
        n_neighbors: int
            Number of neighbors to use for RBF interpolation. Default is 15.
        rbf_function: str
            Type of RBF kernel to use. Default is ``'thin_plate'``.
        """
        if not self.is_fitted:
            raise ValueError('RCA instance has not yet been fitted to observations. Please run\
            the fit method.')
        ntest = test_pos.shape[0]
        test_weights = np.empty((self.n_comp, ntest))
        for j,pos in enumerate(test_pos):
            # determine neighbors
            nbs, pos_nbs = utils.return_neighbors(pos, self.obs_pos, self.weights.T, n_neighbors)
            # train RBF and interpolate for each component
            for i in range(self.n_comp):
                rbfi = Rbf(pos_nbs[:,0], pos_nbs[:,1], nbs[:,i], function=rbf_function)
                test_weights[i,j] = rbfi(pos[0], pos[1])
        return self._transform(test_weights)
        
    def _initialize(self):
        """ Initialization tasks related to noise levels, shifts and flux. Note it includes
        renormalizing observed data, so needs to be ran even if all three are provided."""
        self.init_filters = get_mr_filters(self.shap[:2], opt=self.opt_sig_init, coarse=True)
        # noise levels
        if self.sigs is None:
            transf_data = rca_prox.apply_transform(self.obs_data, self.init_filters)
            sigmads = np.array([1.4826*utils.mad(fs[0]) for fs in transf_data])
            self.sigs = sigmads / np.linalg.norm(self.init_filters[0])
        else:
            self.sigs = np.copy(self.sigs)
        self.sig_min = np.min(self.sigs)
        # intra-pixel shifts
        if self.shifts is None or self.centroids is None:
            thresholds = np.ones(self.shap)
            for i in range(self.shap[2]):
                # don't allow thresholding to be over 80% of maximum observed pixel
                nsig_shifts = min(self.ksig_init,0.8*self.obs_data[:,:,i].max()/self.sigs[i])
                thresholds[:,:,i] *= nsig_shifts*self.sigs[i]
            psf_stack_shift = utils.thresholding_3D(self.obs_data,thresholds,0)
            self.shifts,self.centroids = utils.shift_est(psf_stack_shift)
        self.shift_ker_stack,self.shift_ker_stack_adj = utils.shift_ker_stack(self.shifts,
                                                                              self.upfact)
        # flux levels
        if self.flux is None:
            self.flux = utils.flux_estimate_stack(self.obs_data,rad=4)
        self.flux_ref = np.median(self.flux)
        # Normalize noise levels observed data
        self.sigs /= self.sig_min
        self.obs_data /= self.sigs.reshape(1,1,-1)
    
    def _initialize_graph_constraint(self):
        gber = utils.GraphBuilder(self.obs_data, self.obs_pos, self.n_comp, 
                                  n_eigenvects=self.n_eigenvects, **self.graph_kwargs)
        self.VT, self.alpha, self.distances = gber.VT, gber.alpha, gber.distances
        self.sel_e, self.sel_a = gber.sel_e, gber.sel_a
        self.weights = self.alpha.dot(self.VT)
        
    def _fit(self):
        weights = self.weights
        comp = self.S
        alpha = self.alpha
        #### Source updates set-up ####
        # initialize dual variable and compute Starlet filters for Condat source updates 
        dual_var = np.zeros((self.im_hr_shape))
        self.starlet_filters = get_mr_filters(self.im_hr_shape[:2], opt=self.opt, coarse=True)
        rho_phi = np.sqrt(np.sum(np.sum(np.abs(self.starlet_filters),axis=(1,2))**2))
        
        # Set up source updates, starting with the gradient
        source_grad = grads.SourceGrad(self.obs_data, weights, self.flux, self.sigs, 
                                      self.shift_ker_stack, self.shift_ker_stack_adj, self.upfact)

        # sparsity in Starlet domain prox (this is actually assuming synthesis form)
        sparsity_prox = rca_prox.StarletThreshold(self.starlet_filters, 0)

        # and the linear recombination for the positivity constraint
        lin_recombine = rca_prox.LinRecombine(weights)

        #### Weight updates set-up ####
        # gradient
        weight_grad = grads.CoeffGrad(self.obs_data, comp, self.VT, self.flux, self.sigs, 
                                      self.shift_ker_stack, self.shift_ker_stack_adj, self.upfact)
        
        # cost function
        weight_cost = costObj([weight_grad]) 
        
        # k-thresholding for spatial constraint
        iter_func = lambda x: np.floor(np.sqrt(x))+1
        coeff_prox = rca_prox.KThreshold(iter_func)
        

        for k in range(self.nb_iter):
            " ============================== Sources estimation =============================== "
            # update gradient instance with new weights...
            source_grad.update_A(weights)
            
            # ... update linear recombination weights...
            lin_recombine.update_A(weights)
            
            # ... set optimization parameters...
            beta = source_grad.spec_rad + rho_phi
            tau = 1./beta
            sigma = 1./lin_recombine.norm * beta/2

            # ... update sparsity prox thresholds...
            thresh = utils.reg_format(utils.acc_sig_maps(self.shap,self.shift_ker_stack_adj,self.sigs,
                                                        self.flux,self.flux_ref,self.upfact,weights,
                                                        sig_data=np.ones((self.shap[2],))*self.sig_min))
            thresholds = self.ksig*np.sqrt(np.array([filter_convolve(Sigma_k**2,self.starlet_filters**2) 
                                              for Sigma_k in thresh]))

            sparsity_prox.update_threshold(tau*thresholds)
            
            # and run source update:
            if self.nb_reweight:
                reweighter = cwbReweight(thresholds)
                for _ in range(self.nb_reweight):
                    source_optim = optimalg.Condat(comp, dual_var, source_grad, sparsity_prox,
                                                   Positivity(), linear = lin_recombine,
                                                   max_iter=self.nb_subiter_S, tau=tau, sigma=sigma)
                    comp = source_optim.x_final
                    transf_data = sparsity_prox.starlet_transform(comp)
                    reweighter.reweight(transf_data)
                    thresholds = reweighter.weights 
            else:
                source_optim = optimalg.Condat(comp, dual_var, source_grad, sparsity_prox,
                                               Positivity(), linear = lin_recombine,
                                               max_iter=self.nb_subiter_S, tau=tau, sigma=sigma)
                comp = source_optim.x_final
            
            #TODO: replace line below with Fred's component selection (to be extracted from `low_rank_global_src_est_comb`)
            ind_select = range(comp.shape[2])


            " ============================== Weights estimation =============================== "
            if k < self.nb_iter-1: 
                # update sources and reset iteration counter for K-thresholding
                weight_grad.update_S(comp)
                coeff_prox.reset_iter()
                weight_optim = optimalg.ForwardBackward(alpha, weight_grad, coeff_prox, cost=weight_cost,
                                                beta_param=weight_grad.inv_spec_rad, auto_iterate=False)
                weight_optim.iterate(max_iter=self.nb_subiter_weights)
                alpha = weight_optim.x_final
                weights_k = alpha.dot(self.VT)

                # renormalize to break scale invariance
                weight_norms = np.sqrt(np.sum(weights_k**2,axis=1)) 
                comp *= weight_norms
                weights_k /= weight_norms.reshape(-1,1)
                #TODO: replace line below with Fred's component selection 
                ind_select = range(weights.shape[0])
                weights = weights_k[ind_select,:]
                supports = None #TODO
    
        self.weights = weights
        self.S = comp
        self.alpha = alpha
        source_grad.MX(self.S)
        self.current_rec = source_grad._current_rec

    def _transform(self, weights):
        return self.S.dot(weights)
