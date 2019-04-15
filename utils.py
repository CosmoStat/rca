import scipy.signal as scisig
import gaussfitter
import numpy as np
from modopt.signal.wavelet import filter_convolve

def apply_transform(data, filters):
    """ Transform ``data`` through application of a set of filters.
    
    Parameters
    ----------
    data: np.ndarray
        Data to be transformed. Should be in RCA format (image index is contained
        on last/2nd axis).
    filters: np.ndarray
        Set of filters.
    """
    data = reg_format(np.copy(data))
    return np.array([filter_convolve(im, filters) for im in data])

def acc_sig_maps(shap_im,ker_stack,sig_est,flux_est,flux_ref,upfact,w,sig_data=None):
    shap = w.shape
    map_out = np.zeros((shap_im[0]*upfact,shap_im[1]*upfact,shap[0]))
    for i in range(0,shap[0]):
        map_out[:,:,i] = acc_sig_map(shap_im,ker_stack,sig_est,flux_est,flux_ref,\
        upfact,w[i,:],sig_data=sig_data)
    return map_out
    
def acc_sig_map(shap_im,ker_stack,sig_est,flux_est,flux_ref,upfact,w,sig_data=None):
    """ Computes the square root of :math:`\mathcal{F}^{2*}(\hat\sigma^2)(A^\\top\odot A^\\top)`.
    See equation (27) in RCA paper.
    Note :math:`\mathrm{Var}(B)` has been replaced by the noise level as estimated from the data,
    and here we do not have the term :math:`\mu` (gradient step size in the paper).
    """
    shap = ker_stack.shape
    nb_im = shap[2]
    if sig_data is None:
        sig_data = np.ones((nb_im,))
    var_stack = np.ones((shap_im[0],shap_im[1],nb_im))
    map2 = np.zeros((shap_im[0]*upfact,shap_im[1]*upfact))
    ker_stack_in = np.copy(ker_stack)**2
    for l in range(0,shap[2]):
        var_stack[:,:,l]*=sig_data[l]**2
        map2 += ((w[l]*flux_est[l]/(sig_est[l]*flux_ref))**2)*scisig.convolve(\
        transpose_decim(var_stack[:,:,l],upfact),ker_stack_in[:,:,l],mode='same')
    sigmap =  np.sqrt(map2)
    return sigmap
    
def return_neighbors(new_pos, obs_pos, vals, n_neighbors):
    """ Find the ``n_neighbors`` nearest neighbors."""
    distances = np.linalg.norm(obs_pos-new_pos, axis=1)
    nbs = vals[np.argsort(distances)[:n_neighbors]]
    pos = obs_pos[np.argsort(distances)[:n_neighbors]]
    return nbs, pos

def rca_format(cube):
    """ Switch from "regular" format to "RCA" format (ie. image index is contained
    on last/2nd axis)"""
    return cube.swapaxes(0,1).swapaxes(1,2)

def reg_format(rca_cube):
    """ Switch from "RCA" format to "regular" format (ie. image index is contained
    on 0th axis)."""
    return rca_cube.swapaxes(2,1).swapaxes(1,0)

def decim(im,d,av_en=1,fft=1):
    """ Decimate image to lower resolution."""

    im_filt= np.copy(im)
    im_d = np.copy(im)
    if d>1:
        if av_en==1:
            siz = d+1-(d%2)
            mask = np.ones((siz,siz))/siz**2
            if fft==1:im_filt = scisig.fftconvolve(im, mask, mode='same')
            else:im_filt = scisig.convolve(im, mask, mode='same')
        n1 = int(np.floor(im.shape[0]/d))
        n2 = int(np.floor(im.shape[1]/d))
        im_d = np.zeros((n1,n2))
        i,j=0,0
        for i in range(0,n1):
            for j in range(0,n2):
                im_d[i,j] = im[i*d,j*d]
    if av_en==1:
        return (im_filt,im_d)
    else:
        return im_d

def pairwise_distances(obs_pos):
    """Computes pairwise distances."""
    ones = np.ones(obs_pos.shape[0])
    out0 = np.outer(obs_pos[:,0], ones)
    out1 = np.outer(obs_pos[:,1], ones)
    return np.sqrt((out0 - out0.T)**2 + (out1 - out1.T)**2)

def transpose_decim(im,decim_fact,av_en=0):
    """ Applies the transpose of the decimation matrix."""
    shap = im.shape
    im_out = np.zeros((shap[0]*decim_fact,shap[1]*decim_fact))

    for i in range(0,shap[0]):
        for j in range(0,shap[1]):
            im_out[decim_fact*i,decim_fact*j]=im[i,j]

    if av_en==1:
        siz = decim_fact+1-(decim_fact%2)
        mask = np.ones((siz,siz))/siz**2
        im_out = scisig.fftconvolve(im, mask, mode='same')

    return im_out

def compute_centroid(im,sigw=None,nb_iter=4):
    """ Computes centroid.
    
    Calls:
    
    * gaussfitter.gaussfit
    """
    if sigw is None:
        param=gaussfitter.gaussfit(im,returnfitimage=False)
        #print param
        sigw = (param[3]+param[4])/2
    sigw = float(sigw)
    n1 = im.shape[0]
    n2 = im.shape[1]
    rx = np.array(range(0,n1))
    ry = np.array(range(0,n2))
    Wc = np.ones((n1,n2))
    centroid = np.zeros((1,2))
    # Four iteration loop to compute the centroid
    i=0
    for i in range(0,nb_iter):

        xx = np.ma.outerproduct(rx-centroid[0,0],np.ones(n2))
        yy = np.ma.outerproduct(np.ones(n1),ry-centroid[0,1])
        W = np.exp(-(xx**2+yy**2)/(2*sigw**2))
        centroid = np.zeros((1,2))
        # Estimate Centroid
        Wc = np.copy(W)
        if i == 0:Wc = np.ones((n1,n2))
        totx=0.0
        toty=0.0
        cx=0
        cy=0

        for cx in range(0,n1):
            centroid[0,0] += (im[cx,:]*Wc[cx,:]).sum()*(cx)
            totx += (im[cx,:]*Wc[cx,:]).sum()
        for cy in range(0,n2):
            centroid[0,1] += (im[:,cy]*Wc[:,cy]).sum()*(cy)
            toty += (im[:,cy]*Wc[:,cy]).sum()
        centroid = centroid*np.array([1/totx,1/toty])


    return (centroid,Wc)

def SoftThresholding(data,thresh):
    """ Performs element-wise soft thresholding."""
    thresh_data = np.copy(data)
    belowmask = (np.abs(data) <= thresh)
    abovemask = np.array(1.-belowmask).astype(bool)
    thresh_data[belowmask] = 0.
    thresh_data[abovemask] = (data - np.sign(data)*thresh)[abovemask]
    return thresh_data
    
def HardThresholding(data,thresh):
    """ Performs element-wise hard thresholding."""
    thresh_data = np.copy(data)
    thresh_data[thresh_data < thresh] = 0.
    return thresh_data
    
def kthresholding(x,k):
    """ Applies k-thresholding (keep only ``k`` highest values, set rest to 0).
    """
    k = int(k)
    if k<1:
        print "Warning: wrong k value for k-thresholding"
        k = 1
    if k>len(x):
        return x
    else:
        xout = np.copy(x)*0
        ind = np.argsort(abs(x))
        xout[ind[-k:]] = x[ind[-k:]]
        return xout

def lineskthresholding(mat,k):
    """ Applies k-thresholding to each line of input matrix.
    
    Calls:
    
    * :func:`utils.kthresholding`
    
    """
    mat_out = np.copy(mat)
    shap = mat.shape
    for j in range(0,shap[0]):
        mat_out[j,:] = kthresholding(mat[j,:],k)
    return mat_out

def mad(x):
    """Computes MAD.
    """
    return np.median(np.abs(x-np.median(x)))

def lanczos(U,n=10,n2=None):
    """Generate Lanczos kernel for a given shift.
    """
    if n2 is None:
        n2 = n
    siz = np.size(U)
    H = None
    if (siz == 2):
        U_in = np.copy(U)
        if len(U.shape)==1:
            U_in = np.zeros((1,2))
            U_in[0,0]=U[0]
            U_in[0,1]=U[1]
        H = np.zeros((2*n+1,2*n2+1))
        if (U_in[0,0] == 0) and (U_in[0,1] == 0):
            H[n,n2] = 1
        else:
            i=0
            j=0
            for i in range(0,2*n+1):
                for j in range(0,2*n2+1):
                    H[i,j] = np.sinc(U_in[0,0]-(i-n))*np.sinc((U_in[0,0]-(i-n))/n
                            )*np.sinc(U_in[0,1]-(j-n))*np.sinc((U_in[0,1]-(j-n))/n)

    else :
        H = np.zeros((2*n+1,))
        for i in range(0,2*n):
            H[i] = np.sinc(np.pi*(U-(i-n)))*np.sinc(np.pi*(U-(i-n))/n)
    return H

def shift_est(psf_stack): 
    """Estimates shifts (see SPRITE paper, section 3.4.1., subsection 'Subpixel shifts').
    #TODO? replace this
    
    Calls:
    
    * gaussfitter.gaussfit
    * :func:`utils.compute_centroid`
    """
    shap = psf_stack.shape
    U = np.zeros((shap[2],2))
    param=gaussfitter.gaussfit(psf_stack[:,:,0],returnfitimage=False)
    #(centroid_ref,Wc) = compute_centroid(psf_stack[:,:,0],(param[3]+param[4])/2)
    centroid_out = np.zeros((shap[2],2))
    for i in range(0,shap[2]):
        param=gaussfitter.gaussfit(psf_stack[:,:,i],returnfitimage=False)
        (centroid,Wc) = compute_centroid(psf_stack[:,:,i],(param[3]+param[4])/2)
        U[i,0] = centroid[0,0]-np.double(shap[0])/2
        U[i,1] = centroid[0,1]-np.double(shap[1])/2
        centroid_out[i,0]  = centroid[0,0]
        centroid_out[i,1]  = centroid[0,1]
    return U,centroid_out

def flux_estimate(im,cent=None,rad=4): # Default value for the flux tunned for Euclid PSF at Euclid resolution
    """Estimate flux for one image (see SPRITE paper, section 3.4.1., subsection 'Photometric flux').
    """
    flux = 0
    if cent is None:
        cent = np.array(np.where(im==im.max())).reshape((1,2))
    shap = im.shape
    for i in range(0,shap[0]):
        for j in range(0,shap[1]):
            if np.sqrt((i-cent[0,0])**2+(j-cent[0,1])**2)<=rad:
                flux = flux+im[i,j]
    return flux

def flux_estimate_stack(stack,cent=None,rad=4):
    """Estimate flux for a bunch of images.
    
    Calls:
    
    * :func:`utils.flux_estimate`
    """
    shap = stack.shape
    flux = np.zeros((shap[2],))
    for i in range(0,shap[2]):
        if cent is not None:
            flux[i] = flux_estimate(stack[:,:,i],cent=cent[i,:],rad=rad)
        else:
            flux[i] = flux_estimate(stack[:,:,i],rad=rad)
    return flux

def shift_ker_stack(shifts,upfact,lanc_rad=4):
    """Generate shifting kernels and rotated shifting kernels.
    
    Calls:
    
    * :func:`utils.lanczos`
    """
    shap = shifts.shape
    shift_ker_stack = np.zeros((2*lanc_rad+1,2*lanc_rad+1,shap[0]))
    shift_ker_stack_adj = np.zeros((2*lanc_rad+1,2*lanc_rad+1,shap[0]))

    for i in range(0,shap[0]):

        uin = shifts[i,:].reshape((1,2))*upfact
        shift_ker_stack[:,:,i] = lanczos(uin,n=lanc_rad)
        shift_ker_stack_adj[:,:,i] = np.rot90(shift_ker_stack[:,:,i],2)

    return shift_ker_stack,shift_ker_stack_adj
        
def gen_Pea(distances, e, a):
    """ Computes :math:`P_{e,a}` matrix for given ``e``, ``a`` couple. See Equations (16-17)
    in RCA paper.
    
    Parameters
    ----------
    distances: np.ndarray
        Array of pairwise distances
    
    e: float
        Exponent to which the pairwise distances should be raised.
    
    a: float
        Constant multiplier along Laplacian's diagonal.
    """
    
    Pea = np.copy(distances**e)
    np.fill_diagonal(Pea, 1.)
    Pea = -1./Pea**e
    for i in range(Pea.shape[0]):
        Pea[i,i] = a*(np.sum(-1.*Pea[i]) - 1.)
    return Pea
    
def select_vstar(eigenvects, R):
    """  Pick best eigenvector from a set of :math:`(e,a)`, i.e., solve (35) from RCA paper.
    
    Parameters
    ----------
    eigenvects: np.ndarray
        Array of eigenvects to be tested over.
    
    R: np.ndarray
        :math:`R_i` matrix.
    """
    loss = np.sum(R**2)
    for i,Pea_eigenvects in enumerate(eigenvects):
        for j,vect in enumerate(Pea_eigenvects):
            colvect = np.copy(vect).reshape(1,-1)
            current_loss = np.sum((R - colvect.T.dot(colvect.dot(R)))**2)
            if current_loss < loss:
                loss = current_loss
                eigen_idx = j
                ea_idx = i
                best_VT = np.copy(Pea_eigenvects)
    return ea_idx, eigen_idx, best_VT
    
class GraphBuilder(object):
    """ GraphBuilder.
    
    This class computes the necessary quantities for RCA's graph constraint.
    
    Parameters
    ----------
    obs_data: np.ndarray
        Observed data.
    obs_pos: np.ndarray
        Corresponding positions.
    n_comp: int
        Number of RCA components.
    n_eigenvects: int
        Maximum number of eigenvectors to consider per :math:`(e,a)` couple. Default is ``None``;
        if not provided, *all* eigenvectors will be considered, which can lead to a poor
        selection of graphs, especially when data is undersampled. Ignored if ``VT`` and
        ``alpha`` are provided.
    n_iter: int
        How many alternations should there be when optimizing over :math:`e` and :math:`a`. Default is 3.
    ea_gridsize: int
        How fine should the logscale grid of :math:`(e,a)` values be. Default is 10.
    distances: np.ndarray
        Pairwise distances for all positions. Default is ``None``; if not provided, will be
        computed from given positions.
    auto_run: bool
        Whether to immediately build the graph quantities. Default is ``True``.
    """
    def __init__(self, obs_data, obs_pos, n_comp, n_eigenvects=None, n_iter=3,
                 ea_gridsize=10, distances=None, auto_run=True, verbose=True):
        self.obs_data = obs_data
        self.obs_pos = obs_pos
        self.n_comp = n_comp
        if n_eigenvects is None:
            self.n_eigenvects = self.obs_data.shape[2]
        else:
            self.n_eigenvects = n_eigenvects
        self.n_iter = n_iter
        self.ea_gridsize = ea_gridsize
        self.verbose = verbose
        
        if distances is None:
            self.distances = pairwise_distances(self.obs_pos)
        else:
            self.distances = distances
        if auto_run:
            self._build_graphs()
        
    def _build_graphs(self):
        """Computes graph-constraint related values, see RCA paper sections 5.2 and (especially) 5.5.3.
            
        """
        shap = self.obs_data.shape
        e_max = self.pick_emax()
        if self.verbose:
            print " > power max = ",e_max
        a_range = np.geomspace(0.01, 1.99, self.ea_gridsize)
        e_range = np.geomspace(0.01, e_max, self.ea_gridsize)
        # initialize R matrix with observations
        R = np.copy(np.transpose(self.obs_data.reshape((shap[0]*shap[1],shap[2]))))

        self.sel_a = []
        self.sel_e = []
        idx = []
        list_eigenvects = []
        for _ in range(self.n_comp): 
            e, a, j, best_VT = self.select_params(R, e_range, a_range)
            self.sel_e += [e]
            self.sel_a += [a]
            idx += [j]
            list_eigenvects += [best_VT]
            vect = best_VT[j].reshape(1,-1)
            R -= vect.T.dot(vect.dot(R))
            if self.verbose:
                print " > selected e: {}\tselected a: {}\t chosen index: {}/{}".format(
                                                             e, a, j, self.n_eigenvects)
        self.VT = np.vstack((eigenvect for eigenvect in list_eigenvects))
        self.alpha = np.zeros((self.n_comp, self.VT.shape[0]))
        for i in range(self.n_comp):
            self.alpha[i, i*self.n_eigenvects+idx[i]] = 1
        
    def pick_emax(self, epsilon=1e-15):
        """ Select maximum value of :math:`e` for the greedy search over set of
        :math:`(e,a)` couples, so that the graph is still fully connected.
        """
        nodiag = np.copy(self.distances)
        nodiag[nodiag==0] = 1e20
        dist_ratios = np.min(nodiag,axis=1) / np.max(self.distances, axis=1)
        r_med = np.min(dist_ratios**2)
        return np.log(epsilon)/np.log(r_med)
        
    def select_params(self, R, e_range, a_range):
        """ Selects :math:`(e,a)` parameters and best eigenvector for current :math:`R_i` matrix.
        
        Parameters
        ----------
        R: np.ndarray
            Current :math:`R_i` matrix (as defined in RCA paper, sect. 5.5.3.)
        e_range: np.ndarray
            List of :math:`e` values to be tested.
        a_range: np.ndarray
            List of :math:`a` values to be tested.
        """
        current_a = 0.5
        for i in range(self.n_iter):
            # optimize over e
            Peas = np.array([gen_Pea(self.distances, e, current_a) 
                                                   for e in e_range])
            all_eigenvects = np.array([self.gen_eigenvects(Pea) for Pea in Peas])
            ea_idx, eigen_idx, _ = select_vstar(all_eigenvects, R)
            current_e = e_range[ea_idx]
            
            # optimize over a
            Peas = np.array([gen_Pea(self.distances, current_e, a) 
                                                   for a in a_range])
            all_eigenvects = np.array([self.gen_eigenvects(Pea) for Pea in Peas])
            ea_idx, eigen_idx, best_VT = select_vstar(all_eigenvects, R)
            current_a = a_range[ea_idx]

        return current_e, current_a, eigen_idx, best_VT
        
    def gen_eigenvects(self, mat): 
        """ Computes input matrix's eigenvectors and keep the ``n_eigenvects`` associated with
        the smallest eigenvalues.
        """
        U, s, vT = np.linalg.svd(mat,full_matrices=True)
        vT = vT[-self.n_eigenvects:]
        return vT
            
class CentroidEstimator(object):
    """ Estimate intra-pixel shifts."""
    def __init__(self, im, sig=7.5, n_iter=5, auto_run=True,
                 xc=None, yc=None):
        self.im = im
        self.stamp_size = im.shape
        self.ranges = np.array([np.arange(i) for i in self.stamp_size])
        self.sig = sig
        self.n_iter=5
        self.xc0, self.yc0 = float(self.stamp_size[0])/2, float(self.stamp_size[1])/2
        if xc is None or yc is None:
            self.xc = self.xc0
            self.yc = self.yc0
        else:
            self.xc = xc
            self.yc = yc
        if auto_run:
            self.estimate()
            
    def UpdateGrid(self):
        self.xx = np.outer(self.ranges[0] - self.xc,
                          np.ones(self.stamp_size[1]))
        self.yy = np.outer(np.ones(self.stamp_size[0]),
                          self.ranges[1] - self.yc)
                          
    def EllipticalGaussian(self, e1=0, e2=0):
        """ Computes an elliptical 2D gaussian with arbitrary centroid."""
        # shear it
        gxx = (1-e1)*self.xx - e2*self.yy
        gyy = (1+e1)*self.yy - e2*self.xx
        # compute elliptical gaussian
        return np.exp(-(gxx ** 2 + gyy ** 2) / (2 * self.sig ** 2))
        
    def ComputeMoments(self):
        Q0 = np.sum(self.im*self.window)
        Q1 = np.array([np.sum(np.sum(self.im * self.window, axis=1-i) * self.ranges[i])
                       for i in range(2)])
        #Q2 = np.array([np.sum(self.im*self.window * self.xx**(2-i) * self.yy**i) 
        #               for i in range(3)])
        self.xc = Q1[0]/Q0
        self.yc = Q1[1]/Q0
            
    def estimate(self):
        for _ in range(self.n_iter):
            self.UpdateGrid()
            self.window = self.EllipticalGaussian()
            # Calculate weighted moments.
            self.ComputeMoments()
        return self.xc, self.yc
        
    def return_shifts(self):
        """ Returns intra-pixel shifts, that is, the difference between the estimated centroid
        and the center of the stamp (or pixel grid)."""
        return [self.xc-self.xc0, self.yc-self.yc0]


