from numpy import zeros,size,where,ones,copy,around,double,sinc,random,pi,arange,cos,sin,arccos,transpose,diag,sqrt,arange,floor,exp,array,mean,roots,float64,int,pi,median,rot90,argsort,tile,repeat,squeeze
from numpy.linalg import svd,norm,inv,eigh
import numpy.ma as npma
import scipy.signal as scisig
from pyflann import *
import sys
sys.path.append('../utilities')
import gaussfitter
import datetime,time

from scipy import interpolate

import scipy.linalg as sci_lin

def diagonally_dominated_mat(shap,sig=4.,thresh_en=True,coord_map=None,pol_en=False,pol_mod=False,theta_param=1,cent=None):
    """**[???]**
    
    Calls:
    
    * :func:`optim_utils.dist_map_2`
    * :func:`utils.polar_coord_cloud`"""
    from optim_utils import dist_map_2
    from numpy import sqrt as numsqrt
    coord_cloud = None
    if coord_map is None:
        coord_map = zeros((shap[0],shap[1],2))
        coord_map[:,:,0] = arange(0,shap[0]).reshape((shap[0],1)).dot(ones((1,shap[1])))
        coord_map[:,:,1] = ones((shap[0],1)).dot(arange(0,shap[1]).reshape((1,shap[1])))
        coord_cloud = zeros((2,shap[0]*shap[1]))
        coord_cloud[0,:] = coord_map[:,:,0].reshape((shap[0]*shap[1],))
        coord_cloud[1,:] = coord_map[:,:,1].reshape((shap[0]*shap[1],))
        if pol_en:
            if cent is None:
                cent = array([shap[0]/2,shap[1]/2])
            coord_cloud = polar_coord_cloud(coord_cloud,cent)
            coord_map[:,:,0] = cloud_out[0,:].reshape((shap[0],shap[1]))
            coord_map[:,:,1] = theta_param*cloud_out[1,:].reshape((shap[0],shap[1]))/(2*pi)
            if pol_mod:
                coord_map[:,:,1] *= coord_map[:,:,0]
                coord_cloud[1,:] *= coord_cloud[0,:]
    dist_map = sqrt(dist_map_2(coord_cloud))
    mat = exp(-dist_map**2/sig**2)
    if thresh_en:
        i,j = where(mat>exp(-1.))
        mat*=0
        mat[i,j] = 1.
    mat/=mat.sum()

    return mat

def diagonally_dominated_mat_stack(shap,nb_mat,sig=4.,thresh_en=True,coord_map=None,pol_en=False,\
    pol_mod=False,theta_param=1,cent=None):
    
    """Applies whatever :func:`diagonally_dominated_mat` does to a stack.
    
    Calls:
    
    * :func:`diagonally_dominated_mat`
    """
    mat_ref = diagonally_dominated_mat(shap,sig=sig,thresh_en=thresh_en,coord_map=coord_map,\
    pol_en=pol_en,pol_mod=pol_mod,theta_param=theta_param,cent=cent)
    mat = zeros((shap[0]*shap[1],shap[0]*shap[1],nb_mat))
    for i in range(0,nb_mat):
        mat[:,:,i] = copy(mat_ref)
    return mat

def decim(im,d,av_en=1,fft=1):
    """ Decimate image to lower resolution."""

    im_filt=copy(im)
    im_d = copy(im)
    if d>1:
        if av_en==1:
            siz = d+1-(d%2)
            mask = ones((siz,siz))/siz**2
            if fft==1:im_filt = scisig.fftconvolve(im, mask, mode='same')
            else:im_filt = scisig.convolve(im, mask, mode='same')
        n1 = int(floor(im.shape[0]/d))
        n2 = int(floor(im.shape[1]/d))
        im_d = zeros((n1,n2))
        i,j=0,0
        for i in range(0,n1):
            for j in range(0,n2):
                im_d[i,j] = im[i*d,j*d]
    if av_en==1:
        return (im_filt,im_d)
    else:
        return im_d

def feat_dist_mat(feat_mat):
    """Computes pairwise distances...?
    
    #TODO: maybe some redundancy with :func:`optim_utils.dist_map_2` here?"""
    shap = feat_mat.shape
    mat_out = zeros((shap[0],shap[0]-1))
    a = array(range(0,shap[0]))
    for i in range(0,shap[0]):
        ind_i = where(a!=i)

        for k in range(0,shap[0]-1):
            mat_out[i,k] = sqrt(sum((feat_mat[i,:]-feat_mat[ind_i[0][k],:])**2))
    return mat_out

def transpose_decim(im,decim_fact,av_en=0):
    """ Applies the transpose of the decimation matrix."""
    shap = im.shape
    im_out = zeros((shap[0]*decim_fact,shap[1]*decim_fact))

    for i in range(0,shap[0]):
        for j in range(0,shap[1]):
            im_out[decim_fact*i,decim_fact*j]=im[i,j]

    if av_en==1:
        siz = decim_fact+1-(decim_fact%2)
        mask = ones((siz,siz))/siz**2
        im_out = scisig.fftconvolve(im, mask, mode='same')

    return im_out

def compute_centroid(im,sigw=None,nb_iter=4):
    """ Computes centroid.
    
    #TODO: would be interesting to compare with Sam's moments based computation
    
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
    rx = array(range(0,n1))
    ry = array(range(0,n2))
    Wc = ones((n1,n2))
    centroid = zeros((1,2))
    # Four iteration loop to compute the centroid
    i=0
    for i in range(0,nb_iter):

        xx = npma.outerproduct(rx-centroid[0,0],ones(n2))
        yy = npma.outerproduct(ones(n1),ry-centroid[0,1])
        W = npma.exp(-(xx**2+yy**2)/(2*sigw**2))
        centroid = zeros((1,2))
        # Estimate Centroid
        Wc = copy(W)
        if i == 0:Wc = ones((n1,n2))
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
        centroid = centroid*array([1/totx,1/toty])


    return (centroid,Wc)

def thresholding(x,thresh,thresh_type): 
    """ Performs either soft- (``thresh_type=1``) or hard-thresholding (``thresh_type=0``). Input can be 1D or 2D array.
    """
    xthresh = copy(x)
    n = x.shape

    if len(n)>0:
        n1 = n[0]
    else:
        n1=1
    n2=1
    if len(n)==2:n2 =n[1]
    i,j = 0,0
    if len(n)==2:
        for i in range(0,n1):
            for j in range(0,n2):
                if abs(xthresh[i,j])<thresh[i,j]:xthresh[i,j]=0
                else:
                    if xthresh[i,j]!=0:xthresh[i,j]=(abs(xthresh[i,j])/xthresh[i,j])*(abs(xthresh[i,j])-thresh_type*thresh[i,j])

    elif len(n)==1:
        for i in range(0,n1):
            if abs(xthresh[i])<thresh[i]:xthresh[i]=0
            else:
                if xthresh[i]!=0:xthresh[i]=(abs(xthresh[i])/xthresh[i])*(abs(xthresh[i])-thresh_type*thresh[i])
    elif len(n)==0:
        if abs(xthresh)<thresh:xthresh=0
        else:
            if xthresh!=0:xthresh=(abs(xthresh)/xthresh)*(abs(xthresh)-thresh_type*thresh)

    return xthresh

def kthresholding(x,k):
    """ Applies k-thresholding (keep only k highest values, set rest to 0).
    """
    k = int(k)
    if k<1:
        print "Warning: wrong k value for k-thresholding"
        k = 1
    if k>len(x):
        return x
    else:
        xout = copy(x)*0
        ind = argsort(abs(x))
        xout[ind[-k:]] = x[ind[-k:]]
        return xout

def lineskthresholding(mat,k):
    """ Applies k-thresholding to each line of input matrix.
    
    Calls:
    
    * :func:`utils.kthresholding`
    
    """
    mat_out = copy(mat)
    shap = mat.shape
    for j in range(0,shap[0]):
        mat_out[j,:] = kthresholding(mat[j,:],k)
    return mat_out

def thresholding_3D(x,thresh,thresh_type):
    """Apply thresholding to a set of images (or transport plans I guess).
    
    Calls:
    
    * :func:`utils.thresholding`
    """
    from numpy import copy
    shap = x.shape
    nb_plan = shap[2]
    k=0
    xthresh = copy(x)
    for k in range(0,nb_plan):
        xthresh[:,:,k] = thresholding(copy(x[:,:,k]),thresh[:,:,k],thresh_type)

    return xthresh

def kernel_ext(mat,tol = 0.01): 
    """ Computes input matrix's kernel, defined as the vector space spanned by the eigenvectors corresponding 
    to 1% of the sum of the squared singular values.
    
    #TODO: this is basically just the SVD, all lines between that and the ``return`` are useless.
    """
    U, s, Vt = svd(mat,full_matrices=True)
    e = (s**2).sum()
    eker = 0
    count = 0
    while eker<e*tol:
        count+=1
        eker+=s[-count]**2
    count -=1
    #ker = Vt[-count:,:]
    ker = copy(Vt)

    return ker

def kernel_mat_test_unit(mat,mat_test,tol=0.01):
    """**[???]**
    
    Calls:
    
    * :func:`utils.kernel_ext`
    """
    ker = kernel_ext(mat,tol = tol)

    shap = ker.shape
    nb_vect = shap[0]
    #print "Null space size: ",nb_vect
    loss = (mat_test**2).sum()
    select_ind = -1
    u= None
    uout=None
    for i in range(0,nb_vect):
        u = copy(ker[i,:])
        u = u.reshape((1,shap[1]))
        err = ((mat_test - transpose(u).dot(u.dot(mat_test)))**2).sum()
        if err <loss:
            select_ind = i
            loss = err
            uout = copy(u)
    return loss,uout,ker,select_ind

def kernel_mat_stack_test_unit(mat_stack,mat_test,tol=0):
    """ Computes whatever graph constraint-related quantity :func:`utils.kernel_mat_test_unit` computes
    for a set of matrices.
    
    Calls:
    
    * :func:`utils.kernel_mat_test_unit`
    """
    shap = mat_stack.shape
    nb_mat = shap[2]
    loss = (mat_test**2).sum()
    #print "========== Ref Loss =============",loss
    select_ind = None
    vect_out = None
    ker_out = None
    select_ind_2 = None
    for i in range(0,nb_mat):
        lossi,ui,keri,indi = kernel_mat_test_unit(mat_stack[:,:,i],mat_test,tol=tol)
        #print "========== ith =============",lossi
        if lossi<loss:
            loss = lossi
            vect_out = copy(ui)
            select_ind = i
            select_ind_2 = indi
            ker_out = keri
    return vect_out,select_ind,loss,ker_out,select_ind_2

def log_sampling(val_min,val_max,nb_samp):
    """Literally ``np.logspace`` I think.
    
    #TODO: you know.
    """
    from  numpy import log,double,array,exp
    lval_min = log(val_min)
    lval_max = log(val_max)
    a = double(array(range(0,nb_samp)))/(nb_samp-1)
    a = a*(lval_max-lval_min) + lval_min
    a = exp(a)
    return a

def mad(x):
    """Computes MAD.
    """
    from numpy import *
    return median(abs(x-median(x)))

def cartesian_product(arrays):
    """**[???]***"""
    import numpy as numpy
    broadcastable = numpy.ix_(*arrays)
    broadcasted = numpy.broadcast_arrays(*broadcastable)
    rows, cols = reduce(numpy.multiply, broadcasted[0].shape), len(broadcasted)
    out = numpy.empty(rows * cols, dtype=broadcasted[0].dtype)
    start, end = 0, rows
    for a in broadcasted:
        out[start:end] = a.reshape(-1)
        start, end = end, end + rows
    return out.reshape(cols, rows).T

def get_noise(im,nb_iter=5,k=3):
    """ Estimate noise level for one given image through one soft thresholding iterations.
    See SPRITE paper, appendix... I want to say A.
    
    Calls:
    
    * :func:`utils.mad`
    """
    sig = 1.4826*mad(im)
    for i in range(0,nb_iter):
        im_thresh = im*(abs(im)>k*sig)
        sig = 1.4826*mad(im-im_thresh)
    return sig

def get_noise_arr(arr):
    """Estimate noise for each of a set of image.
    
    Calls:

    * :func:`utils.cartesian_product`
    * :func:`utils.get_noise`
    """
    shap = arr.shape
    ind = list()
    for i in shap[2:]:
        ind.append(arange(0,i))
    coord = cartesian_product(ind)
    noise_map = ones(arr.shape)
    s = slice(None) # equivalent to ':
    for i in range(0,coord.shape[0]):
        sig = get_noise(arr[(s,s)+tuple(coord[i,:])])
        noise_map[(s,s)+tuple(coord[i,:])]*=sig

    return noise_map

def cube_svd(cube,nb_comp=None,ind=None,mean_sub=False):
    """ Performs PCA as initialization.
    
    #TODO: replace with Scikit Learn version? Pretty sure this is where RCA crashes with
    too many inputs
    """
    shap = cube.shape
    if nb_comp is None:
        nb_comp = min(shap[0]*shap[1],shap[2])
    mat = cube.reshape((shap[0]*shap[1],shap[2]))
    data_mean = None
    centered_data = None
    if ind is None:
        ind = range(0,shap[2])
    if mean_sub:
        data_mean = cube.mean(axis=2)
        mat -= data_mean.reshape((shap[0]*shap[1],1)).dot(ones((1,shap[2])))
        centered_data = copy(cube)
        for i in range(0,shap[2]):
            centered_data[:,:,i]-=data_mean
    U, s, Vt = svd(mat[:,ind],full_matrices=False)
    shap_u = U.shape
    coeff = transpose(U[:,0:nb_comp]).dot(mat)
    approx = U[:,0:nb_comp].dot(coeff)
    comp_cube = U[:,0:nb_comp].reshape((shap[0],shap[1],min(nb_comp,shap[2])))
    approx_cube =  approx.reshape((shap[0],shap[1],shap[2]))
    if mean_sub:
        return coeff,comp_cube,approx_cube,data_mean,centered_data
    else:
        return coeff,comp_cube,approx_cube

def lanczos(U,n=10,n2=None):
    """Generate Lanczos kernel for a given shift.
    """
    if n2 is None:
        n2 = n
    siz = size(U)
    H = None
    if (siz == 2):
        U_in = copy(U)
        if len(U.shape)==1:
            U_in = zeros((1,2))
            U_in[0,0]=U[0]
            U_in[0,1]=U[1]
        H = zeros((2*n+1,2*n2+1))
        if (U_in[0,0] == 0) and (U_in[0,1] == 0):
            H[n,n2] = 1
        else:
            i=0
            j=0
            for i in range(0,2*n+1):
                for j in range(0,2*n2+1):
                    H[i,j] = sinc(U_in[0,0]-(i-n))*sinc((U_in[0,0]-(i-n))/n)*sinc(U_in[0,1]-(j-n))*sinc((U_in[0,1]-(j-n))/n)

    else :
        H = zeros((2*n+1,))
        for i in range(0,2*n):
            H[i] = sinc(pi*(U-(i-n)))*sinc(pi*(U-(i-n))/n)
    return H

def mat_to_cube(mat,n1,n2):
    """ Literally ``np.swapaxes`` I think.
    
    #TODO
    """
    shap = mat.shape
    cube = zeros((n1,n2,shap[0]))
    k=0
    for k in range(0,shap[0]):
        cube[:,:,k] = mat[k,:].reshape(n1,n2)
    return cube

def knn_interf(data,nb_neigh,return_index=False):
    """ Computes closest neighbors. 
    
    #TODO: Probably some costly redundancy with :func:`psf_learning_utils.full_displacement` here.
    """
    from numpy import array,float64
    flann = FLANN()
    data_cp = array(data, dtype=float64)
    params = flann.build_index(data_cp)
    result_temp, dists_temp = flann.nn_index(data_cp, nb_neigh+1)
    dists = dists_temp[:,1:nb_neigh+1]
    result = result_temp[:,1:nb_neigh+1]
    if return_index:
        return result,dists,flann
    else:
        return result,dists

def im_gauss_nois_est(im,opt=['-t2','-n2'],filters=None):
    """Compute sigma mad for... What appears to be the first wavelet scale only?
    
    Calls:
    
    * isap.mr_trans_2
    """
    from isap import mr_trans_2
    Result,filters = mr_trans_2(im,filters=filters,opt=opt)
    siz = im.shape
    norm_wav = norm(filters[:,:,0])
    sigma = 1.4826*mad(Result[:,:,0])/norm_wav

    return sigma,filters

def im_gauss_nois_est_cube(cube,opt=None,filters=None,return_map=False):
    """
    Estimate sigma mad for a set of images.
    
    #TODO: Note there is clearly something wrong with ``return_map`` since it fills in a ``map`` but 
    does not return it (just the boolean saying it was filled).
    
    Calls:
    
    * :func:`utils.im_gauss_nois_est`
    """
    shap = cube.shape
    sig = zeros((shap[2],))
    map = None
    if return_map:
        map =ones(shap)

    for i in range(0,shap[2]):
        sig_i,filters = im_gauss_nois_est(cube[:,:,i],opt=opt,filters=filters)
        sig[i] = sig_i
        if return_map:
            map[:,:,i] *= sig[i]
    if return_map:
        return sig,filters,return_map
    else:
        return sig,filters

def polar_coord(coord,cent):
    """Computes polar coordinates for one cartesian position (angles are in radians).
    """
    from numpy import sqrt,arccos,abs,pi
    r = sqrt((coord[0]-cent[0])**2+(coord[1]-cent[1])**2)
    alpha=None
    if r==0:
        alpha=0
    else:
        alpha_ref = arccos(abs(coord[0]-cent[0])/r)
        alpha = alpha_ref
        if (coord[0]-cent[0]<0) and (coord[1]-cent[1]>=0):
            alpha = pi-alpha
        elif (coord[0]-cent[0]<0) and (coord[1]-cent[1]<0):
            alpha = pi+alpha
        elif (coord[0]-cent[0]>=0) and (coord[1]-cent[1]<0):
            alpha = 2*pi-alpha

    return r,alpha

def polar_coord_cloud(coord,cent):
    """Compute polar coordinates for a set of points.
    
    Calls:
    
    * :func:`utils.polar_coord`
    """
    from numpy import zeros
    shap = coord.shape
    out = zeros((2,shap[1]))

    for i in range(0,shap[1]):
        out[0,i],out[1,i] = polar_coord(coord[:,i],cent)

    return out

def shift_est(psf_stack): 
    """Estimates shifts (see SPRITE paper, section 3.4.1., subsection 'Subpixel shifts').
    
    Calls:
    
    * gaussfitter.gaussfit
    * :func:`utils.compute_centroid`
    """
    shap = psf_stack.shape
    U = zeros((shap[2],2))
    param=gaussfitter.gaussfit(psf_stack[:,:,0],returnfitimage=False)
    #(centroid_ref,Wc) = compute_centroid(psf_stack[:,:,0],(param[3]+param[4])/2)
    centroid_out = zeros((shap[2],2))
    for i in range(0,shap[2]):
        param=gaussfitter.gaussfit(psf_stack[:,:,i],returnfitimage=False)
        (centroid,Wc) = compute_centroid(psf_stack[:,:,i],(param[3]+param[4])/2)
        U[i,0] = centroid[0,0]-double(shap[0])/2
        U[i,1] = centroid[0,1]-double(shap[1])/2
        centroid_out[i,0]  = centroid[0,0]
        centroid_out[i,1]  = centroid[0,1]
    return U,centroid_out

def flux_estimate(im,cent=None,rad=4): # Default value for the flux tunned for Euclid PSF at Euclid resolution
    """Estimate flux for one image (see SPRITE paper, section 3.4.1., subsection 'Photometric flux').
    """
    flux = 0
    if cent is None:
        cent = array(where(im==im.max())).reshape((1,2))
    shap = im.shape
    for i in range(0,shap[0]):
        for j in range(0,shap[1]):
            if sqrt((i-cent[0,0])**2+(j-cent[0,1])**2)<=rad:
                flux = flux+im[i,j]
    return flux

def flux_estimate_stack(stack,cent=None,rad=4):
    """Estimate flux for a bunch of images.
    
    Calls:
    
    * :func:`utils.flux_estimate`
    """
    shap = stack.shape
    flux = zeros((shap[2],))
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
    from numpy import rot90
    shap = shifts.shape
    shift_ker_stack = zeros((2*lanc_rad+1,2*lanc_rad+1,shap[0]))
    shift_ker_stack_adj = zeros((2*lanc_rad+1,2*lanc_rad+1,shap[0]))

    for i in range(0,shap[0]):

        uin = shifts[i,:].reshape((1,2))*upfact
        shift_ker_stack[:,:,i] = lanczos(uin,n=lanc_rad)
        shift_ker_stack_adj[:,:,i] = rot90(shift_ker_stack[:,:,i],2)

    return shift_ker_stack,shift_ker_stack_adj
    
def rand_file_name(ext):
    """ Generates random file name. Called by `isap`. Super unsafe!
    #TODO: get rid of it."""
    current_time = datetime.datetime.now().time()
    return 'file'+str(time.clock())+ext
