from numpy import zeros,size,where,ones,copy,around,double,sinc,random,pi,arange,\
cos,sin,arccos,transpose,diag,sqrt,arange,floor,array,mean,roots,float64,int,\
pi,median,rot90,argsort,tile,repeat,squeeze,swapaxes,log
from numpy.linalg import svd,norm,inv,eigh
from numpy.ma import outerproduct,exp
import subprocess
import os
from astropy.io import fits
from scipy.signal import fftconvolve,convolve
import gaussfitter
from pyflann import FLANN
import datetime

def mr_transf_transp(coeff,filters_rot):
    coeff_temp = copy(coeff)
    shap = coeff.shape
    for i in range(0,shap[2]):
        coeff_temp[:,:,i] = fftconvolve(coeff[:,:,i],filters_rot[:,:,i],mode='same')
    out = coeff_temp.sum(axis=2)
    return out

def mr_transf_transp_stack(coeff,filters_rot):
    shap = coeff.shape
    output = zeros((shap[0],shap[1],shap[3]))
    for i in range(0,shap[3]):
        output[:,:,i] = mr_transf_transp(coeff[:,:,:,i],filters_rot)
    return output

def mr_trans(im,opt=None):
    NameImag = rand_file_name('.fits')
    NameResult = rand_file_name('')
    # Writes the input image to a fits file
    fits.writeto(NameImag,im)
    # Performing the system call
    if opt is None:
        subprocess.call(['mr_transform', NameImag, NameResult])
    else:
        subprocess.call(['mr_transform']+opt+[NameImag, NameResult])
    Result = fits.getdata(NameResult+'.mr')

    os.remove(NameImag)
    return swapaxes(swapaxes(Result,0,1),1,2),NameResult+'.mr'

def mr_trans_2(im,filters=None,opt=None):
    shap = im.shape
    if filters is None:
        dirac = zeros((shap[0]-(shap[0]-1)%2,shap[1]-(shap[1]-1)%2)) # Odd dimensions needed
        dirac[int((shap[0]-(shap[0]-1)%2-1)/2),int((shap[1]-(shap[1]-1)%2-1)/2)] = 1.
        filters,file = mr_trans(dirac,opt=opt)
        os.remove(file)
    shap_wav = filters.shape
    output = zeros((shap[0],shap[1],shap_wav[2]))
    for i in range(0,shap_wav[2]):
        output[:,:,i] = fftconvolve(im,filters[:,:,i],mode='same')
    return output,filters

def mr_trans_stack_2(stack,filters=None,opt=None):
    shap = stack.shape
    if filters is None:
        dirac = zeros((shap[0]-(shap[0]-1)%2,shap[1]-(shap[1]-1)%2)) # Odd dimensions needed
        dirac[int((shap[0]-(shap[0]-1)%2-1)/2),int((shap[1]-(shap[1]-1)%2-1)/2)] = 1
        filters,file = mr_trans(dirac,opt=opt)
        os.remove(file)
    shap_wav = filters.shape
    coeff_out = zeros((shap[0],shap[1],shap_wav[2],shap[2]))
    for i in range(0,shap[2]):
        coeff_i,filters = mr_trans_2(stack[:,:,i],filters=filters)
        coeff_out[:,:,:,i] = coeff_i
    return coeff_out,filters

def rot90_stack(cube):
    shap = cube.shape
    cube_out = copy(cube)
    for i in range(0,shap[2]):
        cube_out[:,:,i] = rot90(cube[:,:,i],2)
    return cube_out

def decim(im,d,av_en=1,fft=1):
    im_filt=copy(im)
    im_d = copy(im)
    if d>1:
        if av_en==1:
            siz = d+1-(d%2)
            mask = ones((siz,siz))/siz**2
            if fft==1:im_filt = fftconvolve(im, mask, mode='same')
            else:im_filt = convolve(im, mask, mode='same')
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

def transpose_decim(im,decim_fact,av_en=0):
    shap = im.shape
    im_out = zeros((shap[0]*decim_fact,shap[1]*decim_fact))
    for i in range(0,shap[0]):
        for j in range(0,shap[1]):
            im_out[decim_fact*i,decim_fact*j]=im[i,j]
    if av_en==1:
        siz = decim_fact+1-(decim_fact%2)
        mask = ones((siz,siz))/siz**2
        im_out = fftconvolve(im, mask, mode='same')
    return im_out

def acc_sig_map(shap_im,ker_stack,sig_est,flux_est,flux_ref,upfact,w,sig_data=None):
    shap = ker_stack.shape
    nb_im = shap[2]
    if sig_data is None:
        sig_data = ones((nb_im,))
    var_stack = ones((shap_im[0],shap_im[1],nb_im))
    map2 = zeros((shap_im[0]*upfact,shap_im[1]*upfact))
    ker_stack_in = copy(ker_stack)**2
    for l in range(0,shap[2]):
        var_stack[:,:,l]*=sig_data[l]**2
        map2 += ((w[l]*flux_est[l]/(sig_est[l]*flux_ref))**2)*convolve(\
        transpose_decim(var_stack[:,:,l],upfact),ker_stack_in[:,:,l],mode='same')
    map =  sqrt(map2)
    return map

def acc_sig_maps(shap_im,ker_stack,sig_est,flux_est,flux_ref,upfact,w,sig_data=None):
    shap = w.shape
    map_out = zeros((shap_im[0]*upfact,shap_im[1]*upfact,shap[0]))
    for i in range(0,shap[0]):
        map_out[:,:,i] = acc_sig_map(shap_im,ker_stack,sig_est,flux_est,flux_ref,\
        upfact,w[i,:],sig_data=sig_data)
    return map_out

def mat_to_cube(mat,n1,n2):
    shap = mat.shape
    cube = zeros((n1,n2,shap[0]))
    k=0
    for k in range(0,shap[0]):
        cube[:,:,k] = mat[k,:].reshape(n1,n2)
    return cube

def rect_crop_c(im,n1,n2,cent,dx=0,dy=0):
    nb_im = 1
    if len(im.shape)>2:nb_im = im.shape[2]
    im_crop = None
    if nb_im==1:
        im_crop = im[around(cent[0,0]+dx-n1/2).astype(int):around(cent[0,0]+dx-n1/2).astype(int)+n1,\
        around(cent[0,1]+dy-n2/2).astype(int):around(cent[0,1]+dy-n2/2).astype(int)+n2]
    else:
        im_crop = zeros((n1,n2,nb_im))
        for k in range(0,nb_im):
            imk = im[:,:,k]
            im_crop[:,:,k] = imk[around(cent[k,0]+dx-n1/2).astype(int):around(cent[k,0]+dx-n1/2).astype(int)\
            +n1,around(cent[k,1]+dy-n2/2).astype(int):around(cent[k,1]+dy-n2/2).astype(int)+n2]
    return im_crop

def cube_svd(cube,nb_comp=None,ind=None,mean_sub=False):
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

def log_sampling(val_min,val_max,nb_samp):
    lval_min = log(val_min)
    lval_max = log(val_max)
    a = double(array(range(0,nb_samp)))/(nb_samp-1)
    a = a*(lval_max-lval_min) + lval_min
    a = exp(a)
    return a

def feat_dist_mat(feat_mat):
    shap = feat_mat.shape
    mat_out = zeros((shap[0],shap[0]-1))
    a = array(range(0,shap[0]))
    for i in range(0,shap[0]):
        ind_i = where(a!=i)

        for k in range(0,shap[0]-1):
            mat_out[i,k] = sqrt(((feat_mat[i,:]-feat_mat[ind_i[0][k],:])**2).sum())
    return mat_out

def knn_interf(data,nb_neigh):
    flann = FLANN()
    data_cp = array(data, dtype=float64)
    params = flann.build_index(data_cp)
    result_temp, dists_temp = flann.nn_index(data_cp, nb_neigh+1)
    dists = dists_temp[:,1:nb_neigh+1]
    result = result_temp[:,1:nb_neigh+1]
    return result,dists

def lanczos(U,n=10,n2=None):
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

def shift_ker_stack(shifts,upfact,lanc_rad=4):
    shap = shifts.shape
    shift_ker_stack = zeros((2*lanc_rad+1,2*lanc_rad+1,shap[0]))
    shift_ker_stack_adj = zeros((2*lanc_rad+1,2*lanc_rad+1,shap[0]))
    for i in range(0,shap[0]):
        uin = shifts[i,:].reshape((1,2))*upfact
        shift_ker_stack[:,:,i] = lanczos(uin,n=lanc_rad)
        shift_ker_stack_adj[:,:,i] = rot90(shift_ker_stack[:,:,i],2)
    return shift_ker_stack,shift_ker_stack_adj

def flux_estimate(im,cent=None,rad=4): # Default value for the flux tunned for Euclid PSF at Euclid resolution
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
    shap = stack.shape
    flux = zeros((shap[2],))
    for i in range(0,shap[2]):
        if cent is not None:
            flux[i] = flux_estimate(stack[:,:,i],cent=cent[i,:],rad=rad)
        else:
            flux[i] = flux_estimate(stack[:,:,i],rad=rad)
    return flux

def compute_centroid(im,sigw=None,nb_iter=4):
    if sigw is None:
        param=gaussfitter.gaussfit(im,returnfitimage=False)
        sigw = (param[3]+param[4])/2
    sigw = float(sigw)
    n1 = im.shape[0]
    n2 = im.shape[1]
    rx = array(range(0,n1))
    ry = array(range(0,n2))
    Wc = ones((n1,n2))
    centroid = zeros((1,2))
    # Four iteration loop to compute the centroid
    for i in range(0,nb_iter):
        xx = outerproduct(rx-centroid[0,0],ones(n2))
        yy = outerproduct(ones(n1),ry-centroid[0,1])
        Wc = exp(-(xx**2+yy**2)/(2*sigw**2))
        centroid = zeros((1,2))
        # Estimate Centroid
        if i == 0:Wc = ones((n1,n2))
        totx=0.0
        toty=0.0
        for cx in range(0,n1):
            centroid[0,0] += (im[cx,:]*Wc[cx,:]).sum()*(cx)
            totx += (im[cx,:]*Wc[cx,:]).sum()
        for cy in range(0,n2):
            centroid[0,1] += (im[:,cy]*Wc[:,cy]).sum()*(cy)
            toty += (im[:,cy]*Wc[:,cy]).sum()
        centroid = centroid*array([1/totx,1/toty])
    return centroid

def shift_est(psf_stack): 
    shap = psf_stack.shape
    U = zeros((shap[2],2))
    centroid_out = zeros((shap[2],2))
    for i in range(0,shap[2]):
        param=gaussfitter.gaussfit(psf_stack[:,:,i],returnfitimage=False)
        centroid = compute_centroid(psf_stack[:,:,i],(param[3]+param[4])/2)
        U[i,0] = centroid[0,0]-double(shap[0])/2
        U[i,1] = centroid[0,1]-double(shap[1])/2
        centroid_out[i,0]  = centroid[0,0]
        centroid_out[i,1]  = centroid[0,1]
    return U,centroid_out

def thresholding(x,thresh,thresh_type): # x is a 1D or 2D array, thresh is an array of the same size, thresh_type is 1 or 0, 1 for soft thresholding, 0 for hard thresholding
    xthresh = copy(x)
    n = x.shape
    if len(n)>0:
        n1 = n[0]
    else:
        n1=1
    n2=1
    if len(n)==2:n2 =n[1]
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

def thresholding_3D(x,thresh,thresh_type):
    nb_plan = x.shape[2]
    xthresh = copy(x)
    for k in range(0,nb_plan):
        xthresh[:,:,k] = thresholding(copy(x[:,:,k]),thresh[:,:,k],thresh_type)
    return xthresh

def lineskthresholding(mat,k):
    mat_out = copy(mat)
    shap = mat.shape
    for j in range(0,shap[0]):
        mat_out[j,:] = kthresholding(mat[j,:],k)
    return mat_out

def kthresholding(x,k):
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

def rand_file_name(ext):
    current_time = datetime.datetime.now().time()
    return 'file'+str(current_time)+ext

def mad(x):
    return median(abs(x-median(x)))

def im_gauss_nois_est(im,opt=['-t2','-n2'],filters=None):
    Result,filters = mr_trans_2(im,filters=filters,opt=opt)
    norm_wav = norm(filters[:,:,0])
    sigma = 1.4826*mad(Result[:,:,0])/norm_wav
    return sigma,filters

def im_gauss_nois_est_cube(cube,opt=None,filters=None):
    shap = cube.shape
    sig = zeros((shap[2],))
    for i in range(0,shap[2]):
        sig_i,filters = im_gauss_nois_est(cube[:,:,i],opt=opt,filters=filters)
        sig[i] = sig_i
    return sig

def kernel_ext(mat,tol = 0.01): # The kernel of the matrix mat is defined as the vector space spanned by the eigenvectors corresponding to 1% of the sum of the squared singular values
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
