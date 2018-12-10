import subprocess
import os
from astropy.io import fits
import sys
sys.path.append('../utilities')
import utils
from numpy import *
import scipy.signal as scisig
try:
    import pyct
except ImportError, e:
    pass # module doesn't exist, deal with it.
import numpy

def mr_filter(imag, opt=None):
    """Filters an image by using a multiresolution transform.

        Parameters
        ----------
        imag : 2darray
        Image to filter

        opt : string, optional
        Command line arguments to pass to the executable.

        Returns
        -------
        result : 2darray
        Filtered image

        Examples
        --------
        Filter an image with all default options:
        >>> Result = mr_filter(Imag)

        Same example, but impose the number of scales to be 3, and
        a thresolding at 5 sigma
        >>> Result = mr_filter(Imag, OPT=['-n 3','-s 5'])
        """

    # Declare names for the temporary fits file storing
    # the input and output images
    NameImag = 'xx_imag.fits'
    NameResult = 'xx_result.fits'

    # Writes the input image to a fits file
    fits.writeto(NameImag, imag)

    # Performing the system call
    if opt is None:
        subprocess.call(['mr_filter','-n3 -s5', NameImag, NameResult])
    else:
        subprocess.call(['mr_filter']+opt+[NameImag, NameResult])

    # Recovering the output of mr_filter
    Result = fits.getdata(NameResult)

    # Erasing temporary files
    os.remove(NameImag)
    os.remove(NameResult)

    # Returning the processed image
    return Result

def mr_trans(im,opt=None,path='',exe_path=''):
    NameImag = path+utils.rand_file_name('.fits')
    NameResult = path+utils.rand_file_name('')
    # Writes the input image to a fits file
    fits.writeto(NameImag,im)
    # Performing the system call
    mr_exe = exe_path+'mr_transform'
    if opt is None:
        subprocess.call([mr_exe, NameImag, NameResult])
    else:
        subprocess.call([mr_exe]+opt+[NameImag, NameResult])
    Result = fits.getdata(NameResult+'.mr')

    os.remove(NameImag)
    return swapaxes(swapaxes(Result,0,1),1,2),NameResult+'.mr'

def mr_trans_1d(x,opt=None,path='../data/',exe_path='',info_filename=None,save_info_en=True):
    NameSignal = path+utils.rand_file_name('.fits')
    NameResult = path+utils.rand_file_name('')
    # Writes the input signal to a fits file
    fits.writeto(NameSignal,x)

    # Performing the system call
    mr_exe = exe_path+'mr1d_trans'
    if opt is None:
        subprocess.call([mr_exe, NameSignal, NameResult])
    else:
        if save_info_en:
            if info_filename is None:
                info_filename = path+utils.rand_file_name('.mr')
            opt.append('-w '+info_filename)
        subprocess.call([mr_exe]+opt+[NameSignal, NameResult])
    Result = fits.getdata(NameResult+'.fits')
    os.remove(NameSignal)
    os.remove(NameResult+'.fits')
    return Result,info_filename

def mr_recons_1d(x,info_filename,path='../data/',exe_path=''):
    trans_file = path+utils.rand_file_name('.fits')
    # Writes the input coefficients to a fits file
    fits.writeto(trans_file,x)
    NameResult = path+utils.rand_file_name('')
    # Performing the system call
    mr_exe = exe_path+'mr1d_recons'
    subprocess.call([mr_exe, trans_file, info_filename,NameResult])
    Result = fits.getdata(NameResult+'.fits')
    os.remove(NameResult+'.fits')
    os.remove(trans_file)

    return Result

def mr_trans_2(im,filters=None,opt=None,exe_path=''):
    shap = im.shape
    if filters is None:
        dirac = zeros((shap[0]-(shap[0]-1)%2,shap[1]-(shap[1]-1)%2)) # Odd dimensions needed
        dirac[int((shap[0]-(shap[0]-1)%2-1)/2),int((shap[1]-(shap[1]-1)%2-1)/2)] = 1
        filters,file = mr_trans(dirac,opt=opt,exe_path=exe_path)
        os.remove(file)
    shap_wav = filters.shape
    output = zeros((shap[0],shap[1],shap_wav[2]))

    for i in range(0,shap_wav[2]):
        temp = scisig.fftconvolve(im,filters[:,:,i],mode='same')
        output[:,:,i]  =temp
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

def mr_transf_transp(coeff,filters_rot):
    coeff_temp = copy(coeff)
    shap = coeff.shape
    for i in range(0,shap[2]):
        coeff_temp[:,:,i] = scisig.fftconvolve(coeff[:,:,i],filters_rot[:,:,i],mode='same')
    out = coeff_temp.sum(axis=2)

    return out

def mr_transf_transp_stack(coeff,filters_rot):

    shap = coeff.shape
    output = zeros((shap[0],shap[1],shap[3]))

    for i in range(0,shap[3]):
        output[:,:,i] = mr_transf_transp(coeff[:,:,:,i],filters_rot)

    return output

def mr_trans_stack(stack,opt=None,path='../data/',clean_en=False):
    shap = stack.shape
    coeff_list  = list()
    mr_list = list()
    for i in range(0,shap[2]):
        Result,File = mr_trans(stack[:,:,i],opt=opt,path=path)
        coeff_list.append(Result)
        if clean_en:
            os.remove(File)
        else:
            mr_list.append(File)

    if clean_en:
        return coeff_list
    else:
        return coeff_list,mr_list

def mr_read(mr_filename,exe_path='../../CPP/sprite/build/'):
    NameImag = '../data/'+utils.rand_file_name('.fits')
    subprocess.call([exe_path+'mr_read_write', mr_filename, NameImag,'0'])
    Result = fits.getdata(NameImag)
    os.remove(NameImag)
    return Result

def mr_write(fits_filename,mr_filename,exe_path='../../CPP/sprite/build/'):
    subprocess.call([exe_path+'mr_read_write',fits_filename,mr_filename,'1'])

def mr_recons(mr_file):
    NameImag = '../data/'+utils.rand_file_name('.fits')
    subprocess.call(['mr_recons', mr_file,NameImag])
    Result = fits.getdata(NameImag)
    os.remove(NameImag)

    return Result
def mr_recons_stack(mr_files,shap):
    output = zeros(shap)
    for i in range(0,len(mr_files)):
        output[:,:,i] =  mr_recons(mr_files[i])

    return output


def mr_recons_coeff(coeff,mr_filename,exe_path='../../CPP/sprite/bin/'):
    NameImag = '../data/'+utils.rand_file_name('.fits')
    fits.writeto(NameImag,coeff)
    mr_write(NameImag,mr_filename,exe_path)
    Result = mr_recons(mr_filename)
    os.remove(NameImag)
    return Result

def mr1d_trans(signal,opt=None):
    """Apply a multiresolution transform to a signal.

        Parameters
        ----------
        imag : 1darray
        Signal transform

        opt : string, optional
        Command line arguments to pass to the executable.

        Returns
        -------
        result : 1darray
        Transformed signal

        Examples
        --------
        Filter a signal with all default options:
        >>> Result = mr1d_filter(signal)

        Same example, but impose the number of scales to be 3, and
        a thresolding at 5 sigma
        >>> Result = mr_filter(Imag, OPT='-n 3 -s 5')
        """

    # Declare names for the temporary fits file storing
    # the input and output images
    NameSignal = 'xx_imag.fits'
    NameResult = 'xx_result.fits'

    # Writes the input image to a fits file
    fits.writeto(NameSignal, signal)

    # Performing the system call
    if opt is None:
        subprocess.call(['mr1d_trans', NameSignal, NameResult])
    else:
        subprocess.call(['mr1d_trans', opt, NameSignal, NameResult])

    # Recovering the output of mr_filter
    Result = fits.getdata(NameResult)

    # Erasing temporary files
    os.remove(NameSignal)
    os.remove(NameResult)

    # Returning the processed image
    return Result

def mr1d_trans(signal,opt=None):
    """Apply a multiresolution transform to a signal.

        Parameters
        ----------
        imag : 1darray
        Signal transform

        opt : string, optional
        Command line arguments to pass to the executable.

        Returns
        -------
        result : 1darray
        Transformed signal

        Examples
        --------
        Filter a signal with all default options:
        >>> Result = mr1d_filter(signal)

        Same example, but impose the number of scales to be 3, and
        a thresolding at 5 sigma
        >>> Result = mr_filter(Imag, OPT='-n 3 -s 5')
        """

    # Declare names for the temporary fits file storing
    # the input and output images
    NameSignal = 'xx_imag.fits'
    NameResult = 'xx_result.fits'

    # Writes the input image to a fits file
    fits.writeto(NameSignal, signal)

    # Performing the system call
    if opt is None:
        subprocess.call(['mr1d_trans', NameSignal, NameResult])
    else:
        subprocess.call(['mr1d_trans', opt, NameSignal, NameResult])

    # Recovering the output of mr_filter
    Result = fits.getdata(NameResult)

    # Erasing temporary files
    os.remove(NameSignal)
    os.remove(NameResult)

    # Returning the processed image
    return Result


def mr1d_recons(wav_coeff, info_filename):
    """Reconstruct a signal from multiresolution transform coefficients.

        Parameters
        ----------"""

    # Declare names for the temporary fits file storing
    # the input and output images
    NameSignal = 'xx_imag.fits'
    NameResult = 'xx_result.fits'

    # Writes the input image to a fits file
    fits.writeto(NameSignal, wav_coeff)

    # Performing the system call
    subprocess.call(['mr1d_recons', NameSignal,info_filename,NameResult])

    # Recovering the output of mr_filter
    Result = fits.getdata(NameResult)

    # Erasing temporary files
    os.remove(NameSignal)
    os.remove(NameResult)

    # Returning the processed image
    return Result

def mr1d_filter_beta(signal, info_filename,nsigma,thresh_type,opt=None,sig=1):
    """Filters a signal by using a multiresolution transform. The transform is assumed to be undecimated.

        Parameters
        ----------
        imag : 1darray
        Signal to filter

        opt : string, optional
        Command line arguments to pass to the executable.

        Returns
        -------
        result : 1darray
        Filtered signal

        Examples
        --------
        Filter a signal with all default options:
        >>> Result = mr1d_filter(signal)

        Same example, but impose the number of scales to be 3, and
        a thresolding at 5 sigma
        >>> Result = mr_filter(Imag, OPT='-n 3 -s 5')
        """

    # Declare names for the temporary fits file storing
    # the input and output images
    NameSignal = 'xx_imag.fits'
    NameResult = 'xx_result.fits'
    print nsigma
    # Writes the input image to a fits file
    fits.writeto(NameSignal, signal)

    # Performing the system call
    if opt is None:
        opt = "-w "+info_filename
    else:
        opt = "-w "+info_filename+" "+opt
    subprocess.call(['mr1d_trans', opt, NameSignal, NameResult])


    # Recovering the output of mr_filter
    sig_trans = fits.getdata(NameResult)
    N = sig_trans.shape[0]
    l = sig_trans.shape[1]
    i=0
    dirac = zeros((l,))
    dirac[int(l/2 -1)]=1
    # Erasing temporary files
    os.remove(NameSignal)
    os.remove(NameResult)
    fits.writeto(NameSignal, dirac)
    subprocess.call(['mr1d_trans', opt, NameSignal, NameResult])
    filt = fits.getdata(NameResult)

    for i in range(0,N-1):
        scale_i = sig_trans[i,:]
        filt_i = filt[i,:]
        norm_i = sqrt((filt_i**2).sum())
        thresh_i = ones((l,))*norm_i*sig*nsigma
        sig_trans[i,:] = util.thresholding(scale_i,thresh_i,thresh_type)

    # Erasing temporary files
    os.remove(NameSignal)
    os.remove(NameResult)

    Result = mr1d_recons(sig_trans,info_filename)

    os.remove(info_filename)

    # Returning the processed image
    return Result



def stack_pyct_fwd(im,curv_obj=None,nb_sc=4,nb_dir=8,corr_en=False):
    shap_im = im.shape
    if curv_obj is None:
        curv_obj = pyct.fdct2([shap_im[0],shap_im[1]],nb_sc,nb_dir,True,norm=True)
    corr_coeff=None # MT = corr_coeff*M^-1
    if corr_en:
        a = numpy.random.randn(shap_im[0],shap_im[1])
        b = curv_obj.fwd(a)
        corr_coeff = (b*numpy.conj(b)).sum()/(a**2).sum()

    curv1 = curv_obj.fwd(im[:,:,0])
    nb_pts = len(curv1)
    curv_out = zeros((shap_im[2],nb_pts))
    curv_out[0,:] = copy(curv1)
    for i in range(1,shap_im[2]):
        curv_out[i,:] = curv_obj.fwd(im[:,:,i])
    if corr_en:
        return curv_out,curv_obj,corr_coeff
    else:
        return curv_out,curv_obj
def stack_pyct_inv(trans,curv_obj=None,nb_sc=4,nb_dir=8,shap = None):
    shap_trans = trans.shape
    if curv_obj is None:
        curv_obj = pyct.fdct2([shap[0],shap[1]],nb_sc,nb_dir,True,norm=True)

    im1 = curv_obj.inv(trans[0,:])
    if shap is None:
        shap = im1.shape
    im_out = zeros((shap[0],shap[1],shap_trans[0]))
    im_out[:,:,0] = im1
    for i in range(1,shap_trans[0]):
        im_out[:,:,i] = curv_obj.inv(trans[i,:])

    return im_out
