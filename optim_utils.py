from numpy import *
import numpy as np
import sys
sys.path.append('../utilities')
import utils
from pyflann import *
import psf_learning_utils
from modopt.opt.cost import costObj
import grads as grad
import operators as lambdaops
import modopt.opt.proximity as prox
import proxs as lambdaprox
import modopt.opt.algorithms as optimalg
from modopt.opt.linear import Identity

def dist_map_2(arr): # The samples are in the columns
    """Pairwise distances...?"""
    from numpy.linalg import norm
    from numpy import ones,transpose,fill_diagonal

    nb_samp = arr.shape[1]
    norm_mat = ((norm(arr,axis=0).reshape((nb_samp,1)))**2).dot(ones((1,nb_samp)))
    scalar_prod = transpose(arr).dot(arr)
    res = norm_mat+transpose(norm_mat)-2*scalar_prod
    i,j = where(res<0)
    res[i,j] = 0
    return res

def non_uniform_smoothing_bas_mat(weights):
    """ **[???]**
    
    """
    shap = weights.shape
    A = zeros((shap[0],shap[0]))
    B = zeros((shap[0],shap[0]))

    range_ref = array(range(0,shap[0]))

    for i in range(0,shap[0]):
        ik = where(range_ref!=i)
        A[ik[0],i] = -weights[i,:]
        B[i,i] = sum(weights[i,:])

    return A,B

def non_uniform_smoothing_mat_2(weights,e): # e>=0
    """ **[???]**
    
    Calls:
    
    * :func:`optim_utils.non_uniform_smoothing_bas_mat`
    """
    A,B = non_uniform_smoothing_bas_mat(weights)
    mat_out = A.dot(transpose(A))+ e*(A.dot(B)+ B.dot(transpose(A))) + (e**2)*B.dot(B) # B is diagonal matrix
    return mat_out

def non_uniform_smoothing_mat_dist_1(dist,expo_range,e):
    """ **[???] Computes some distance matrix, *again*.**
    
    Calls:
    
    * :func:`optim_utils.non_uniform_smoothing_mat_2`
    """
    dist_med = np.median(dist)
    nb_samp = len(expo_range)
    nb_im = dist.shape[0]
    mat_stack = zeros((nb_im,nb_im,nb_samp))
    for i in range(0,nb_samp):
        dist_weights = (dist_med/dist)**expo_range[i]
        dist_weigths = dist_weights/dist_weights.max()
        mat_stack[:,:,i] = non_uniform_smoothing_mat_2(dist_weights,e)
    return mat_stack

def non_uniform_smoothing_mat_dist_2(dist,expo,e_range):
    """ Same as :func:`optim_utils.non_uniform_smoothing_mat_dist_1`, but with
    constant ``expo`` and varying ``e`` (as opposed to constant ``e`` and 
    varying ``expo``.
    
    Calls:
    
    * :func:`optim_utils.non_uniform_smoothing_mat_2`
    """
    dist_med = np.median(dist)
    nb_samp = len(e_range)
    nb_im = dist.shape[0]
    mat_stack = zeros((nb_im,nb_im,nb_samp))
    for i in range(0,nb_samp):
        dist_weights = (dist_med/dist)**expo
        dist_weigths = dist_weights/dist_weights.max()
        mat_stack[:,:,i] = non_uniform_smoothing_mat_2(dist_weights,e_range[i])
    return mat_stack

def notch_filt_optim_2(test_mat,dist,expo_range,e_range,nb_iter=2,tol=0.01):
    """**[???]** Finds notch filter hyperparameters I guess?
    
    Calls:
    
    * :func:`optim_utils.non_uniform_smoothing_mat_dist_1`
    * :func:`utils.kernel_mat_stack_test_unit`
    * :func:`optim_utils.non_uniform_smoothing_mat_dist_2`
    """
    expo_out = None
    e_out = 0.5
    loss = None
    vect = None
    ker = None
    j2 = None
    for i in range(0,nb_iter):
        mat_stack = non_uniform_smoothing_mat_dist_1(dist,expo_range,e_out)
        vect,j,loss,ker,j2 = utils.kernel_mat_stack_test_unit(mat_stack,test_mat,tol=tol)
        print "=========== Loss out ==========: ",loss
        expo_out = expo_range[j]
        mat_stack = non_uniform_smoothing_mat_dist_2(dist,expo_out,e_range)
        vect,j,loss,ker,j2 = utils.kernel_mat_stack_test_unit(mat_stack,test_mat,tol=tol)
        print "=========== Loss out ==========: ",loss
        e_out = e_range[j]

    return expo_out,e_out,loss,vect,ker,j2

def analysis(cube,sig,field_dist,p_min = 0.01,e_min=0.01,e_max=1.99,nb_max=30,tol=0.01):
    """Computes graph-constraint related values, see RCA paper sections 5.2 and (especially) 5.5.3.
    
    
    Calls:
    
    * :func:`utils.knn_interf`
    * :func:`optim_utils.pow_law_select`
    * :func:`utils.feat_dist_mat`
    * :func:`utils.log_sampling`
    * :func:`optim_utils.notch_filt_optim_2`
    * :func:`utils.mat_to_cube`
    """
    nb_samp_opt = 10
    shap = cube.shape
    nb_neighs = shap[2]-1
    neigh,dists = utils.knn_interf(field_dist,nb_neighs)
    p_max = pow_law_select(dists,nb_neighs)
    print "power max = ",p_max

    print "Done..."
    dists_unsorted = utils.feat_dist_mat(field_dist)
    e_range = utils.log_sampling(e_min,e_max,nb_samp_opt)
    p_range = utils.log_sampling(p_min,p_max,nb_samp_opt)
    res_mat = copy(transpose(cube.reshape((shap[0]*shap[1],shap[2]))))

    list_comp = list()
    list_e = list()
    list_p = list()
    list_ind = list()
    list_ker = list()
    err = 1e20
    nb_iter = 0
    while nb_iter<nb_max: # err>sig and
        expo_out,e_out,loss,vect,ker,j = notch_filt_optim_2(res_mat,dists_unsorted,p_range,e_range,nb_iter=3,tol=tol)
        list_e.append(e_out)
        list_p.append(expo_out)
        list_comp.append(vect)
        list_ind.append(j)
        list_ker.append(ker)
        nb_iter+=1
        res_mat = res_mat-transpose(vect).dot(vect.dot(res_mat))
        print "nb_comp: ",nb_iter," residual: ",loss," e: ",e_out," p: ",expo_out,"chosen index: ",j,"/",shap[2]
        err = sum(res_mat**2)

    e_vect = zeros((nb_iter,))
    p_vect = zeros((nb_iter,))
    weights = zeros((nb_iter,shap[2]))
    ker = zeros((nb_iter*shap[2],shap[2]))
    ind = zeros((nb_iter,nb_iter*shap[2]))
    for i in range(0,nb_iter):
        e_vect[i] = list_e[i]
        p_vect[i] = list_p[i]
        weights[i,:] = list_comp[i].reshape((shap[2],))
        ker[i*shap[2]:(i+1)*shap[2],:] = list_ker[i]
        ind[i,i*shap[2]+list_ind[i]] = 1


    res_mat = copy(transpose(cube.reshape((shap[0]*shap[1],shap[2]))))
    proj_coeff = weights.dot(res_mat)
    comp = utils.mat_to_cube(proj_coeff,shap[0],shap[1])

    proj_data = transpose(weights).dot(proj_coeff)
    proj_data = utils.mat_to_cube(proj_data,shap[0],shap[1])

    return e_vect,p_vect,weights,comp,proj_data,ker,ind

def pow_law_select(dist_weights,nb_neigh,min_val=10**(-15)):
    """ **[???] but related to proximity constrains hyperparameters**
    """

    a = dist_weights[:,0]/dist_weights[:,nb_neigh-1]
    r_med = a.min()
    print "r_med: ",r_med,nb_neigh
    p = log(min_val)/log(r_med)
    return p

def polychromatic_psf_field_est_2(im_stack_in,spectrums,wvl,D,opt_shift_est,nb_comp,field_pos=None,nb_iter=4,nb_subiter=100,mu=0.3,\
                        tol = 0.1,sig_supp = 3,sig=None,shifts=None,flux=None,nsig_shift_est=4,pos_en = True,simplex_en=False,\
                        wvl_en=True,wvl_opt=None,nsig=3,graph_cons_en=False):
    """ Main LambdaRCA function.
    
    Calls:
    
    * :func:`utils.get_noise_arr`
    * :func:`utils.diagonally_dominated_mat_stack` 
    * :func:`psf_learning_utils.full_displacement` 
    * :func:`utils.im_gauss_nois_est_cube` 
    * :func:`utils.thresholding_3D` 
    * :func:`utils.shift_est` 
    * :func:`utils.shift_ker_stack` 
    * :func:`utils.flux_estimate_stack` 
    * :func:`optim_utils.analysis` 
    * :func:`utils.cube_svd`
    * :func:`grads.polychrom_eigen_psf`
    * :func:`grads.polychrom_eigen_psf_coeff_graph`
    * :func:`grads.polychrom_eigen_psf_coeff`
    * :func:`psf_learning_utils.field_reconstruction`
    * :func:`operators.transport_plan_lin_comb_wavelet`
    * :func:`operators.transport_plan_marg_wavelet`
    * :func:`operators.transport_plan_lin_comb`
    * :func:`operators.transport_plan_lin_comb_coeff`
    * :func:`proxs.simplex_threshold`
    * :func:`proxs.Simplex`
    * :func:`proxs.KThreshold`
    """

    im_stack = copy(im_stack_in)
    if wvl_en:
        from utils import get_noise_arr

    print "--------------- Transport architecture setting ------------------"
    nb_im = im_stack.shape[-1]
    shap_obs = im_stack.shape
    shap = (shap_obs[0]*D,shap_obs[1]*D)
    P_stack = utils.diagonally_dominated_mat_stack(shap,nb_comp,sig=sig_supp,thresh_en=True)
    i,j = where(P_stack[:,:,0]>0)
    supp = transpose(array([i,j]))
    t = (wvl-wvl.min()).astype(float)/(wvl.max()-wvl.min())

    neighbors_graph,weights_neighbors,cent,coord_map,knn = psf_learning_utils.full_displacement(shap,supp,t,\
    pol_en=True,cent=None,theta_param=1,pol_mod=True,coord_map=None,knn=None)

    print "------------------- Forward operator parameters estimation ------------------------"
    centroids = None
    if sig is None:
        sig,filters = utils.im_gauss_nois_est_cube(copy(im_stack),opt=opt_shift_est)

    if shifts is None:
        map = ones(im_stack.shape)
        for i in range(0,shap_obs[2]):
            map[:,:,i] *= nsig_shift_est*sig[i]
        print 'Shifts estimation...'
        psf_stack_shift = utils.thresholding_3D(copy(im_stack),map,0)
        shifts,centroids = utils.shift_est(psf_stack_shift)
        print 'Done...'
    else:
        print "---------- /!\ Warning: shifts provided /!\ ---------"
    ker,ker_rot = utils.shift_ker_stack(shifts,D)
    sig /=sig.min()
    for k in range(0,shap_obs[2]):
        im_stack[:,:,k] = im_stack[:,:,k]/sig[k]
    print " ------ ref energy: ",(im_stack**2).sum()," ------- "
    if flux is None:
        flux = utils.flux_estimate_stack(copy(im_stack),rad=4)

    if graph_cons_en:
        print "-------------------- Spatial constraint setting -----------------------"
        e_opt,p_opt,weights,comp_temp,data,basis,alph  = analysis(im_stack,0.1*prod(shap_obs)*sig.min()**2,field_pos,nb_max=nb_comp)

    print "------------- Coeff init ------------"
    A,comp,cube_est = utils.cube_svd(im_stack,nb_comp=nb_comp)

    i=0
    print " --------- Optimization instances setting ---------- "

    # Data fidelity related instances
    polychrom_grad = grad.polychrom_eigen_psf(im_stack, supp, neighbors_graph, \
                weights_neighbors, spectrums, A, flux, sig, ker, ker_rot, D)

    if graph_cons_en:
        polychrom_grad_coeff = grad.polychrom_eigen_psf_coeff_graph(im_stack, supp, neighbors_graph, \
                weights_neighbors, spectrums, P_stack, flux, sig, ker, ker_rot, D, basis)
    else:
        polychrom_grad_coeff = grad.polychrom_eigen_psf_coeff(im_stack, supp, neighbors_graph, \
                weights_neighbors, spectrums, P_stack, flux, sig, ker, ker_rot, D)


    # Dual variable related linear operators instances
    dual_var_coeff = zeros((supp.shape[0],nb_im))
    if wvl_en and pos_en:
        lin_com = lambdaops.transport_plan_lin_comb_wavelet(A,supp,weights_neighbors,neighbors_graph,shap,wavelet_opt=wvl_opt)
    else:
        if wvl_en:
            lin_com = lambdaops.transport_plan_marg_wavelet(supp,weights_neighbors,neighbors_graph,shap,wavelet_opt=wvl_opt)
        else:
            lin_com = lambdaops.transport_plan_lin_comb(A, supp,shap)

    if not graph_cons_en:
        lin_com_coeff = lambdaops.transport_plan_lin_comb_coeff(P_stack, supp)

    # Proximity operators related instances
    id_prox = Identity()
    if wvl_en and pos_en:
        noise_map = get_noise_arr(lin_com.op(polychrom_grad.MtX(im_stack))[1])
        dual_var_plan = np.array([zeros((supp.shape[0],nb_im)),zeros(noise_map.shape)])
        dual_prox_plan = lambdaprox.simplex_threshold(lin_com, nsig*noise_map,pos_en=(not simplex_en))
    else:
        if wvl_en:
            # Noise estimation
            noise_map = get_noise_arr(lin_com.op(polychrom_grad.MtX(im_stack)))
            dual_var_plan = zeros(noise_map.shape)
            dual_prox_plan = prox.SparseThreshold(lin_com, nsig*noise_map)
        else:
            dual_var_plan = zeros((supp.shape[0],nb_im))
            if simplex_en:
                dual_prox_plan = lambdaprox.Simplex()
            else:
                dual_prox_plan = prox.Positivity()

    if graph_cons_en:
        iter_func = lambda x: floor(sqrt(x))
        prox_coeff = lambdaprox.KThreshold(iter_func)
    else:
        if simplex_en:
            dual_prox_coeff = lambdaprox.Simplex()
        else:
            dual_prox_coeff = prox.Positivity()

    # ---- (Re)Setting hyperparameters
    delta  = (polychrom_grad.inv_spec_rad**(-1)/2)**2 + 4*lin_com.mat_norm**2
    w = 0.9
    sigma_P = w*(np.sqrt(delta)-polychrom_grad.inv_spec_rad**(-1)/2)/(2*lin_com.mat_norm**2)
    tau_P = sigma_P
    rho_P = 1

    # Cost function instance
    cost_op = costObj([polychrom_grad])

    condat_min = optimalg.Condat(P_stack, dual_var_plan, polychrom_grad, id_prox, dual_prox_plan, lin_com, cost=cost_op,\
                 rho=rho_P,  sigma=sigma_P, tau=tau_P, rho_update=None, sigma_update=None,
                 tau_update=None, auto_iterate=False)
    print "------------------- Transport plans estimation ------------------"

    condat_min.iterate(max_iter=nb_subiter) # ! actually runs optimisation
    P_stack = condat_min.x_final
    dual_var_plan = condat_min.y_final

    obs_est = polychrom_grad.MX(P_stack)
    res = im_stack - obs_est

    for i in range(0,nb_iter):
        print "----------------Iter ",i+1,"/",nb_iter,"-------------------"

        # Parameters update
        polychrom_grad_coeff.set_P(P_stack)
        if not graph_cons_en:
            lin_com_coeff.set_P_stack(P_stack)
            # ---- (Re)Setting hyperparameters
            delta  = (polychrom_grad_coeff.inv_spec_rad**(-1)/2)**2 + 4*lin_com_coeff.mat_norm**2
            w = 0.9
            sigma_coeff = w*(np.sqrt(delta)-polychrom_grad_coeff.inv_spec_rad**(-1)/2)/(2*lin_com_coeff.mat_norm**2)
            tau_coeff = sigma_coeff
            rho_coeff = 1

        # Coefficients cost function instance
        cost_op_coeff = costObj([polychrom_grad_coeff])

        if graph_cons_en:
            beta_param = polychrom_grad_coeff.inv_spec_rad# set stepsize to inverse spectral radius of coefficient gradient
            min_coeff = optimalg.ForwardBackward(alph, polychrom_grad_coeff, prox_coeff, beta_param=beta_param, 
                                                 cost=cost_op_coeff,auto_iterate=False)
        else:
            min_coeff = optimalg.Condat(A, dual_var_coeff, polychrom_grad_coeff, id_prox, dual_prox_coeff, lin_com_coeff, cost=cost_op_coeff,\
                                            rho=rho_coeff,  sigma=sigma_coeff, tau=tau_coeff, rho_update=None, sigma_update=None,\
                                            tau_update=None, auto_iterate=False)

        print "------------------- Coefficients estimation ----------------------"
        min_coeff.iterate(max_iter=nb_subiter) # ! actually runs optimisation
        if graph_cons_en:
            prox_coeff.reset_iter()
            alph = min_coeff.x_final
            A = alph.dot(basis)
        else:
            A = min_coeff.x_final
            dual_var_coeff = min_coeff.y_final

        # Parameters update
        polychrom_grad.set_A(A)
        if not wvl_en:
            lin_com.set_A(A)
        if wvl_en:
            # Noise estimate update
            noise_map = get_noise_arr(lin_com.op(polychrom_grad.MtX(im_stack))[1])
            dual_prox_plan.update_weights(noise_map)

        # ---- (Re)Setting hyperparameters
        delta  = (polychrom_grad.inv_spec_rad**(-1)/2)**2 + 4*lin_com.mat_norm**2
        w = 0.9
        sigma_P = w*(np.sqrt(delta)-polychrom_grad.inv_spec_rad**(-1)/2)/(2*lin_com.mat_norm**2)
        tau_P = sigma_P
        rho_P = 1

        # Cost function instance
        condat_min = optimalg.Condat(P_stack, dual_var_plan, polychrom_grad, id_prox, dual_prox_plan, lin_com, cost=cost_op,\
                     rho=rho_P,  sigma=sigma_P, tau=tau_P, rho_update=None, sigma_update=None,
                     tau_update=None, auto_iterate=False)
        print "------------------- Transport plans estimation ------------------"

        condat_min.iterate(max_iter=nb_subiter) # ! actually runs optimisation
        P_stack = condat_min.x_final
        dual_var_plan = condat_min.y_final

        # Normalization
        for j in range(0,nb_comp):
            l1_P = sum(abs(P_stack[:,:,j]))
            P_stack[:,:,j]/= l1_P
            A[j,:] *= l1_P
            if graph_cons_en:
                alph[j,:] *= l1_P
        polychrom_grad.set_A(A)
        # Flux update
        obs_est = polychrom_grad.MX(P_stack)
        err_ref = 0.5*sum((obs_est-im_stack)**2)
        flux_new = (obs_est*im_stack).sum(axis=(0,1))/(obs_est**2).sum(axis=(0,1))
        print "Flux correction: ",flux_new
        polychrom_grad.set_flux(polychrom_grad.get_flux()*flux_new)
        polychrom_grad_coeff.set_flux(polychrom_grad_coeff.get_flux()*flux_new)

        obs_est = polychrom_grad.MX(P_stack)
        res = im_stack - obs_est
        err_rec = 0.5*sum(res**2)
        print "err_ref : ",err_ref," ; err_rec : ", err_rec
        # Computing residual


    psf_est = psf_learning_utils.field_reconstruction(P_stack,shap,supp,neighbors_graph,weights_neighbors,A)

    return psf_est,P_stack,A,res

