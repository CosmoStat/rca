from numpy import zeros,size,where,ones,copy,around,double,sinc,random,pi,arange,\
cos,sin,arccos,transpose,diag,sqrt,arange,floor,exp,array,mean,roots,float64,int,\
pi,median,rot90,argsort,tile,repeat,squeeze,log,int
from numpy.linalg import svd,norm,inv,eigh
from scipy.signal import fftconvolve,convolve
import utils_for_rca
import copy as cp
import os

def pos_proj(z,tol=0): # Puts negative entries of z to zero
    u = copy(z)
    u = u.reshape(size(z))
    i = where(u<0)
    u[i[0]] = 0
    shap = z.shape
    u = u.reshape(shap)
    return u

def pos_proj_mat(m,tol=0):
    u = copy(m)
    shap = m.shape
    j=0
    for j in range(0,shap[1]):
        u[:,j]=pos_proj(squeeze(m[:,j]),tol=tol)
    return u

def pos_proj_mat_2(m1,m2):
    m = m1+m2
    I1 = (m>=0)
    I2 = (m<0)
    m1_proj = m1*I1 + I2*(m1-m2)/2
    m2_proj = m2*I1 + I2*(m2-m1)/2
    return m1_proj,m2_proj

def pos_proj_cube(m,tol=0):
    u = copy(m)
    shap = m.shape
    k=0
    for k in range(0,shap[2]):
        u[:,:,k]=pos_proj_mat(squeeze(m[:,:,k]),tol=0)
    return u

def sr_stack_op(input):
    S = input[0]
    upfact = input[1]
    ker = input[2]
    ker_adj = input[3]
    shap = ker_adj.shape
    A = input[4]
    sig = input[5]
    flux = input[6]
    flux_ref = median(flux)
    nb_im = shap[2]
    shapS = S.shape
    output = zeros((shapS[0]/upfact,shapS[1]/upfact,nb_im))
    for i in range(0,nb_im):
        im_i = zeros((shapS[0],shapS[1]))
        for j in range(0,shapS[2]):
            im_i+=S[:,:,j]*A[j,i]
        output[:,:,i] = (flux[i]/(flux_ref*sig[i]))*utils_for_rca.decim(fftconvolve(im_i,ker[:,:,i],mode='same'),upfact,av_en=0)
    return output

def sr_stack_trans(input):
    y = input[0]
    upfact = input[1]
    ker = input[2]
    ker_adj = input[3]
    A = input[4]
    shapA = A.shape
    sig = input[5]
    flux = input[6]
    flux_ref = median(flux)
    shapy = y.shape
    output = zeros((shapy[0]*upfact,shapy[1]*upfact,shapA[0]))
    nb_im = shapy[2]

    for i in range(0,nb_im):
        im_i = (flux[i]/(flux_ref*sig[i]))*convolve(utils_for_rca.transpose_decim(y[:,:,i],upfact),ker_adj[:,:,i],mode='same')
        for j in range(0,shapA[0]):
            output[:,:,j] += im_i*A[j,i]
    return output

def sr_stack_trans_op_src(input): # M^TM
    output1 = sr_stack_op(input)
    input2 = cp.deepcopy(input)
    input2[0] = output1
    output = sr_stack_trans(input2)
    return output

def pow_meth(opname,op_param,siz,tol=0.5,ainit=None,nb_iter_max=30,opt_vect=None):

    a = None
    if ainit is not None:
        a = ainit
    else:
        if size(siz)==1:
            a = random.randn(siz)
        elif size(siz)==2:
            a = random.randn(siz[0],siz[1])
        elif size(siz)==3:
            a = random.randn(siz[0],siz[1],siz[2])
        a = a/sqrt(((a)**2).sum())

    L_old = 10
    L=1
    i=0
    print "----------- spec_rad est ---------- : "

    while i<nb_iter_max and (100*abs(L-L_old))/L>tol:
        print L
        op_param[0] = a
        b = opname(op_param)
        if opt_vect is not None:
            b-=(b*opt_vect).sum()*opt_vect
        L_old = L
        L = sqrt(((b)**2).sum())
        a = b/L
        i+=1

    if i==nb_iter_max:
        print "Warning max number if iterations reached in pow_meth"
    return a,L

def non_uniform_smoothing_bas_mat(weights):
    shap = weights.shape
    A = zeros((shap[0],shap[0]))
    B = zeros((shap[0],shap[0]))
    range_ref = array(range(0,shap[0]))
    for i in range(0,shap[0]):
        ik = where(range_ref!=i)
        A[ik[0],i] = -weights[i,:]
        B[i,i] = (weights[i,:]).sum()
    return A,B

def non_uniform_smoothing_mat_2(weights,e): # e>=0
    A,B = non_uniform_smoothing_bas_mat(weights)
    mat_out = A.dot(transpose(A))+ e*(A.dot(B)+ B.dot(transpose(A))) + (e**2)*B.dot(B) # B is diagonal matrix
    return mat_out

def non_uniform_smoothing_mat_dist_1(dist,expo_range,e):
    dist_med = median(dist)
    nb_samp = len(expo_range)
    nb_im = dist.shape[0]
    mat_stack = zeros((nb_im,nb_im,nb_samp))
    for i in range(0,nb_samp):
        dist_weights = (dist_med/dist)**expo_range[i]
        dist_weigths = dist_weights/dist_weights.max()
        mat_stack[:,:,i] = non_uniform_smoothing_mat_2(dist_weights,e)
    return mat_stack

def non_uniform_smoothing_mat_dist_2(dist,expo,e_range):
    dist_med = median(dist)
    nb_samp = len(e_range)
    nb_im = dist.shape[0]
    mat_stack = zeros((nb_im,nb_im,nb_samp))
    for i in range(0,nb_samp):
        dist_weights = (dist_med/dist)**expo
        dist_weigths = dist_weights/dist_weights.max()
        mat_stack[:,:,i] = non_uniform_smoothing_mat_2(dist_weights,e_range[i])
    return mat_stack

def notch_filt_optim_2(test_mat,dist,expo_range,e_range,nb_iter=2,tol=0.01):
    expo_out = None
    e_out = 0.5
    loss = None
    vect = None
    ker = None
    j2 = None
    for i in range(0,nb_iter):
        mat_stack = non_uniform_smoothing_mat_dist_1(dist,expo_range,e_out)
        vect,j,loss,ker,j2 = utils_for_rca.kernel_mat_stack_test_unit(mat_stack,test_mat,tol=tol)
        expo_out = expo_range[j]
        mat_stack = non_uniform_smoothing_mat_dist_2(dist,expo_out,e_range)
        vect,j,loss,ker,j2 = utils_for_rca.kernel_mat_stack_test_unit(mat_stack,test_mat,tol=tol)
        e_out = e_range[j]

    return expo_out,e_out,loss,vect,ker,j2

def analysis(cube,sig,field_pos,pmin = 0.01,emin=0.01,emax=1.99,nb_max=30,tol=0.01):
    nb_samp_opt = 10
    shap = cube.shape
    nb_neighs = shap[2]-1
    neigh,dists = utils_for_rca.knn_interf(field_pos,nb_neighs)
    p_max = pow_law_select(dists,nb_neighs)
    p_min = 0.01
    e_min = 0.01
    e_max = 1.99
    dists_unsorted = utils_for_rca.feat_dist_mat(field_pos)
    e_range = utils_for_rca.log_sampling(e_min,e_max,nb_samp_opt)
    p_range = utils_for_rca.log_sampling(p_min,p_max,nb_samp_opt)
    res_mat = copy(transpose(cube.reshape((shap[0]*shap[1],shap[2]))))
    list_comp = list()
    list_e = list()
    list_p = list()
    list_ind = list()
    list_ker = list()
    err = (cube**2).sum()
    nb_iter = 0
    while err>sig and nb_iter<nb_max:
        expo_out,e_out,loss,vect,ker,j = notch_filt_optim_2(res_mat,dists_unsorted,p_range,e_range,nb_iter=3,tol=tol)
        list_e.append(e_out)
        list_p.append(expo_out)
        list_comp.append(vect)
        list_ind.append(j)
        list_ker.append(ker)
        nb_iter+=1
        res_mat = res_mat-transpose(vect).dot(vect.dot(res_mat))
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
    comp = utils_for_rca.mat_to_cube(proj_coeff,shap[0],shap[1])
    proj_data = transpose(weights).dot(proj_coeff)
    proj_data = utils_for_rca.mat_to_cube(proj_data,shap[0],shap[1])
    return e_vect,p_vect,weights,comp,proj_data,ker,ind

def pow_law_select(dist_weights,nb_neigh,min_val=10**(-15)):
    a = dist_weights[:,0]/dist_weights[:,nb_neigh-1]
    r_med = a.min()
    p = log(min_val)/log(r_med)
    return p

def man_inv(mat,cond=None):
    U, s, Vt = svd(mat,full_matrices=False)
    eig = zeros((s.shape[0],s.shape[0]))
    for i in range(0,s.shape[0]):
        if cond is not None:
            if s[i]>s[0]/cond:
                eig[i,i] = (s[i])**(-1)
            else:
                eig[i,i] = cond/s[0]
        else:
            if s[i]>0:
                eig[i,i] = (s[i])**(-1)
    inv = U.dot(eig.dot(Vt))
    return inv

def lsq_mult_coeff(im,src,man=True):
    shap = src.shape
    mat = zeros((shap[2],shap[2]))
    for i in range(0,shap[2]):
        for j in range(0,shap[2]):
            mat[i,j] = (src[:,:,i]*src[:,:,j]).sum()
    v = zeros((shap[2],1))
    for k in range(0,shap[2]):
        v[k,0] = (im*src[:,:,k]).sum()
    output = None
    if man:
        output = man_inv(mat).dot(v)
    else:
        output = inv(mat).dot(v)
    return output,mat,v

def lsq_mult_coeff_stack(im,src,man=True):
    shap1 = src.shape
    shap2 = im.shape
    coeff_out = zeros((shap1[2],shap2[2]))
    v = zeros((shap1[2],shap2[2]))
    mat = zeros((shap1[2],shap1[2],shap1[3]))
    for i in range(0,shap2[2]):
        outi,mati,vi = lsq_mult_coeff(im[:,:,i],src[:,:,:,i],man=man)
        coeff_out[:,i] = outi.reshape((shap1[2],))
        mat[:,:,i] = copy(mati)
        v[:,i] = vi.reshape((shap1[2],))
    return coeff_out,mat,v

def non_unif_smoothing_mult_coeff_pos_cp_6(im,src,nb_iter=100,verbose=True):
    # to is the weight related to the primal variable;  basis is a concatenation
    # of the optimal notch filter operator eigenvectors
    Ainit,mat,v = lsq_mult_coeff_stack(im,src)
    shap = src.shape
    spec_rad = zeros((shap[3]))
    for k in range(0,shap[3]):
        U, s, Vt = svd(mat[:,:,k],full_matrices=False)
        spec_rad[k] = s.max()#*(basis[:,k]**2).sum()
    spec_norm = spec_rad.sum()
    A = copy(Ainit)
    t = 1
    Ax = copy(A)
    rad = zeros((shap[3],))
    for k in range(0,shap[3]):
        res = copy(im[:,:,k])
        for l in range(0,shap[2]):
            res-=Ainit[l,k]*src[:,:,l,k]
        rad[k] = (res**2).sum()
    if verbose:
        print "--->> ref res: <<---",rad.sum()
    ref_res = rad.sum()
    cost = 1
    cost_old=0
    i=0
    while (i < nb_iter) and (100*abs((cost-cost_old)/cost)>0.01 or cost>1.1*ref_res) :
        res = copy(im)
        for k in range(0,shap[3]):
            for l in range(0,shap[2]):
                res[:,:,k]-=A[l,k]*src[:,:,l,k]
        if verbose:
            print " -------- mse: ",(res**2).sum(),"-----------"
        cost_old = cost
        cost = (res**2).sum()
        temp = Ainit*0
        for k in range(0,shap[3]):
            for l in range(0,shap[2]):
                temp[l,k]+=(src[:,:,l,k]*res[:,:,k]).sum()
        grad = -temp
        Ay = A - grad/spec_norm
        Ax_old = copy(Ax)
        Ax = copy(Ay)
        told = t
        t = (1+sqrt(4*t**2 +1))/2
        lambd = 1 + (told-1)/t
        A  = Ax_old + lambd*(Ax-Ax_old)
        i+=1
    return A

def non_unif_smoothing_mult_coeff_pos_cp_5(im,src,src_hr,tree,basis,alpha_init,\
                                            theta=0.1,p_smth_mat_inv=None,to=None,eps=0.01,nb_iter=100,tol=0.01,Ainit=None,\
                                            pos_en=False,reg_param=1000,spars_en=True, verbose=True):
    # to is the weight related to the primal variable;
    # basis is a concatenation of the optimal notch filter
    # operator eigenvectors
    shap = src.shape
    shap1 = src_hr.shape
    src_mat = zeros((shap[2],shap[2]))
    for i in range(0,shap[2]):
        for j in range(i+1,shap[2]):
            src_mat[i,j] = (src_hr[:,:,i]*src_hr[:,:,j]).sum()
    src_mat = src_mat+transpose(src_mat)
    for i in range(0,shap[2]):
        src_mat[i,i] = (src[:,:,i]**2).sum()
    U, s, Vt = svd(src_mat,full_matrices=False)
    U, s2, Vt = svd(basis.dot(transpose(basis)),full_matrices=False)
    spec_rad_pos = s.max()*s2.max()
    nb_neigh = tree.shape[1]
    Atemp,mat,v = lsq_mult_coeff_stack(im,src)
    spec_rad = zeros((shap[3]))
    rad = zeros((shap[3],))
    for k in range(0,shap[3]):
        res = copy(im[:,:,k])
        for l in range(0,shap[2]):
            res-=Ainit[l,k]*src[:,:,l,k]
        rad[k] = (res**2).sum()
    if verbose:
        print "--->> ref res: <<---",rad.sum()
    ref_res = rad.sum()
    cost = 1
    cost_old=0
    for k in range(0,shap[3]):
        U, s, Vt = svd(mat[:,:,k],full_matrices=False)
        spec_rad[k] = s.max()#*(basis[:,k]**2).sum()
    spec_norm=spec_rad.sum()*s2.max()
    shapb = basis.shape
    alpha = copy(alpha_init)*0
    i=0
    t = 1
    alphax = copy(alpha)
    shap_alpha = alpha.shape
    supports = zeros((shap_alpha[0],shap_alpha[1],min(nb_iter,shapb[0])))
    while (i < min(nb_iter,shapb[0])) and (100*abs((cost-cost_old)/cost)>0.01 or cost>1.1*ref_res) :
        A = alpha.dot(basis)
        res = copy(im)
        for k in range(0,shap[3]):
            for l in range(0,shap[2]):
                res[:,:,k]-=A[l,k]*src[:,:,l,k]
        if verbose:
            print " -------- mse: ",(res**2).sum(),"-----------"
        cost_old = cost
        cost = (res**2).sum()
        temp = Ainit*0
        for k in range(0,shap[3]):
            for l in range(0,shap[2]):
                temp[l,k]+=(src[:,:,l,k]*res[:,:,k]).sum()
        grad = -temp.dot(transpose(basis))
        alphay = alpha - grad/spec_norm
        alphax_old = copy(alphax)
        if spars_en:
            alphax = utils_for_rca.lineskthresholding(alphay,int(floor(sqrt(i)))+1)
            supports[:,:,i] = copy(alphax)
        else:
            alphax = copy(alphay)
        told = t
        t = (1+sqrt(4*t**2 +1))/2
        lambd = 1 + (told-1)/t
        alpha  = alphax_old + lambd*(alphax-alphax_old)
        supp = where(abs(alpha[0,:])>0)
        i+=1
    mat_out = alpha.dot(basis)
    return mat_out,alpha,supports


def low_rank_global_src_est_comb(input,weights,y,ksig=4,eps=0.9,ainit=None,nb_iter=100,tol=0.01,nb_rw=0,\
                                pos_en=True,opt='coif5',nb_scale=None,wav_en=False,wavr_en=False,optr = None,\
                                Y2=None,Y3=None,V=None,cY3=None,only_noise_est_en=False,\
                                select_en=True,thresh_perc = 0.01,rad=None,filters=None,filters_rot=None,iter_min=10,\
                                mu=0.1,nb_sc = 4,nb_dir=8,verbose=True):

    S = copy(input[0]) # (Main variable)
    ref_mse = ((y-sr_stack_op(input))**2).sum()
    l1_norm = abs(S*weights).sum()
    if verbose:
        print "--->>ref res: <<---",ref_mse
        print "--->>ref l1 norm: <<---",l1_norm
    ker_adj = input[3]
    shap = ker_adj.shape
    A = input[4]
    nb_im = shap[2]
    shapS = S.shape
    ind_select = range(0,shapS[2])
    spec_rad = 0
    # Spectral radius setting
    eig_max,spec_rad1 = pow_meth(sr_stack_trans_op_src,input,shapS,ainit=ainit)
    U, s, Vt = svd(A.dot(transpose(A)),full_matrices=False)
    spec_rad2 = sqrt(s[0])
    if Y2 is None:
        Y2 = zeros((shapS[0],shapS[1],nb_im))
    rweights_an = None
    weights_an = None
    spec_rad3 = 0
    if wavr_en:
        if Y3 is None or filters is None:
            Y3,filters = utils_for_rca.mr_trans_stack_2(S*0,opt=optr)
            filters_rot = utils_for_rca.rot90_stack(filters)
            Y3[:,:,-1,:]*=0 # Puts the coarse scales to 0
        for i in range(0,filters.shape[2]):
            spec_rad3 +=(abs(filters[:,:,i]).sum())**2
        spec_rad3 = sqrt(spec_rad3)
        weights_an_temp,filters_temp = utils_for_rca.mr_trans_stack_2(weights,filters=filters**2)
        weights_an = 4*sqrt(weights_an_temp[:,:,:-1,:])
        rweights_an = copy(Y3[:,:,:-1,:])
    if wavr_en and pos_en:
        spec_rad = spec_rad1+spec_rad2+spec_rad3
    elif pos_en:
        spec_rad = spec_rad1+spec_rad2
    elif wavr_en:
        spec_rad = spec_rad1+spec_rad3
    else:
        spec_rad = spec_rad1
    rweights = ones((shapS[0],shapS[1],shapS[2]))
    cost=0
    cost_old=1
    iter=0
    shapy = y.shape
    ones_mat = ones((shapS[0],shapS[1],shapS[2]))
    input1 = cp.deepcopy(input)
    input2 = cp.deepcopy(input)
    for l in range(0,nb_rw+1):
        if verbose:
            print l+1,"th pass/",nb_rw+1
        while (iter<nb_iter) and (100*abs(cost-cost_old)/abs(cost_old)>0.001 or iter<iter_min):
            if verbose:
                print "tol: ",100*abs(cost-cost_old)/abs(cost_old)
            iter+=1
            input1[0] = copy(S)
            est = sr_stack_op(input1)
            res = est - y
            input2[0] = res
            gradS = sr_stack_trans(input2)
            temp30 = None
            if pos_en:
                temp30 = S*0
                for i in range(0,shapS[2]):
                    for j in range(0,nb_im):
                        temp30[:,:,i]+=Y2[:,:,j]*A[i,j]
            temp31 = None
            ctemp31 = None
            if wavr_en:
                temp31 = utils_for_rca.mr_transf_transp_stack(Y3,filters_rot)
            temp3 = None
            ctemp3 = None
            if wavr_en and pos_en:
                temp3 = temp30+temp31
            elif pos_en:
                temp3 = temp30
            elif wavr_en:
                temp3 = temp31
            else:
                temp3 = S*0
            P = S - mu*(temp3+gradS)/spec_rad
            if wavr_en is not True:
                P = utils_for_rca.thresholding_3D(P,mu*weights*rweights/spec_rad,1)
            Y = 2*P-S
            S = S+(1-1.0/(iter+1))*(P-S)
            if pos_en:
                temp2 = Y2*0
                for i in range(0,nb_im):
                    for j in range(0,shapS[2]):
                        temp2[:,:,i]+=Y[:,:,j]*A[j,i]
                tY2 = -pos_proj_cube(-Y2-temp2/spec_rad)
                Y2 = Y2+(1-1.0/(iter+1))*(tY2-Y2)
            if wavr_en:
                temp3,filters = utils_for_rca.mr_trans_stack_2(Y,filters=filters)
                for k in range(0,shapS[2]):
                    temp4 = None
                    if l==0:
                        temp4 = utils_for_rca.thresholding_3D(Y3[:,:,:-1,k]+temp3[:,:,:-1,k]*spec_rad,mu*weights_an[:,:,:,k],1)
                    else:
                        temp4 = utils_for_rca.thresholding_3D(Y3[:,:,:-1,k]+temp3[:,:,:-1,k]*spec_rad,mu*rweights_an[:,:,:,k]*weights_an[:,:,:,k],1)
                    tY3 = Y3[:,:,:-1,k]+temp3[:,:,:-1,k] - temp4/spec_rad
                    Y3[:,:,:-1,k] = Y3[:,:,:-1,k]+(1-1.0/(iter+1))*(tY3-Y3[:,:,:-1,k])

            cost_old = cost
            cost = (res**2).sum()
            # ----------- Sanity check ----------- #
            if verbose:
                print "Mini val: ",est.min(),"; Residual: ",cost," Gradient's norm: ",(gradS**2).sum()," spectral norm: ",spec_rad
            if wavr_en is not True:
                if verbose:
                    print "Weighted l1 norm direct domain: ",(abs(weights*rweights*S)).sum()
            else:
                trans_data,filters = utils_for_rca.mr_trans_stack_2(S,filters=filters)
                if verbose:
                    if l==0:
                        print "Weighted l1 norm analysis: ",(abs(trans_data[:,:,:-1,:])).sum()
                    else:
                        print "Weighted l1 norm analysis: ",(abs(trans_data[:,:,:-1,:]*rweights_an)).sum()
        if wavr_en:
            rweights_an = (1+(abs(trans_data[:,:,:-1,:])/weights_an))**(-1)
            if verbose:
                print "Weight max: ",rweights_an.max()," Weight min: ",rweights_an.min()
        else:
            rweights  = (1+(abs(S)/weights))**(-1)
        iter = 0
        cost=2
        cost_old=1
    list_ind = list()
    for id in range(0,shapS[2]):
        a = ksig*sqrt((S[:,:,id]**2).sum())/sqrt((weights[:,:,id]**2).sum())
        if verbose:
            print "Src ",id," PNR: ",a
        if a>=1:
            list_ind.append(id)
    if select_en and len(list_ind)>0:
        ind_select = tuple(list_ind)
        S = S[:,:,ind_select]
    return filters,filters_rot,Y2,Y3,cY3,S,ind_select

def rca_main_routine(psf_stack_in,field_pos,upfact,opt,nsig,sparsity_en=True,\
                                                        pix_sparsity=True,dist_weight_deg=1,shifts=None,\
                                                        opt_shift_est=['-t2','-n2'],nsig_shift_est=5,sig_est=None,flux_est=None,nb_iter=2,\
                                                        nb_subiter=300,nb_comp_max=10,tol=0.1,\
                                                        positivity_en=False,nb_rw=1,lsq_en=False,\
                                                        shifts_regist = True,wavr_en=False,verbose=True):

    psf_stack = copy(psf_stack_in)
    shap = psf_stack.shape
    if nb_comp_max>shap[2]:
        print "/!\ Warning: number of components higher than the number of images; reduced to ",shap[2]
    siz_in = upfact*array(shap[0:-1])
    if sparsity_en is False:
        print "----- Sparsity disable -----"
    " ============================== Degradation operator parameters estimation =============================== "
    centroids = None
    if sig_est is None:
        print 'Noise level estimation...'
        sig_est,filters_lr = utils_for_rca.im_gauss_nois_est_cube(psf_stack_in,opt=opt_shift_est)
        #print sig_est
        print 'Done.'
    if shifts_regist:
        if shifts is None:
            map = ones((shap[0],shap[1],shap[2]))
            for i in range(0,shap[2]):
                nsig_shifts = min(nsig_shift_est,0.8*psf_stack_in[:,:,i].max()/sig_est[i])
                map[:,:,i] *= nsig_shifts*sig_est[i]
            print 'Shifts estimation...'
            psf_stack_shift = utils_for_rca.thresholding_3D(psf_stack_in,map,0)
            shifts,centroids = utils_for_rca.shift_est(psf_stack_shift)
            print 'Done.'

        else:
            print "------------ /!\ Warning: shifts provided /!\ ---------"
    else:
        print "------------ /!\ Warning: no registration /!\ ---------"
        shifts = zeros((shap[2],2))
    if flux_est is None:
        flux_est = utils_for_rca.flux_estimate_stack(psf_stack,rad=4)
    flux_ref = median(flux_est)
    shift_ker_stack,shift_ker_stack_adj = utils_for_rca.shift_ker_stack(shifts,upfact)
    sig_min = sig_est.min()
    sig_min_vect = ones((shap[2],))*sig_min
    sig_est = sig_est/sig_min
    nb_im = shap[2]
    for k in range(0,shap[2]):
        psf_stack[:,:,k] = psf_stack[:,:,k]/sig_est[k]
    if verbose:
        print " ------ ref energy: ",(psf_stack**2).sum()," ------- "
        print "flux min: ",flux_est.min()," flux max: ",flux_est.max(),"sig min: ",\
                sig_est.min()," sig max: ",sig_est.max()
    weights = None
    ref_weights = None
    w = ones((shap[2],))
    im_hr = zeros((upfact*shap[0],upfact*shap[1],shap[2]))
    " ================================ FOV distances settings ================================ "
    print "Contructing PSF tree..."
    nb_neighs = shap[2]-1
    neigh,dists = utils_for_rca.knn_interf(field_pos,nb_neighs)
    p_max = pow_law_select(dists,nb_neighs)
    if verbose:
        print "power max = ",p_max
    p_min = 0.01
    print "Done..."
    dists_unsorted = utils_for_rca.feat_dist_mat(field_pos)
    dist_med = median(dists)
    dist_weights = (dist_med/dists_unsorted)**dist_weight_deg
    dist_weigths = dist_weights/dist_weights.max()
    spec_rad_smooth = None # Spatial smoothing operator gradient lip constant estimation
    siz_in = upfact*array(shap[0:-1])
    eig_max = None
    " ============================== Inputs variables ============================== "
    input = list()
    x = zeros((upfact*shap[0],upfact*shap[1]))
    input.append(x)
    input.append(upfact)
    input.append(shift_ker_stack)
    input.append(shift_ker_stack_adj)
    w = ones((shap[2],))
    w_old = ones((shap[2],))
    im_hr = zeros((upfact*shap[0],upfact*shap[1],shap[2]))
    input.append(w)
    input.append(sig_est)
    input.append(flux_est)
    comp = zeros((upfact*shap[0],upfact*shap[1],nb_comp_max)) # Main variable
    u,mr_file = utils_for_rca.mr_trans(comp[:,:,0],opt=opt)
    os.remove(mr_file)
    shap2 = u.shape
    comp_lr = zeros((shap[0],shap[1],nb_comp_max,shap[2]))
    i = 0
    Y1 = None
    Y2 = None
    Y3 = None
    cY3 = None
    V = None
    filters=None
    filters_rot = None
    rad = None
    select_en = False
    to = None
    if verbose:
        print " ------ ref energy: ",(psf_stack**2).sum()," ------- "
    if shap[2]==2:
        if verbose:
            print "Number of sources insufficient to use the spatial constraint; activating the lsq"
        lsq_en = True

    " ============================== Weights init ============================== "
    res = copy(psf_stack)
    if lsq_en:
        if shifts_regist:
            weights,coeff_res,cube_est = utils_for_rca.cube_svd(utils_for_rca.rect_crop_c(res,\
            int(0.9*shap[0]),int(0.9*shap[1]),centroids),nb_comp=nb_comp_max)
        else:
            weights,coeff_res,cube_est = utils_for_rca.cube_svd(res,nb_comp=nb_comp_max)
        weights = weights[0:nb_comp_max,:]
        for l in range(0,nb_comp_max):
            a = sqrt((weights[l,:]**2).sum())
            if a>0:
                weights[l,:] /= a
    else:
        if shifts_regist:
            e_opt,p_opt,weights,comp_temp,data,ker,alph_ref  = analysis(utils_for_rca.rect_crop_c(res,int(0.9*shap[0])\
            ,int(0.9*shap[1]),centroids),int(0.9*shap[0])*int(0.9*shap[1])*sig_min**2,field_pos,tol=0,nb_max=nb_comp_max)
        else:
            e_opt,p_opt,weights,comp_temp,data,ker,alph_ref  = analysis(res,int(0.9*shap[0])*int(0.9*shap[1])*sig_min**2,\
                                                                       field_pos,tol=0,nb_max=nb_comp_max)
        alph = alph_ref
        for l in range(0,nb_comp_max):
            a = sqrt((weights[l,:]**2).sum())
            if a>0:
                weights[l,:] /= a
    weights_temp,coeff_res_temp,cube_est = utils_for_rca.cube_svd(res,nb_comp=nb_comp_max)
    if verbose:
        print "============================>>>> min res: ",((cube_est-res)**2).sum()," <<<<============================"
    ainit = None
    input_ref = list()
    input_ref.append(copy(comp))
    input_ref.append(upfact)
    input_ref.append(shift_ker_stack)
    input_ref.append(shift_ker_stack_adj)
    input_ref.append(weights)
    ref_weights = copy(weights)
    input_ref.append(sig_est)
    input_ref.append(flux_est)
    survivors = ones((nb_comp_max,))



    for k in range(0,nb_iter):
        " ============================== Sources estimation =============================== "
        thresh = nsig*utils_for_rca.acc_sig_maps(shap,shift_ker_stack_adj,sig_est,flux_est,\
        flux_ref,upfact,weights,sig_data=sig_min_vect)
        if k==1:
            select_en = True
        filters,filters_rot,Y2,Y3,cY3,comp,ind_select = \
        low_rank_global_src_est_comb(input_ref,thresh,psf_stack,ksig=nsig,eps=0.8\
        ,ainit=ainit,nb_iter=nb_subiter,tol=1,nb_rw=nb_rw,Y2=Y2,V=V,rad=None,\
        select_en=select_en,wavr_en=wavr_en,optr=opt,pos_en=positivity_en,\
        filters=filters,filters_rot=filters_rot,verbose=verbose)
        comp_lr = zeros((shap[0],shap[1],comp.shape[2],shap[2]))
        survivors = zeros((nb_comp_max,))
        for l in range(0,comp.shape[2]):
            for p in range(0,shap[2]):
                comp_lr[:,:,l,p] = (flux_est[p]/(sig_est[p]*flux_ref))*\
                utils_for_rca.decim(fftconvolve(comp[:,:,l],\
                shift_ker_stack[:,:,p],mode='same'),upfact,av_en=0)
        id0 = where((survivors==1))
        id = id0[0]
        " ============================== Weights estimation =============================== "
        n_max = nb_iter-1
        if k < n_max:
            weights_k = None
            if lsq_en:
                weights_k = non_unif_smoothing_mult_coeff_pos_cp_6(psf_stack,comp_lr,nb_iter=nb_subiter*2)
            else:
                weights_k,alph,supports = non_unif_smoothing_mult_coeff_pos_cp_5\
                (psf_stack,comp_lr,comp,neigh,ker,alph[ind_select,:],\
                to=to,nb_iter=nb_subiter*2,tol=0.1,Ainit=weights[ind_select,:],\
                pos_en=positivity_en)
            list_surv = list()
            for l in range(0,comp.shape[2]):
                a = sqrt((weights_k[l,:]**2).sum())
                if a>0:
                    list_surv.append(l)
                    if wavr_en:
                        comp[:,:,l] *= a
                        weights_k[l,:] /= a
            ind_select = tuple(list_surv)
            weights = weights_k[ind_select,:]
            input_ref[0] = comp[:,:,ind_select]
            input_ref[4] = weights
    for l in range(0,shap[2]):
        for p in range(0,comp.shape[2]):
            im_hr[:,:,l] = im_hr[:,:,l]+weights[p,l]*comp[:,:,p]
            res[:,:,l] = psf_stack[:,:,l]-(flux_est[l]/(sig_est[l]*flux_ref))*\
            utils_for_rca.decim(fftconvolve(im_hr[:,:,l],shift_ker_stack[:,:,l],mode='same')\
            ,upfact,av_en=0)

    return im_hr,comp,weights,res,sig_est*sig_min,flux_est,\
           shifts,alph,alph_ref,e_opt,p_opt,ker,supports
