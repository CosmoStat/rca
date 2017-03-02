#! /usr/bin/env Python
# -*- coding: utf-8 -*-

"""RCA

Resolved Component Analysis

http://www.cosmostat.org/software/rca/

:Author: Fred Ngole Mboula

:Version: 1.0

:Date: 16/02/2017

"""

import warnings
from rca_agrs import get_args
from rca_lib import rca_main_routine
from astropy.io import fits


def run_script():
    '''Run the script

    This method reads the input FITS files, calls the RCA main routine from
    rca_lib.py and outputs the results in FITS format to the specified file
    name.

    '''

    print args.use_positivity

    psf_data = fits.getdata(args.psf_file)
    pos_data = fits.getdata(args.pos_file)
    mr_trans_opt = ['-t' + args.wavelet_type, '-n' + args.wavelet_scales]

    res = rca_main_routine(psf_data, pos_data, args.upsample, mr_trans_opt,
                           args.sparse_strength, nb_iter=args.n_altmin,
                           nb_subiter=args.n_iter, nb_comp_max=args.n_atoms,
                           nb_rw=args.n_reweight,
                           sparsity_en=args.use_sparsity,
                           wavr_en=args.use_wavelets,
                           positivity_en=args.use_positivity,
                           sig_est=args.noise_estimate,
                           flux_est=args.flux_estimate)

    fits.writeto(args.output_file, res[0])


def main():

    warnings.simplefilter("ignore")

    try:
        global args
        args = get_args()
        run_script()

    except Exception, err:
        print err
        return 1


if __name__ == "__main__":
    main()
