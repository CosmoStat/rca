# -*- coding: utf-8 -*-

'''RCA ARGUMENTS

Command line arguments for rca.py

'''

import sys
import argparse as ap
import ConfigParser as cp
from os.path import isfile


def get_args():

    '''Get options

    This method defines the command line agruments for RCA and reads the
    specified configuration file.

    Returns
    -------
    The agrument parser

    '''

    formatter = ap.ArgumentDefaultsHelpFormatter
    config_parser = ap.ArgumentParser(description=__doc__,
                                      formatter_class=formatter,
                                      add_help=False)
    config_parser.add_argument('-c', '--config_file',
                               default='rca_config.ini',
                               required=False,
                               help='Specify configuration file',
                               metavar='FILE')
    args, remaining_argv = config_parser.parse_known_args()

    defaults = {'upsample': 2,
                'wavelet_scales': '3',
                'wavelet_type': '2',
                'sparse_strength': 4,
                'n_atoms': 10,
                'n_altmin': 2,
                'n_iter': 200,
                'n_reweight': 1,
                'use_sparsity': True,
                'use_wavelets': True,
                'use_positivity': True}

    if isfile(args.config_file):
        config = cp.SafeConfigParser()
        config.read([args.config_file])
        defaults.update(dict(config.items('Defaults')))

    ints = ['upsample', 'sparse_strength', 'n_atoms', 'n_altmin',
            'n_iter', 'n_reweight']
    bools = ['use_sparsity', 'use_wavelets', 'use_positivity']
    defaults = {k: int(v) if k in ints else eval(v) if k in bools else v
                for k, v in defaults.items()}

    parser = ap.ArgumentParser(parents=[config_parser],
                               formatter_class=formatter)
    parser.set_defaults(**defaults)
    parser.add_argument('--psf_file', help='Input PSFs file')
    parser.add_argument('--pos_file', help='PSF positions file')
    parser.add_argument('--output_file', help='Output file name')
    parser.add_argument('--upsample', help='Upsampling factor')
    parser.add_argument('--wavelet_scales', help='Number of wavelet sacales')
    parser.add_argument('--wavelet_type', help='Type of wavelet to use')
    parser.add_argument('--sparse_strength',
                        help='Strength of sparsity constraint')
    parser.add_argument('--n_atoms', help='Number of atoms to learn')
    parser.add_argument('--n_altmin', help='Number of alternate minimizations')
    parser.add_argument('--n_iter', help='Number of iterations')
    parser.add_argument('--n_reweight', help='Number of reweightings')
    parser.add_argument('--use_sparsity', help='Option to use sparsity')
    parser.add_argument('--use_wavelets', help='Option to use wavelets')
    parser.add_argument('--use_positivity', help='Option to use positivity')
    parser.add_argument('--noise_estimate', help='Estimate of the noise')
    parser.add_argument('--flux_estimate', help='Estimate of the flux')

    return parser.parse_args(remaining_argv)
