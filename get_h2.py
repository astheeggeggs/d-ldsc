#!/usr/bin/env python

from __future__ import division
import src.ldscoregxe as ld
import src.parse as ps
import src.sumstats as sumstats
import src.regressions as reg
import src.printing as pr
import numpy as np
import pandas as pd
from subprocess import call
from itertools import product
import time, sys, traceback, argparse

try:
    x = pd.DataFrame({'A': [1, 2, 3]})
    x.drop_duplicates(subset='A')
except TypeError:
    raise ImportError('LDSC requires pandas version > 0.15.2')

parser = argparse.ArgumentParser()
parser.add_argument('--out', default='ldsc', type=str,
    help='Output filename prefix. If --out is not set, LDSC will use ldsc as the '
    'defualt output filename prefix.')
parser.add_argument('--h2', default=None, type=str,
    help='Filename prefix for a .chisq file for one-phenotype LD Score regression. '
    'LDSC will automatically append .chisq or .chisq.gz to the filename prefix.'
    '--h2 requires at minimum also setting the --ref-ld and --w-ld flags.')
# Basic Flags for Working with Variance Components
parser.add_argument('--ref-ld', default=None, type=str,
    help='Use --ref-ld to tell LDSC which LD Scores to use as the predictors in the LD '
    'Score regression. '
    'LDSC will automatically append .ldscore/.ldscore.gz to the filename prefix.')
parser.add_argument('--ref-ld-chr', default=None, type=str,
    help='Same as --ref-ld, but will automatically concatenate ldscore files split '
    'across 22 chromosomes. LDSC will automatically append ldscore/ldscore.gz '
    'to the filename prefix. If the filename prefix contains the symbol @, LDSC will '
    'replace the @ symbol with chromosome numbers. Otherwise, LDSC will append chromosome '
    'numbers to the end of the filename prefix.'
    'Example 1: --ref-ld-chr ld/ will read ld/1.ldscore.gz ... ld/22.ldscore.gz'
    'Example 2: --ref-ld-chr ld/@_kg will read ld/1_kg.ldscore.gz ... ld/22_kg.ldscore.gz')
parser.add_argument('--w-ld', default=None, type=str,
    help='Filename prefix for file with LD Scores with sum r^2 taken over SNPs included '
    'in the regression. LDSC will automatically append .ldscore/.ldscore.gz.')
parser.add_argument('--w-ld-chr', default=None, type=str,
    help='Same as --w-ld, but will read files split into 22 chromosomes in the same '
    'manner as --ref-ld-chr.')
parser.add_argument('--no-intercept', action='store_true',
    help = 'This constrains the LD Score regression intercept to equal 1.')
parser.add_argument('--M', default=None, type=str,
    help='# of SNPs (if you don\'t want to use the .M files that came with your .ldscore.gz files)')
parser.add_argument('--two-step', default=None, type=float,
    help='Test statistic bound for use with the two-step estimator. Not compatible with --no-intercept and --constrain-intercept.')
parser.add_argument('--chisq-max', default=None, type=float,
    help='Max chi^2.')
# Flags you should almost never use
parser.add_argument('--n-blocks', default=200, type=int,
    help='Number of block jackknife blocks.')
parser.add_argument('--not-M-5-50', default=False, action='store_true',
    help='This flag tells LDSC to use the .M file instead of the .M_5_50 file.')
# Transform to liability scale
parser.add_argument('--samp-prev',default=None,
    help='Sample prevalence of binary phenotype (for conversion to liability scale).')
parser.add_argument('--pop-prev',default=None,
    help='Population prevalence of binary phenotype (for conversion to liability scale).')
parser.add_argument('--additive', default=False, action='store_true', 
    help='Include the additive component in the LD score regression: this option assumes that the column labelled Z_A is the additive Z scores.')
parser.add_argument('--additive-orig', default=False, action='store_true', 
    help='Include the additive component in the LD score regression: this option assumes that the column labelled Z is the additive Z scores.')
parser.add_argument('--dominance', default=False, action='store_true', 
    help='Include the dominance component in the LD score regression.')
parser.add_argument('--gxe', default=False, action='store_true', 
    help='Include the gene x environment component in the LD score regression.')
parser.add_argument('--write-h2', default=False, action='store_true', 
    help='Write heritability estimates the .h2 file.')
parser.add_argument('--N-CHR', default=22, type=int,
    help='Primarily for simulation studies and debugging, if the number of'
    ' `chromosomes` is not 22')
parser.add_argument('--average-cov-kurtosis', default=None, type=float,
    help='Average Kurtosis of the covariate used in the GxE interaction')
parser.add_argument('--cov-case-control', default=False, action='store_true',
    help='Is the GxE covariate being considered boolean?')
parser.add_argument('--cov-case-prev', default=None, type=float,
    help='What is the prevalence of being a `case` for the GxE covariate?')
parser.add_argument('--pheno-name', default=None, type=str,
    help='What is the name of the phenotype?')

if __name__ == '__main__':

    args = parser.parse_args()
    if args.out is None:
        raise ValueError('--out is required.')

    log = pr.Logger(args.out +'.log')
    try:
        defaults = vars(parser.parse_args(''))
        opts = vars(args)
        non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]
        header = 'Call: \n'
        header += './get_h2.py \\\n'
        options = ['--'+x.replace('_','-')+' '+str(opts[x])+' \\' for x in non_defaults]
        header += '\n'.join(options).replace('True','').replace('False','')
        header = header[0:-1]+'\n'
        log.log(header)
        log.log('Beginning analysis at {T}'.format(T=time.ctime()))
        start_time = time.time()

        # Determine heritability estimate.
        if (args.h2) and (args.ref_ld or args.ref_ld_chr) and (args.w_ld or args.w_ld_chr):
            if args.ref_ld and args.ref_ld_chr:
                raise ValueError('Cannot set both --ref-ld and --ref-ld-chr.')
            if args.w_ld and args.w_ld_chr:
                raise ValueError('Cannot set both --w-ld and --w-ld-chr.')
            if (args.samp_prev is not None) != (args.pop_prev is not None):
                raise ValueError('Must set both or neither of --samp-prev and --pop-prev.')

            args.frqfile = None
            args.frqfile_chr = None

            sumstats.estimate_h2(args, log)

        # Bad flags.
        else:
            print(header)
            print('Error: no analysis selected.')
            print('get_h2.py -h describes options.')
    except Exception:
        ex_type, ex, tb = sys.exc_info()
        log.log(traceback.format_exc(ex))
        raise
    finally:
        log.log('Analysis finished at {T}'.format(T=time.ctime()))
        time_elapsed = round(time.time()-start_time,2)
        log.log('Total time elapsed: {T}'.format(T=pr.sec_to_str(time_elapsed)))
