#!/usr/bin/env python

from __future__ import division
import src.ldscoregxe as ld
import src.parse as ps
import src.printing as pr
import numpy as np
import pandas as pd
from subprocess import call
from itertools import product
import time, sys, traceback, argparse
# import matplotlib.pyplot as plt
import scipy.stats as ss

try:
    x = pd.DataFrame({'A': [1, 2, 3]})
    x.drop_duplicates(subset='A')
except TypeError:
    raise ImportError('LDSC requires pandas version > 0.15.2')

def get_y(args, log):

    if args.bfile:
        args.N_CHR = 1

    h2 = args.h2_A + args.h2_D + args.h2_AC
    
    for sim in range(args.n_sims):
        if args.bfile:
            snp_file, snp_obj = args.bfile+'.bim', ps.PlinkBIMFile
            ind_file, ind_obj = args.bfile+'.fam', ps.PlinkFAMFile
            array_file, array_obj = args.bfile+'.bed', ld.PlinkBEDFile
            
        else:
            raise ValueError("Must pass sequence data.")

        # Read bim/snp
        array_snps = snp_obj(snp_file)
        m = len(array_snps.IDList)
        log.log('Read list of {m} SNPs from {f}'.format(m=m, f=snp_file))

        keep_snps = None

        # Read fam
        array_indivs = ind_obj(ind_file)
        n = len(array_indivs.IDList)
        log.log('Read list of {n} individuals from {f}'.format(n=n, f=ind_file))
        
        keep_indivs = None

        # Read genotype array
        log.log('Reading genotypes from {fname}'.format(fname=array_file))
        geno_array = array_obj(array_file, n, array_snps, keep_snps=keep_snps,
            keep_indivs=keep_indivs, mafMin=args.maf)

        if args.bfile: m_tot = geno_array.m
        
        if sim==0:
            # Initiate the matrices to be filled with phenotypes.
            phi_sims = np.zeros(((geno_array.n, args.n_sims)))
            if args.h2_AC > 0:
                C_sims = np.zeros((geno_array.n, args.n_sims))

        phi, beta_A, beta_D, beta_AC, C = geno_array.get_pheno(args, log)
        print(phi[0:10])

        # Add noise term and normalise
        phi_sims[:,sim] = geno_array.add_noise(phi, h2)
        if args.h2_AC > 0:
            C_sims[:,sim] = C.reshape(n)

    # Write the phis to file.
    sim_num = range(args.n_sims)

    phi_df = pd.DataFrame.from_records(np.c_[array_indivs.df, phi_sims])
    phi = np.tile('phi', args.n_sims)
    phi_sim_names = [''.join(str(i) for i in z) for z in zip(phi,sim_num)]
    phi_df.columns = array_indivs.__colnames__ + phi_sim_names
    phi_outname = args.out + '.phi'
    phi_df.to_csv(phi_outname, sep="\t", header=False, index=False)

    if args.h2_AC > 0:
        C_df = pd.DataFrame.from_records(np.c_[array_indivs.df, C_sims])
        C = np.tile('C', args.n_sims)
        C_sim_names = [''.join(str(i) for i in z) for z in zip(C,sim_num)]
        C_df.columns = array_indivs.__colnames__ + C_sim_names
        C_outname = args.out + '.covariate'
        C_df.to_csv(C_outname, sep="\t", header=False, index=False)

    np.seterr(divide='raise', invalid='raise')

parser = argparse.ArgumentParser()
parser.add_argument('--out', default='ldsc', type=str,
    help='Output filename prefix. If --out is not set, LDSC will use ldsc as the '
    'defualt output filename prefix.')
# Basic LD Score Estimation Flags'
parser.add_argument('--bfile', default=None, type=str,
    help='Prefix for Plink .bed/.bim/.fam file')
parser.add_argument('--maf', default=0, type=float,
    help='Minor allele frequency lower bound. Default is MAF > 0.')
# DEV: add in a flag for a covariate - for environment interaction.
parser.add_argument('--gxefile', default=None, type=str,
    help='Name of covariate file in - in the same format as the phenotype file.')
parser.add_argument('--beta-prop_A', default=1, type=float,
    help='The proportion of SNPs that contribute to the phenotype. [Default: 1].')
parser.add_argument('--beta-prop_D', default=1, type=float,
    help='The proportion of SNPs that contribute to the phenotype. [Default: 1].')
parser.add_argument('--beta-prop_AC', default=1, type=float,
    help='The proportion of SNPs that contribute to the phenotype. [Default: 1].')
parser.add_argument('--gxe-non-independent', default=False, action='store_true',
    help='Make the gxe and additive betas have dependency at a site.'
    'Currently I just have them as being the same, scaled by the input heritability')
parser.add_argument('--h2_A', default=0.4, type=float,
    help='The heritability for the simulation study. [Default: 0.4].')
parser.add_argument('--h2_D', default=0.4, type=float,
    help='The heritability for the simulation study. [Default: 0.4].')
parser.add_argument('--h2_AC', default=0.4, type=float,
    help='The heritability for the simulation study. [Default: 0.4].')
parser.add_argument('--C-bool', default=False, action='store_true',
    help='Is the environmental covariate that you consider boolean?')
parser.add_argument('--C-bool-p', default=0.5, type=float,
    help='Probability of "success" for a boolean covariate which is independent of the genetic data.')
parser.add_argument('--chunk-size', default=50, type=int,
    help='Chunk size for LD Score calculation. Use the default.')
parser.add_argument('--n-sims', default=1, type=int,
    help='Number of simulated sets of phenotypes for the same bfile.')
parser.add_argument('--N-CHR', default=1, type=int,
    help='Primarily for simulation studies and debugging, if the number of'
    ' `chromosomes` is not 22')

if __name__ == '__main__':

    args = parser.parse_args()
    if args.out is None:
        raise ValueError('--out is required.')

    log = pr.Logger(args.out+'.log')
    try:
        defaults = vars(parser.parse_args(''))
        opts = vars(args)
        non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]
        header = 'Call: \n'
        header += './get_y.py \\\n'
        options = ['--'+x.replace('_','-')+' '+str(opts[x])+' \\' for x in non_defaults]
        header += '\n'.join(options).replace('True','').replace('False','')
        header = header[0:-1]+'\n'
        log.log(header)
        log.log('Beginning analysis at {T}'.format(T=time.ctime()))
        start_time = time.time()

        if args.bfile is not None:
            # Get the phenotypes.
            get_y(args, log)

        # Bad flags.
        else:
            print(header)
            print('Error: no analysis selected.')
            print('get_y.py -h describes options.')
    except Exception:
        ex_type, ex, tb = sys.exc_info()
        log.log( traceback.format_exc(ex) )
        raise
    finally:
        log.log('Analysis finished at {T}'.format(T=time.ctime()) )
        time_elapsed = round(time.time()-start_time,2)
        log.log('Total time elapsed: {T}'.format(T=pr.sec_to_str(time_elapsed)))
