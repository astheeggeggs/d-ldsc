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

try:
    x = pd.DataFrame({'A': [1, 2, 3]})
    x.drop_duplicates(subset='A')
except TypeError:
    raise ImportError('LDSC requires pandas version > 0.15.2')

def get_sumstats(args, log):

    if args.bfile:
        snp_file, snp_obj = args.bfile + '.bim', ps.PlinkBIMFile
        ind_file, ind_obj = args.bfile + '.fam', ps.PlinkFAMFile
        array_file, array_obj = args.bfile + '.bed', ld.PlinkBEDFile

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

    if args.phifile:
        y_cols = list(range(2,args.ny+2))
        y_names = [''.join(str(i) for i in z) for z in zip(np.tile('y', args.ny), range(args.ny))]
        PhenoFile = ps.__ID_List_Factory__(['IID', 'FID'] + y_names, 0, '.phi',
            usecols=[0, 1] + y_cols)
        pheno_file, pheno_obj = args.phifile, PhenoFile

        y = pheno_obj(pheno_file)
        
        # Check that this matches the .fam file.
        if any(y.df['FID'] != array_indivs.df['FID']) or any(y.df['IID'] != array_indivs.df['IID']):
            raise ValueError('The provided .phi file does not match the provided .fam file.')
       
        # Create a vector of phenotypes that is then passed to the function.
        y = np.array(y.df[y_names])
        # Ensure that the phenotype is normalised.
        y = (y - np.mean(y, axis=0)) / np.std(y, axis=0)

        # These are the additive betas.
        X_A, X_D = geno_array.nextSNPs(geno_array.m)

        sumstats_df = pd.DataFrame.from_records(np.c_[geno_array.df])
        sumstats_df.columns = geno_array.colnames

        if args.additive:
            # Determine the additive betas.
            beta_A = np.dot(y.T, X_A) / geno_array.n
            Z_A = np.multiply(np.sqrt(geno_array.n), beta_A)
            if args.ny > 1:
                Z_A_index = [''.join(str(i) for i in z) for z in zip(np.tile('Z_A', args.ny), range(args.ny))]
            else:
                Z_A_index = ['Z_A']

            sumstats_df = sumstats_df.reindex(columns=list(sumstats_df.columns.values) + Z_A_index)
            sumstats_df[Z_A_index] = Z_A.T

        if args.dominance:
            # Determine the dominance betas.
            beta_D = np.dot(y.T, X_D) / geno_array.n
            Z_D = np.multiply(np.sqrt(geno_array.n), beta_D)
            if args.ny > 1:
                Z_D_index = [''.join(str(i) for i in z) for z in zip(np.tile('Z_D', args.ny), range(args.ny))]
            else:
                Z_D_index = ['Z_D']
            sumstats_df = sumstats_df.reindex(columns=list(sumstats_df.columns.values) + Z_D_index)
            sumstats_df[Z_D_index] = Z_D.T

        if args.gxe:
            # First, read in the covariate file.
            C_cols = list(range(2, args.nC + 2))
            C_names = [''.join(str(i) for i in z) for z in zip(np.tile('C', args.nC), range(args.nC))]
            CovFile = ps.__ID_List_Factory__(['IID', 'FID'] + C_names, 0, '.covariate',
                usecols=[0, 1] + C_cols)
            if args.gxefile is None:
                raise Exception("Covariate file not passed.")
            cov_file, cov_obj = args.gxefile, CovFile

            C = cov_obj(cov_file)

            # Check that this matches the .fam file
            if any(C.df['FID'] != array_indivs.df['FID']) or any(C.df['IID'] != array_indivs.df['IID']):
                raise ValueError('The provided .cov file does not match the provided .fam file.')

            C = np.array(C.df[C_names])
            # Ensure that the covariate is normalised.
            C = (C - np.mean(C, axis=0)) / np.std(C, axis=0)
            
            # Determine the gene by environment interaction betas.
            beta_AC = np.dot((C*y).T, X_A) / geno_array.n
            Z_AC = np.multiply(np.sqrt(geno_array.n), beta_AC)
            if args.ny > 1:
                Z_AC_index = [''.join(str(i) for i in z) for z in zip(np.tile('Z_AC', args.ny), range(args.ny))]
            else:
                Z_AC_index = ['Z_AC']

            sumstats_df = sumstats_df.reindex(columns=list(sumstats_df.columns.values) + Z_AC_index)
            sumstats_df[Z_AC_index] = Z_AC.T
            
        # DEV: For non simulation studies, this will have to be general.
        N = np.tile(geno_array.n, geno_array.m)
        sumstats_df['N'] = N
        sumstats_outname = args.out + '.sumstats'
        # Save the summary statistics file.
        sumstats_df.drop(['CHR', 'CM', 'BP', 'MAF'], axis=1).to_csv(sumstats_outname, sep='\t', header=True, index=False)

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
parser.add_argument('--phifile', default=None, type=str,
    help='Name of the phenotype file used to generate the beta summary statistics.')
parser.add_argument('--additive', default=False, action='store_true',
    help='Determine the additive beta summary statistics from the passed genotype and '
    'phenotype files.')
parser.add_argument('--dominance', default=False, action='store_true',
    help='Determine the dominance beta summary statistics from the passed genotype and '
    'phenotype files.')
parser.add_argument('--gxe', default=False, action='store_true',
    help='Determine the gene by environment summary statistics from the passed genotype '
    'and phenotype files.')
parser.add_argument('--gxefile', default=None, type=str,
    help='Name of covariate file in - in the same format as the phenotype file.')
parser.add_argument('--ny', default=1, type=int,
    help='Number of measured phenotypes in the .phi file.')
parser.add_argument('--nC', default=1, type=int,
    help='Number of covariates in the .cov file.')

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
        header += './get_sumstats.py \\\n'
        options = ['--'+x.replace('_','-')+' '+str(opts[x])+' \\' for x in non_defaults]
        header += '\n'.join(options).replace('True','').replace('False','')
        header = header[0:-1]+'\n'
        log.log(header)
        log.log('Beginning analysis at {T}'.format(T=time.ctime()))
        start_time = time.time()

        if args.bfile is not None:
            # Get the LD scores.
            get_sumstats(args, log)

        # Bad flags.
        else:
            print(header)
            print('Error: no analysis selected.')
            print('get_sumstats.py -h describes options.')
    except Exception:
        ex_type, ex, tb = sys.exc_info()
        log.log( traceback.format_exc(ex) )
        raise
    finally:
        log.log('Analysis finished at {T}'.format(T=time.ctime()) )
        time_elapsed = round(time.time()-start_time,2)
        log.log('Total time elapsed: {T}'.format(T=pr.sec_to_str(time_elapsed)))
