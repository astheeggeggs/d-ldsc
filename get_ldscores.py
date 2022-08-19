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

def __filter__(fname, noun, verb, merge_obj):
    merged_list = None
    if fname:
        f = lambda x,n: x.format(noun=noun, verb=verb, fname=fname, num=n)
        x = ps.FilterFile(fname)
        c = 'Read list of {num} {noun} to {verb} from {fname}'
        print(f(c, len(x.IDList)))
        merged_list = merge_obj.loj(x.IDList)
        len_merged_list = len(merged_list)
        if len_merged_list > 0:
            c = 'After merging, {num} {noun} remain'
            print(f(c, len_merged_list))
        else:
            error_msg = 'No {noun} retained for analysis'
            raise ValueError(f(error_msg, 0))

        return merged_list

def get_ldscores(args, log):

    if args.bfile:
        snp_file, snp_obj = args.bfile+'.bim', ps.PlinkBIMFile
        ind_file, ind_obj = args.bfile+'.fam', ps.PlinkFAMFile
        array_file, array_obj = args.bfile+'.bed', ld.PlinkBEDFile

    # Read .bim/snp.
    array_snps = snp_obj(snp_file)
    m = len(array_snps.IDList)
    log.log('Read list of {m} SNPs from {f}'.format(m=m, f=snp_file))

    if args.extract is not None:
        keep_snps = __filter__(args.extract, 'SNPs', 'include', array_snps)
    else:
        keep_snps = None

    # Read .fam.
    array_indivs = ind_obj(ind_file)
    n = len(array_indivs.IDList)
    log.log('Read list of {n} individuals from {f}'.format(n=n, f=ind_file))
    
    keep_indivs = None

    # Read genotype array.
    log.log('Reading genotypes from {fname}'.format(fname=array_file))
    geno_array = array_obj(array_file, n, array_snps, keep_snps=keep_snps,
        keep_indivs=keep_indivs, mafMin=args.maf)

    # Determine block widths.
    x = np.array((args.ld_wind_snps, args.ld_wind_kb, args.ld_wind_cm), dtype=bool)
    if np.sum(x) != 1:
        raise ValueError('Must specify exactly one --ld-wind option')

    if args.ld_wind_snps:
        max_dist = args.ld_wind_snps
        coords = np.array(range(geno_array.m))
    elif args.ld_wind_kb:
        max_dist = args.ld_wind_kb * 1000
        coords = np.array(array_snps.df['BP'])[geno_array.kept_snps]
    elif args.ld_wind_cm:
        max_dist = args.ld_wind_cm
        coords = np.array(array_snps.df['CM'])[geno_array.kept_snps]

    block_left = ld.getBlockLefts(coords, max_dist)
    if block_left[len(block_left)-1] == 0 and not args.yes_really:
        error_msg = 'Do you really want to compute whole-chromosome LD Score? If so, set the '
        error_msg += '--yes-really flag (warning: it will use a lot of time / memory)'
        raise ValueError(error_msg)

    scale_suffix = ''
    if args.pq_exp is not None:
        log.log('Computing LD with pq^{S}.'.format(S=args.pq_exp))
        msg = 'Note that LD Scores with pq raised to a nonzero power are '
        msg += 'not directly comparable to normal LD Scores.'
        log.log(msg)
        scale_suffix = '_S{S}'.format(S=args.pq_exp)
        pq = np.matrix(geno_array.maf*(1-geno_array.maf)).reshape((geno_array.m, 1))
        pq = np.power(pq, args.pq_exp)
        annot_matrix = pq # Note that in the dominance case, we multiply by (pq)^2.
    else:
        annot_matrix = None

    df = pd.DataFrame.from_records(np.c_[geno_array.df])
    df.columns = geno_array.colnames
    out_fname = args.out + '.ldscore'

    log.log('Estimating LD Scores.')

    if args.additive:
        log.log('Additive LD scores...')
        # Note, in the below, annot_matrix is None, unless we specify pq_exp.
        lN_A = geno_array.ldScoreVarBlocks_add(block_left, args.chunk_size, annot=annot_matrix)
        df['L2_A{scale_suffix}'.format(scale_suffix=scale_suffix)] = lN_A.reshape(geno_array.m,1)

    if args.dominance:
        log.log('Dominance LD scores...')
        geno_array._currentSNP = 0
        # Note, in the below, annot_matrix is None, unless we specify pq_exp.
        lN_D = geno_array.ldScoreVarBlocks_dom(block_left, args.chunk_size, annot=annot_matrix)
        df['L2_D{scale_suffix}'.format(scale_suffix=scale_suffix)] = lN_D.reshape(geno_array.m,1)

    log.log('Writing LD Scores for {N} SNPs to {f}'.format(f=out_fname, N=len(df)))
    df.drop(['CM','MAF'], axis=1).to_csv(out_fname, sep='\t', header=True, index=False,
        float_format='%.3f')
    
    M = []
    M_5_50 = []

    if args.additive:
        M += [geno_array.m]
        M_5_50 += [np.sum(geno_array.maf > 0.05)]
    if args.dominance:
        M += [geno_array.m]
        M_5_50 += [np.sum(geno_array.maf > 0.05)]

    # print .M
    fout_M = open(args.out +'.M','w')
    # print('\t'.join(map(str,M)), file=fout_M)
    print >>fout_M, '\t'.join(map(str,M))
    fout_M.close()

    # print .M_5_50
    fout_M_5_50 = open(args.out +'.M_5_50','w')
    # print('\t'.join(map(str,M_5_50)), file=fout_M_5_50)
    print >>fout_M_5_50, '\t'.join(map(str,M_5_50))
    fout_M_5_50.close()

    # Print LD Score summary.
    pd.set_option('display.max_rows', 200)
    log.log('\nSummary of LD Scores in {F}'.format(F=out_fname))
    t = df.iloc[:,4:].describe()
    print(t)
    print(str(t))
    log.log(str(t.iloc[1:,:]))

    # Print NaN instead of weird errors.
    np.seterr(divide='ignore', invalid='ignore')  
    
    # Print correlation matrix including all LD Scores and sample MAF.
    log.log('')
    log.log('MAF/LD Score Correlation Matrix')
    log.log(str(df.iloc[:,4:].corr()))

    np.seterr(divide='raise', invalid='raise')

parser = argparse.ArgumentParser()
parser.add_argument('--out', default='ldsc', type=str,
    help='Output filename prefix. If --out is not set, LDSC will use ldsc as the '
    'defualt output filename prefix.')
# Basic LD Score Estimation Flags'
parser.add_argument('--bfile', default=None, type=str,
    help='Prefix for Plink .bed/.bim/.fam file')
parser.add_argument('--extract', default=None, type=str,
    help='File with SNPs to include in LD Score estimation. '
    'The file should contain one SNP ID per row.')
parser.add_argument('--additive', default=False, action='store_true',
    help='Generate additive LD scores.')
parser.add_argument('--dominance', default=False, action='store_true',
    help='Generate dominance LD scores.')
parser.add_argument('--maf', default=0, type=float,
    help='Minor allele frequency lower bound. Default is MAF > 0.')
# Filtering / Data Management for LD Score
parser.add_argument('--ld-wind-snps', default=None, type=int,
    help='Specify the window size to be used for estimating LD Scores in units of '
    '# of SNPs. You can only specify one --ld-wind-* option.')
parser.add_argument('--ld-wind-kb', default=None, type=float,
    help='Specify the window size to be used for estimating LD Scores in units of '
    'kilobase-pairs (kb). You can only specify one --ld-wind-* option.')
parser.add_argument('--ld-wind-cm', default=None, type=float,
    help='Specify the window size to be used for estimating LD Scores in units of '
    'centiMorgans (cM). You can only specify one --ld-wind-* option.')
# Flags you should almost never use
parser.add_argument('--chunk-size', default=50, type=int,
    help='Chunk size for LD Score calculation. Use the default.')
parser.add_argument('--yes-really', default=False, action='store_true',
    help='Yes, I really want to compute whole-chromosome LD Score.')
parser.add_argument('--pq-exp', default=None, type=float,
    help='Setting this flag causes LDSC to compute LD Scores with the given scale factor, '
    'in the additive case this results in \ell^A_j := \sum_k (p_k(1-p_k))^a r_A^2_{jk}'
    'in the dominance case this results in \ell^D_j := \sum k (p_k(1-p_k))^{2a} r_D^2_{jk}, where '
    'p_k denotes the MAF of SNP j and a is the argument to --pq-exp.')

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
        header += './get_ldscores.py \\\n'
        options = ['--'+x.replace('_','-')+' '+str(opts[x])+' \\' for x in non_defaults]
        header += '\n'.join(options).replace('True','').replace('False','')
        header = header[0:-1]+'\n'
        log.log(header)
        log.log('Beginning analysis at {T}'.format(T=time.ctime()))
        start_time = time.time()
        
        if args.bfile is not None:
            get_ldscores(args, log)

        # Bad flags.
        else:
            print(header)
            print('Error: no analysis selected.')
            print('get_ldscores.py -h describes options.')
    except Exception:
        ex_type, ex, tb = sys.exc_info()
        log.log( traceback.format_exc(ex) )
        raise
    finally:
        log.log('Analysis finished at {T}'.format(T=time.ctime()) )
        time_elapsed = round(time.time()-start_time,2)
        log.log('Total time elapsed: {T}'.format(T=pr.sec_to_str(time_elapsed)))
