# This module deals with getting all the data needed for LD Score regression from files
# into memory and checking that the input makes sense. There is no math here. LD Score
# regression is implemented in the regressions module.

from __future__ import division
import numpy as np
import pandas as pd
from scipy.stats import norm
import itertools as it
import src.parse as ps
import src.regressions as reg
import sys
import traceback
import copy
import os
import re

# Complementary bases.
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

# Bases.
BASES = COMPLEMENT.keys()

# True iff strand ambiguous.
STRAND_AMBIGUOUS = {''.join(x): x[0] == COMPLEMENT[x[1]]
                    for x in it.product(BASES, BASES)
                    if x[0] != x[1]}

# SNPS we want to keep (pairs of alleles).
VALID_SNPS = {x for x in map(lambda y: ''.join(y), it.product(BASES, BASES))
              if x[0] != x[1] and not STRAND_AMBIGUOUS[x]}

# T iff SNP 1 has the same alleles as SNP 2 (allowing for strand or ref allele flip).
MATCH_ALLELES = {x for x in map(lambda y: ''.join(y), it.product(VALID_SNPS, VALID_SNPS))
                 # strand and ref match
                 if ((x[0] == x[2]) and (x[1] == x[3])) or
                 # ref match, strand flip
                 ((x[0] == COMPLEMENT[x[2]]) and (x[1] == COMPLEMENT[x[3]])) or
                 # ref flip, strand match
                 ((x[0] == x[3]) and (x[1] == x[2])) or
                 ((x[0] == COMPLEMENT[x[3]]) and (x[1] == COMPLEMENT[x[2]]))}  # strand and ref flip

# T iff SNP 1 has the same alleles as SNP 2 w/ ref allele flip.
FLIP_ALLELES = {''.join(x):
                ((x[0] == x[3]) and (x[1] == x[2])) or  # strand match
                # strand flip
                ((x[0] == COMPLEMENT[x[3]]) and (x[1] == COMPLEMENT[x[2]]))
                for x in MATCH_ALLELES}

def _splitp(fstr):
    flist = fstr.split(',')
    flist = [os.path.expanduser(os.path.expandvars(x)) for x in flist]
    return flist

def smart_merge(x, y):
    # Check if SNP columns are equal. If so, save time by using concat instead of merge.
    if len(x) == len(y) and (x.index == y.index).all() and (x.SNP == y.SNP).all():
        x = x.reset_index(drop=True)
        y = y.reset_index(drop=True).drop('SNP', 1)
        out = pd.concat([x, y], axis=1)
    else:
        out = pd.merge(x, y, how='inner', on='SNP')
    return out

def _read_ref_ld(args, log):
    # Read reference LD Scores.
    ref_ld = _read_chr_split_files(args.N_CHR, args.ref_ld_chr, args.ref_ld, log,
                                   'reference panel LD Score', ps.ldscore_fromlist)
    log.log('Read reference panel LD Scores for {N} SNPs.'.format(N=len(ref_ld)))
    return ref_ld

def _read_M(args, log, n_annot):
    # Read M (--M, --M-file, etc).
    if args.M:
        try:
            M_annot = [float(x) for x in _splitp(args.M)]
        except ValueError as e:
            raise ValueError('Could not cast --M to float: ' + str(e.args))
    else:
        if args.ref_ld:
            M_annot = ps.M_fromlist(
                _splitp(args.ref_ld), common=(not args.not_M_5_50))
        elif args.ref_ld_chr:
            M_annot = ps.M_fromlist(
                _splitp(args.ref_ld_chr), args.N_CHR, common=(not args.not_M_5_50))

    try:
        M_annot = np.array(M_annot).reshape((1, n_annot))
    except ValueError as e:
        print(M_annot)
        print(n_annot)
        raise ValueError(
            '# terms in --M must match # of LD Scores in --ref-ld.\n' + str(e.args))

    return M_annot

def _read_w_ld(args, log):
    # Read regression SNP LD.
    if (args.w_ld and ',' in args.w_ld) or (args.w_ld_chr and ',' in args.w_ld_chr):
        raise ValueError(
            '--w-ld must point to a single fileset (no commas allowed).')
    w_ld = _read_chr_split_files(args.N_CHR, args.w_ld_chr, args.w_ld, log,
                                 'regression weight LD Score', ps.ldscore_fromlist)
    col_names = ['SNP']
    
    if not (args.additive or args.additive_orig):
        w_ld = w_ld.drop(['L2_A_0'], axis=1)
    if not args.dominance:
        w_ld = w_ld.drop(['L2_D_0'], axis=1)

    for i in range(len(w_ld.columns)):
        if w_ld.columns[i].find('L2_') != -1:
            col_names += [w_ld.columns[i] + '_weights']

    w_ld.columns = col_names  # prevent colname conflicts w/ ref ld

    log.log('Read regression weight LD Scores for {N} SNPs.'.format(N=len(w_ld)))
    return w_ld

def _read_chr_split_files(N_CHR, chr_arg, not_chr_arg, log, noun, parsefunc, **kwargs):
    # Read files split across 22 chromosomes (annot, ref_ld, w_ld).
    try:
        if not_chr_arg:
            log.log('Reading {N} from {F} ...'.format(F=not_chr_arg, N=noun))
            out = parsefunc(_splitp(not_chr_arg), **kwargs)
        elif chr_arg:
            f = ps.sub_chr(chr_arg, '[1-22]')
            log.log('Reading {N} from {F} ...'.format(F=f, N=noun))
            out = parsefunc(_splitp(chr_arg), N_CHR, **kwargs)
    except ValueError as e:
        log.log('Error parsing {N}.'.format(N=noun))
        raise e

    return out

def _read_sumstats(args, log, fh, alleles=False, dropna=False):
    # Parse summary statistics.
    log.log('Reading summary statistics from {S} ...'.format(S=fh))
    sumstats, n_codings = ps.sumstats(fh,
        args.additive_orig, args.additive, args.dominance, args.gxe,
        alleles=alleles, dropna=dropna)
    log_msg = 'Read summary statistics for {N} SNPs.'
    log.log(log_msg.format(N=len(sumstats)))
    m = len(sumstats)
    sumstats = sumstats.drop_duplicates(subset='SNP')
    if m > len(sumstats):
        log.log('Dropped {M} SNPs with duplicated rs numbers.'.format(M=m - len(sumstats)))

    return sumstats, n_codings

def _check_variance(log, M_annot, ref_ld):
    # Remove zero-variance LD Scores.
    ii = ref_ld.ix[:, 1:].var() == 0 # NB there is a SNP column here
    if ii.all():
        raise ValueError('All LD Scores have zero variance.')
    else:
        log.log('Removing partitioned LD Scores with zero variance.')
        ii_snp = np.array([True] + list(~ii))
        ii_m = np.array(~ii)
        ref_ld = ref_ld.ix[:, ii_snp]
        M_annot = M_annot[:, ii_m]

    return M_annot, ref_ld, ii

def _merge_and_log(ld, sumstats, noun, log):
    # Wrap smart merge with log messages about # of SNPs.
    sumstats = smart_merge(ld, sumstats)
    msg = 'After merging with {F}, {N} SNPs remain.'
    if len(sumstats) == 0:
        raise ValueError(msg.format(N=len(sumstats), F=noun))
    else:
        log.log(msg.format(N=len(sumstats), F=noun))

    return sumstats

def _read_ld_sumstats(args, log, fh, alleles=False, dropna=True):
    sumstats, n_codings = _read_sumstats(args, log, fh, alleles=alleles, dropna=dropna)
    ref_ld = _read_ref_ld(args, log)
    n_annot = len(ref_ld.columns) - 1
    M_annot = _read_M(args, log, n_annot)
    M_annot, ref_ld, novar_cols = _check_variance(log, M_annot, ref_ld)
    w_ld = _read_w_ld(args, log)
    sumstats = _merge_and_log(ref_ld, sumstats, 'reference panel LD', log)
    sumstats = _merge_and_log(sumstats, w_ld, 'regression SNP LD', log)
    w_ld_cnames = sumstats.columns[-(len(w_ld.columns)-1):]
    ref_ld_cnames = ref_ld.columns[1:len(ref_ld.columns)]
    return M_annot, w_ld_cnames, ref_ld_cnames, sumstats, novar_cols

def _check_ld_condnum(args, log, ref_ld):
    # Check condition number of LD Score matrix.
    if len(ref_ld.shape) >= 2:
        cond_num = int(np.linalg.cond(ref_ld))
        if cond_num > 100000:
            warn = "WARNING: LD Score matrix condition number is {C}. "
            warn += "Remove collinear LD Scores. "
            raise ValueError(warn.format(C=cond_num))

def _warn_length(log, sumstats):
    if len(sumstats) < 200000:
        log.log(
            'WARNING: number of SNPs less than 200k; this is almost always bad.')

def estimate_h2(args, log):
    # Estimate h2.
    args = copy.deepcopy(args)
    if args.samp_prev is not None and args.pop_prev is not None:
        args.samp_prev, args.pop_prev = map(
            float, [args.samp_prev, args.pop_prev])
    
    intercept_h2 = None
    if args.no_intercept: intercept_h2 = 1

    M_annot, w_ld_cname, ref_ld_cnames, sumstats, novar_cols = _read_ld_sumstats(
        args, log, args.h2)
    ref_ld = np.array(sumstats[ref_ld_cnames])
    
    # Perform some checks and throw some warnings.
    _check_ld_condnum(args, log, ref_ld_cnames)
    _warn_length(log, sumstats)
    
    n_snp = len(sumstats)
    n_blocks = min(n_snp, args.n_blocks)
    n_annot = len(ref_ld_cnames)
    chisq_max = args.chisq_max

    old_weights = False

    if args.two_step is None and intercept_h2 is None:
        args.two_step = 30
    else:
        old_weights = True
        if args.chisq_max is None:
            chisq_max = max(0.001*sumstats.N.max(), 80)

    old_weights = True
    args.two_step = None

    transpose = lambda x: np.array(x).reshape((n_snp, 1))

    if args.additive or args.additive_orig:
        chisq_A = transpose(sumstats.Z_A**2)
    if args.dominance:
        chisq_D = transpose(sumstats.Z_D**2)
    if args.gxe:
        chisq_AC = transpose(sumstats.Z_AC**2)

    if args.additive or args.additive_orig:
        where_additive_weights = []
        for i in range(len(w_ld_cname)):
            if w_ld_cname[i].find('_A_') != -1:
                where_additive_weights += [i]
        where_additive_refs = []
        for i in range(len(ref_ld_cnames)):
            if ref_ld_cnames[i].find('_A_') != -1:
                where_additive_refs += [i]

    if args.dominance:
        where_dominance_weights = []
        for i in range(len(w_ld_cname)):
            if w_ld_cname[i].find('_D_') != -1:
                where_dominance_weights += [i]

        where_dominance_refs = []
        for i in range(len(ref_ld_cnames)):
            if ref_ld_cnames[i].find('_D_') != -1:
                where_dominance_refs += [i]

    # First let's do the additive component.
    if args.additive or args.additive_orig:

        if chisq_max is not None:
            ii_A = np.ravel(chisq_A < chisq_max)
            sumstats_A = sumstats.ix[ii_A, :]
            log.log('Removed {M} SNPs with chi^2 > {C} ({N} SNPs remain)'.format(
                    C=chisq_max, N=np.sum(ii_A), M=n_snp-np.sum(ii_A)))
            n_snp = np.sum(ii_A)  # lambdas are late-binding, so this works
            ref_ld = np.array(sumstats_A[ref_ld_cnames])
            chisq_A = chisq_A[ii_A].reshape((n_snp, 1))
        else:
            sumstats_A = sumstats

        if args.two_step is not None:
            log.log('Using two-step estimator with cutoff at {M}.'.format(M=args.two_step))
        
        log.log('Sanity check: running weighted regression of chi^2_A against {R}, using weights defined by {W}.'.format(
            R=ref_ld_cnames[where_additive_refs], W=w_ld_cname[where_additive_weights]))
        hsqhat_A = reg.Hsq(chisq_A,
                         transpose(ref_ld[:,where_additive_refs]),
                         transpose(sumstats_A[w_ld_cname[where_additive_weights]]), transpose(sumstats_A.N),
                         M_annot[0,0].reshape((1,1)), n_blocks=n_blocks, intercept=intercept_h2,
                         twostep=args.two_step, old_weights=old_weights)

        log.log(hsqhat_A.summary(ref_ld_cnames, P=args.samp_prev, K=args.pop_prev))
    
    # Then let's do the dominance component:
    if args.dominance:
        if chisq_max is not None:
            ii_D = np.ravel(chisq_D < chisq_max)
            sumstats_D = sumstats.ix[ii_D, :]
            log.log('Removed {M} SNPs with chi^2 > {C} ({N} SNPs remain)'.format(
                    C=chisq_max, N=np.sum(ii_D), M=n_snp-np.sum(ii_D)))
            n_snp = np.sum(ii_D)  # lambdas are late-binding, so this works
            ref_ld = np.array(sumstats_D[ref_ld_cnames])
            chisq_D = chisq_D[ii_D].reshape((n_snp, 1))
        else:
            sumstats_D = sumstats

        log.log('Sanity check: running weighted regression of chi^2_D against {R}, using weights defined by {W}.'.format(
            R=ref_ld_cnames[where_dominance_refs], W=w_ld_cname[where_dominance_weights]))
        hsqhat_D = reg.Hsq(chisq_D,
                         transpose(ref_ld[:,where_dominance_refs]),
                         transpose(sumstats_D[w_ld_cname[where_dominance_weights]]), transpose(sumstats_D.N),
                         M_annot[0,1].reshape((1,1)), n_blocks=n_blocks, intercept=intercept_h2,
                         twostep=args.two_step, old_weights=old_weights)

        log.log(hsqhat_D.summary(ref_ld_cnames, P=args.samp_prev, K=args.pop_prev))

    # Finally let's do a gene by environment component:
    if args.gxe:
        if chisq_max is not None:
            ii_AC = np.ravel(chisq_AC < chisq_max)
            sumstats_AC = sumstats.ix[ii_AC, :]
            log.log('Removed {M} SNPs with chi^2 > {C} ({N} SNPs remain)'.format(
                    C=chisq_max, N=np.sum(ii_AC), M=n_snp-np.sum(ii_AC)))
            n_snp = np.sum(ii_AC)  # lambdas are late-binding, so this works
            ref_ld = np.array(sumstats_AC[ref_ld_cnames])
            chisq_AC = chisq_AC[ii_AC].reshape((n_snp, 1))
        else:
            sumstats_AC = sumstats

        # Old version, that doesn't account for kurtosis.
        hsqhat_AC_old = reg.Hsq(chisq_AC,
                         transpose(ref_ld[:,where_additive_refs]),
                         transpose(sumstats_AC[w_ld_cname[where_additive_weights]]), transpose(sumstats_AC.N),
                         M_annot[0,0].reshape((1,1)), n_blocks=n_blocks, intercept=intercept_h2,
                         twostep=args.two_step, old_weights=old_weights)

        log.log(hsqhat_AC_old.summary(ref_ld_cnames, P=args.samp_prev, K=args.pop_prev))

        # Need to get the kurtosis in here.
        # Should provide option to include the kurtosis.
        # If the user doesn't have the average kurtosis, then can set the default to be that of 
        # a normal distribution if the phenotype is continuous, or based on the case control
        # proportions if it's a boolean variable.

        # Need to figure out what the weights should be.

        if args.average_cov_kurtosis is None:
            if args.cov_case_control is False and args.cov_case_prev is None:
                # Assume that the variable is normally distributed.
                log.log("Assuming that this covariate is N(0,1) distributed.")
                args.average_cov_kurtosis = 3
            else:
                if args.cov_case_prev is not None:
                    # Work out what we expect that the kurtosis will be 
                    # under a Bernoulli variable.
                    log.log("Evaluating expected average kurtosis assuming the supplied case/control proportions for the covariate")
                    args.average_cov_kurtosis = (1 - 3 * args.cov_case_prev * (1 - args.cov_case_prev)) / (args.cov_case_prev * (1 - args.cov_case_prev))
                    log.log("Expected Kurtosis is {K}".format(K=args.average_cov_kurtosis))
                else:
                    log.log("Assuming a normalised case/control covariate with equal cases and controls")
                    args.average_cov_kurtosis = 1

        hsqhat_AC = reg.Hsq(chisq_AC,
                         (transpose(ref_ld[:,where_additive_refs]) + args.average_cov_kurtosis - 1),
                         transpose(sumstats_AC[w_ld_cname[where_additive_weights]]), transpose(sumstats_AC.N),
                         M_annot[0,0].reshape((1,1)), n_blocks=n_blocks, intercept=intercept_h2,
                         twostep=args.two_step, old_weights=old_weights)

        log.log(hsqhat_AC.summary(ref_ld_cnames, P=args.samp_prev, K=args.pop_prev))

    if args.write_h2:
        if args.samp_prev is not None and args.pop_prev is not None:
            c = reg.h2_obs_to_liab(1, args.samp_prev, args.pop_prev)
        else:
            c = 1

        h2 = []

        if args.pheno_name is not None:
            h2.extend([args.pheno_name])

        if args.additive or args.additive_orig:
            h2.extend(hsqhat_A.to_print_to_file(c))

        if args.dominance:
            h2.extend(hsqhat_D.to_print_to_file(c))

        if args.gxe:
            h2.extend(hsqhat_AC.to_print_to_file(c))

        print(h2)
        # Output columns order.
        # [ 
        # 'phenotype',
        # 'mean_chi2',
        # 'lambdaGC',
        # 'intercept',
        # 'intercept_se',
        # 'intercept_z',
        # 'intercept_p',
        # 'ratio',
        # 'ratio_se',
        # 'h2_observed',
        # 'h2_observed_se',
        # 'h2_liability',
        # 'h2_liability_se',
        # 'h2_z',
        # 'h2_p'
        # ]

        fout_h2 = open(args.out +'.h2','w')
        fout_h2.write('\t'.join(map(str,h2)) + '\n')
        fout_h2.close()

    return 1
