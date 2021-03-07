from __future__ import division
import numpy as np
import bitarray as ba
# import matplotlib.pyplot as plt
import scipy.stats as ss

def getBlockLefts(coords, max_dist):
    # Converts coordinates + max block length to the a list of coordinates of the leftmost
    # SNPs to be included in blocks.

    # Parameters
    # ----------
    # coords : array
    #     Array of coordinates. Must be sorted.
    # max_dist : float
    #     Maximum distance between SNPs included in the same window.

    # Returns
    # -------
    # block_left : 1D np.ndarray with same length as block_left
    #     block_left[j] :=  min{k | dist(j, k) < max_dist}.

    M = len(coords)
    j = 0
    block_left = np.zeros(M)
    for i in range(M):
        while j < M and abs(coords[j] - coords[i]) > max_dist:
            j += 1

        block_left[i] = j

    return block_left

class __GenotypeArrayInMemory__(object):
    # Parent class for various classes containing interfaces for files with genotype
    # matrices, e.g., plink .bed files, etc
    def __init__(self, fname, n, snp_list, keep_snps=None, keep_indivs=None, mafMin=None):
        self.m = len(snp_list.IDList)
        self.n = n
        self.keep_snps = keep_snps
        self.keep_indivs = keep_indivs
        self.df = np.array(snp_list.df[['CHR', 'SNP', 'BP', 'CM']])
        self.colnames = ['CHR', 'SNP', 'BP', 'CM']
        self.mafMin = mafMin if mafMin is not None else 0
        self._currentSNP = 0
        (self.nru, self.geno) = self.__read__(fname, self.m, n)
        # Filter individuals
        if keep_indivs is not None:
            keep_indivs = np.array(keep_indivs, dtype='int')
            if np.any(keep_indivs > self.n):
                raise ValueError('keep_indivs indices out of bounds')

            (self.geno, self.m, self.n) = self.__filter_indivs__(self.geno, keep_indivs, self.m,
                self.n)

            if self.n > 0:
                print('After filtering, {n} individuals remain'.format(n=self.n))
            else:
                raise ValueError('After filtering, no individuals remain')

        # Filter SNPs
        if keep_snps is not None:
            keep_snps = np.array(keep_snps, dtype='int')
            if np.any(keep_snps > self.m):  # if keep_snps is None, this returns False
                raise ValueError('keep_snps indices out of bounds')

        (self.geno, self.m, self.n, self.kept_snps, self.freq) = self.__filter_snps_maf__(
            self.geno, self.m, self.n, self.mafMin, keep_snps)

        if self.m > 0:
            print('After filtering, {m} SNPs remain'.format(m=self.m))
        else:
            raise ValueError('After filtering, no SNPs remain')

        self.df = self.df[self.kept_snps, :]
        self.maf = np.minimum(self.freq, np.ones(self.m)-self.freq)
        self.sqrtpq = np.sqrt(self.freq*(np.ones(self.m)-self.freq))
        self.df = np.c_[self.df, self.maf]
        self.colnames.append('MAF')

    def __read__(self, fname, m, n):
        raise NotImplementedError

    def __filter_indivs__(geno, keep_indivs, m, n):
        raise NotImplementedError

    def __filter_maf_(geno, m, n, maf):
        raise NotImplementedError

    def __l2_unbiased__(self, x, n):
        sq = np.square(x)
        return sq - 1/n

    def ldScoreVarBlocks_add(self, block_left, c):
        l2 = lambda x: self.__l2_unbiased__(x, self.n)
        return self.__corSumVarBlocks__(block_left, c, l2, self.nextSNPs_add)

    def ldScoreVarBlocks_dom(self, block_left, c):
        l2 = lambda x: self.__l2_unbiased__(x, self.n)
        return self.__corSumVarBlocks__(block_left, c, l2, self.nextSNPs_dom)

    # general methods for calculating sums of Pearson correlation coefficients
    def __corSumVarBlocks__(self, block_left, c, l2, snp_getter):
        # Parameters
        # ----------
        # block_left : np.ndarray with shape (M, )
        #     block_left[i] = index of leftmost SNP included in LD Score of SNP i.
        #     if c > 1, then only entries that are multiples of c are examined, and it is
        #     assumed that block_left[a*c+i] = block_left[a*c], except at
        #     the beginning of the chromosome where the 0th SNP is included in the window.

        # c : int
        #     Chunk size.
        # func : function
        #     Function to be applied to the genotype correlation matrix.
        #     lambda x: np.square(np.square(x)). For L1, lambda x: x.
        # snp_getter : function(int)
        #     The method to be used to get the next SNPs (normalized genotypes? Normalized
        #     genotypes with the minor allele as reference allele? etc)

        # Returns
        # -------
        # cor_sum : np.array with shape (M, 1)
        #     Estimates.

        m, n = self.m, self.n
        block_sizes = np.array(np.arange(m) - block_left)
        block_sizes = (np.ceil(block_sizes / c) * c).astype(int)

        cor_sum = np.zeros((m, 1))
        # b = index of first SNP for which SNP 0 is not included in LD Score
        b = np.nonzero(block_left > 0)
        if np.any(b):
            b = b[0][0]
        else:
            b = m
        b = int(np.ceil(b/c)*c)  # Round up to a multiple of c

        if b > m:
            c = 1
            b = m
        l_A = 0  # l_A := index of leftmost SNP in matrix A

        A = snp_getter(b)
        rfuncAB = np.zeros((b, c))
        rfuncBB = np.zeros((c, c))
        # Chunk inside of block
        for l_B in range(0, b, c):  # l_B := index of leftmost SNP in matrix B
            B = A[:, l_B:l_B+c]
            np.dot(A.T, B / n, out=rfuncAB)
            rfuncAB = l2(rfuncAB)
            cor_sum[l_A:l_A+b, :] += np.sum(rfuncAB, axis=1).reshape(b,1)

        # Chunk to right of block
        b0 = int(b)
        md = int(c*np.floor(m/c))
        end = md + 1 if md != m else md

        for l_B in range(b0, end, c):

            # Update the block
            old_b = b

            b = block_sizes[l_B]
            if l_B > b0 and b > 0:
                # block_size can't increase more than c
                # block_size can't be less than c unless it is zero
                # Both of these things make sense
                A = np.hstack((A[:, old_b-b+c:old_b], B))
                l_A += old_b-b+c
            elif l_B == b0 and b > 0:
                A = A[:, b0-b:b0]
                l_A = b0-b
            elif b == 0:  # no SNPs to left in window, e.g., after a sequence gap
                A = np.array(()).reshape((n, 0))
                l_A = l_B
            if l_B == md:
                c = m - md
                rfuncAB = np.zeros((b, c))
                rfuncBB = np.zeros((c, c))
            if b != old_b:
                rfuncAB = np.zeros((b, c))

            B = snp_getter(c)

            np.dot(A.T, B / n, out=rfuncAB)
            rfuncAB = l2(rfuncAB)
            cor_sum[l_A:l_A+b, :] += np.sum(rfuncAB, axis=1).reshape(b,1)
            cor_sum[l_B:l_B+c, :] += np.sum(rfuncAB, axis=0).reshape(c,1)
            np.dot(B.T, B / n, out=rfuncBB)
            rfuncBB = l2(rfuncBB)
            cor_sum[l_B:l_B+c, :] += np.sum(rfuncBB, axis=0).reshape(c,1)

        print(cor_sum)
        return cor_sum

class PlinkBEDFile(__GenotypeArrayInMemory__):
    # Interface for Plink .bed format
    def __init__(self, fname, n, snp_list, keep_snps=None, keep_indivs=None, mafMin=None):
        self._bedcode = {
            2: ba.bitarray('11'),
            9: ba.bitarray('10'),
            1: ba.bitarray('01'),
            0: ba.bitarray('00')
            }

        __GenotypeArrayInMemory__.__init__(self, fname, n, snp_list, keep_snps=keep_snps,
            keep_indivs=keep_indivs, mafMin=mafMin)

    def __read__(self, fname, m, n):
        if not fname.endswith('.bed'):
            raise ValueError('.bed filename must end in .bed')

        fh = open(fname, 'rb')
        magicNumber = ba.bitarray(endian="little")
        magicNumber.fromfile(fh, 2)
        bedMode = ba.bitarray(endian="little")
        bedMode.fromfile(fh, 1)
        e = (4 - n % 4) if n % 4 != 0 else 0
        nru = n + e
        self.nru = nru
        # Check magic number.
        if magicNumber != ba.bitarray('0011011011011000'):
            raise IOError("Magic number from Plink .bed file not recognized")

        if bedMode != ba.bitarray('10000000'):
            raise IOError("Plink .bed file must be in default SNP-major mode")

        # Check file length.
        self.geno = ba.bitarray(endian="little")
        self.geno.fromfile(fh)
        self.__test_length__(self.geno, self.m, self.nru)
        return (self.nru, self.geno)

    def __test_length__(self, geno, m, nru):
        exp_len = 2*m*nru
        real_len = len(geno)
        if real_len != exp_len:
            s = "Plink .bed file has {n1} bits, expected {n2}"
            raise IOError(s.format(n1=real_len, n2=exp_len))

    def __filter_indivs__(self, geno, keep_indivs, m, n):
        n_new = len(keep_indivs)
        e = (4 - n_new % 4) if n_new % 4 != 0 else 0
        nru_new = n_new + e
        nru = self.nru
        z = ba.bitarray(m*2*nru_new, endian="little")
        for e, i in enumerate(keep_indivs):
            z[2*e::2*nru_new] = geno[2*i::2*nru]
            z[2*e+1::2*nru_new] = geno[2*i+1::2*nru]

        self.nru = nru_new
        return (z, m, n_new)

    def __filter_snps_maf__(self, geno, m, n, mafMin, keep_snps):
        
        # Credit to Chris Chang and the Plink2 developers for this algorithm
        # Modified from plink_filter.c
        # https://github.com/chrchang/plink-ng/blob/master/plink_filter.c

        # Genotypes are read forwards (since we are cheating and using endian="little")

        # A := (genotype) & 1010...
        # B := (genotype) & 0101...
        # C := (A >> 1) & B

        # Then

        # a := A.count() = missing ct + hom major ct
        # b := B.count() = het ct + hom major ct
        # c := C.count() = hom major ct

        # Which implies that

        # missing ct = a - c
        # # of indivs with nonmissing genotype = n - a + c
        # major allele ct = b + c
        # major allele frequency = (b+c)/(2*(n-a+c))
        # het ct + missing ct = a + b - 2*c

        # Why does bitarray not have >> ????

        nru = self.nru
        m_poly = 0
        y = ba.bitarray()
        if keep_snps is None:
            keep_snps = range(m)
        kept_snps = []
        freq = []
        for e, j in enumerate(keep_snps):
            z = geno[2*nru*j:2*nru*(j+1)]
            A = z[0::2]
            a = A.count()
            B = z[1::2]
            b = B.count()
            c = (A & B).count()
            major_ct = b + c  # number of copies of the major allele
            n_nomiss = n - a + c  # number of individuals with nonmissing genotypes
            f = major_ct / (2*n_nomiss) if n_nomiss > 0 else 0
            het_miss_ct = a+b-2*c  # remove SNPs that are only either het or missing
            if np.minimum(f, 1-f) > mafMin and het_miss_ct < n:
                freq.append(f)
                y += z
                m_poly += 1
                kept_snps.append(j)

        return (y, m_poly, n, kept_snps, freq)

    def __causal_SNPs__(self, p_causal, h2):
        m, n = self.m, self.n
        beta = np.zeros(m)

        if p_causal < 1:
            n_SNPs = np.round(m*p_causal).astype(int)
            SNP = np.random.choice(m, size=n_SNPs, replace=False)
        else:
            SNP = np.arange(m)
            n_SNPs = m

        beta[SNP] = np.array(np.random.normal(0, np.sqrt(h2/n_SNPs), n_SNPs))

        return beta

    def get_pheno(self, args, log, add_noise=True):
        # This will always be 'nextSNPs'
        n, m = self.n, self.m
        snp_getter = self.nextSNPs
        C = np.zeros((n,1))
        phi = np.zeros(n)

        if args.h2_A > 0:
            beta_A = self.__causal_SNPs__(args.beta_prop_A, args.h2_A)
        else:
            beta_A = np.zeros(m)
        
        if args.h2_D > 0:
            beta_D = self.__causal_SNPs__(args.beta_prop_D, args.h2_D)
        else:
            beta_D = np.zeros(m)

        if args.h2_AC > 0:
            if args.gxe_non_independent:
                beta_AC = np.sqrt(args.h2_AC / args.h2_A) * beta_A
                print(beta_A[0:10])
                print(beta_AC[0:10])
            else:
                beta_AC = self.__causal_SNPs__(args.beta_prop_AC, args.h2_AC)
        else:
            beta_AC = np.zeros(m)
            # Also, create the covariate.
            # DEV: May want to pass the covariate in the creation of the phenos.
        if args.h2_AC > 0:
            if args.C_bool:
                C = np.random.binomial(1, args.C_bool_p, size=n)
                C = (C - np.mean(C)) / np.std(C)
                C = C.reshape((n,1))
                log.log('Boolean covariate.')
                log.log('Average kurtosis for this covariate: {K}.'.format(K=np.sum(C**4)/n))
            else:
                C = np.random.normal(0, 1, n)
                C = (C - np.mean(C)) / np.std(C)
                C = C.reshape((n,1))
                log.log('Normally distributed covariate. Kurtosis should be around 3.')
                log.log('Average kurtosis for this covariate: {K}.'.format(K=np.sum(C**4)/n))
            
        else:
            beta_AC = np.zeros(m)

        while self._currentSNP < m:
            if (self._currentSNP % 10000 == 0):
                print(self._currentSNP)
            start = self._currentSNP
            current_c = min(args.chunk_size, (m - self._currentSNP))
            X_A, X_D = snp_getter(current_c)
            stop = self._currentSNP
            phi += np.dot(X_A, beta_A[start:stop]) + np.dot(X_D, beta_D[start:stop]) + np.dot(X_A * C, beta_AC[start:stop])
        
        self._currentSNP = 0

        return phi, beta_A, beta_D, beta_AC, C

    def add_noise(self, phi, h2):
        n = self.n
        noise = np.zeros(n)

        if h2 < 1:
            noise = np.array(np.random.normal(0, np.sqrt(1-h2), n))
            
        phi += noise

        # Renormalise the phenotype.
        phi = (phi - np.mean(phi)) / np.std(phi)
        # Reset the current SNP in case we want to do multiple simulations.
        self._currentSNP = 0
        return phi

    def nextSNPs_add(self, b):
        X_A, X_D = self.nextSNPs(b)
        return X_A

    def nextSNPs_dom(self, b):
        X_A, X_D = self.nextSNPs(b)
        return X_D

    def nextSNPs(self, b):
        try:
            b = int(b)
            if b <= 0:
                raise ValueError("b must be > 0")
        except TypeError:
            raise TypeError("b must be an integer")

        if self._currentSNP + b > self.m:
            s = '{b} SNPs requested, {k} SNPs remain'
            raise ValueError(s.format(b=b, k=(self.m-self._currentSNP)))

        c = self._currentSNP
        n = self.n
        nru = self.nru
        slice = self.geno[2*c*nru:2*(c+b)*nru]

        X_A = np.array(slice.decode(self._bedcode), dtype="float64").reshape((b, nru)).T
        # DEV - Need to be very careful here...
        # X_A = (self.freq > 0.5)*(-X_A+2) + (self.freq <= 0.5)*X_A
        # DEV: SUPER IMPORTANT WE FIGURE OUT WHAT'S GOING ON HERE.
        X_D = np.copy(X_A)

        X_A = X_A[0:n, :]
        X_D = X_D[0:n, :]

        Y_A = np.zeros(X_A.shape)
        Y_D = np.copy(Y_A)

        for j in range(0, b):

            newsnp_A = X_A[:, j]
            newsnp_D = X_D[:, j]
            ii = newsnp_A != 9

            avg_A = np.mean(newsnp_A[ii])
            newsnp_A[np.logical_not(ii)] = avg_A

            denom_A = np.std(newsnp_A)

            if denom_A == 0:
                denom_A = 1

            Y_A[:, j] = (newsnp_A - avg_A) / denom_A

            # Now the dominance effects.
            p = avg_A / 2
            
            newsnp_D[newsnp_D == 0] = - p / (1 - p)
            newsnp_D[newsnp_D == 2] = - (1 - p) / p

            # Average genotype imputation (everything that's 9 (unknown), gets
            # sent to the average of the dominance encoded column).
            avg_D = np.mean(newsnp_D[ii])
            newsnp_D[np.logical_not(ii)] = avg_D

            denom_D = np.std(newsnp_D)

            if denom_D == 0:
                denom_D = 1

            # and renormalise
            Y_D[:,j] = (newsnp_D - avg_D) / denom_D

        self._currentSNP += b

        return Y_A, Y_D

