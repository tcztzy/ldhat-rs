#ifndef LDHAT_H
#define LDHAT_H
/*! \file ldhat.h
    \brief LDhat
*/
/*! \def DEBUG 0
    \brief Debug option
*/
#define DEBUG 0 /*Debug option*/
/*! Constants for composite likelihood part
 */
/*! \def NSHUFF 1000
    \brief Number of permutations in tests for recombination
*/
#define NSHUFF 1000
/*! \def NRUN 1000000
    \brief Number of proposals in IS estimation of coalescent likelihoods
*/
#define NRUN 1000000
/*! \def ADD 10000
    \brief Size of extra PTs to be added when more memory needed
*/
#define ADD 10000
/*! \def SEQ_MAX 1000
    \brief Max number of sequences*/
#define SEQ_MAX 1000
/*! \def MAXNAME 65535
    \brief Max length of sequences names
*/
#define MAXNAME 65535
/*! \def MAXLINE 65535
    \brief Max length of line
*/
#define MAXLINE 65535 /*Max length of line*/
/*! \def MAXW 50
    \brief MAXW*2 = Max number of SNPs to consider for likelihood - i.e. ignore
   SNPS > MAXW apart
*/
#define MAXW 50
#define BURNIN 100000

/*! Constants for block routines
 */

/*! Site type */
struct site_type {
  int pt[16];        /*!< Haplotype pair type*/
  int nt;            /*!< Number of such type in data*/
  double ld_stat[3]; /*!< LD statistics for pair type*/
  int miss;          /*!< Integer indicating whether PT contains missing data*/
  double lkptmx;     /*!< Max Likelihood for PT*/
  double rmpt;       /*!< Rho_max for PT*/
  int rm;            /*!< Min number of recombination events for pair - 0 or 1*/
};

struct data_sum {

  int nseq;         /*!<Number of sequences*/
  int lseq;         /*!<No. segregating sites*/
  double tlseq;     /*!<Total length of sequence (may be in kb)*/
  int w;            /*!<Size of window to be used in analysis*/
  int hd;           /*!<Haploid (1) or diploid (2) data*/
  char lc;          /*!<Crossing-over (L) or gene conversion (C) model*/
  int ptt;          /*!<Number of pair types*/
  double avpwd;     /*!<Average pairwise differences*/
  double varpwd;    /*!<Sample variance in pairwise differences*/
  int rmin;         /*!<Lower bound on minimum number of recombination events*/
  double rwak;      /*!<4Ner estimated by Wakeley 1997*/
  double avc;       /*!<Average conversion tract length*/
  double th;        /*!<Theta per site*/
  double rho;       /*!<Rho for whole gene (or gamma for conversion model)*/
  int rho_i;        /*!<position of maximum rho*/
  double rho_drive; /*!<Rho to be used in driving simulations*/
  double fit_obs;   /*!<Observed fit*/
  double *rmap;     /*!<Array for recombination map*/
  double lkmax;     /*!<Maximum composite likelihood*/
  double **lksurf;  /*!<Array for likelihood surface*/
  double ld[4];     /*!<LD statistics*/
  double rme; /*!<Number of points for rho in coalescent likelihood estimation*/
  double rmax;  /*!<Max rho in coalescent likelihood estimation*/
  int rce;      /*!<Maximum rho for estimation - can be >> rmax*/
  int rcat;     /*!<Number of categories for estimating rho - can be >> rme*/
  double fit;   /*!<Fit for simulations*/
  double clr;   /*!<Composite likelihood ratio - from simulations*/
  int ng[2];    /*!<Counters for P values in simulations*/
  int n_update; /*!<Number of updates in MCMC*/
  int r_update; /*!<Number of updates between samples from MCMC*/
  double bpen;  /*!<Block penalty for MCMC*/
  int exact; /*!<Switch to speed up events if exact set of likelihoods will be
                inputed*/
  char prefix[MAXNAME + 1]; /*!< Prefix for output filenames */
};

struct block {
  int num;           /*!< position in array of active blocks*/
  double rate;       /*!< Recombination rate (per kb) in blocks*/
  int pos;           /*!< Which SNP the block starts at*/
  int size;          /*!< Length of block in SNPs*/
  struct block *bpr; /*!< Pointer to RH block*/
  struct block *bpl; /*!< Pointer to LH block*/
};
#endif
