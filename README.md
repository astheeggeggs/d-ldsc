# d-ldsc
Dominance extension to the population LD-score regression framework.

`d-ldsc` is a command line tool to determine the non-additive contribution of site specific effects to heritability of a trait. This code started as a fork of the `ldsc` [codebase](https://github.com/bulik/ldsc), but has subsequently diverged. However, much of their excellent tutorials and gotchas still apply, so be sure to check out the orginal `ldsc` [wiki](https://github.com/bulik/ldsc/wiki). `d-ldsc` builds on first LD-score regression paper, see the citation below. Also, much of this README is based on the original `ldsc` README.

`d-ldsc` also computes additive and dominance LD Scores, and has functionality to create phenotypes (with a specified polygenic underlying architecture comprising an additive and dominance contributions), and summary statistics (additive and dominance marginal associations).

The code is slightly different to use than `ldsc`. Rather than a single command with a collection of flags to do most of the work, we have split the executables into a collection of files to perform distinct tasks. They are:

+ `get_y.py` 
+ `get_sumstats.py`
+ `get_ldscores.py`
+ `get_h2.py`

These tools allow us to simulate polygenic continuous phenotypes with infinitesimal and spike/slab genetic architecture for additive and dominance heritability, evaluate the associated summary statistcs, determine additive and dominance LD-score in a reference panel, and determine estimates of the additive and dominance heritability, respectively. 

Run the `-h` for with each of the commands for detailed information about the various option flags. 


## Getting Started

In order to download `d-ldsc`, you should clone this repository via the commands
```  
git clone https://github.com/astheeggeggs/d-ldsc.git
cd d-ldsc
```

In order to install the Python dependencies, you will need the [Anaconda](https://store.continuum.io/cshop/anaconda/) Python distribution and package manager. After installing Anaconda, run the following commands to create an environment with LDSC's dependencies:

```
conda env create --file environment.yml
source activate d-ldsc
```

Once the above has completed, you can run:

```
./get_y.py -h
./get_sumstats.py -h
./get_ldscores.py -h
./get_h2.py -h
```
to print a list of all command-line options. If these commands fail with an error, then something as gone wrong during the installation process. 

## Updating d-ldsc

You can update to the newest version of `d-ldsc` using `git` in the usual way. Navigate to your `d-ldsc/` directory (e.g., `cd ldsc`), then run
```
git pull
```
If `d-ldsc` is up to date, you will see 
```
Already up-to-date.
```
otherwise, you will see `git` output similar to 
```
remote: Counting objects: 3, done.
remote: Compressing objects: 100% (3/3), done.
remote: Total 3 (delta 0), reused 0 (delta 0), pack-reused 0
Unpacking objects: 100% (3/3), done.
From https://github.com/astheeggeggs/d-ldsc
   95f4db3..a6a6b18  master     -> origin/master
Updating 95f4db3..a6a6b18
Fast-forward
 README.md | 15 +++++++++++++++
 1 file changed, 15 insertions(+)
 ```
which tells you which files were changed. If you have modified the `d-ldsc` source code, `git pull` may fail with an error such as `error: Your local changes to the following files would be overwritten by merge:`. 

In case the Python dependencies have changed, you can update the `d-ldsc` environment using

```
conda env update --file environment.yml
```

## Support

Before contacting us, please try the following:

1. The original `ldsc` github has an excellent [wiki](https://github.com/bulik/ldsc/wiki) which contains tutorials on [estimating LD Score](https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial), [heritability]. Much of these tutorials are applicable to usage of `d-ldsc`, you can also consult the orginal [FAQ](https://github.com/bulik/ldsc/wiki/FAQ).
2. The method to estimate additive heritability is described in the citation below. We have preprint incoming detailing estimation of dominance heritability in the UK-biobank, which also includes the method to estimate dominance heritability. A working version of the methods portion of supplement to our preprint is available in this repository and gives a detailed derivation of the approach. Find it at https://github.com/astheeggeggs/d-ldsc/methods/methods.pdf.

If that doesn't work, you can get in touch with us by submitting an issue to this repository.

## Citation

If you use the dominance LD-score regression software please cite this github repository and

[Bulik-Sullivan, et al. LD Score Regression Distinguishes Confounding from Polygenicity in Genome-Wide Association Studies.
Nature Genetics, 2015.](http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3211.html)

## License

This project is licensed under GNU GPL v3.

## Authors

Duncan Palmer (MGH and Broad Institute of MIT and Harvard)

Original `ldsc` codebase authors:
Brendan Bulik-Sullivan (Broad Institute of MIT and Harvard)
Hilary Finucane (MIT Department of Mathematics)
