RCA
===
Resolved Component Analysis

v2.0.2

Before using RCA from this repo, you might want to check out [MCCD](https://github.com/CosmoStat/mccd), a more recent PSF modelling approach. MCCD includes all aspects of RCA, but extends it in several ways. In particular, it can simultaneously fit all detectors of a multi-CCD mosaic camera. It has been shown to outperform RCA. Reference: [Liaudat et al., 2020](https://arxiv.org/abs/2011.09835)

## Description
RCA is a PSF modelling python package. It enforces several constraints, notably some related to sparsity and spatial structure, to build a spatially-varying PSF model from observed, noisy and possibly undersampled star stamps. Some modicum of documentation can be found [here](https://cosmostat.github.io/rca/) - see also [quick start](#quick-start) below.

## Requirements
The following python packages are required:

  - numpy
  - scipy
  - [ModOpt](https://github.com/CEA-COSMIC/ModOpt)
  
You will also need a compiled version of the sparse2d module of [ISAP](http://www.cosmostat.org/software/isap); alternatively, you should be able to install [PySAP](https://github.com/CEA-COSMIC/pysap) and let it handle the compilation and installation of sparse2d.

## Installation
After installing all dependencies, RCA just needs to be cloned and python-installed:

```
git clone https://github.com/CosmoStat/rca.git
cd rca
python setup.py install
```

## References
  - [Ngolè et al., 2016](https://arxiv.org/abs/1608.08104) - Inverse Problems, 32(12)
  - [Schmitz et al., 2020](https://arxiv.org/abs/1906.07676) - A&A, 636, A78
  
## Quick start
The basic syntax to run RCA is as follows:

```
from rca import RCA

# initialize RCA instance:
rca_fitter = RCA(4)

# fit it to stars
rca_fitter.fit(stars, star_positions)

# return PSF model at positions of interest
psfs = rca_fitter.estimate_psf(galaxy_positions)
```
A complete list of the parameters for `RCA` and its `fit` and `estimate_psf` methods can be found in [the documentation](https://cosmostat.github.io/rca/rca.html#module-rca). The main ones to take into account are:

  - RCA initialization:
    - `n_comp`, the number of eigenPSFs to learn ("r" in the papers)
    - `upfact`, the upsampling factor if superresolution is required ("m_d" or "D" in the papers)
  - `fit`:
    - `obs_data` should contain your observed stars (see note below for formatting conventions)
    - `obs_pos`, their respective positions
    - either `shifts` (with their respective centroid shifts wrt. a common arbitrary grid) or, if they are to be estimated from the data, a rough estimation of the `psf_size` (for the window function - can be given in FWHM, R^2 or Gaussian sigma)
  - `estimate_psf`:
    - `test_pos`, the positions at which the PSF should be estimated

The rest can largely be left to default values for basic usage.

Note `RCA.fit` expects the data to be stored in a `(p, p, n_stars)` array, that is, with the indexing (over objects) in the _last_ axis. You can use `rca.utils.rca_format` to convert to the proper format from the more conventional `(n_stars, p, p)`.

An example with stars from an HST ACS image can also be found in the `example` folder.

## Changelog
This is "v2" of RCA, which has been very largely overhauled. "v1" can still be accessed [here](https://github.com/CosmoStat/rca/commit/60845d44de56a9df58bed724ff2a1fbdae288c04). The most significant changes are:

  - Speed of the fitting step has been dramatically increased.
  - Interpolation to galaxy (or test star, or any other) positions has been added.
  - The source update step is now performed using the [Condat (2013 - pdf warning)](http://www.gipsa-lab.fr/~laurent.condat/publis/Condat-optim-JOTA-2013.pdf) algorithm, as implemented in ModOpt. In particular, this means the problem solved is now actually in synthesis form.
  - A lot of flexibility was added to the graph constraint, and can be accessed through the `rca.utils.GraphBuilder` arguments (which can be passed on to `RCA.fit`).
  - Simplicity and ease-of-use were favoured, but came at the cost of some flexibility. In particular:
    - "v2" always uses both sparsity in the Starlet domain and the graph constraint (both of which could in principle be turned off in "v1"); 
    - the package should now be ran from a Python session (whereas it could in principle be launched as a standalone executable before);
    - if shifts are not provided, we now expect a rough estimate of the PSF size to be given as input.
    
Either of the first two features could easily be re-added to the new version in the future. The latter does not seem too unreasonable, since you likely have a pretty decent idea of the size of your PSF - in fact, in the likely scenario where you obtained star stamps from [SExtractor](https://github.com/astromatic/sextractor), then you should already have both the shifts _and_ a pretty good estimation of your FWHM.
