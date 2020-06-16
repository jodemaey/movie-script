# Movie-script

## About

A python script to generate movie with the output of the [MAOOAM](https://github.com/Climdyn/MAOOAM) model.

(c) 2016-2020 Jonathan Demaeyer.
See [LICENSE.txt](./LICENSE.txt) for license information.  

## Requirements

This code needs: 
 - Python 2
 - FFmpeg
 - Matplotlib
 - Numpy
 
See [environment.yml](./environment.yml) for a list of requirements.
This file can also be used to create an Anaconda environment:

    conda env create -f environment.yml
    conda activate movie-script

## Usage

    ./movie-script.py <data-filename> <ageom> <ogeom>

Example:
  
    ./movie-script.py ./data/test.dat 2x4 2x4


* **WARNING 1/2** : Assume that the spectral "geometry" of the model is contiguous.
               E.g. a "2x4" geometry means that all the mode with x- and
               y-wavenumber <= 2 and 4 respectively are included in the model.

* **WARNING 2/2** : Assume the mode indexing convention proposed in 

  De Cruz, L., Demaeyer, J. and Vannitsem, S.: The Modular Arbitrary-Order
  Ocean-Atmosphere Model: MAOOAM v1.0, Geosci. Model Dev., 9, 2793-2808,
  [doi:10.5194/gmd-9-2793-2016](http://dx.doi.org/10.5194/gmd-9-2793-2016), 2016.

  By default, Lua and python implementation of MAOOAM use this indexing convention,
  but the the fortran one is more flexible, and care must taken when configuring it.
  
## Example

* Representation of the atmospheric streamfunction modes % of total variance and of the meridional profile of the zonally averaged 
U wind at 500 hPa:

  ![](./misc/atm_example.gif)

* Representation of the oceanic streamfunction modes % of total variance and 3D phase space projection:

  ![](./misc/oc_example.gif)
