SHARPs
======

The Solar Dynamics Observatory (SDO) takes about a terabyte and a half of data a day, which is more data than any other satellite in the NASA Heliophysics Division. One of the data products released by the Solar Dynamics Observatory science team is called [Space-weather HMI Active Region Patches](http://link.springer.com/article/10.1007%2Fs11207-014-0529-3), or SHARPs. SHARPs include patches of vector magnetic field data taken by the Helioseismic and Magnetic Imager (HMI) instrument aboard SDO. These patches encapsulate automatically-detected active regions. SHARP data also include spaceweather keywords describing these active regions. [Bobra & Couvidat (2015)](http://arxiv.org/abs/1411.1405), [Bobra & Ilonidis (2016)](https://arxiv.org/abs/1603.03775), and [Jonas et al. (2018)](http://adsabs.harvard.edu/abs/2018SoPh..293...48J) used machine-learning algorithms to show that these spaceweather keywords are useful for predicting solar activity. 

Users can access the SHARP data with a [SunPy](https://sunpy.org/) affiliated package called [drms](https://drms.readthedocs.io/en/stable/). If you use `drms` in your research, please cite [The SunPy Community et al. 2020](https://dx.doi.org/10.3847/1538-4357/ab4f7a) and [Glogowski et al. 2019](https://joss.theoj.org/papers/10.21105/joss.01614).

### Contents

This repository contains several notebooks and functions designed to interact with and understand SHARP data. The `requirements.txt` file lists all the packages necessary to run the notebooks and functions in this repository.

**Getting Started**
    
* The `plot_swx_d3.ipynb` notebook is a good place to get started. This notebook queries data using the SunPy affiliated package called [`drms`](https://joss.theoj.org/papers/10.21105/joss.01614), generates plots of keywords and images, demonstrates how to query the SHARP data, and exports data in a variety of formats.

**Space-weather Keywords**

* The `calculate_sharpkeys.py` file contains all the functions to calculate spaceweather keywords from vector magnetic field data. Sample data are included in this repository under the `files` directory. For an explanation of the variable `cdelt1_arcsec`, in the function `get_data()`, see `cdelt1_arcsec.pdf`. See `SHARP_Issue_Tracker.md` for a list of known issues with the SHARP data.

**Coordinates**

* The `feature_extraction.ipynb` notebook identifies which pixels in an image taken by the Atmospheric Imaging Assembly (AIA) instrument on SDO fall within the SHARP bounding box by applying coordinate transformations. 
* The `active_region_distances.ipynb` notebook calculates the distance between two SHARP regions.

**Visualizations**

* The `hedgehog.ipynb` notebook develops an aesthetically pleasing way to visualize a vector magnetic field using SHARP data.
* The `movie.ipynb` notebook generates movies of SHARP data.

**Disambiguation**

* The `disambiguation.py` file contains several functions that disambiguate the azimuthal component of the vector magnetic field data and construct the field vector in spherical coordinate components on a CCD grid. This works on both the SHARP data and full-disk data. See `disambiguate_data.py` for some examples.

**Parallelization**

* See [`sharp-features.ipynb`](https://gitlab.com/wtbarnes/aia-on-pleiades/-/blob/master/notebooks/tidy/sharp-features.ipynb) for a demonstration on how to efficiently calculate SHARP keywords using significant computational resources (in this case, the [NASA Ames Pleiades supercomputer](https://www.nas.nasa.gov/hecc/resources/pleiades.html)) and a Python package for parallelization called [Dask](https://dask.org/).

### Citation

If you use the [Space-weather HMI Active Region Patch](http://link.springer.com/article/10.1007%2Fs11207-014-0529-3) data in your research, please consider citing our paper. Here is the bibtex entry for the paper:

```
@ARTICLE{2014SoPh..289.3549B,
   author = {{Bobra}, M.~G. and {Sun}, X. and {Hoeksema}, J.~T. and {Turmon}, M. and 
	{Liu}, Y. and {Hayashi}, K. and {Barnes}, G. and {Leka}, K.~D.
	},
    title = "{The Helioseismic and Magnetic Imager (HMI) Vector Magnetic Field Pipeline: SHARPs - Space-Weather HMI Active Region Patches}",
  journal = {\solphys},
archivePrefix = "arXiv",
   eprint = {1404.1879},
 primaryClass = "astro-ph.SR",
 keywords = {Active regions, magnetic fields, Flares, relation to magnetic field, Instrumentation and data management},
     year = 2014,
    month = sep,
   volume = 289,
    pages = {3549-3578},
      doi = {10.1007/s11207-014-0529-3},
   adsurl = {http://adsabs.harvard.edu/abs/2014SoPh..289.3549B},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
