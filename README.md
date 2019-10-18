calculating-spaceweather-keywords
=================================

The Solar Dynamics Observatory (SDO) takes about a terabyte and a half of data a day, which is more data than any other satellite in NASA history. SDO data are stored in a [publicly-available, web-accessible pSQL database](http://jsoc.stanford.edu/ajax/lookdata.html) at Stanford University. These data are also accessible via a JSON API called [jsoc_info](http://jsoc.stanford.edu/jsocwiki/AjaxJsocConnect) and a python library called [drms](https://drms.readthedocs.io/en/stable/).

One of the data products released by the Solar Dynamics Observatory is called [Space-weather HMI Active Region Patches](http://link.springer.com/article/10.1007%2Fs11207-014-0529-3), or SHARPs. SHARPs include patches of vector magnetic field data taken by the Helioseismic and Magnetic Imager (HMI) instrument aboard SDO. These patches encapsulate automatically-detected active regions. SHARP data also include spaceweather keywords describing these active regions. [Bobra & Couvidat (2015)](http://arxiv.org/abs/1411.1405), [Bobra & Ilonidis (2016)](https://arxiv.org/abs/1603.03775), and [Jonas et al. (2018)](http://adsabs.harvard.edu/abs/2018SoPh..293...48J) used machine-learning algorithms to show that these spaceweather keywords are useful for predicting solar activity. 

### Contents

This repository contains several codes designed to show you how to interact with and understand SHARP data (view ipython notebooks on the [ipython notebook viewer](http://nbviewer.ipython.org/)).

* `plot_swx_d3.ipynb` generates interactive [d3](https://d3js.org/) plots of keywords and images, and movies using image data.
* `movie.ipynb` generates movies of SHARP data.
* `hedgehog.ipynb` provides a way to visualize the vector magnetic field in SHARP data.
* `feature_extraction.ipynb` takes images from another instrument on SDO, called the Atmospheric Imaging Assembly (AIA), to determine which AIA pixels fall within the SHARP bounding boxes; this code also contains examples of how to automatically extract features from AIA data.
* `calculate_swx_fits.py` contains all the functions to calculate spaceweather keywords from vector magnetic field data. `calculate_swx_workflow.ipynb` provides a workflow to calculate these keywords by fetching the vector magnetic field data from the JSOC database using the [`drms` package](https://joss.theoj.org/papers/10.21105/joss.01614) and using [Dask](https://dask.org/) to parallelize the calculations.
* `disambiguation.py` contains several functions that disambiguate the azimuthal component of the vector magnetic field data and construct the field vector in spherical coordinate components on a CCD grid. See `disambiguate_data.py` for some examples.

Sample data are included in this repository under the test_fits_files directory. All SDO data are publicly available. 

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
