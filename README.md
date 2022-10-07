## Documentation and latest stable build
[![javadoc](https://javadoc.io/badge2/io.github.rocsg/fijirelax/javadoc.svg)](https://javadoc.io/doc/io.github.rocsg/fijirelax)
[![Maven Central](https://img.shields.io/maven-central/v/io.github.rocsg/fijirelax.svg?label=Maven%20Central)](https://search.maven.org/search?q=g:%22io.github.rocsg%22%20AND%20a:%22fijirelax%22)


## Summary

FijiRelax is a generic tool for 3D+t MRI analysis and exploration using multi-echo spin-echo sequences. This work was supported by the French Ministry of Agriculture, France AgriMer, CNIV and IFV, within VITIMAGE and VITIMAGE-2024 projects (program Plan National Dépérissement du Vignoble). This tools is developed in the context of :
- the [Vitimage 1](https://www.plan-deperissement-vigne.fr/recherches/programmes-de-recherche/vitimage) and [Vitimage 2024](https://www.plan-deperissement-vigne.fr/recherches/programmes-de-recherche/vitimage-2024) projects.
- the [Aplim](https://umr-agap.cirad.fr/recherche/projets-de-recherche/aplim) flagship project.


## Plugin features

- Proton density, T1 and T2 maps computation from multi-echo spin-echo sequences (multiple TR and/or TE)
- Parameters estimation by fitting noise-corrected mono- and biexponential decay models
- Automatic correction of spatial drift and deformations for long T1 or T2 sequences
- Exploration of T1/T2 distribution in ROI over time
- Operable through a GUI, or scriptable for batch processing of large datasets
<img src="https://github.com/Rocsg/FijiRelax/blob/master/images/fijirelax-snap-glob-explorer.png" width="800" caption="Time-lapse exploration of parameters in a plant under drought stress">
 
## Dataset for testing purpose

A comprehensive dataset can be found on Zenodo at [https://doi.org/10.5281/zenodo.4518730](https://doi.org/10.5281/zenodo.4518730) 



## Installation

The following video guide you throughout the installation process, and take a first tour of FijiRelax functions.

[![Installation](https://i.ytimg.com/vi/8jEVQjRbFcU/hqdefault.jpg)](https://www.youtube.com/watch/8jEVQjRbFcU)


In order to install FijiRelax on your computer, please follow these steps:

1\. *(if needed) *Download and install Fiji from https://fiji.sc/ ; start Fiji, and let it automatically update. Then restart Fiji.

2\. Open Fiji, run the **Update manager** (Help > Update). Click on "OK" to close the first popup windows, then click on the button **Manage update sites...**.

3\. In this list, activate **ImageJ-ITK** by checking the corresponding checkboxes. Don't close the window, or reopen it if you read this too late.

4\. Add the **Fijiyama** repository (by clicking on the button **Add update site**, and filling the fields : name = "/plugins/fijiyama", site = https://sites.imagej.net/Fijiyama), then check the associated checkbox. Now you can click on **Close** and apply the modifications.

5\. Restart Fiji: a new **FijiRelax** entry should be available in the menu (Plugins > Analyze"). If not, go back to the Update Manager, and check that the repositories **ImageJ-ITK** and **Fijiyama** are correctly selected.


## Preparing your data

FijiRleax needs properly formatted dataset:
- Nifti 4D images, or a set of Nifti 3D images
- Dicom dirs with 3D images, or a set of dirs with 2D images


## The interface

FijiRelax interface have four main panels :
- With the first panel, you can import / open / export data.
- The second panel holds the processing routines.
- The third panel contains the explorer button.
- The fourth panel has additional helper functions.

<img src="https://github.com/Rocsg/FijiRelax/blob/master/images/fijirelax-snap-main-window.png" width="300" caption="FijiRelax main window">


 

## Tutorials

**Tutorial part 1: proton density, T1 and T2 time-series from 3D dicom data of a sorgho plant**

[![Tutorial part 1](https://i.ytimg.com/vi/nhWRZN9puFg/hqdefault.jpg)](https://www.youtube.com/watch/nhWRZN9puFg)

**Tutorial part 2: from 4D HyperMaps to time-lapse plant physiology monitoring**
[![Tutorial part 2](https://i.ytimg.com/vi/tiJnq_xN-dY/hqdefault.jpg)](https://www.youtube.com/watch/tiJnq_xN-dY)

## HyperMap data structure

The output image is a 4D MR hyperimage. The "channels" slicer helps you to explore the 4th dimension, that is the images computed, and the input spin echo images. In detail :

-   Channels 1,2 3 are respectively the M0 map, T1 map, T2 map (see this information in the slice title, just upside the image pixels)
-   Channels 4,5, ..... NR-3 are the successive NR repetition times of the "T1 sequence", in increasing order.
-   Channels NR-2,..... NR-2+NE are the successive NE echo times of the "T2 sequence", in increasing order.

  
Unit for the channels 2 and 3 are milliseconds, what mean you can use it like it, without any additional conversion.  
For time-lapse experiments, one can compute such a 4D MR hyperimage at successive timepoints, and register and combine them in a 5D MR hyperimage (the same, with an additional slicer to walk through time). Registration and data combining can be done using the series registration mode of the [Fijiyama](/plugins/fijiyama) plugin.

  
## The science behind

This plugin compute M0, T1 and T2 maps pixelwise from a given set of spin-echo sequences, acquired with different repetition times and/or different echo times.

First a 3d registration is computed to align precisely the successive images, using libraries of the [Fijiyama](/plugins/fijiyama)  plugin. Then the rice noise level is estimated, and the M0, T1 and T2 parameters are estimated, fitting mono or bi-exponential curves, corrected with the measured rice noise. For more information, see the paper in next section.

## Citing this work

- Romain Fernandez, Cédric Moisy, Christophe Goze-Bac, Maïda Cardoso, Rahima Sidi-Boulenouar, Jean-Luc Verdeil, 2021  «FijiRelax: Fast and noise-corrected estimation of MRI relaxation maps in 3D + t» *under review*

## Software dependencies acknowledgements

- Johannes Schindelin et al for [Fiji](/software/fiji) (Schindelin et al., 2012)
- Karl Schmidt for MRI Analysis Calculator and CurveFitters

## License

This program is an open-source **free software**: it can be redistributed and/or modified under the terms of the **GNU General Public License** as published by the Free Software Foundation ([http://www.gnu.org/licenses/gpl.txt](http://www.gnu.org/licenses/gpl.txt)).

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

