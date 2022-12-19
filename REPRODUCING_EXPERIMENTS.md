# Overview

This document contains instructions for reproducing the experiments described in the FijiRelax paper.

# Dataset

The dataset associated with FijiRelax test scripts is available at Zenodo :
https://zenodo.org/record/4518731
Download the latest version : https://doi.org/10.5281/zenodo.4518731

# Reproducing figures of the paper

First, follow the instructions given in the Readme file to proceed to installation, and download the dataset. Then use the beanshell scripts to reproduce the experiments of the paper. You can find them in the FijiRelax repository, in test/Beanshell_Scripts/*
Each figure is associated with a corresponding script, the correspondence is made with an explicit naming of the scripts. Drag-slide a script into Fiji, then click run.

# Reproducing performance benchmarking

The figure 1-f is a performance benchmarking using open data provided with QMRLab. Reproducing this experiment involves installing other software to make the benchmarking. The performance tests only measure the computation time for the fitting process and do not include any loading time for the application or data. They were run on a laptop Dell Latitude with an Intel® Xeon(R) E-2186M CPU @ 2.90GHz × 8 cores/16 threads and a 32 GB RAM, running Ubuntu 20.04. 
QMRLab last release indicated was 2.4.1, and was run with octave. MyRelax last relase indicated was 1.0.0, and was run with Python 3.7. FijiRelax was v4.0.4, executed with Java-SE 1.8. In both cases, power mode was activated in "Performance mode" (no sleep for the processors)


## Computation time measurement with FijiRelax

Load the script 1-f and click run. The corresponding computation time should appear in the log window of Fiji (in our case : 1.646 seconds, 2.369 seconds, 2.519 seconds)


## Computation time measurement with MyRelax

MyRelax can be found at https://github.com/fragrussu/MyRelax
A comprehensive guide help to install it.
Install the five dependancies, then clone the repository, and navigate to ./MyRelax/myrelax
We recommend to make a clean processing arborescence including an "Input" and an "Output" directory, that we denote $INPUT and $OUTPUT thereafter. Before starting, the $INPUT directorty should contain the input data required. This input data can be found in the Zenodo downloaded archive, in the subdirectory $ZENODO_ARCHIVE/Experiments_for_figures_and_tables/Input_data/Images/Human_brain_nifti/*
Please note that you can maximize the performance by identifying what is the optimal number of threads that can be handled for parallel computing. This number will be denoted $NCPU

The script used is getT2T2star.py . In order to make it limit performance measurements to curve fitting, and to not include loading times, we change a bit the last lines of the scripts. That way, the script print the exact fitting time. The modified version of this Python script can be found in the FijiRelax repository in test/Python_Scripts/

To run the experiment, use the following command :
python getT2T2star.py $INPUT/SEdata.nii $INPUT/TES.txt $OUTPUT/  --mask $INPUT/Mask.nii.gz   --ncpu $NCPU

The result should be displayed in the Python command-line interface (7.479 seconds in our case)



## Computation time measurement with QMRLab
QMRlab can be found at https://github.com/qMRLab/qMRLab you can clone it to a directory we'll denote $QMRLAB_DIR
While the documentation is comprehensive, reproducing the experiment with open tools (octave) can be tricky. 
We applied the documented installation manuel until a point. We had to install octave, but also octave-dev, and finally install manually the needed octave packages from the shell in order to go through a dead-lock within the octave version manager:
$bash-shell$> sudo apt-get install octave octave-dev
$bash-shell$> sudo apt-get install octave-control octave-image octave-io octave-optim octave-signal octave-statistics

Then : 
$bash-shell$> cd $QMRLAB_DIR
$bash-shell$> octave
$octave-shell$> startup
$octave-shell$> startup
$octave-shell$> model=mono_t2;
$octave-shell$> qMRgenBatch(model)

The two last lines generate the matlab instructions and download the corresponding dataset, in a directory that you have to select (we denote it $DEMO_DIR). Within the octave shell, navigate to $DEMO_DIR/mono_t2_demo/. The matlab script to be run is mono_t2_batch.m

In order to run the performances tests, we insert some code to display time elapsed during the fit.
The modified versions of the matlab script (mono_t2_batch_with_offset.m and mono_t2_batch_without_offset) can be found in the FijiRelax repository in test/Matlab_Scripts/

To run the experiment, use the following command :
$octave-shell$> mono_t2_batch_without_offset
The result should be displayed in the octave command-line interface

This gives you the computation time for the first case of QMRLab (436.822 seconds in our case). Repeat it with the second script to get the second case (451.18 seconds in our case).



 
