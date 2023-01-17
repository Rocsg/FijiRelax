---
title: 'FijiRelax: Fast and noise-corrected estimation of MRI relaxation maps in 3D + t'
tags:
  - Java
  - Fiji
  - Image analysis
  - MRI
  - Relaxometry
  - Phenotyping
authors:
  - name: Romain Fernandez
    orcid: 0000-0003-3670-044X
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Cedric Moisy
    orcid: 0000-0003-4637-4713
    affiliation: 3
affiliations:
 - name: CIRAD, UMR AGAP Institut, F-34398 Montpellier, France
   index: 1
 - name: AGAP Institut, Univ Montpellier, CIRAD, INRAE, Institut Agro, Montpellier, France.
   index: 2
 - name: IFV, French Institute of Vine and Wine, UMT Géno-Vigne, UMR AGAP, F-34398 Montpellier, France.
   index: 3
date: 17 October 2022
bibliography: paper.bib

---

# Summary
Quantitative Magnetic Resonance Imaging (MRI) is a relevant technique to measure water availability and binding in animal or vegetal tissues. Among others, MRI is useful to decipher the impacts of stresses on living samples [@sidi2019multiscale]. Reproducible measurements of anatomical and physiological changes can be obtained by quantitative mapping of biomarkers (see \autoref{fig:figure1}-a) such as proton density (PD) and magnetization relaxation times (T1, T2) in 2D images [@borisjuk2012surveying ; @van2016quantitative]. However, the application of quantitative MRI to 3D + t microscopic imaging is still challenging. Detection of weak signals associated with small voxels requires longer acquisition times, which favors spatial drift and tissue deformations [@han2015drift]. In addition, low signal-to-noise ratio (SNR) environments dramatically reduce the estimation accuracy of the parameters [@raya2010t2]. Non-specialists facing these issues with multi-echo spin-echo data may be discouraged when using available open implementations of T1-T2 relaxometry: there is a lack of a user-friendly tools capable of processing large data in a reasonable time. We propose FijiRelax, an efficient Fiji plugin for the estimation of relaxation parameter maps from 3D + t multi-echo spin-echo sequences, integrating Rician noise and spatial drift corrections. The plugin is distributed and documented as an ImageJ plugin ([Plugin page](https://imagej.net/plugins/fijirelax)), available in a graphical interface and scriptable in BeanShell to allow batch processing for high-throughput experiments. 

# Statement of need
FijiRelax is a generic tool for multi-echo spin-echo T1-T2 relaxometry capable of processing a wide variety of MRI images ranging from a plant stem to a human brain (see \autoref{fig:figure1}-b). It has been designed for three types of scientists: i) end-users using a GUI, ii) advanced users able to use a scripting language to process large number of images, and iii) developers able to adapt and extend the application with new functionalities.

* **End-users using a GUI**: this mode is recommended for scientists who are not specialists in image processing or programming. Download FijiRelax through the official Fiji release, and follow the step-by-step installation instructions, as well as the hands-on tutorials built on the test dataset hosted at Zenodo [@fijirelaxDataset]. Then, use the graphical user interface to import and process your own Bruker/NIFTI/Custom data, explore the relaxation maps in space and time using the graphical relaxation curve explorer and export your results as 2D/3D/4D TIFF images. This mode is also recommended for studying new datasets or new biological questions. Among the interface features, the plugin provides a graphical explorer to visualize the relaxation curves, and the estimated PD-weighted T1 and T2 distributions over customizable areas of interest. In 5D hypermaps, the distributions at each timepoint can be displayed simultaneously, giving access to valuable information on water distribution in tissues and its evolution during the monitoring period. 

* **Advanced users**: this mode can be used by scientists with programming skills. Load the sample BeanShell scripts provided in the test directory of the Github repository ([https://github.com/rocsg/fijirelax](https://github.com/rocsg/fijirelax)) by dragging them into the Fiji interface and run the scripts to reproduce the results shown in \autoref{fig:figure1}: import a dataset, convert it to a HyperMap (see \autoref{fig:figure1}-e), compute the parameter maps. Then, adapt these scripts to your needs to process your own data and batch-process multiple experiments.

* **Developers**: this mode is for programmers fluent with Java and Maven. Start by exploring the FijiRelax API: [API Overview](https://javadoc.io/doc/io.github.rocsg/fijirelax/latest/index.html). Build your own tools on top of the FijiRelax library, provided as a jar file hosted at Maven central repository ([Artifact](https://search.maven.org/artifact/io.github.rocsg/fijirelax)), by indicating FijiRelax as a dependency in your POM file and run the unit tests. FijiRelax is hosted on a Github public repository ([https://github.com/rocsg/fijirelax](https://github.com/rocsg/fijirelax)) and developers can offer to contribute to its development and extend it by requesting features, or proposing new features.

# State of the field
Open-access implementations of T1-T2 relaxometry curve fitting have been released with scripting capabilities in Python [@GRUSSU2020116884] and with a graphical interface in MATLAB [@karakuzu2020qmrlab]. FijiRelax offers both scripting capabilities and a graphical interface, and is 4 times faster than the Python implementation and 100 times faster than the MATLAB implementation (see computation time benchmark in \autoref{fig:figure1}-f), while including noise-corrected fitting and spatial drift correction. Main features of FijiRelax are:

* Rice noise integration in the curve fitting algorithm, leading to a more accurate estimation in low SNR situations (see \autoref{fig:figure1}-d)

* Spatial drift correction by automatic registration using Fijiyama libraries [@fernandez2021fijiyama] before map computation (see \autoref{fig:figure1}-c). 

* Handling a 2D/3D time-series, capture data and acquisition metadata and storing it as a single TIFF file ("HyperMap"), which can be visualized with Fiji and investigated with the FijiRelax curve explorer (see the 2D/3D to 5D workflow in \autoref{fig:figure1}-e).

With these capabilities, this plugin is suitable for studies in low signal-to-noise situations such as microscopic MRI and time-lapse MRI. FijiRelax works efficiently with a wide range of data (see \autoref{fig:figure1}-b) and is suitable for large 3D datasets and time-lapse experiments. We believe that FijiRelax will facilitate the implementation of quantitative MRI approaches to open new avenues in MRI-based tissue phenotyping.

# Acknowledgements
We acknowledge the contributions of Jean-Luc Verdeil, Christophe Goze-Back, Anne-Sophie Spilmont, Maïda Cardoso, Dimiter Prodanov, Janne Holopainen and Christophe Pradal. 
This work was funded by the Plan Deperissement de la Vigne (France Agrimer) and APLIM Project (Agropolis Foundation), and realized within the MaCS4Plants network.

# Figures
![Figure 1: Quantitative MRI using FijiRelax. a) Computation of PD, T1, and T2 maps (bottom line) from spin-echo sequences collected on a living Sorghum stem (upper line). Middle: visualization of spin-echo values with the GUI. b) Result of the map computation from a benchmarking dataset of a human brain (left) and a grapevine stem (right). c) Results of the drift artifacts correction with the registration feature. Maps shown are computed before (left) and after (right) registration. d) Comparison of the Rice fit model (FijiRelax feature) with fitting models used in other open-source software. From left to right: expected results, exponential fit results, offset fit results, and Rice fit results. e) FijiRelax workflow. Yellow boxes: input/output data; grey boxes: processing operations; white stars: compulsory steps. f) Benchmark of FijiRelax against other open-source equivalent software.\label{fig:figure1}](images/figure.png){ width=100% }


# References

