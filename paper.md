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
 - name: CIRAD, UMR AGAP, F-34398 Montpellier, France
   index: 1
 - name: AGAP, Univ Montpellier, CIRAD, INRAE, Institut Agro, Montpellier, France.
   index: 2
 - name: IFV, French Institute of Vine and Wine, UMT Géno-Vigne, UMR AGAP, F-34398 Montpellier, France.
   index: 3
date: 17 October 2022
bibliography: paper.bib

---

# Summary

Quantitative Magnetic Resonance Imaging (MRI) is a relevant phenotyping technique for deciphering the impacts of abiotic and biotic stresses on living samples [@sidi2019multiscale]. Reproducible measurements of anatomical and physiological changes can be obtained by quantitative mapping of biomarkers such as proton density (PD) and magnetization relaxation times (T1, T2) in 3D + t [@borisjuk2012surveying ; @van2016quantitative]. However, the application of quantitative MRI to microscopic imaging is still challenging. Detection of weak signals requires longer acquisition times, which favors spatial drift and tissue deformations [@han2015drift]. In addition, low signal-to-noise ratio (SNR) environments dramatically reduce parameters estimation accuracy [@raya2010t2]. Non-specialists facing these issues with multi-echo spin-echo data may be discouraged when using available open implementations of T1-T2 relaxometry: there is a lack of a user-friendly tool capable of processing a large amount of data in a reasonable time. We propose FijiRelax, an efficient plugin for the estimation of relaxation parameter maps for 3D + t multi-echo spin-echo sequences, integrating spatial drift correction and Rician noise correction. The plugin is distributed and documented as an ImageJ plugin ([Plugin page](https://imagej.net/plugins/fijirelax)), available as a graphical interface and scriptable in BeanShell to allow batch processing for high-throughput experiments. 

# Statement of need
FijiRelax is a generic tool capable of processing a wide variety of MRI images ranging from a plant stem to a human brain (see \autoref{fig:figure1}). It is 100 times faster than open-access implementations of T1-T2 relaxometry curve fitting in MatLab [@karakuzu2020qmrlab] and Python [@GRUSSU2020116884] (see computation time benchmark in \autoref{fig:figure1}). It provides unique features such as:

* Rice noise integration in the curve fitting algorithm, leading to a more accurate estimation in low SNR situations (see \autoref{fig:figure1})

* Spatial drift correction by automatic registration before map computation (see \autoref{fig:figure1}). 

With these capabilities, this plugin is suitable for low signal-to-noise situations such as microscopic MRI and time-lapse MRI. 

FijiRelax has been designed for three types of scientists: i) end-users using a GUI, ii) advanced users able to use a scripting languages to process large number of images, and iii) developers able to adapt and extend the application with new functionalities.

* **End-users using a GUI**: scientists who are not image processing or programming specialists can download FijiRelax through the official Fiji release, and follow the step-by-step installation instructions, as well as the hands-on tutorials built on the test dataset hosted at Zenodo [@fijirelaxDataset]. Then, they can use the graphical user interface to import and process their own Bruker/NIFTI/Custom data, explore the relaxation maps in space and time using the graphical relaxation curve explorer and export their results as 2D/3D/4D TIFF images. This mode is preferred for studying new datasets or new biological questions. Among the interface features, the plugin provides a graphical explorer to visualize the relaxation curves, and the estimated PD-weighted T1 and T2 distribution over customizable areas of interest. In 5D hypermaps, the distributions at each time-point can be displayed simultaneously, giving access to valuable information about the evolution of tissue water distribution during the monitoring period.

* **Advanced users**: scientists with programming skills can load the sample BeanShell scripts by dragging them into the Fiji interface and run the scripts to obtain the results shown in \autoref{fig:figure1}. Then they can adapt these scripts to their needs, including processing their own data and batch-processing multiple experiments.

* **Developers**: programmers fluent with Java and Maven can get started by exploring the FijiRelax API: [API Overview](https://javadoc.io/doc/io.github.rocsg/fijirelax/latest/index.html). They can build their own tools on top of the FijiRelax library, provided as a jar hosted at maven central repository ([Artifact](https://search.maven.org/artifact/io.github.rocsg/fijirelax), by indicating FijiRelax as a dependency in their POM file. FijiRelax is hosted on a public repository on github ([https://github.com/rocsg/fijirelax](https://github.com/rocsg/fijirelax)) and developers can offer to contribute to the development of the the plugin.

FijiRelax works efficiently with a wide range of data (see \autoref{fig:figure1}) and is suitable for large 3D datasets and time-lapse experiments. We believe that FijiRelax will facilitate the implementation of quantitative MRI approaches to open new avenues in MRI-based tissue phenotyping.

# Acknowledgements

We acknowledge the contributions of Jean-Luc Verdeil, Christophe Goze-Back, Anne-Sophie Spilmont, Maïda Cardoso and Dimiter Prodanov. 
This work was funded by the Plan Deperissement de la Vigne (France Agrimer) and APLIM Project (Agropolis Foundation).

# Figures

![Figure 1: Quantitative MRI using FijiRelax. a) Computation of PD, T1, and T2 maps (bottom line) from spin-echo sequences collected on a living Sorghum stem (upper line). Middle: visualization of spin-echo values with the GUI. b) Result of map computation from a benchmarking dataset of a human brain (left) and a grapevine stem (right). c) Results of drift artifacts correction with the registration feature. Maps shown are computed before (left) and after (right) registration. d) Comparison of the Rice fit model (FijiRelax feature) with fitting model used in other open-source software. From left to right: expected results, exponential fit results, offset fit results, and Rice fit results. e) FijiRelax workflow. Yellow boxes: input/output data; grey boxes: processing operations; white stars: compulsory steps. f) Benchmark of FijiRelax against other open-source equivalent software\label{fig:figure1}](images/figure.png){ width=100% }


# References

