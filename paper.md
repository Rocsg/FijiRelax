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
    affiliation: 2
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

Quantitative Magnetic Resonance Imaging (MRI) is a relevant phenotyping technique for deciphering impacts of abiotic and biotic stresses on living samples. Reproducible measurements of anatomical and physiological modifications can be obtained by quantitative mapping of biomarkers such as proton density (PD) and magnetization relaxation times (T1, T2) in 3D + t. However, applying quantitative MRI to plant tissues at a microscopic scale is still challenging. Detection of weak signals requires longer acquisition times, which favors spatial drift and tissue deformations. In addition, low signal-to-noise ratio (SNR) environments dramatically reduce parameters estimation accuracy. Finally, non-specialists might be deterred by long computation times and the lack of user-friendly solutions. 

# Statmement of need

FijiRelax is a Fiji plugin allowing calculation and exploration of 3D + t relaxation parameter maps from a series of spin-echo sequences. It is a generic tool capable of processing a wide variety of MRI images ranging from a plant stem to a human brain. Its performance was compared to other open-source solutions. FijiRelax provides an efficient and generic tool to facilitate the time-lapse analysis of resonance magnetic images and their application to plant and animal biology, as well as medicine. FijiRelax was designed to be used by both researchers and by students in courses on magnetical resonance imaging. 


# Acknowledgements

We acknowledge contributions from Jean-Luc Verdeil, Christophe Goze-Back, Anne-Sophie Spilmont and Maïda Cardoso.
This work was funded by the Plan Deperissement de la Vigne (France Agrimer) and APLIM Project (Agropolis Foundation).


# Figures
![Figure 1: Quantitative MRI using FijiRelax. \newline
a) Computation of PD, T1, and T2 maps (bottom line) from spin-echo sequences collected on a living Sorghum stem (upper line). Middle: visualization of spin-echo values with the GUI. \newline
 \newline
b) Result of map computation from a benchmarking dataset of a human brain (left) and a grapevine stem (right).  \newline
c) Results of drift artifacts correction with the registration feature. Maps shown are computed before (left) and after (right) registration.  \newline
d) Comparison of the Rice fit model (FijiRelax feature) with fitting model used in other open-source software. From left to right: expected results, exponential fit results, offset fit results, and Rice fit results. \newline
e) FijiRelax workflow. Yellow boxes: input/output data; grey boxes: processing operations; white stars: compulsory steps. \newline
f) Benchmark of FijiRelax against other open-source equivalent software



# References

