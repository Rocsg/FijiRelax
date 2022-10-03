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
![Figure 1: Quantitative MRI using FijiRelax.<br>
a) PD, T1, and T2 maps computed from spin-echo sequences collected on a living Sorghum stem.<br>
Top: Spin-echo sequences acquired through micro-RMI. Recovery times = {600, 1200, 1800, 10000 ms}, Echo time spacing= 11.15 ms. <br>
Middle: Spin-echo data modeling versus MRI measurements. The magnitude of the spin-echo signal (blue crosses) corresponds to points collected on the relaxation decay curves after TR seconds of longitudinal relaxation (red curve) and TE seconds of transversal relaxation (green curves).
Bottom: Computed PD, T1, and T2 maps.<br>
b) Maps computed from a grapevine stem (left) and human brain (right) datasets. <br>
c) Correction of drift artifacts using the registration feature. T1 values estimated from simulated spin-echo images showing a translation drift (250 µm total). Maps are computed before (left) and after (right) registration. <br>
d) Comparison of different fitting models for T2 values estimation.<br>
e) FijiRelax workflow. Yellow boxes: input/output data; grey boxes: processing operations; white stars: compulsory steps.<br>
f) Investigation of T1/T2 distribution in a specific ROI using the graphical explorer. <br>
g) Computation times measured during a classical phenotyping experiment.<br>
](images/figure.png){ width=100% }



# References

