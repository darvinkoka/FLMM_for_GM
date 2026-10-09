# Functional Ground Motion Models (FGMM)

## Overview
This repository aims to introduce and implement **Functional Linear Mixed Models (FLMM)** 
within the framework of Ground Motion Modeling. 

Traditionally, empirical Ground Motion Models (GMMs) are developed to predict seismic 
intensity measures (e.g., spectral accelerations) at a discrete set of oscillator periods.
While effective, this scalar approach often neglects the intrinsic continuous nature of 
the seismic response and requires secondary empirical models to define cross-period correlations. 
This project seeks to overcome these limitations by treating the response spectra as 
continuous mathematical functions.

## Scope of the Project
The primary objective of this project is to develop a robust functional statistical
framework for ground motion prediction. By leveraging Functional Linear Mixed Models,
this research focuses on:

* **Continuous Spectral Modeling:** Shifting the paradigm from discrete scalar 
predictions to the modeling of the entire response spectrum as a continuous functional data object.
* **Hierarchical Variance Partitioning:** Utilizing the mixed-effects framework to 
rigorously account for the nested and hierarchical data structures typical of strong-motion databases 
(i.e., properly separating between-event, between-site, and within-event residual components).
* **Intrinsic Covariance Structure:** Naturally capturing and preserving the correlation
structure across the entire frequency domain without relying on external, 
empirically derived correlation matrices.