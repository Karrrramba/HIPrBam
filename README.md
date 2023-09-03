Hydrophobic Interaction chromatography Protein porfiling with Bayesian Additive Models 

The HIPrBAM project aims at identifying differences in protein elution profiles from 2D LC-MS/MS experiments.
This is achieved by fitting two BAMs for each protein in the data set and quantyfying differences between the models by the likelihood ratio test (LRT).

This is an ongoing project.
An extensive Readme file will be provided shortly.

Some things still need some work:
- For some proteins, the simple model fits much better than the nested model.
- Simulated data is not quite comparable to Ctrl 1 in terms of distribution. 
- There is an error when fitting BRAF data from size-exclusion chromatography experiments.

Things to follow:
- Comparison of LRT and two-sample t-test results.
- Profiling of protein complexes.
