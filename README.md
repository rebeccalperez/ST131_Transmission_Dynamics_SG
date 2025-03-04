This repository contains all codes used in the analysis and figure generation for Perez et al.:
## Transmission Dynamics of Escherichia Coli Sequence Type 131 Amongst Households – a “One Health” Prospective Cohort Study

In this prospective cohort study, we adopted a longitudinal, One Health sampling approach to identify the transmission pathways of 
E. coli ST131 as a gut colonising strain in households. We quantified the carriage duration and acquisition risks of E. coli ST131 
via Markov models and mapped transmission events using both epidemiological and genomic data. Further, we conducted univariate and
multivariate analysis of risk factors associated with *E. coli* ST131 carriage and with carrier status (persistent, intermittent). 


## Main analysis: 
### [Regression analysis](https://github.com/rebeccalperez/ST131_Transmission_Dynamics_SG/blob/main/ST131_Multivariate_Univariate_Regression_Analysis.Rmd) 
### [MCMC models: data formatting and runcode](https://github.com/rebeccalperez/ST131_Transmission_Dynamics_SG/blob/main/2-State_MCMC_Runcode.Rmd)
### [MCMC stan models for all participants](https://github.com/rebeccalperez/ST131_Transmission_Dynamics_SG/blob/main/Two_State_MCMC_T10.stan)
### [MCMC stan model for persistent carrier subset](https://github.com/rebeccalperez/ST131_Transmission_Dynamics_SG/blob/main/Two_State_MCMC_T5.stan)

## Figure generation: 
### [Fig 2: ST131 density](https://github.com/rebeccalperez/ST131_Transmission_Dynamics_SG/blob/main/Figures/ST131_Isolate_Density_By_Sample_Plots.Rmd) 
### [Fig 3: MCMC results](https://github.com/rebeccalperez/ST131_Transmission_Dynamics_SG/tree/main/Plots/MCMC_Results) 

## Supplementary analyses and figures:
### [S1: Sensitivity analysis](https://github.com/rebeccalperez/ST131_Transmission_Dynamics_SG/blob/main/Figures/Supplementary/Persistent_Carrier_Sens_Analysis.Rmd)
### [S2: Hidden Markov model](https://github.com/rebeccalperez/ST131_Transmission_Dynamics_SG/blob/main/Figures/Supplementary/HMM_False_Neg_Carriage_Estimation.Rmd)
### [S3: Isolate sources](https://github.com/rebeccalperez/ST131_Transmission_Dynamics_SG/blob/main/Figures/Supplementary/ST131_Isolate_Sources_Plot.Rmd) 
### [S4: Antimicrobial susceptibility](https://github.com/rebeccalperez/ST131_Transmission_Dynamics_SG/blob/main/Figures/Supplementary/AMR_Phenotype_Heatmap.Rmd)


