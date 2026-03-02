# BIOS6624 Project 1: Impact of Hard Drug Use on HAART Treatment Response

This project investigates:
	• how baseline hard drug use (heroin, cocaine, etc.) influences HIV treatment response two years after initiating Highly Active Antiretroviral Therapy (HAART).

This project employs both frequentist and bayesian framework approach to ensure robust results:
	1. Frequentist Framework used Multiple Linear Regression models adjusted for baseline values and demographic covariates.
	2. Bayesian Framework: Bayesian Linear Regression models were fitted using brms.
	  • Inference include 95% HPDI, probability of direction (pd), and Pr(Clinically Meaningful Harm)
	3. Adherence: Mediation analysis to determine if drug use effects are driven by poor medication adherence (ADH).
	  
Data
	•	The cleaned dataset being used for analysis can be found from "Project1/Data/hiv_merge.csv", which includes 401 subjects with 371 non-users and 30 hard drug users. 

Code
	•	The R Markdown (.Rmd) file "Project0/Code/Project1_word.Rmd" fully reproduces the analysis.
	•	The R Markdown script:
	•	Performs all data cleaning and variable construction
	•	Generates tables summarizing coefficients of hard drug use
	•	Fits linear and bayesian models for four treatment response outcomes
	•	Produces diagnostic figures (Trace plots, ACF plots, posterior density plots)

The final output is a Word document containing all results for this project.