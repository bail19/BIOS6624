# BIOS6624 Project 0: Salivary Cortisol and DHEA Patterns of Changes over Time

Project Overview
This project examines the patterns of changes over time for salivary cortisol and DHEA levels collected in an at-home setting using a novel saliva collection device (SPIT booklet). 
The study evaluates:
	1.	Agreement between participant-recorded sampling times (booklet) and electronically recorded times (MEMs cap).
	2.	Adherence to protocol-specified sampling times (+30 minutes and +10 hours post-waking).
	3.	Patterns of changes in cortisol and DHEA levels over the day.

Mixed-effects models are used to characterize hormone trajectories while accounting for repeated measurements within subjects.

Data
	•	The cleaned analysis dataset can be found from "Project0/Data/Project0_Clean_v2.csv".
	•	The dataset includes sampling times, wake times, electronic monitoring data, and hormone measurements.

Code
	•	The R Markdown (.Rmd) file "Project0/Code/Project0_word.Rmd" fully reproduces the analysis.
	•	The R Markdown script:
	•	Performs all data cleaning and variable construction
	•	Generates tables summarizing agreement, adherence, and hormone levels
	•	Fits mixed-effects models for cortisol and DHEA
	•	Produces all figures (agreement plots, density plots, spaghetti plots with fitted curves)

The final output is a Word document containing all results for this project.