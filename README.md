# Preeclampsia Early Assessment of Risk from Liquid Biopsy (PEARL)

This repository contains code for the study:  

Mohamed Adil, Teodora R. Kolarova, Anna-Lisa Doebley, Leah A. Chen, Cara Tobey, Patricia Galipeau, Sam Rosen, Michael Yang, Brice Colbert, Robert D. Patton, Thomas W. Persse, Erin Kawelo, Jonathan B. Reichel, Colin C. Pritchard, Shreeram Akilesh, Christina M. Lockwood, Gavin Ha†, Raj Shree†.  
<b>Preeclampsia risk prediction from non-invasive prenatal cell-free DNA screening.</b>   
Nature Medicine. 2025 Feb 12. doi: [10.1038/s41591-025-03509-w](https://doi.org/10.1038/s41591-025-03509-w) Online ahead of print.


## System requirements

### Dependencies
Python/3.9.6-GCCcore-11.2.0

| Package	|	Version |
| ------- | ------- |
Jupyterlab |	3.1.6  
Joblib	|	1.0.1  
Matplotlib	| 3.7.1  
Pandas	|	2.0.2  
scikit-learn	| 1.1.1  
scipy	|	1.12.0  
seaborn	|	0.11.2  
sklearn	|	0.0.post1  
xgboost	|	1.6.1  

## Input files  
Meta data - Supplementary_Tables.xlsx  
Features data - Raw_feature_tables.xlsx  

## Installation guide  
Load jupyter notebook in jupyterlab  

## Demo/ Instructions for use  
Update path to meta and feature data.  
Update path for output files.  
`Run > Run All Cells` 

## Expected Outputs  
Trained Griffin-FF model  
Trained Griffin-PE models   
Figure 2G & 2H  




