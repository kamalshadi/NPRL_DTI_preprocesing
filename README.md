# Introduction
This script runs HCP-like minimal preprocessing pipeline for diffusion weighted images.

# Input Requirement
The input requirement is based on NPRL group at Emory acquisition protocol. For each subject following files must be available:

* b0 image right to left acquisition
* b0 image left to right acquisition
* Diffusion image right to left phase-encoded acquisition
* Diffusion image left to right phase-encoded acquisition
* bvector (.bvec) for each diffusion image
* b-value (.bval) for each diffusion image

## Directory structure
Put all subjects' data in the folder **\<study>/Inputs** . **\<study>** is an optional name for the
study and the parent folder of Inputs. In the Inputs directory each subject have her own
subfolder. Picture below is a sample directory structure for a study with 3 participants.

![directory structure for the script](https://github.com/kamalshadi/NPRL_DTI_preprocesing/blob/master/sample.png)

# Environment Requirement
The script is written for python2.7 and only runs on linux-like terminals.

# Usage
1- open your Terminal

2 - type (option -t is a readout time of the acquisition)
```
python preprocess.py -i path_to_Study -t <fsl_readout_time>
```
After the code finishes, you can find all preprocessed data in **\<study>/Outputs** directory
