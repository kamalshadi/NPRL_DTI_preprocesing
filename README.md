# Introduction
This script runs HCP-like minimal preprocessing pipeline for diffusion weighted images.

# Input Requirement
The input requirement is based on NPRL group at Emory acquisition protocol. For each subject following files must be available:

* b0 image right to left acquisition

* b0 image left to right acquisition

* diffusion image right to left acquisition

* diffusion image left to right acquisition

* bvector (.bvec) for each diffusion image

* b-value (.bval) for each diffusion image


3- open Mac Terminal
4 - type (option -t is a readout time of the acquisition, it was 0.0508 sec for the TMS-EEG
dataset but it is 0.059 for the PAS dataset- these values are reported by MRIconvert tool)
python preprocess.py -i path_to_Study -t 0.059
5 - After the code finishes, you can find all preprocessed data in Outputs directory

1 - Put all subjects data in the folder <study>/Inputs . <study> is an optional name for the
study and the parent folder of Inputs. In the Inputs directory each subject have her own
subfolder. Picture below is a sample directory structure for a study with 3 participants.

2 - Within each subject folder, following files must exists:
