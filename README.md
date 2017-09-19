# Introduction
This script runs HCP-like minimal preprocessing pipeline for diffusion weighted images.

# Input Requirement
The input requirement is based on NPRL group at Emory acquisition protocol. For each subject following files must be available:

* b0 image right to left acquisition
* b0 image left to right acquisition
* Diffusion image right to left phase-encoded acquisition
* Diffusion image left to right phase-encoded acquisition
* bvector (.bvec) for each diffusion image
* b-value (.bval) for each diffusion image<br/>


The naming scheme for the files must be set in the _**config.json**_ file placed in the root of the stydy folder. Sample config.json file is included in this repository.
```
{
  "filenames": {
    "b0_rl": [ #should contain two strings identifying b0 image right to left images
      "b0",  
      "RL"
    ],
    "b0_lr": [ #should contain two strings identifying b0 image left to right images
      "b0",
      "LR"
    ],
    "dwi_rl": [  #should contain two strings identifying weighted image right to left images
      "D",
      "RL"
    ],
    "dwi_lr": [ #should contain two strings identifying weighted image left to right images
      "D",
      "LR"
    ]
  },
  "volumes": {
    "b0_rl": 6, #number of b0 right to left volumes
    "b0_lr": 6, #number of b0 left to right volumes
    "dwi_rl": [
      6, #number of b0 imaged included in weighted acquisions (right to left)
      64 #number of weighted right to left volumes
    ],
    "dwi_lr": [
      6, #number of b0 imaged included in weighted acquisions (left to right)
      64 #number of weighted left to right volumes
    ]
  }
}
```
## Directory structure
Put all subjects' data in the folder **\<study>/Inputs** . **\<study>** is an optional name for the
study and the parent folder of Inputs. In the **Inputs** directory each subject have her own
subfolder. Picture below is a sample directory structure for a study with 3 participants.

<img src="https://github.com/kamalshadi/NPRL_DTI_preprocesing/blob/master/sample.png" alt="directory structure for the script" style="width: 200px;" style="text-align: center, horizontal-align: middle;"/>
# Environment Requirement
The script is written for python2.7 and only runs on linux-like terminals.

# Usage
1- open your Terminal

2 - type (option -t is a readout time of the acquisition)
```
python preprocess.py -i path_to_Study -t <fsl_readout_time>
```
After the code finishes, you can find all preprocessed data in **\<study>/Outputs** directory
