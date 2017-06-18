1 - Put all subjects data in the folder <study>/Inputs . <study> is an optional name for the
study and the parent folder of Inputs. In the Inputs directory each subject have her own
subfolder. Picture below is a sample directory structure for a study with 3 participants.

2 - Within each subject folder, following files must exists:
1. b0 image right to left acquisition - (b0 up)
2. b0 image left to right acquisition - (b0 dn)
3. diffusion image right to left acquisition - (bitc up)
4. diffusion image left to right acquisition - (bitc dn)
5. bvector (.bvec)
6. b-value (.bval)
3- open Mac Terminal
4 - type (option -t is a readout time of the acquisition, it was 0.0508 sec for the TMS-EEG
dataset but it is 0.059 for the PAS dataset- these values are reported by MRIconvert tool)
python preprocess.py -i path_to_Study -t 0.059
5 - After the code finishes, you can find all preprocessed data in Outputs directory
