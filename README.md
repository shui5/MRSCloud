MRSCloud Documentation

MRSCloud is written in MATLAB to simulate metabolite basis functions for MRS quantifications. This work has only been focused and tested on brain MRS data using 3T scanners.

Source code is published on Github https://github.com/shui5/MRSCloud

Features: MRSCloud supports simulation of metabolite basissets in 3T using one-dimensional projection method, coherence pathways filters, and precalculation of propagators to accelerate the overall process. This tool is built based on FID-A (https://github.com/CIC-methods/FID-A).

Usage: 

Modify simMRS.json:
- Specify the metabolites of interest. Following are the full list of metabolites available for simulation.
- Common metabolites for the healthy brain: ["Asc","Asp","Cr","EA","GABA","GPC","GSH","Gln","Glu","Gly","H2O","Lac","mI","NAA","NAAG","PCh","PCr","PE","Ser","sI","Tau"],
- Metabolites for specific interests: ["Ala","Ace","AcO","AcAc","Cit","Cystat","HCar","Lys","Thr","bHG","Tyros","Val","Phenyl","bHB","Gua","iLe","Pyr","Suc","Tryp"],
- Exogenous compounds: ["EtOH","MSM"],

- FieldStr: The exiting version of MRSCloud only supports generation of basis sets in 3T.
- Options for vendor: GE/Philips/Siemens/Universal_Philips/Universal_Siemens
- Options for mega_or_hadam: UnEdited/MEGA/HERMES/HERCULES
- Options for localization: PRESS/sLASER
- TE: specific the echo time in ms.They are fixed at 68 ms and 80 ms for HERMES and HERCULES, respectively
- editOn and editOff are flexible for MEGA-PRESS. They are fixed for HERCULES.
- editTp is the duration of the editing pulses. It is normally 14 ms for MEGA-PRESS/HERMES and 20 ms for HERCULES.
- spatial_points: 41 is acceptable and 101 is ideal. The higher the number of spatial points, no longer it takes for the simulation.
- Keep the private parameters unchanged.
- Set the output directories. work_dir is not necessary for MATLAB.
- Run the run_simulations_cloud.m script.

Remark 1: Cystat, HCar, iLe, and Lys have complex J coupling spin systems and will significantly increase the simulation time)
Remark 2: Glc (glucose) is available in the spin system. Simulation will take extremely long for Glc due to its complex spin system.
Remark 3: Product sequence and rf waveform are not shared in the GitHub repo.

Should you publish material that made use of MRSCloud, please cite the following publication:

Hui SCN, Saleh MG, Zöllner HJ, Oeltzschner G, Fan H, Li Y, Song Y, Jiang H, Near J, Lu H, Mori S, Edden RAE. MRSCloud: A cloud-based MRS tool for basis set simulation. Magn Reson Med, 2022; 88: 1994-2004.

Acknowledgements: This work has been supported by the National Institutes of Health, Grant/Award Numbers: K99 DA051315, P41 EB031771, R00 AG062230,
R01 EB016089, R01 EB023963, R21 AG060245
