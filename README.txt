MRSCloud Documentation

MRSCloud is a cloud-based MATLAB executable for performing metabolites simulation to generate basis sets for 3T magnetic resonance spectroscopy (MRS) data. Visit https://braingps.mricloud.org/mrs-cloud for the online platform. Registration is required to use.

Source code is published on Github https://github.com/shui5/MRSCloud

Features:
MRSCloud supports simulations for up to 35 metabolites (29 commonly seen brain metabolites and 6 other metabolites for specific interests (see following list for details). Metabolites with complex J coupling spin system including (Cystat, HCar and Lys) may significantly increase simulation time.

1.  Alanine(Ala)
2.  Ascorbic Acid(Asc)
3.  Aspartic Acid(Asp)
4.  Citrate(Cit)
5.  Creatine(Cr)
6.  Ethanolamine(EA)
7.  Ethanol(EtOH)
8.  Gamma‐Aminobutyric Acid(GABA)
9.  Glycerophosphocholine(GPC)
10. Glutathione(GSH)
11. Glutamine(Gln)
12. Glutamate(Glu)
13. Glycine(Gly)
14. Water(H2O)
15. Lactate(Lac)
16. myo‐Inositol(mI)
17. N‐Acetyl Aspartate(NAA)
18. N‐Acetyl Aspartyl Glutamate(NAAG)
19. Choline‐containing Compounds(PCh)
20. Phosphocreatine(PCr)
21. Phosphorylethanolamine(PE)
22. Phenylalanine(Phenyl)
23. scyllo‐Inositol(sI)
24. Serine(Ser)
25. Taurine(Tau)
26. Threonine(Thr)
27. Tyrosine(Tyros)
28. Valine(Val)
29. β‐Hydoxybutyrate(bHB)

% Additional metabolites that are not common in normal brain.

30. Acetate(AcO)
31. Acetone(Ace)
32. Acetoacetate(AcAc)
33. Cystathionine(Cystat)
34. Homocarnosine(HCar)
35. Lysine(Lys)
36. 2‐Hydroxyglutarate(2HG)

Usage:
Visit https://braingps.mricloud.org/mrs-cloud 
1.  On the user interface, select metabolites that are interested to be included in the basis set.
2.  Select localization method (PRESS/sLASER).
3.  Select vendor (GE/Philips/GE/Universal_Philips/Universal_Siemens).
3.1 Visit https://pubmed.ncbi.nlm.nih.gov/30682536/ for more details of the Universal_Philips and Universal_Siemens sequences. 
4.  Select sequence (UnEdited/MEGA/HERMES/HERCULES).
4.1 If MEGA is selected, users have the option to change echo time (TE), target of the editing pulse (Edit On/ Edit Off), duration of the editing pulse and the edited target. Default setting is for GABA at TE 68 ms editing at 1.9 ppm for the "On" pulse and 7.5 ppm for the "Off" pulse.
5.  Submit the request
6.  Go to the 'My job status' tap to check the status of the simulation. 
7.  Download the basis set (in .zip format) when it is ready.
8.  Unzip the file and the basis set is saved in two formats (.BASIS and .mat). The basis set in .BASIS format is compatible with LCModel and the one in .mat format is compatible with Osprey.
9.  The .mat basis set contains all FIDs and spectral data which can be converted to other formats for other fitting tools.

10. For advanced users, visit https://github.com/shui5/MRSCloud for the source code.
11. Setup parameters in simMRS.json
12. Run run_simulations_cloud.m
13. Replace own tau time, TE1, waveforms etc. in load_parameters.m
14. Make sure to remove FID-A and Gannet from your MATLAB path.

For any questions, feedback, suggestions, or critique, please contact Steve Hui <stevehui@jhu.edu> or Richard Edden <edden@jhu.edu>.

Should you publish material that made use of MRSCloud, please cite the following publications:

Hui SCN, Saleh MG, Zöllner HJ, Oeltzschner G, Fan H, Li Y, Song Y, Jiang H, Near J, Lu H, Mori S, Edden RAE. MRSCloud: A cloud-based MRS tool for basis set simulation. Magn Reson Med. 2022 Jul 1. doi: 10.1002/mrm.29370. Epub ahead of print. PMID: 35775808.

Simpson R, Devenyi GA, Jezzard P, Hennessy TJ, Near J. Advanced processing and simulation of MRS data using the FID appliance (FID-A)-An open source, MATLAB-based toolkit. Magn Reson Med. 2017 Jan;77(1):23-33. doi: 10.1002/mrm.26091. Epub 2015 Dec 30. PMID: 26715192.

Acknowledgements
This work has been supported by NIH grants K99 DA051315, P41 EB031771, R00 AG062230,
R01 EB016089, R01 EB023963, R21AG060245
