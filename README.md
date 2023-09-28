# ml-hea-corrosion-code-data
Data, and code for "machine learning accelerated discovery of corrosion-resistant high-entropy alloys".

- "single-phase-HEA" includes the original and processed experimental dataset, and files for the machine learning models used to predict single-phase formability of high-entropy alloys.

- "eam-sampling" includes scripts to reproduce the fast sampling of high-entropy alloy configurations using embedded atom method (EAM). The generated structure files are also included.

- "first-principles-calculations" includes example scripts to perform high-throughput DFT calculations and how to generate simple representative bulk and surface structures with one to five elements.
The data for each single point calculation is also included.

- "mtp-training" includes the training data and a script to train MTP for AlCrFeNiCo.

- "mtp-simulations" includes scripts that use trained MTP to perform simulations targeting the three corrosion metrics and supporting results shown in SI. The relevant data are also included.

- "generate-sqs-structure" includes the input and output files for generating special quasi-random structures using [atat](https://www.sciencedirect.com/science/article/pii/S0364591613000540).


