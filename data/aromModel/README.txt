
This folder contains the source files used to train SVM models for aromaticity
prediction that are based on the NSPDK graph kernel.

For each model an archive with the molecules including aromaticity information
in SMILES format is provided. This data was used as input to the GGL tool
'aromModelNSPDK' to train and test the according SVM models.

The tool 'aromModel2impl' was used to transform the sparse vector data into the
C++ implementation format needed for integrating the model into the library.

