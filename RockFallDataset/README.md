# RockFall Reduced Size Dataset

## Requirements
This script requires at least a **MATLAB R2019a**

## Files Description
+ **RockFallFull.mat**: File containing all samples available (original dataset) 
+ **prepareDatasetRockFall.m**: script to generate dataset

## Setting Parameters

The following parameters should be set prior to the execution. All of them have a default value included. If both files (*RockFallFull.mat* and *prepareDatasetRockFall.m*) are located in the same folder, the code is ready to run.

Parameters description:
+ **rawDatasetFile**: Original dataset filepath (if it's located in the same folder, just include the file name) 
+ **preName**: Prefix for output file name
+ **header**: Variables name string from the original dataset (not necessary to modify)
+ **nHeader**: Choose variables from the original dataset to include. Example: [1 2 3 4 5 6] includes all variables available.
+ **extraInput**: Set to 1 if you want to include meta-variables / set to 0 if you want to keep just the variables from the original dataset
+ **dataPerc**: Size of the reduced dataset (percentage of the original dataset). Example: 0.1 -> 10% of the size
+ **pTrain**: Size of the training subset (percentage of the reduced dataset). Example: 0.8 -> 80% of the reduced dataset used for training
+ **maxDif**: maximum acceptable difference for the variance between original and reduced dataset is set to. Example: 0.01 -> 1% maximum acceptable difference

## Program Execution

After setting all variables according to desired output, execute the code. The output files are *-.CSV* files, their name starts with the *preName* previously defined followed by 'Train' or 'Test' sufix.

ATTENTION: Remove the output files before creating new files with the same name. New data will be appended to the previously inserted data.


