
clear;
clc;
close all;
format long;

% Load original dataset
rawDatasetFile = 'RockFallFull.mat';
load(rawDatasetFile);

% Define prefix for output file name
preName = 'RockFallExt10';


%% 


% Variables name string (for output file)
header = {'y','VertDrop','Alpha','Kn','Kt','Theta'}; % 'AvgImpact','Rough'

% Choose variables to include:
% [1 2 3 4 5 6] includes all variables
nHeader = [1 2 3 4 5 6];


% Include meta-variables:
% 1: YES
% 0: NO
extraInput = 1;

% Size of the reduced dataset (percentage of the original dataset)
% Example: 0.1 -> 10% of the size
dataPerc = 0.1;

% Size of the training subset (percentage of the reduced dataset)
% Example: 0.8 -> 80% of the reduced dataset used for training
pTrain = 0.8;

% The maximum acceptable difference for the variance between original
% and reduced dataset is set to: 
% Example: 0.01 -> 1% maximum acceptable difference
maxDif = 0.01;



%figure,plot(RockFallFull(:,1));

%% Create data set
nTotal = size(RockFallFull,1);
nDataRed = round(nTotal*dataPerc);      % Calculate number of samples in the reduced size dataset
nTrain = round(nDataRed*pTrain);        % Calculate number of samples for training subset
nTest = nDataRed - nTrain;              % Number of samples for the testing subset

varDataFull = var(RockFallFull(:,1));   % Variance of the original dataset

cntTemp = 0;

scoreVar = 1;                           
% Verify if the new smaller data set has the same variance as the original
% data set. 
while scoreVar > maxDif
    randFull = randperm(nTotal);
    RockFallReduc = RockFallFull(randFull(1:nDataRed),nHeader);
    
    
    varDataReduc = var(RockFallReduc(:,1)); % Variance of the new reduced size dataset
    
    scoreVar = abs(1 - varDataReduc/varDataFull);   % Percentage of difference between original and new datasets
    cntTemp = cntTemp + 1;                          % Counting failures until success
end
    
% Include meta-variables if requested
if extraInput == 1
    aux = [RockFallReduc(:,2:end).^2 ...
        RockFallReduc(:,2:end).^3 ...
        RockFallReduc(:,2:end).^(1/2) ...
        RockFallReduc(:,2:end).^(1/3) ...
        RockFallReduc(:,2:end).^(-1/2) ...
        RockFallReduc(:,2:end).^(-1/3) ...
        RockFallReduc(:,2:end).^(-1)...
        RockFallReduc(:,2:end).^(-2)...
        RockFallReduc(:,2:end).^(-3)];
    
    RockFallReduc = [RockFallReduc aux];
    
    % Include meta-variables names in the header (string)
    nFeat = size(aux,2);
    nFeatOrg = size(nHeader,2);
    newHeader = {'y'};
    for i = 1:nFeat+nFeatOrg-1
        newHeader{i+1} = strcat('x',num2str(i));
    end
    header = newHeader;  
    

end

% Split subsets
dataTrain = RockFallReduc(1:nTrain,:);
dataTest = RockFallReduc(nTrain+1:end,:);


%% Save data

strTrain = strcat(preName,'Train.csv');
writecell(header,strTrain);
writematrix(dataTrain,strTrain,'WriteMode','append');

strTest = strcat(preName,'Test.csv');
writecell(header,strTest);
writematrix(dataTest,strTest,'WriteMode','append');