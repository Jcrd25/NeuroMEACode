function ExtractedName = M_AnalyzedNamer(FileName)
%This function is meant to be used to extract key information from the
%FileName. This version work for MEA files recorded from the Hippocampus.
%
%Input
%   FileName - The name of the file obtained from the h5 file
%Output
%   ExtractedName - The name meant to be used to label the dataset



%Example:
%FileName      = H:\Documents\Analyzed_MATLAB\MEAData_2018-04-19_Slice1.mat
%ExtractedName = MEAData_2018-04-19_Slice1
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN: June 28, 2017
%LAST MODIFIED: June 28, 2017
%v1.0

%This is used to remove the text before the M in MATLAB
[~,R] = strtok(FileName, 'M');
%This is used to remove the text before the M in MEA
[~,R] = strtok(R, 'M');
%Remove the file name extension
[ExtractedName] = strtok(R, '.');