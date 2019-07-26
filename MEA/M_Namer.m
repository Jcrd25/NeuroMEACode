function ExtractedName = M_Namer(FileName)
%This function is meant to be used to extract key information from the
%FileName. This version work for MEA files recorded from the Hippocampus.
%
%Input
%   FileName - The name of the file obtained from the h5 file
%Output
%   ExtractedName - The name meant to be used to label the dataset



%Example:
%FileName      = H:\Documents\DATA\2017_06_13\Slice 1\2017-06-13T13-07-19Hippocampus_Slice1_Baseline.h5
%ExtractedName = Hippocampus_Slice1_Baseline
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN: June 28, 2017
%LAST MODIFIED: June 28, 2017
%v1.0

%This is used to remove the text before the H in Hippocampus
[~,R] = strtok(FileName, 'H');
%Remove the file name extension
[ExtractedName] = strtok(R, '.');



