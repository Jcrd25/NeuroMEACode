function [NAME] = M_MATLABSaver(FileName)
%This code is meant to be used in order to determine the name of the file
%you want to save
%MEAData

SplitFilePathUsed = strsplit(FileName,'\');

%Convert the last cell array that contains the File Name
CharSFPU = char(SplitFilePathUsed(end));

SplitFileName = strsplit(CharSFPU,'_');

%Obtain the date and the Slice number from the file name
Date     = char(SplitFileName(1));
SliceNum = char(SplitFileName(3));

NAME = strcat('MEAData_', Date,'_',SliceNum);