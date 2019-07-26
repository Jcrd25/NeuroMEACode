function [DownSampledData, UpdatedSamplingRate] = MEA_downsampler(MEAData, SamplingRate, DesiredSamplingRate)
%Use this function to downsample the data recroded from a MEA. It will the
%recorded frequency to a desired frequency. Also known as decimating.
%
%Inputs
%   MEAlData            - The voltage trace of the desired channel in µV
%   SamplingRate        - Sampling rate used in the data adquisition
%   DesiredSamplingRate - Desired end sampling rate
%Output
%   DownSampledData     - The Voltage trace downsampled to the desired
%                         rate
%   UpdatedSamplingRate - The new sampling rate at which the downsampled
%                         data was downsampled to
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN: May 13, 2017
%LAST MODIFIED: May 13, 2017
%v1.0


DownSampledValue = floor(SamplingRate/DesiredSamplingRate); 
% DownSampledData = zeros(
%Determine the amount of channels
SizeVariable = size(MEAData.RawData.VoltageTrace);
NumChannels = SizeVariable(2);
TDownSampledData = [];

%Downsample the data to the desired frequency for all Channels
for i = 1:NumChannels
    VoltageTrace =  MEAData.RawData.VoltageTrace(:,i);             
    DownSampledChannelData = VoltageTrace(1:DownSampledValue:end,:);    
    %Add the ith channel data to the downsampled data set.
        %OPTIMIZATION:It works through concatanation, this is not the fastest 
        %nor RAM saving way to do it
    TDownSampledData = [TDownSampledData DownSampledChannelData];
end

%Set the new sampling rate
DownSampledData     = TDownSampledData;
UpdatedSamplingRate = SamplingRate/DownSampledValue;             
 
