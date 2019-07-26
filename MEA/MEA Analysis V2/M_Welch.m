function [Power, Frequency] = M_Welch(ChannelData, SamplingRate)
%This funciton is meant to be used in order to calcualate the power at the
%different corresponding frequencies of the inputed data by using the Fast
%Fourier Transform.
%
%Inputs
%   ChannelData  - The voltage trace of the desired channel in µV 
%   SamplingRate - Sampling rate used in the data adquisition
%Output
%   Power   
%   Frequency
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN: January 9, 2018
%LAST MODIFIED: January 9, 2018
%v1.0

%Determine the length of the recording
N = length(ChannelData);
%Calculate the Power using welch's power spectral density estimate
[ppx,freq] = pwelch(ChannelData,[],[],[], SamplingRate);


Power     = ppx;
Frequency = freq;