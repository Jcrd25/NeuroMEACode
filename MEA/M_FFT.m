function  [Power, Frequency] = M_FFT(ChannelData, SamplingRate)
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
%WRITTEN: Spetember 15, 2017
%LAST MODIFIED: September 15, 2017
%v1.0

%Determine the length of the recording
N = length(ChannelData);
%Calculate the Fast Fourier Transform
ChannelData_dft = fft(ChannelData);
%Take the positive values only
ChannelData_dft = ChannelData_dft(1:N/2+1);
%Caculate the power
psd_ChannelData = (1/(SamplingRate*N)) * abs(ChannelData_dft).^2;
psd_ChannelData(2:end-1) = 2*psd_ChannelData(2:end-1);
%define the frequency for x-axis
freq = 0:SamplingRate/length(ChannelData):SamplingRate/2;

Power     = psd_ChannelData;
Frequency = freq;