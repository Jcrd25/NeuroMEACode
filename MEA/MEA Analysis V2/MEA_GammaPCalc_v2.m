function [GammaPower,HighGamma,LowGamma] = MEA_GammaPCalc_v2(MEAData,ChannelNumber)
%This function is used in order to calculate the power in the frequencies
%in the gamma band range (30-100 Hz)
%Inputs
%   
%Output
%   GammaPower
%
%NOTE:
%Standard frequency ranges based on literature and multichannel systems 
%Delta : 0.5 - 4 Hz
%Theta : 4   - 8 Hz
%Alpha : 8   - 13 Hz
%Beta  : 13  - 30 Hz
%Gamma : 30  - 80 Hz
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN: September 15, 2017
%LAST MODIFIED: Septermber 15, 2017
%v1.0

%Set the default settings if they have not been specified
if ~isfield(MEAData.RawData, 'SamplingRate')
     %if Setting does not exist, so default it 20 kHz
      SamplingRate = 20000;
      fprintf('Data not downsampled.\n Using the default setting of 20 kHz\n')
else
    SamplingRate = MEAData.RawData.SamplingRate;
end

[Power, Frequency ] = M_AverageFFT( MEAData.RawData.VoltageTrace(:,ChannelNumber), MEAData.Settings.InfoChannel.Label(ChannelNumber), MEAData.Settings.FileName, SamplingRate, 'N');


%Define the gamma frequency range
Freq1  = 30;
Freq2  = 100;
Freq1H = 60;
Freq2H = 100;
Freq1L = 30;
Freq2L = 60;

%Calculate the Area under the curve as an estimate of the total power in
%the gamma band
[GammaPower] = M_AUC(Power, Frequency, Freq1 , Freq2 );
[HighGamma]  = M_AUC(Power, Frequency, Freq1H, Freq2H);
[LowGamma]   = M_AUC(Power, Frequency, Freq1L, Freq2L);