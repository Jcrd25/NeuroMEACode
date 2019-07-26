function MEAPowerSpectrumCalc(MEAData, ChannelNumber, Setting, Plot, InterFileInterval)
%This function is meant to be used to calculate the power spectrum of the
%data adquired with the MEA multichannel systems. 
%
%Inputs
%   MEAData       -Data
%   ChannelNumber - Desired Channel Number to get the power spectrum
%   Setting       - Used to determine what function to use
%       Setting 1 - Uses the FFT and plots a Power Spectrum with log
%                   scale (in dB) also known as a Periodogram
%       Setting 2 - Uses the FFT and plots a Power Spectrum without the log
%                   scale( in µV^2) also known as a Periodogram
%       Setting 3 - Generates a spectrogram
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN: April 13, 2017
%LAST MODIFIED: May 23, 2017
%v1.0



% import the required information
% ChannelDataMEA = MEAData.ChannelDataVoltageTrace; %Imports the Voltage trace in µV
% define length and sampleling frequency
% lenghtData  = size(ChannelDataMEA,2);
% NumChannels = size(ChannelDataMEA,1);
%
if ~isfield(MEAData.RawData, 'DownSamplingRate')
     %if Setting does not exist, so default it 20 kHz
      SamplingRate = 20000;
      fprintf('Data not downsampled.\n Using the default setting of 20 kHz\n')
else
    SamplingRate = MEAData.RawData.DownSamplingRate;
end
% RecordingTime = lenghtData/SamplingRate; %seconds

%Check if setting is set
if ~exist('Setting','var')
     %if Setting does not exist, so default it all
      Setting = 1;
      fprintf('No Setting specified.\n Using the default setting\n')
end

%Check if setting is set
if ~exist('Plot','var')
     %if Setting does not exist, so default it all
      Plot = 'Y';
      fprintf('No Setting for Plotting specified.\n Using the default setting\n')
end

% Check the setting to determine which MEA_FFT function is going to be used

if Setting == 1
    for ii = 1:length(MEAData)
        MEA_FFT( MEAData(ii).RawData.DownSampledData(:,ChannelNumber), MEAData(ii).Settings.InfoChannel.Label(ChannelNumber), MEAData(ii).Settings.FileName, SamplingRate)
    end
    
elseif Setting == 2
    for ii = 1:length(MEAData)
       MEA_FFT_NdB( MEAData(ii).RawData.DownSampledData(:,ChannelNumber), MEAData(ii).Settings.InfoChannel.Label(ChannelNumber), MEAData(ii).Settings.FileName, SamplingRate)
    end
    
elseif Setting == 3
    for ii = 1:length(MEAData)
        MEA_FFT_Spectrogram( MEAData(ii).RawData.DownSampledData(:,ChannelNumber), MEAData(ii).Settings.InfoChannel.Label(ChannelNumber), MEAData(ii).Settings.FileName, SamplingRate)
    end
elseif Setting == 4
    for ii = 1:length(MEAData)
        M_AverageFFT( MEAData(ii).RawData.DownSampledData(:,ChannelNumber), MEAData(ii).Settings.InfoChannel.Label(ChannelNumber), MEAData(ii).Settings.FileName, SamplingRate, Plot);
    end
elseif Setting == 5
    for ii = 1:length(MEAData)
        MEA_FFT_SpectCont( MEAData(ii).RawData.DownSampledData(:,ChannelNumber), MEAData(ii).Settings.InfoChannel.Label(ChannelNumber), MEAData(ii).Settings.FileName, SamplingRate, InterFileInterval)
    end
end