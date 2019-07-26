 function PeakData = MEA_Spike_Detector_v2(MEAData, ChannelNumber,ShowDetectedPeaks)
%This function is meant to be used to detect Spikes in recordings using the
%multichannel systems 60 channel MEA. 
%
%Inputs
%   MEAData           - Data
%   ChannelNumber     - Desired Channel Number to get the Spikes
%   ShowDetectedPeaks - Y or N input used to choose if the graphs with the
%                       detected peaks should be shown
%Outputs
%   PeakData          - Data about the Spike peaks
%       SpikePeakTime    = Time of Peaks;
%       SpikePeakVoltage = Voltage value at the peak (Amplitude);
%       SpikePeakLoc     = Index number of the peak;
%       Time             = Time;
%       SpikeNum         = Number of the Spikes detected;
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN:  November 3, 2016
%LAST MODIFIED: February 6, 2018
%v1.0
%Modified from SpikeDetection
%WORK IN PROGRESS

%If variables and settings are not specified use default settings
if ~exist('ChannelNumber','var')
     %if ChannelNumber does not exist, default set it to 1st channel
      ChannelNumber = 1;
end

if ~exist('ShowDetectedPeaks','var')
     %if ChannelNumber does not exist, default set it to Y
      ShowDetectedPeaks = 'Y';
end

SamplingRate = MEAData.RawData.SamplingRate;
VoltageTrace = MEAData.RawData.VoltageTrace(:,ChannelNumber);


%May want to include a step to filter the signal
%**May want to remove lower frequency

%Define the parameter to be used for the High Pass Filter
%High Pass
Fstop = 75;
Fpass = 100;
Astop = 65;
Apass = 0.5;
Fs    = SamplingRate;


d = designfilt('highpassiir','StopbandFrequency',Fstop ,...
  'PassbandFrequency',Fpass,'StopbandAttenuation',Astop, ...
  'PassbandRipple',Apass,'SampleRate',Fs,'DesignMethod','butter');

%Filter the Data
VoltageTrace = filtfilt(d, VoltageTrace);

%Calculate the time values
RecordingTime = (length(VoltageTrace))/SamplingRate; %in seconds
Time = (1/SamplingRate):1/SamplingRate:RecordingTime;

%Calculate the Root Mean Square of the signal
RMS_Signal = rms(VoltageTrace);

%Set threshold to be 8 times the Root Mean Square
SetThreshold = 6*RMS_Signal;


%Used to determine how much to use as minPeakProminence
   %which is how many units (of Voltage) does it have to drop off
   %before detecting the next sweep.
%    MinPeakProminence = 20;
   
   %MIGHT HAVE TO ADJUST THIS TO THE POSIBILITY OF MULTIPLE UNITS
   %Also called Dead Time set to 3 ms
   %MinPeakDistance is the minumun time between peaks.
   
   MinPeakDistance=  3 * (1/1000) * SamplingRate ;

   %Detect the location of all peaks in the data that are above the
   %threshold. 
   
   [SpikePeaksDetected, LocOfPeaks]=findpeaks(abs(VoltageTrace),'MinPeakHeight', SetThreshold,'MinPeakDistance',MinPeakDistance);
   %Determine the amount of peaks determined
   SizeofPeaksDetected = size(SpikePeaksDetected);
   
   if ShowDetectedPeaks == 'Y'
       
       %To check code confirm data with Clampfit Analysis
       fprintf('Detected %d spikes\n', SizeofPeaksDetected(1))
       
       %Plot the data with the peaks identified and visually inspect it.
       PlotPeakDetection(Time,VoltageTrace,LocOfPeaks)
       ThresholdPlotLine = ones(1,length(Time));
       ThresholdPlotLine = ThresholdPlotLine*SetThreshold;
       plot(Time, ThresholdPlotLine,'--r');
       plot(Time, ThresholdPlotLine*(-1),'--r');
   elseif ShowDetectedPeaks == 'N'
   else
       fprintf('ERROR\n **REVISE CODE**')
   end

  
  
   SpikePeakTime    = Time(LocOfPeaks);
   SpikePeakVoltage = VoltageTrace(LocOfPeaks);
   SpikePeakLoc     = LocOfPeaks;
   
   %Define the PeakData Struct
   %**This step is a bit redundant**
   PeakData.SpikePeakTime    = SpikePeakTime;
   PeakData.SpikePeakVoltage = SpikePeakVoltage;
   PeakData.SpikePeakLoc     = SpikePeakLoc;
   PeakData.Time             = Time;
   PeakData.SpikeNum         = SizeofPeaksDetected;
     
   