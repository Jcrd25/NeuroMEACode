function MEADataPlotter(MEAData,ChannelToPlot,TimeSegmentDesired)
%This function will plot all channels of thedata adquired from the 60
%pMEA200 30 iR from multichannel systems. It plots the channels in a 8x8
%plot. The y axis is in Potential in µV and the x axis is time in s.
%
%Inputs
%    MEAData- this is struct generated using the
%Optional inputs
%    TimeSegmentDesired - it is index of the desired tiume segment to plot.
%                         Must take into consideration the sampling rate.  
%                         Default is to plot the entire trace.
%Outputs
%   Plot of the voltage trace over time
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN: April 13, 2017
%LAST MODIFIED: May 23, 2017
%v1.0

%Set the default parameters if the parameters are not defined 
if ~isfield(MEAData.RawData, 'DownSamplingRate')
     %if Setting does not exist, so default it 20 kHz
      SamplingRate = 20000;
      fprintf('ATTENTION: Data not downsampled.\n Using the default setting of 20 kHz\n')
else
    SamplingRate = MEAData.RawData.DownSamplingRate; 
end

if ~isfield(MEAData.RawData, 'DownSampledData')
     %if Down Sampled Data not present use the raw trace
      VoltageTrace = MEAData.RawData.VoltageTrace;
      fprintf('ATTENTION: Data not downsampled.\n Using the raw voltage trace instead\n')
else
    VoltageTrace = MEAData.RawData.DownSampledData;
end

%Define the MEA channel order. It is required in order to plot the channel
%in the correct subplot spot in the final graph.
if MEAData.Settings.MeaName == '60pMEA200/30iR'
    MEAChannelOrder = [0 ,21,31,41,51,61,71, 0,...
                       12,22,32,42,52,62,72,82,...
                       13,23,33,43,53,63,73,83,...
                       14,24,34,44,54,64,74,84,...
                       15,25,35,45,55,65,75,85,...
                       16,26,36,46,56,66,76,86,...
                       17,27,37,47,57,67,77,87,...
                       0 ,28,38,48,58,68,78, 0]';
                   
elseif MEAData.Settings.MeaName == '60pMEA100/30iR'
    MEAChannelOrder = [11,21,31,41,51,61,71,81,91,101 ...
                       12,22,32,42,52,62,72,82,92,102 ...
                       13,23,33,43,53,63,73,83,93,103 ...
                       14,24,34,44,54,64,74,84,94,104 ...
                       15,25,35,45,55,65,75,85,95,105 ...
                       16,26,36,46,56,66,76,86,96,106 ]';
else
    fprintf('WARNING: Not a known MEA layout')
end


% %If you want to plot the filtered data uncomment the next lines
% %Filter the data using a low pass Filter

% %Design the low pass filter
% Fpass = 500;
% Fstop = 700;
% Apass = 0.5;
% Astop = 65;
% Fs = SamplingRate;
% 
% d = designfilt('lowpassiir', ...
%   'PassbandFrequency',Fpass,'StopbandFrequency',Fstop, ...
%   'PassbandRipple',Apass,'StopbandAttenuation',Astop, ...
%   'DesignMethod','butter','SampleRate',Fs);

% %Bandpass filter
% Forder = 20;
% Freq1  = 30;
% Freq2  = 120;
% Fs     = SamplingRate;
% 
% d = designfilt('bandpassiir','FilterOrder',Forder, ...
%          'HalfPowerFrequency1',Freq1,'HalfPowerFrequency2',Freq2, ...
%          'SampleRate',Fs);
%      
% %Filter the data
% VoltageTrace = filtfilt(d,VoltageTrace);


lenghtData  = size(VoltageTrace,1);
RecordingTime = 0:1/SamplingRate:(lenghtData-1)/SamplingRate;
LabelsChan = zeros(length(MEAData.Settings.InfoChannel.Label),1);

if ~exist('TimeSegmentDesired','var')
     %if TimeSegmentDesired does not exist, so default it all
      TimeSegmentDesired = 1:size(VoltageTrace,1);
end

if ~exist('ChannelToPlot','var')
     %if ChannelToPlot does not exist, default set it to all channels
      ChannelToPlot = 1:60;
end

%This is the for the MEA200/30iR layout
if MEAData.Settings.MeaName == '60pMEA200/30iR'
    for ii =1:length(MEAData(1).Settings.InfoChannel.Label)
        %Change the label of the reference channel to 15
        if ii == 15
            LabelsChan(ii) = 15;
        else
            LabelsChan(ii) = str2double(MEAData(1).Settings.InfoChannel.Label(ii));
        end
    end

    LabelsChan = LabelsChan';
    figure
    if length(ChannelToPlot) > 1
        hold on
        %For each channel find the data set and plot it in the correct subplot
        for iii = 1:length(LabelsChan)
            %this part searched for the data of the current channel being plotted
            for iv = 1:length(MEAChannelOrder)
                %if the label is in the correct spot, plot the data in the
                %subplot
                if LabelsChan(iii) == MEAChannelOrder(iv)
                    subplot(8,8,iv)
                    plot(RecordingTime(TimeSegmentDesired), VoltageTrace(TimeSegmentDesired,iii))
                    %change title
                    P1 = sprintf('%d', LabelsChan(iii));
                    title(P1)
                    %set x and y limits to be the same for all channels. Can be
                    %changed if needed.
                    ylim([-100, 100])
                    xlim([RecordingTime(TimeSegmentDesired(1)), RecordingTime(TimeSegmentDesired(end))])
                end
            end
        end
        hold off
    else
        %if you only want to plot one channel the code will do the following
        for iii = 1:length(LabelsChan)
            %search the index of the channel you want to plot
                if LabelsChan(iii) == ChannelToPlot
                    plot(RecordingTime(TimeSegmentDesired), VoltageTrace(TimeSegmentDesired,iii))
                    P1 = sprintf('Channel %d', LabelsChan(iii));
                    title(P1)            
                    xlim([RecordingTime(TimeSegmentDesired(1)), RecordingTime(TimeSegmentDesired(end))])
                    xlabel('Time (s)')
                    ylabel('Voltage (\muV)')
                end
        end
    end

%This is the for the 60 channel perfused MEA100/30iR layout
elseif MEAData.Settings.MeaName == '60pMEA100/30iR'
    for ii =1:length(MEAData(1).Settings.InfoChannel.Label)
        %Change the label of the reference channel to 14
        if ii == 15
        LabelsChan(ii) = 14;
        else
        LabelsChan(ii) = str2double(MEAData(1).Settings.InfoChannel.Label(ii));
        end
    end

    LabelsChan = LabelsChan';
    figure
    if length(ChannelToPlot) > 1
        hold on
        %For each channel find the data set and plot it in the correct subplot
        for iii = 1:length(LabelsChan)
            %this part searched for the data of the current channel being plotted
            for iv = 1:length(MEAChannelOrder)
                %if the label is in the correct spot, plot the data in the
                %subplot
                if LabelsChan(iii) == MEAChannelOrder(iv)
                    subplot(6,10,iv)
                    plot(RecordingTime(TimeSegmentDesired), VoltageTrace(TimeSegmentDesired,iii))
                    %change title
                    P1 = sprintf('%d', LabelsChan(iii));
                    title(P1)
                    %set x and y limits to be the same for all channels. Can be
                    %changed if needed.
                    ylim([-100, 100])
                    xlim([RecordingTime(TimeSegmentDesired(1)), RecordingTime(TimeSegmentDesired(end))])
                end
            end
        end
        hold off
    else
        %if you only want to plot one channel the code will do the following
        for iii = 1:length(LabelsChan)
            %search the index of the channel you want to plot
                if LabelsChan(iii) == ChannelToPlot
                    plot(RecordingTime(TimeSegmentDesired), VoltageTrace(TimeSegmentDesired,iii))
                    P1 = sprintf('Channel %d', LabelsChan(iii));
                    title(P1)            
                    xlim([RecordingTime(TimeSegmentDesired(1)), RecordingTime(TimeSegmentDesired(end))])
                    xlabel('Time (s)')
                    ylabel('Voltage (\muV)')
                end
        end
    end
end


    

