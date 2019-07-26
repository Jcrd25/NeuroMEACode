%MEA Analysis Master Script
%Import Data from HDF5 File
%Select HDF5 File and create a struct that contains the files to be
%analyzed
%You may select mulitple files but then you need to adjust the following
%script segments to analyze the desired channels at a time.
  
HDF5DataFile = uipickfiles;    

MEAData = HDF5_MEA_data_structuremaker_v2(HDF5DataFile);
%%
% for ChannelNumber = [1:60] 
%     MEA_SpectogramPlotter(MEAData, ChannelNumber, 1000,'Kainate')
%     set(findall(gca, '-property', 'FontSize'),'FontSize',24);
% end
%%
%Plot the raw trace

%Comment in the desired filter to be used

%Bandpass filter
SamplingRate = MEAData(1).RawData.SamplingRate;
Forder = 20;
Freq1  = 20;

Freq2  = 60;
Fs     = SamplingRate;

d = designfilt('bandpassiir','FilterOrder',Forder, ...
         'HalfPowerFrequency1',Freq1,'HalfPowerFrequency2',Freq2, ...
         'SampleRate',Fs);
     
%Lowpass filter      
% Fpass = 100;
% Fstop = 120;
% Apass = 0.5;
% Astop = 65;
% Fs = SamplingRate;
% 
% d = designfilt('lowpassiir', ...
%   'PassbandFrequency',Fpass,'StopbandFrequency',Fstop, ...
%   'PassbandRipple',Apass,'StopbandAttenuation',Astop, ...
%   'DesignMethod','butter','SampleRate',Fs);

for ChannelNumber = [25]
    
    for stage = [LastBaseline, LastKainate, LastDrug]
        figure
        SamplingRate = MEAData(stage).RawData.SamplingRate;
        y = MEAData(stage).RawData.VoltageTrace(:,ChannelNumber);
        x =  0:1/SamplingRate:(length(y)-1)/SamplingRate;
        %Comment out if not using a filter
        y= filtfilt(d,y);
        plot(x,y)
        ylim([-20 20])
        xlim([29 29.5])
    end
end
%%
%Power Spectrum at the desired file

Plot = 'S';
%Select the file index of the file you want to look at

%Define the three file indexes that are going to be used in order to
%calculate gamma power.
tic
LastBaseline = 1;
LastKainate  = 2;
LastDrug     = 3; 
FreqUpper = 100;
FreqLower = 15; 
LineNum = 0;
Color = {'k','b','r'};
for ChannelNumber = 1:60
    if Plot == 'N'
%         close all
        
    elseif Plot == 'S'
        figure
        hold on
    end
    
    for ii = [LastBaseline,LastKainate,LastDrug]
        LengthofRecording = length(MEAData(ii).RawData.VoltageTrace(:,ChannelNumber));
        TimeWindow = 2; %in seconds
        SamplingRate = MEAData(ii).RawData.SamplingRate;
        NumWindows = floor((LengthofRecording/SamplingRate)/TimeWindow);
        Indexs = 1:TimeWindow*SamplingRate:LengthofRecording;
        for DataWindows = 1:NumWindows
            [Complete_Power, Frequency] = M_AverageFFT( MEAData(ii).RawData.VoltageTrace(Indexs(DataWindows):Indexs(DataWindows+1),ChannelNumber), MEAData(ii).Settings.InfoChannel.Label(ChannelNumber), MEAData(ii).Settings.FileName, MEAData(ii).RawData.SamplingRate, Plot,Color{ii});            
        end
%         set(findall(gca, 'Type', 'Line'),'LineWidth',6);
        set(findall(gca, '-property', 'FontSize'),'FontSize',24);
        
        %Calculate the area under the curve for the data in the 
        [GammaPower,HighGamma,LowGamma] = MEA_GammaPCalc_v2(MEAData(ii),ChannelNumber);
        
        MEAData(ii).AnalyzedData(ChannelNumber).Gamma.TotalPower = GammaPower;
        MEAData(ii).AnalyzedData(ChannelNumber).Gamma.LowGamma   = LowGamma;
        MEAData(ii).AnalyzedData(ChannelNumber).Gamma.HighGamma  = HighGamma;
        
        MEAData(ii).AnalyzedData(ChannelNumber).AllPower         = (Complete_Power);
        MEAData(ii).AnalyzedData(ChannelNumber).Frequency        = Frequency;
        
        
        
    end
    if Plot == 'N'
    elseif Plot == 'S'
%         legend('Baseline','Kainate','Kainate + MK801')
        hold off
    end
    
end
toc
%%
%Sonogram
TimeCounter = 0;
for ChannelNumber = 1:60
    
    z=[];
     for ii = [LastBaseline,LastKainate,LastDrug]
        LengthofRecording = length(MEAData(ii).RawData.VoltageTrace(:,ChannelNumber));
        TimeWindow = 2; %in seconds
        SamplingRate = MEAData(ii).RawData.SamplingRate;
        NumWindows = floor((LengthofRecording/SamplingRate)/TimeWindow);
        Indexs = 1:TimeWindow*SamplingRate:LengthofRecording;
        for DataWindows = 1:NumWindows
            TimeCounter = TimeCounter+1;
            [Complete_Power, Frequency] = M_AverageFFT( MEAData(ii).RawData.VoltageTrace(Indexs(DataWindows):Indexs(DataWindows+1),ChannelNumber), MEAData(ii).Settings.InfoChannel.Label(ChannelNumber), MEAData(ii).Settings.FileName, MEAData(ii).RawData.SamplingRate, Plot,Color{ii});            
            for FreqI = 1:100
            z(TimeCounter,FreqI) = Complete_Power(FreqI);
            end
        end
        x = 0:TimeWindow:TimeCounter*TimeWindow;
        y = 1:100;
    imagesc(x,y,z)
    ylabel('Frequency (Hz)')
    xlabel('Time (min)')
    c = colorbar;
    ylabel(c, 'Fold Change')
    colormap jet
    ylim([0 100])
%     caxis([0 7]);
    ham = gca;
    set(ham, 'YDir', 'normal')
    set(findall(gca, '-property', 'FontSize'),'FontSize',24); 
    end
end
%%
%Normalize the data to the baseline power
for ChannelNumber = 1:60
    for ii = [LastKainate LastDrug]
        MEAData(ii).AnalyzedData(ChannelNumber).NormData = M_NormGamma(MEAData(LastBaseline).AnalyzedData(ChannelNumber).AllPower, MEAData(ii).AnalyzedData(ChannelNumber).AllPower);
    end
end
%%
%**Not working properly**
%Do the preliminary automatic Gamma detection
GammaPresence     = zeros(1, 60);
LowGammaPresence  = zeros(1, 60);
HighGammaPresence = zeros(1, 60);

for ChanNum = 1:60
%     Determine if there is Gamma
    if MEAData(LastKainate).AnalyzedData(ChanNum).NormalizedGamma.TotalPower
        GammaPresence(ii) = 1;
    end
%     Determine if there is low gamma
    if MEAData(LastKainate).AnalyzedData(ChanNum).Gamma.LowGamma
        LowGammaPresence(ii) = 1;
    end
%     Determine if there is high gamma
    if MEAData(LastKainate).AnalyzedData(ChanNum).Gamma.HighGamma
        HighGammaPresence(ii) = 1;
    end
    
end
    %%
    %plot the data
    for ChannelNumber = [1:60];
%     close allcloseclose 
    figure
    hold on
    for ii = [LastBaseline, LastKainate]
        x = MEAData(ii).AnalyzedData(ChannelNumber).Frequency(1:101);
        y = 10*log10(MEAData(ii).AnalyzedData(ChannelNumber).AllPower(1:101));
        
        plot(x,y)
        xlabel('Frequency (Hz)')
        ylabel('Fold increase from baseline')
        
        ChannelNumberDouble = str2double(MEAData(1).Settings.InfoChannel.Label(ChannelNumber));
        P1 = sprintf('Channel %d', ChannelNumberDouble);  
        title(P1)

        set(findall(gca, 'Type', 'Line'),'LineWidth',6);
        set(findall(gca, '-property', 'FontSize'),'FontSize',24);
    end
    legend('Baseline','Kainate','Kainate + Bicuculline')
    hold off
end
%
%%
% Obtain The labels of the corresponding channels for the brain regiojn
% interested in
LastBaseline = 1;
LastKainate  = 2;
LastDrug     = 3; 

FreqUpper = 60;
FreqLower = 15; 
%Get the channel labels 
LabelsofInterest = [];
for ii = 1:length(MEAData(1).Settings.InfoChannel.Label)
    for  n = 1:length(ChannelsofInterest)
        %Change the label of the reference channel to 15
   
        
        if ChannelsofInterest(n) == str2double(MEAData(1).Settings.InfoChannel.Label(ii))
          LabelsofInterest(n) = ii;  
        end
    end
end
%%
%Plot normalized power spectrum of the desired channels of interest
c=1;
X2 =[];
Y2 =[];
for ChannelNumber = [LabelsofInterest]
%     close all
    figure
    hold on
    for ii = [2,3]
        x = MEAData(ii).AnalyzedData(ChannelNumber).Frequency(1:100);
        y = MEAData(ii).AnalyzedData(ChannelNumber).NormData(1:100);
        plot(x,y)
        xlabel('Frequency (Hz)')
        ylabel('Fold increase from baseline')
        
        ChannelNumberDouble = str2double(MEAData(1).Settings.InfoChannel.Label(ChannelNumber));
        P1 = sprintf('Channel %d', ChannelNumberDouble);  
        title(P1)
        set(findall(gca, 'Type', 'Line'),'LineWidth',6);
        set(findall(gca, '-property', 'FontSize'),'FontSize',24);
    end
    ylim([0,10])
    hold off
end
%%
%Calculate the different oscillation parameters for the kainate induced
%oscillations
c=1;
X2 =[];
Y2 =[];
condition = LastKainate;

for ChannelNumber = [LabelsofInterest];
    for ii = [2,3]
        x = MEAData(ii).AnalyzedData(ChannelNumber).Frequency(1:101);
        y = MEAData(ii).AnalyzedData(ChannelNumber).NormData(1:101);
        X2(c,:) = x(11:101);
        Y2(c,:) = y(11:101);
        c=c+1;

    end
end
%for exporting the individual changes in power for each frequency
Y2 = Y2';
X2 = X2';

%Calculate the Q value and determine the peak frequency
counter = 1;
EENLG = [];
EEQ   = [];
EEPF  = [];
for ChannelNumber = [LabelsofInterest]
    for ii = [condition]
        [QValue, PeakFreq] = M_QvalueCalculator(MEAData(ii).AnalyzedData(ChannelNumber).NormData, MEAData(ii).AnalyzedData(ChannelNumber).Frequency, FreqUpper, FreqLower);
        
        MEAData(ii).AnalyzedData(ChannelNumber).QValue        = QValue;
        MEAData(ii).AnalyzedData(ChannelNumber).PeakFrequency = PeakFreq;
        %Easy Export Q value, Peak Freq and the power at desired frequency
        EEQ(counter)  = QValue ;
        EEPF(counter) = PeakFreq;
        EENLG(counter,1) = MEAData(ii).AnalyzedData(ChannelNumber).Gamma.LowGamma / MEAData(LastBaseline).AnalyzedData(ChannelNumber).Gamma.LowGamma;
        counter =counter +1;
    end
end

EEQ = EEQ';
EEPF = EEPF';

%Due to bicuculline getting rid of the oscillations, peak power and
%frequency are not calculated, this might change if instead the drug alters
%the oscillations 
counter = 1;
condition = LastDrug;

for ChannelNumber = [LabelsofInterest]
    for ii = [condition]
        EENLG(counter,2) = MEAData(ii).AnalyzedData(ChannelNumber).Gamma.LowGamma / MEAData(LastBaseline).AnalyzedData(ChannelNumber).Gamma.LowGamma;
        counter =counter +1;
    end
end

%%
%Spike detection algorithm
for ii = [LastBaseline,LastKainate, LastDrug]
    for ChannelNumber = 1:60
        PeakData(ChannelNumber) = MEA_Spike_Detector_v2(MEAData(ii), ChannelNumber,'N');
        
        %Save the desired data
        MEAData(ii).PeakData(ChannelNumber).SpikePeakTime    = PeakData(ChannelNumber).SpikePeakTime;
        MEAData(ii).PeakData(ChannelNumber).SpikePeakVoltage = PeakData(ChannelNumber).SpikePeakVoltage;
        MEAData(ii).PeakData(ChannelNumber).SpikePeakLoc     = PeakData(ChannelNumber).SpikePeakLoc;            
        MEAData(ii).PeakData(ChannelNumber).SpikeNum         = PeakData(ChannelNumber).SpikeNum;
    end
end


%Spike Interspike calculator

for ii = [LastBaseline, LastKainate, LastDrug]
    for ChannelNumber = 1:60
        MEAData(ii).PeakData(ChannelNumber).ISI = MEA_ISICalculator(MEAData(ii).PeakData, ChannelNumber) ;
    end
end

%Calculate the spike frequecies
for ii = [LastBaseline, LastKainate, LastDrug]
    for ChannelNumber = 1:60
        Time = length(MEAData(ii).RawData.VoltageTrace(ChannelNumber))/...
            MEAData(ii).RawData.SamplingRate;
        [MEAData(ii).PeakData(ChannelNumber).InstFreq, ...
            MEAData(ii).PeakData(ChannelNumber).GlobFreq]= ...
            MEA_SpikeFreqCalculator(MEAData(ii).PeakData, ChannelNumber,Time) ;
    end
end

%Calculate the different burst parameters
for ii = [LastBaseline, LastKainate, LastDrug]
    for ChannelNumber = 1:60
        BurstData = MEA_BurstParamCalc(MEAData(ii), ChannelNumber);
        
        %Save desired data
        MEAData(ii).BurstData(ChannelNumber).NumBurst         = BurstData.NumBurst;
        MEAData(ii).BurstData(ChannelNumber).FrequencyofBurst = BurstData.FrequencyofBurst;
        MEAData(ii).BurstData(ChannelNumber).AveBL            = BurstData.AveBurstLength; 
        MEAData(ii).BurstData(ChannelNumber).AveSPB           = BurstData.AveSPB; 
    end
end

%%
%Easy Export the Burst data and spike frequency data from Channels of 
%Interest for kainate
counter =1;
condition = LastKainate;

%Set ii as the desired condition
EEABL  = [];
EEASPB = [];
EENBF  = []';
for ChannelNumber = [LabelsofInterest]
    %Calculate fold change
    for ii  = [condition]
        ChangeBL = (MEAData(ii).BurstData(ChannelNumber).AveBL -...
            MEAData(LastBaseline).BurstData(ChannelNumber).AveBL)/...
            MEAData(LastBaseline).BurstData(ChannelNumber).AveBL;
        ChangeSPB = (MEAData(ii).BurstData(ChannelNumber).AveSPB - ...
            MEAData(LastBaseline).BurstData(ChannelNumber).AveSPB)/...
            MEAData(LastBaseline).BurstData(ChannelNumber).AveSPB;
        ChangeBF = (MEAData(ii).BurstData(ChannelNumber).FrequencyofBurst - ...
            MEAData(LastBaseline).BurstData(ChannelNumber).FrequencyofBurst)/...
            MEAData(LastBaseline).BurstData(ChannelNumber).FrequencyofBurst;
        ChangeISF = (MEAData(ii).PeakData(ChannelNumber).InstFreq -...
            MEAData(LastBaseline).PeakData(ChannelNumber).InstFreq)/...
            MEAData(LastBaseline).PeakData(ChannelNumber).InstFreq;
        ChangeGSF = (MEAData(ii).PeakData(ChannelNumber).GlobFreq -...
            MEAData(LastBaseline).PeakData(ChannelNumber).GlobFreq)/...
            MEAData(LastBaseline).PeakData(ChannelNumber).GlobFreq;
    end
    %Easy Export
    for ii = [condition]
        EEABL(counter,1)  = ChangeBL;
        EEASPB(counter,1) = ChangeSPB;
        EENBF(counter,1)  = ChangeBF;
        EEISF(counter,1)  = ChangeISF;
        EEGSF(counter,1)  = ChangeGSF;
        counter =counter +1;
    end
end


%Export the data for the bicuculline 
%Easy Export the Burst data from Channels of Interest for Bicuculline
counter =1;
condition = LastDrug;

for ChannelNumber = [LabelsofInterest]
    %Calculate fold change
    for ii  = [condition]
        ChangeBL = (MEAData(ii).BurstData(ChannelNumber).AveBL-...
            MEAData(LastBaseline).BurstData(ChannelNumber).AveBL)/...
            MEAData(LastBaseline).BurstData(ChannelNumber).AveBL;
        ChangeSPB = (MEAData(ii).BurstData(ChannelNumber).AveSPB - ...
            MEAData(LastBaseline).BurstData(ChannelNumber).AveSPB)/...
            MEAData(LastBaseline).BurstData(ChannelNumber).AveSPB;
        ChangeBF = (MEAData(ii).BurstData(ChannelNumber).FrequencyofBurst - ...
            MEAData(LastBaseline).BurstData(ChannelNumber).FrequencyofBurst)/...
            MEAData(LastBaseline).BurstData(ChannelNumber).FrequencyofBurst;
        ChangeISF = (MEAData(ii).PeakData(ChannelNumber).InstFreq -...
            MEAData(LastBaseline).PeakData(ChannelNumber).InstFreq)/...
            MEAData(LastBaseline).PeakData(ChannelNumber).InstFreq;
        ChangeGSF = (MEAData(ii).PeakData(ChannelNumber).GlobFreq -...
            MEAData(LastBaseline).PeakData(ChannelNumber).GlobFreq)/...
            MEAData(LastBaseline).PeakData(ChannelNumber).GlobFreq;
    end
    %Easy Export
    for ii = [condition]
        EEABL(counter,2)  = ChangeBL;
        EEASPB(counter,2) = ChangeSPB;
        EENBF(counter,2)  = ChangeBF;
        EEISF(counter,2)  = ChangeISF;
        EEGSF(counter,2)  = ChangeGSF;
        counter =counter +1;
    end
end

%%
%Save the MEAData Struct
[NAME] = M_MATLABSaver(MEAData(1).Settings.FileName);
NAME = char(NAME);
Filepath = 'H:\Documents\Analyzed_MATLAB\';
Name = strcat(Filepath, NAME);
Name = char(Name);
save(Name, 'MEAData'); 

%%
%Load the previously saved data
x = load(NAME)