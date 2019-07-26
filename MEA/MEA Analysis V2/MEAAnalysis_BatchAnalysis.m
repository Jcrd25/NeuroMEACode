%MEA Analysis Master Script
%Import Data from HDF5 File
%Select HDF5 File and create a struct that contains the files to be
%analyzed
%You may select mulitple files but then you need to adjust the following
%script segments to analyze the desired channels at a time.

NumFilesToStudy = 2;

for ii = 1:NumFilesToStudy
    FileG(ii).HDF5DataFile = uipickfiles;  
end
%%
FileG(1).CA1Channels = VarName1';
FileG(1).CA3Channels = VarName2';

FileG(2).CA1Channels = VarName3';
FileG(2).CA3Channels = VarName4';

% FileG(3).CA1Channels = VarName5';
% FileG(3).CA3Channels = VarName6';
% 
% FileG(4).CA1Channels = VarName7';
% FileG(4).CA3Channels = VarName8';
% 
% FileG(5).CA1Channels = VarName9';
% FileG(5).CA3Channels = VarName10';
% 
% FileG(6).CA1Channels = VarName11';
% FileG(6).CA3Channels = VarName12';

%%
tic
for FG = 1:NumFilesToStudy

HDF5DataFile = FileG(FG).HDF5DataFile;
CA1Channels = FileG(FG).CA1Channels;
CA3Channels = FileG(FG).CA3Channels;
    

FilesToLoad = length(HDF5DataFile);

LastBaseline = 30;
LastKainate  = 60;
LastDrug     = length(HDF5DataFile); 
    
%load the desired HDF5 file adn calculate all of the desired parameters to
%be saved in the data structure
for FileNumber = 1:FilesToLoad
    
    MEAData = HDF5_MEA_data_extractor(HDF5DataFile{FileNumber});
    
    %Calculate the power at the frequencies for the currently loaded HDF5
    %file

    Plot = 'N';
       
    FreqUpper = 100;
    FreqLower = 15; 
    
    %Obtain the relevant data for all of the channels
    for ChannelNumber = [1:60]
        
        [Complete_Power, Frequency] = M_AverageFFT( MEAData.RawData.VoltageTrace(:,ChannelNumber), MEAData.Settings.InfoChannel.Label(ChannelNumber), MEAData.Settings.FileName, MEAData.RawData.SamplingRate, Plot);
        %Calculate the area under the curve for the data in the 
        [GammaPower,HighGamma,LowGamma] = MEA_GammaPCalc_v2(MEAData,ChannelNumber);
        
        %Organize the data in the structure
        MEAData.AnalyzedData(ChannelNumber).Gamma.TotalPower = GammaPower;
        MEAData.AnalyzedData(ChannelNumber).Gamma.LowGamma   = LowGamma;
        MEAData.AnalyzedData(ChannelNumber).Gamma.HighGamma  = HighGamma;
        MEAData.AnalyzedData(ChannelNumber).AllPower         = (Complete_Power);
        MEAData.AnalyzedData(ChannelNumber).Frequency        = Frequency;
        
        %Calculate the Peak data
        PeakData(ChannelNumber) = MEA_Spike_Detector_v2(MEAData, ChannelNumber,'N');
        
        %Save the desired data
        MEAData.PeakData(ChannelNumber).SpikePeakTime    = PeakData(ChannelNumber).SpikePeakTime;
        MEAData.PeakData(ChannelNumber).SpikePeakVoltage = PeakData(ChannelNumber).SpikePeakVoltage;
        MEAData.PeakData(ChannelNumber).SpikePeakLoc     = PeakData(ChannelNumber).SpikePeakLoc;            
        MEAData.PeakData(ChannelNumber).SpikeNum         = PeakData(ChannelNumber).SpikeNum;
        
        %calculate the interspike intervals
        MEAData.PeakData(ChannelNumber).ISI = MEA_ISICalculator(MEAData.PeakData, ChannelNumber) ;
        
        Time = length(MEAData.RawData.VoltageTrace(ChannelNumber))/...
            MEAData.RawData.SamplingRate;
        
        %Calculate spike frequency usind both instantaneous and global
        %definitions
        [MEAData.PeakData(ChannelNumber).InstFreq, ...
            MEAData.PeakData(ChannelNumber).GlobFreq]= ...
            MEA_SpikeFreqCalculator(MEAData.PeakData, ChannelNumber,Time) ;
        
        %Calculate the desired burst data
        BurstData = MEA_BurstParamCalc(MEAData, ChannelNumber);
        
        %Save desired data
        MEAData.BurstData(ChannelNumber).NumBurst         = BurstData.NumBurst;
        MEAData.BurstData(ChannelNumber).FrequencyofBurst = BurstData.FrequencyofBurst;
        MEAData.BurstData(ChannelNumber).AveBL            = BurstData.AveBurstLength; 
        MEAData.BurstData(ChannelNumber).AveSPB           = BurstData.AveSPB; 
        
        %Normalize the gamma data to the baseline
        if FileNumber > LastBaseline
            MEAData.AnalyzedData(ChannelNumber).NormData = ...
                M_NormGamma(MEADataS(LastBaseline).AnalyzedData(ChannelNumber).AllPower, ...
                MEAData.AnalyzedData(ChannelNumber).AllPower);
        end
    end
    %Save the calculated data in the data stucture
    MEADataS(FileNumber).AnalyzedData = MEAData.AnalyzedData;
    MEADataS(FileNumber).Settings     = MEAData.Settings;
    MEADataS(FileNumber).PeakData     = MEAData.PeakData; 
    MEADataS(FileNumber).BurstData    = MEAData.BurstData; 
    
end


%Use predetermined label for the desired brain regions to be analyzed.



LabelsCA3 = [];
LabelsCA1 = [];

%Use the channel numbers to obtain the channel labels or ID
for ii = 1:length(MEAData(1).Settings.InfoChannel.Label)
    CA1Chan = length(CA1Channels);
    CA3Chan = length(CA3Channels);
    
   for  n = 1:CA1Chan
        %Get the labels for CA1
        if CA1Channels(n) == str2double(MEAData(1).Settings.InfoChannel.Label(ii))
          LabelsCA1(n) = ii;  
        end
   end
   for  n = 1:CA3Chan
        %Get the labels for CA1
        if CA3Channels(n) == str2double(MEAData(1).Settings.InfoChannel.Label(ii))
          LabelsCA3(n) = ii;  
        end
   end 
   
end

%Save the channel labels to the data struct
MEADataS(1).Channels.CA1Channels = CA1Channels; 
MEADataS(1).Channels.CA3Channels = CA3Channels;
MEADataS(1).Channels.CA1Labels   = LabelsCA1;
MEADataS(1).Channels.CA3Labels   = LabelsCA3;
toc

%Save the MEAData Struct
[NAME] = M_MATLABSaver(MEADataS(1).Settings.FileName);
NAME = char(NAME);
Filepath = 'H:\Documents\Analyzed_MATLAB\';
Name = strcat(Filepath, NAME);
Name = char(Name);
save(Name, 'MEADataS'); 

clearvars -except FileG NumFilesToStudy
end