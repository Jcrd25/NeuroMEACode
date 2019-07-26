
%Neuronal Culture Master Script for spike analysis

%Identify the desired files to be analyzed.
HDF5DataFile = uipickfiles;   

%%
Filepath = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\NeuronalCulture\Shigeki Collab\March 8 Culture\';
Filenum = 1;
%%
tic
for ii = Filenum:length(HDF5DataFile)
    
    
    MEAData = HDF5_NC_MEA_data_structuremaker_v2(HDF5DataFile{ii}); 

    %Spike Analyzis extractor
    %Spike detection algorithm

        for ChannelNumber = 1:60


            PeakData(ChannelNumber) = MEA_Spike_Detector_v2(MEAData, ChannelNumber,'N');

            %Save the desired data
            MEANCData.PeakData(ChannelNumber).SpikePeakTime    = PeakData(ChannelNumber).SpikePeakTime;
            MEANCData.PeakData(ChannelNumber).SpikePeakVoltage = PeakData(ChannelNumber).SpikePeakVoltage;
            MEANCData.PeakData(ChannelNumber).SpikePeakLoc     = PeakData(ChannelNumber).SpikePeakLoc;            
            MEANCData.PeakData(ChannelNumber).SpikeNum         = PeakData(ChannelNumber).SpikeNum;
            %Obtain the data that will be stored from the original HDF5 file
            MEANCData.Settings = MEAData.Settings;
            %Save the time varibale
            Time            = length(MEAData.RawData.VoltageTrace(:,ChannelNumber))/...
                             MEAData.RawData.SamplingRate;
            MEANCData.Time(ChannelNumber).RecordingTime  = Time;

        end
%     clear MEAData    

    %Spike Interspike calculator


        for ChannelNumber = 1:60
            MEANCData.PeakData(ChannelNumber).ISI = MEA_ISICalculator(MEANCData.PeakData, ChannelNumber) ;
        end

    %Calculate the spike frequecies

        for ChannelNumber = 1:60

            [MEANCData.PeakData(ChannelNumber).InstFreq, ...
                MEANCData.PeakData(ChannelNumber).GlobFreq]= ...
                MEA_SpikeFreqCalculator(MEANCData.PeakData, ChannelNumber,Time) ;
        end


    %Calculate the different burst parameters

        for ChannelNumber = 1:60
            BurstData = MEA_BurstParamCalc(MEANCData, ChannelNumber);

            %Save desired data
            MEANCData.BurstData(ChannelNumber).NumBurst         = BurstData.NumBurst;
            MEANCData.BurstData(ChannelNumber).FrequencyofBurst = BurstData.FrequencyofBurst;
            MEANCData.BurstData(ChannelNumber).AveBL            = BurstData.AveBurstLength; 
            MEANCData.BurstData(ChannelNumber).AveSPB           = BurstData.AveSPB; 
        end

        ArrayBurstData = NC_ArrayBurstDetector(MEANCData);
        
    %Save the MEAData Struct
    [NAME] = NC_MATLABSaver(MEANCData(1).Settings.FileName);
    NAME = char(NAME);
    Name = strcat(Filepath, NAME);
    Name = char(Name);
    save(Name, 'MEANCData');
    
    [NAME] = NC_MATLABSaver(MEANCData(1).Settings.FileName);
    NAME = char(NAME);
    Name = strcat(Filepath,'RAW_', NAME);
    Name = char(Name);
    save(Name, 'MEAData','-v7.3');
    Filenum = Filenum+1;
    clearvars -except HDF5DataFile Filepath Filenum
    toc
end
toc