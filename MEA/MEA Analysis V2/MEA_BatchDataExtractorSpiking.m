function FileData = MEA_BatchDataExtractorSpiking(FileName)
%This function is meant to extract all the File Data and give out a very
%summarized data set
%Input
%   FileName - The name of the file to be loaded
%Output
%   FileData - Summarized version of all the time recordings of a single
%   slice
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN: April 13, 2018
%LAST MODIFIED: November 28, 2018
%v1.0


    load(FileName)
    
    NumOfFiles = length(MEADataS);
    %Get the power for each time point
    CA1Channels = MEADataS(1).Channels.CA1Labels;
    CA3Channels = MEADataS(1).Channels.CA3Labels;
    
    %check all files
    for FileNum = 1:NumOfFiles
        %CA1
        %Check if there were channels in CA1 that are being analyzed
        if ~isempty(CA1Channels) == 1
            CA1Counter   = 1;
            %extract the data form each channel in CA1
            for ii = 1:length(CA1Channels)
                %Shift obtained the desired freuqencies
               CA1NumberofSpikes(CA1Counter) = MEADataS(FileNum).PeakData(CA1Channels(ii)).SpikeNum;
               CA1InstFreq(CA1Counter) = MEADataS(FileNum).PeakData(CA1Channels(ii)).InstFreq;
               CA1NumerofBursts(CA1Counter) = MEADataS(FileNum).BurstData(CA1Channels(ii)).NumBurst;
               CA1Counter = CA1Counter + 1;
            end
            %save the data in the structure
        FileData.Time(FileNum).CA1NumberofSpikes  = CA1NumberofSpikes;
        FileData.Time(FileNum).CA1InstFreq        = CA1InstFreq;
        FileData.Time(FileNum).CA1NumerofBursts   = CA1NumerofBursts;
        
        end
        %CA3
        %Check if there were channels in CA3 that are being analyzed
        if ~isempty(CA3Channels) == 1
            CA3Counter   = 1; 
            %extract the data form each channel in CA3
            for ii = 1:length(CA3Channels)
               CA3NumberofSpikes(CA3Counter) = MEADataS(FileNum).PeakData(CA3Channels(ii)).SpikeNum;
               CA3InstFreq(CA3Counter) = MEADataS(FileNum).PeakData(CA3Channels(ii)).InstFreq;
               CA3NumerofBursts(CA3Counter) = MEADataS(FileNum).BurstData(CA3Channels(ii)).NumBurst;
               CA3Counter = CA3Counter + 1;
            end
            %save the raw power into the sructure
        FileData.Time(FileNum).CA3NumberofSpikes  = CA3NumberofSpikes;
        FileData.Time(FileNum).CA3InstFreq        = CA3InstFreq;
        FileData.Time(FileNum).CA3NumerofBursts   = CA3NumerofBursts;         
        end 
    end