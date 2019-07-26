function FileData = MEA_BatchDataExtractor(FileName,LastBaseline, FirstFreq, SecondFreq)
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
%LAST MODIFIED: July 20, 2018
%v1.0


    load(FileName)
    
    NumOfFiles = length(MEADataS);
    %Get the power for each time point
    CA1Channels = MEADataS(1).Channels.CA1Labels;
    CA3Channels = MEADataS(1).Channels.CA3Labels;
    
    %check all files
    for FileNum = 1:LastBaseline
        %CA1
        %Check if there were channels in CA1 that are being analyzed
        if ~isempty(CA1Channels) == 1
            CA1Counter   = 1;
            %extract the data form each channel in CA1
            for ii = 1:length(CA1Channels)
                %Shift obtained the desired freuqencies
                RAWPowerCA1(CA1Counter,FirstFreq:SecondFreq) = MEADataS(FileNum).AnalyzedData(CA1Channels(ii)).AllPower(FirstFreq+1:SecondFreq+1);
                CA1Counter = CA1Counter + 1;
            end
            %save the data in the structure
        FileData.Time(FileNum).RawPowerCA1 = RAWPowerCA1;  
        
        end
        %CA3
        %Check if there were channels in CA3 that are being analyzed
        if ~isempty(CA3Channels) == 1
            CA3Counter   = 1; 
            %extract the data form each channel in CA3
            for ii = 1:length(CA3Channels)
                RAWPowerCA3(CA3Counter,FirstFreq:SecondFreq) = MEADataS(FileNum).AnalyzedData(CA3Channels(ii)).AllPower(FirstFreq+1:SecondFreq+1);
                CA3Counter = CA3Counter + 1;
            end
            %save the raw power into the sructure
        FileData.Time(FileNum).RawPowerCA3 = RAWPowerCA3;    
        end
    end
    %Extract the normalized data
    for FileNum = (LastBaseline+1):NumOfFiles
        if ~isempty(CA1Channels) == 1 
            %CA1
            CA1Counter   = 1;
            %Get all of the normalized data from the CA1 channels
            for ii = 1:length(CA1Channels)
                %Shift obtained the desired freuqencies
                PowerAllCA1(CA1Counter,FirstFreq:SecondFreq) = MEADataS(FileNum).AnalyzedData(CA1Channels(ii)).NormData(FirstFreq+1:SecondFreq+1);
                CA1LowGamma(CA1Counter) = MEADataS(FileNum).AnalyzedData(CA1Channels(ii)).Gamma.LowGamma;
                RAWPowerCA1(CA1Counter,FirstFreq:SecondFreq) = MEADataS(FileNum).AnalyzedData(CA1Channels(ii)).AllPower(FirstFreq+1:SecondFreq+1);
                CA1Counter = CA1Counter + 1;
            end
            %Calculate the mean of each frequency
            for ii = 1:length(PowerAllCA1)
                MeanCA1(ii) = mean(PowerAllCA1(:,ii));
                MeanCA1LG   = mean(CA1LowGamma);
            end
            %Save data in structure
            FileData.Time(FileNum).PowerAllCA1 = PowerAllCA1;
            FileData.Time(FileNum).RawPowerCA1 = RAWPowerCA1;
            FileData.Time(FileNum).CA1LG       = CA1LowGamma;
            FileData.Time(FileNum).CA1MeanLG   = MeanCA1LG;
            FileData.Time(FileNum).MeanCA1     = MeanCA1;
            FileData.ChannelID.CA1Channels     = CA1Channels;
        end
        %CA3
        if ~isempty(CA3Channels) == 1
            CA3Counter   = 1; 
            %extract the normalized data for all CA3 channels
            for ii = 1:length(CA3Channels)
                PowerAllCA3(CA3Counter,FirstFreq:SecondFreq) = MEADataS(FileNum).AnalyzedData(CA3Channels(ii)).NormData(FirstFreq+1:SecondFreq+1);
                CA3LowGamma(CA3Counter) = MEADataS(FileNum).AnalyzedData(CA3Channels(ii)).Gamma.LowGamma;
                RAWPowerCA3(CA3Counter,FirstFreq:SecondFreq) = MEADataS(FileNum).AnalyzedData(CA3Channels(ii)).AllPower(FirstFreq+1:SecondFreq+1);
                CA3Counter = CA3Counter + 1;
            end
            %Calculate the mean of each frequency
            for ii = 1:length(PowerAllCA3)
                MeanCA3(ii) = mean(PowerAllCA3(:,ii));
                MeanCA3LG   = mean(CA3LowGamma);
            end
            %save data in the structure
            FileData.Time(FileNum).PowerAllCA3 = PowerAllCA3;
            FileData.Time(FileNum).RawPowerCA3 = RAWPowerCA3;
            FileData.Time(FileNum).CA3LG       = CA3LowGamma;
            FileData.Time(FileNum).CA3MeanLG   = MeanCA3LG;
            FileData.Time(FileNum).MeanCA3     = MeanCA3;
            FileData.ChannelID.CA3Channels     = CA3Channels;
        end     
        %save the frequencies        
        FileData.Time(FileNum).Frequencies = FirstFreq:SecondFreq;       
    end 