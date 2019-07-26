function NC_BurstPlotter(CultureData,DIV,SavingInfo)
%This function will export the Burst Data for use in excel or other
%software that can generate figures.
%
%Inputs
%   CultureData - Data contaning all the cultures and the different days
%   they were in culture to be plotted
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN:  October 24, 2018
%LAST MODIFIED: October 24, 2018

%Define the three key varibales that will determine the size of all of the
%matrices to be created for the variables of interest
NumberOfCultures = length(CultureData);
NumberOfDays     = length(DIV);
NumChannels      = 60;

SpikeFreqThreshold  = 0.2*CultureData(1).Time;

%predefine the active channel variables
NumberOfBursts        = zeros(NumberOfDays, NumberOfCultures);
FrequencyOfBursts     = zeros(NumberOfDays, NumberOfCultures);
AverageBurstLength    = zeros(NumberOfDays, NumberOfCultures);
AverageSpikesPerBurst = zeros(NumberOfDays, NumberOfCultures);
AverageGlobSpikeFreq  = zeros(NumberOfDays, NumberOfCultures);
AverageInstSpikeFreq  = zeros(NumberOfDays, NumberOfCultures);

for CultureID = 1:NumberOfCultures
    for DayID = 1:NumberOfDays
        ActiveChanNum = 0;
        ActiveChan         = zeros(NumChannels,1);
        NumBurst           = zeros(NumChannels,1);
        FreqBurst          = zeros(NumChannels,1);
        AveBurstLength     = zeros(NumChannels,1);
        AveSpikesPerBurst  = zeros(NumChannels,1);
        AveInstSpikeFreq   = zeros(NumChannels,1);
        AvePercSpikesBurst = zeros(NumChannels,1);
        for Chan = 1:NumChannels
            %calculate the number of active channels
            if CultureData(CultureID).DIV(DayID).PeakData(Chan).SpikeNum >= SpikeFreqThreshold
                ActiveChan(Chan) = true;
                ActiveChanNum = ActiveChanNum+1;
                NumBurst(Chan)...
                    = CultureData(CultureID).DIV(DayID).BurstData(Chan).NumBurst;
                FreqBurst(Chan)...
                    = CultureData(CultureID).DIV(DayID).BurstData(Chan).FrequencyofBurst;
                AveBurstLength(Chan)...
                    = CultureData(CultureID).DIV(DayID).BurstData(Chan).AveBL;
                AveSpikesPerBurst(Chan)...
                    = CultureData(CultureID).DIV(DayID).BurstData(Chan).AveSPB;
                AveGlobSpikeFreq(Chan)...
                    = CultureData(CultureID).DIV(DayID).PeakData(Chan).GlobFreq;
                AveInstSpikeFreq(Chan)...
                    = CultureData(CultureID).DIV(DayID).PeakData(Chan).InstFreq;
                AvePercSpikesBurst(Chan)...
                    = (NumBurst(Chan) * AveSpikesPerBurst(Chan) * 100)/...
                    CultureData(CultureID).DIV(DayID).PeakData(Chan).SpikeNum;
                
            else
                ActiveChan(Chan) = false;
                NumBurst(Chan)           = NaN;
                FreqBurst(Chan)          = NaN;
                AveBurstLength(Chan)     = NaN;
                AveSpikesPerBurst(Chan)  = NaN;
                AveGlobSpikeFreq(Chan)   = NaN;
                AveInstSpikeFreq(Chan)   = NaN;
                AvePercSpikesBurst(Chan) = NaN;
            end
        end
        NumberOfBursts(DayID, CultureID)           = nanmean(NumBurst);
        FrequencyOfBursts(DayID,CultureID)         = nanmean(FreqBurst);
        AverageBurstLength(DayID, CultureID)       = nanmean(AveBurstLength);
        AverageSpikesPerBurst(DayID, CultureID)    = nanmean(AveSpikesPerBurst);
        AverageGlobSpikeFreq(DayID, CultureID)     = nanmean(AveGlobSpikeFreq);
        AverageInstSpikeFreq(DayID, CultureID)     = nanmean(AveInstSpikeFreq);
        NumofActiveChannels(DayID, CultureID)      = ActiveChanNum;
        AveragePercSpikesInBurst(DayID, CultureID) = nanmean(AvePercSpikesBurst);
        
        if NumofActiveChannels == 0
            ArrayBurstLength(DayID,CultureID) = NaN;
            ArrayFreqBurst(DayID,CultureID)   = NaN;
            ArraySPB(DayID,CultureID)         = NaN;
            ArrayNumBursts(DayID,CultureID)   = NaN;
        else
            CultureData(CultureID).DIV(DayID).Time(1).RecordingTime = CultureData(1).Time;
            ArrayBurstData = NC_ArrayBurstDetector(CultureData(CultureID).DIV(DayID));
            ArrayBurstLength(DayID,CultureID) = ArrayBurstData.AveBurstLength;
            ArrayFreqBurst(DayID,CultureID)   = ArrayBurstData.FrequencyofBurst;
            ArraySPB(DayID,CultureID)         = ArrayBurstData.AveSPB;
            ArrayNumBursts(DayID,CultureID)   = ArrayBurstData.NumBurst; 
        end
    end
    
end




%Save the data
ExcelFileName = strcat(SavingInfo.SaveDirectory, ...       %This is the file path
    'NeuronalCulture_', SavingInfo.Treatment,'_DATA','.xlsx'); %this is the actual file name

%Set the column and row headers

for ii = 1:NumberOfCultures
    RowHeader{ii,1} = sprintf('Culture %d',ii);
end
%determine the title of each column
ColumnHeader{1} = 'Replicate';
for ii = 1:NumberOfDays
    ColumnHeader{ii+1} = sprintf('DIV %d',DIV(ii));
end

%Prepare the output cell arrays witht he appropriate headings
OutNumberOfBursts        = [ColumnHeader; RowHeader num2cell(NumberOfBursts')];
OutFrequencyOfBursts     = [ColumnHeader; RowHeader num2cell(FrequencyOfBursts')];
OutAverageBurstLength    = [ColumnHeader; RowHeader num2cell(AverageBurstLength')];
OutAverageSpikesPerBurst = [ColumnHeader; RowHeader num2cell(AverageSpikesPerBurst')];
OutAverageGlobSpikeFreq  = [ColumnHeader; RowHeader num2cell(AverageGlobSpikeFreq')];
OutAverageInstSpikeFreq  = [ColumnHeader; RowHeader num2cell(AverageInstSpikeFreq')];
OutNumberofActiveChan    = [ColumnHeader; RowHeader num2cell(NumofActiveChannels')];
OutPercSpikesInBurst     = [ColumnHeader; RowHeader num2cell(AveragePercSpikesInBurst')];
OutArrayBurstLength  = [ColumnHeader; RowHeader num2cell(ArrayBurstLength')];
OutArrayFreqBurst    = [ColumnHeader; RowHeader num2cell(ArrayFreqBurst')];
OutArraySPB          = [ColumnHeader; RowHeader num2cell(ArraySPB')];
OutArrayNumBursts    = [ColumnHeader; RowHeader num2cell(ArrayNumBursts')];




%write to excel spreadsheet 
xlswrite(ExcelFileName, OutNumberOfBursts, 'NumberOfBursts')
xlswrite(ExcelFileName, OutFrequencyOfBursts, 'FrequencyOfBursts')
xlswrite(ExcelFileName, OutAverageBurstLength, 'AverageBurstLength')
xlswrite(ExcelFileName, OutAverageSpikesPerBurst, 'AverageSpikesPerBurst')
xlswrite(ExcelFileName, OutAverageGlobSpikeFreq, 'AverageGlobalSpikesFreq')
xlswrite(ExcelFileName, OutAverageInstSpikeFreq, 'AverageInstantaneousSpikeFreq')        
xlswrite(ExcelFileName, OutNumberofActiveChan, 'NumberofActiveChannels')
xlswrite(ExcelFileName, OutPercSpikesInBurst, 'PercentSpikesinBurst')   
xlswrite(ExcelFileName, OutArrayBurstLength, 'ArrayBurstLenght') 
xlswrite(ExcelFileName, OutArrayFreqBurst, 'FreqArrayBurst') 
xlswrite(ExcelFileName, OutArraySPB, 'ArraySPB') 
xlswrite(ExcelFileName, OutArrayNumBursts, 'NumArrayBursts') 


