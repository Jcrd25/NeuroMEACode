function NC_ExcelExport(CultureData,DIV,SavingInfo)
%This function will all of the data for use in excel or other
%software that can generate figures.
%
%Inputs
%   CultureData - Data contaning all the cultures and the different days
%   they were in culture to be plotted
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN:  October 24, 2018
%LAST MODIFIED: March 18, 2019

%Define the three key varibales that will determine the size of all of the
%matrices to be created for the variables of interest
NumberOfCultures = length(CultureData);
NumberOfDays     = length(DIV);
NumChannels      = 60;

SpikeFreqThreshold  = 0.2*CultureData(1).Time;

AverageGlobSpikeFreq  = zeros(NumberOfDays, NumberOfCultures);
AverageInstSpikeFreq  = zeros(NumberOfDays, NumberOfCultures);

for CultureID = 1:NumberOfCultures
    for DayID = 1:NumberOfDays
        for Chan = 1:NumChannels
            if CultureData(CultureID).DIV(DayID).PeakData(Chan).SpikeNum > SpikeFreqThreshold
                CultureData(CultureID).DIV(
            end
        end
        for ChanNum = 1:NumChannels
            Channel(ChanNum).AverageGlobSpikeFreq(DayID)  = AveGlobSpikeFreq(ChanNum);
            Channel(ChanNum).AverageInstSpikeFreq(DayID)  = AveInstSpikeFreq(ChanNum);
        end 
    end
    figure
    hold on
    for ChanNum = 1:NumChannels
        x = DIV;
        y = Channel(ChanNum).AverageGlobSpikeFreq;
        yall(ChanNum,:) = y;
        plot(x,y,'k','LineWidth', 0.5)
    end
    
    plot(x,mean(yall),'k','LineWidth', 6)
    xlabel('Days in culture')
    ylabel('Spike Frequency (Hz)')
    if strcmp(SavingInfo.Treatment, 'Vehicle')
        P1 = sprintf('Spike Frequency in %s Treated Culture %s', SavingInfo.Treatment, SavingInfo.VehicleName{CultureID});
    elseif strcmp(SavingInfo.Treatment, 'MK801')
       P1 = sprintf('Spike Frequency in %s Treated Culture %s', SavingInfo.Treatment, SavingInfo.TreatedName{CultureID}); 
    end
    title(P1)
    hold off
    FigureName = strcat(SavingInfo.SaveDirectory, ...       %This is the file path
    'NeuronalCulture_',SavingInfo.VehicleName{CultureID},'Treatment_', SavingInfo.Treatment,'_DATA'); %this is the actual file name
    saveas(gcf,FigureName,'tiffn')
    yaverage(CultureID,:) = mean(yall);
end

figure
hold on
for ii = 1:NumberOfCultures
    plot(x,yaverage(ii,:),'k','LineWidth', 0.5)
end
plot(x,mean(yaverage),'k','LineWidth', 6)
xlabel('Days in culture')
ylabel('Spike Frequency (Hz)')
P1 = sprintf('Average Spike Frequency in %s Treated Cultures', SavingInfo.Treatment);
title(P1)
hold off
FigureName = strcat(SavingInfo.SaveDirectory, ...       %This is the file path
'NeuronalCultureAverage_','Treatment_', SavingInfo.Treatment,'_DATA'); %this is the actual file name
saveas(gcf,FigureName,'tiffn')
hold off




% %Save the data
% ExcelFileName = strcat(SavingInfo.SaveDirectory, ...       %This is the file path
%     'NeuronalCulture_', SavingInfo.Treatment,'_DATA','.xlsx'); %this is the actual file name
% 
% %Set the column and row headers
% 
% for ii = 1:NumberOfCultures
%     RowHeader{ii,1} = sprintf('Culture %d',ii);
% end
% %determine the title of each column
% ColumnHeader{1} = 'Replicate';
% for ii = 1:NumberOfDays
%     ColumnHeader{ii+1} = sprintf('DIV %d',DIV(ii));
% end
% 
% %Prepare the output cell arrays witht he appropriate headings
% OutAverageGlobSpikeFreq = [ColumnHeader; RowHeader num2cell(AverageGlobSpikeFreq')];
% OutAverageInstSpikeFreq = [ColumnHeader; RowHeader num2cell(AverageInstSpikeFreq')];
% 
% %write to excel spreadsheet 
% xlswrite(ExcelFileName, OutAverageGlobSpikeFreq, 'AverageGlobalSpikesFreq')
% xlswrite(ExcelFileName, OutAverageInstSpikeFreq, 'AverageInstantaneousSpikeFreq')        
            