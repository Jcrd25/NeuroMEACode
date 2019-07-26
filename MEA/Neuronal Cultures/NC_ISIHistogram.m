function [ISIExport,CumProbabilityExport] = NC_ISIHistogram(CultureData,DIV,SavingInfo, MEANumber)
%This function will generate histograms of the culutures and have the ISI 
%in an exportable manner.
%
%Inputs
%   CultureData - Data contaning all the cultures and the different days
%   they were in culture to be plotted
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN:  May 6, 2019
%LAST MODIFIED: May 6, 2019

%Define the three key varibales that will determine the size of all of the
%matrices to be created for the variables of interest
NumberOfCultures = length(CultureData);
NumberOfDays     = length(DIV);
NumChannels      = 60;

SpikeFreqThreshold  = 0.2*CultureData(1).Time;

% edges = 0:10:350;
figure(1)

for DayID = 1:NumberOfDays 
    xall= [];
    for CultureID = 1:NumberOfCultures
        x =[];
        for Chan = 1:NumChannels
            if CultureData(CultureID).DIV(DayID).PeakData(Chan).SpikeNum >= SpikeFreqThreshold
            x = [x CultureData(CultureID).DIV(DayID).PeakData(Chan).ISI'];
            end
        end
        if ~isempty(x) == 1
%             figure
%             histogram(x)
%             xlabel('ISI (ms')
%             ylabel('counts')
%             P1 = sprintf('ISI of MEAID %d DIV %d',CultureID,DIV(DayID));
%             title(P1)
%             ISIExport(CultureID).Day(DayID).ExportCumulValues = ecdf(x);
            ISIExport(CultureID).Day(DayID).ISIs = x;
            xall = [xall x];
            figure
            cdfplot(x)
            P1 = sprintf('Cum Prob of %s DIV %d',MEANumber{CultureID},DIV(DayID));
            title(P1)
            FigureName = strcat(SavingInfo.SaveDirectory, ...       %This is the file path
            'NeuronalCulture_',MEANumber{CultureID},'Treatment_', SavingInfo.Treatment,'_Figure'); %this is the actual file name
            saveas(gcf,FigureName,'tiffn')
        end
    end
%     figure
%     histogram(xall)
%     xlabel('ISI (ms)')
%     ylabel('counts')
%     P1 = sprintf('ISI of DIV %d',DIV(DayID));
%     title(P1)
    figure(1)
    [Exporty,Exportx] = ecdf(xall);
    hold on
    cdfplot(xall)
    xlabel('ISI (ms)')
    P1 = sprintf('Comulative Probability of %s treated', SavingInfo.Treatment);
    title(P1)
    CumProbabilityExport(DayID).CumulValues = [Exporty,Exportx];
end

