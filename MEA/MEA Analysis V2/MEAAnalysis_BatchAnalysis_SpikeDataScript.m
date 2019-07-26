%This script is meant to be used to generate figures containing all
%generated from files used using one condition.
%MK801 injection
FileName{1} = 'H:\Documents\Analyzed_MATLAB\MEAData_2018-04-19_Slice1.mat';
FileName{2} = 'H:\Documents\Analyzed_MATLAB\MEAData_2018-04-19_Slice2.mat';
FileName{3} = 'H:\Documents\Analyzed_MATLAB\MEAData_2018-04-10_Slice1.mat';
%%
%PBS injection
FileName{1} = 'H:\Documents\Analyzed_MATLAB\MEAData_2018-04-11_Slice2.mat';
FileName{2} = 'H:\Documents\Analyzed_MATLAB\MEAData_2018-04-03_Slice1.mat';
FileName{3} = 'H:\Documents\Analyzed_MATLAB\MEAData_2018-05-03_Slice2.mat';
%%
%Jenkins WT
FileName{1} = 'H:\Documents\Analyzed_MATLAB\Jenkins\MEAData_2017-12-01_Slice2.mat';
FileName{2} = 'H:\Documents\Analyzed_MATLAB\Jenkins\MEAData_2017-12-06_Slice3.mat';
FileName{3} = 'H:\Documents\Analyzed_MATLAB\Jenkins\MEAData_2018-03-02_Slice2.mat';
FileName{5} = 'H:\Documents\Analyzed_MATLAB\Jenkins\MEAData_2018-01-17_Slice2.mat';
FileName{4} = 'H:\Documents\Analyzed_MATLAB\Jenkins\MEAData_2018-01-11_Slice2.mat';
%%
%Jenkins Mutant
FileName{1} = 'H:\Documents\Analyzed_MATLAB\Jenkins\MEAData_2017-12-08_Slice3.mat';
FileName{2} = 'H:\Documents\Analyzed_MATLAB\Jenkins\MEAData_2017-12-14_Slice1.mat';
FileName{3} = 'H:\Documents\Analyzed_MATLAB\Jenkins\MEAData_2017-12-14_Slice2.mat';
FileName{4} = 'H:\Documents\Analyzed_MATLAB\Jenkins\MEAData_2018-02-22_Slice1.mat';
FileName{5} = 'H:\Documents\Analyzed_MATLAB\Jenkins\MEAData_2018-02-22_Slice3.mat';


%%
%Bicuculline block
FileName{1} = 'H:\Documents\Analyzed_MATLAB\Bicuculline Recrodings\MEAData_2018-08-27_Slice2.mat';
% FileName{2} = 'H:\Documents\Analyzed_MATLAB\Bicuculline Recrodings\MEAData_2018-08-28_Slice1.mat';
FileName{2} = 'H:\Documents\Analyzed_MATLAB\Bicuculline Recrodings\MEAData_2018-09-12_Slice1.mat';
FileName{3} = 'H:\Documents\Analyzed_MATLAB\Bicuculline Recrodings\MEAData_2018-09-12_Slice2.mat';
FileName{4} = 'H:\Documents\Analyzed_MATLAB\Bicuculline Recrodings\MEAData_2018-09-14_Slice1.mat';
FileName{5} = 'H:\Documents\Analyzed_MATLAB\Bicuculline Recrodings\MEAData_2018-09-14_Slice2.mat';
%%
%MK801 Block

FileName{1} = 'H:\Documents\Analyzed_MATLAB\MK801 Bath application\MEAData_2018-09-18_Slice1.mat';
FileName{2} = 'H:\Documents\Analyzed_MATLAB\MK801 Bath application\MEAData_2018-09-18_Slice2.mat';

%%
%Caculate the average power spectrum for each brain region during the last
%kainate application for each slice
LastBaseline = 30;
LastKainate  = 60;
% LastBaseline = 1;
% LastKainate  = 31;
for File = 1:length(FileName)

    Files(File) = MEA_BatchDataExtractorSpiking(FileName{File});
         
 end
%%
%Firing activity plot section in CA1
NumofFiles = 68;
BinSize = 2;
PlotFields = {'CA1NumberofSpikes', 'CA1NumerofBursts'};
PlotUnitsY = {'Number of Spikes/min' ,'Number of Bursts/min'};
PlotTitles = {'Number of Spikes in CA1', 'Channel Bursting in CA1'};
for PlotFieldIx = 1:length(PlotFields)
    for File = 1:length(Files)
        figure
        hold on
        [indData] = M_SliceSpikePlotter(Files(File), PlotFields{PlotFieldIx});
        [TotalChannels, TotalTimePoints] = size(indData);
        for ChannelNum = 1:TotalChannels   
            plot((1:TotalTimePoints)*BinSize,indData(ChannelNum,:),'k','LineWidth',0.5)   
        end
        plot((1:TotalTimePoints)*BinSize,mean(indData),'k','LineWidth',6);
        AllFileData(File,:) = mean(indData);
        ylabel(PlotUnitsY{PlotFieldIx})
        xlabel('Time (min)')
        set(findall(gca, '-property', 'FontSize'),'FontSize',24);
        title(PlotTitles{PlotFieldIx},'Interpreter', 'none')
    end
    
    figure
    hold on
    for File = 1:length(Files)
        plot((1:TotalTimePoints)*BinSize,AllFileData(File,:)/BinSize,'k','LineWidth',0.5)   
    end
    plot((1:TotalTimePoints)*BinSize,mean(AllFileData)/BinSize,'k','LineWidth',6);
    ylabel(PlotUnitsY{PlotFieldIx})
    xlabel('Time (min)')
    set(findall(gca, '-property', 'FontSize'),'FontSize',24);
    P1 = ['Mean ' PlotTitles{PlotFieldIx}];
    title(P1,'Interpreter', 'none')
end
clearvars -except FileName Files AllFileData
%%
%Firing activity plot section in CA3
NumofFiles = 68;
BinSize = 2;
PlotFields = {'CA3NumberofSpikes', 'CA3NumerofBursts'};
PlotUnitsY = {'Number of Spikes/min' ,'Number of Bursts/min'};
PlotTitles = {'Number of Spikes in CA3', 'Channel Bursting in CA3'};
for PlotFieldIx = 1:length(PlotFields)
    for File = 1:length(Files)
        figure
        hold on
        [indData] = M_SliceSpikePlotter(Files(File), PlotFields{PlotFieldIx});
        [TotalChannels, TotalTimePoints] = size(indData);
        for ChannelNum = 1:TotalChannels   
            plot((1:TotalTimePoints)*BinSize,indData(ChannelNum,:),'k','LineWidth',0.5)   
        end
        plot((1:TotalTimePoints)*BinSize,mean(indData),'k','LineWidth',6);
        AllFileData(File,:) = mean(indData);
        ylabel(PlotUnitsY{PlotFieldIx})
        xlabel('Time (min)')
        set(findall(gca, '-property', 'FontSize'),'FontSize',24);
        title(PlotTitles{PlotFieldIx},'Interpreter', 'none')
    end
    
    figure
    hold on
    for File = 1:length(Files)
        plot((1:TotalTimePoints)*BinSize,AllFileData(File,:)/BinSize,'k','LineWidth',0.5)   
    end
    plot((1:TotalTimePoints)*BinSize,mean(AllFileData)/BinSize,'k','LineWidth',6);
    ylabel(PlotUnitsY{PlotFieldIx})
    xlabel('Time (min)')
    set(findall(gca, '-property', 'FontSize'),'FontSize',24);
    P1 = ['Mean ' PlotTitles{PlotFieldIx}];
    title(P1,'Interpreter', 'none')
end
clearvars -except FileName Files AllFileData
