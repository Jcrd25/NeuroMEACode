%Neuronal Culture Master Script for spike analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Edit these variables
%Identify the desired files to be analyzed.
SavedDataDirectory = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\NeuronalCulture\Shigeki Collab\Data files';
SaveToDirectory = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\NeuronalCulture\Shigeki Collab\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MatlabFiles = BatchLoader(SavedDataDirectory);   
CultureNames = fieldnames(MatlabFiles);

%determine the number of DIV and the number of MEAs used
for CultureInd = 1:length(CultureNames)
    for FileInd = 1:length(MatlabFiles.(CultureNames{CultureInd}))
        SplitFileName = strsplit(MatlabFiles.(CultureNames{CultureInd})(FileInd).name,'_');
        if strcmp(SplitFileName{1}, 'RAW') == 1
            MatlabFiles.(CultureNames{CultureInd})(FileInd).DataType = 'RAW';
        elseif strcmp(SplitFileName{1}, 'NC') == 1
            MatlabFiles.(CultureNames{CultureInd})(FileInd).DataType = 'NC';
        end
        for ii = 1:length(SplitFileName)
            if ~isempty(strfind(SplitFileName{ii}, 'DIV'))
                MatlabFiles.(CultureNames{CultureInd})(FileInd).DIV = SplitFileName{ii};
                MatlabFiles.(CultureNames{CultureInd})(FileInd).Treatment = SplitFileName{ii+1};
            elseif ~isempty(strfind(SplitFileName{ii}, 'MEA-'))
                MatlabFiles.(CultureNames{CultureInd})(FileInd).MEA = SplitFileName{ii};
            end
        end
    end
end
clear ii CultureInd FileInd

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Data Viewer
% Save2 = [SaveToDirectory 'Figures\'];
% %This part will generate graphs for all the cultures and 
% for CultureInd = 1:length(CultureNames)
%    for NCFileInd = 1:(length(MatlabFiles.(CultureNames{CultureInd})))
%        if MatlabFiles.(CultureNames{CultureInd})(NCFileInd).DataType == 'NC'
%         DATAVIEWERNC(MatlabFiles.(CultureNames{CultureInd})(NCFileInd),...
%             Save2, CultureNames{CultureInd})
%        end
%    end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Excel Exporter
%This part will generate excel files with the data
Save2 = [SaveToDirectory 'ExcelData\'];
Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';NumChan = 60;DIVs = {'DIV10','DIV14'};
ExportParams = {'NumActiveChannels';'NumSpikes';'GlobalFreq';'InstFreq';'NumBurst';'FreqBurst';...
    'AveBurstLength';'AveSpikesPerBurst'};
ExpermCond = {'scramble', 'EGFP','RaiCV','Rai141'};
for Conditions = 1:length(ExpermCond)
    for DIV = 1:length(DIVs)
        Counter = 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Get the data and organize it in a script
        for CultureInd = 1:length(CultureNames)
            for NCFileInd = 1:(length(MatlabFiles.(CultureNames{CultureInd})))
                if strcmp(MatlabFiles.(CultureNames{CultureInd})(NCFileInd).DataType, 'NC') == 1
                    if strcmp(ExpermCond{Conditions}, MatlabFiles.(CultureNames{CultureInd})(NCFileInd).Treatment) == 1
                        if strcmp( DIVs{DIV}, MatlabFiles.(CultureNames{CultureInd})(NCFileInd).DIV) == 1
                            Data = DataExcel(MatlabFiles.(CultureNames{CultureInd})(NCFileInd));
                            TempData.(ExpermCond{Conditions}).(DIVs{DIV}).MEA(Counter).Data = Data;
                            TempData.(ExpermCond{Conditions}).(DIVs{DIV}).MEA(Counter).ID =...
                                MatlabFiles.(CultureNames{CultureInd})(NCFileInd).MEA;
                            Counter = Counter+1;
                        end
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Export Data to excel
        for ParamInt = 1:length(ExportParams) 
            if strcmp(ExportParams{ParamInt}, 'NumActiveChannels') == 1               
                FirstColumn{1} = '';
                FirstColumn{2} = 'Number of Active Channels';
                FirstColumn = FirstColumn';
                Array = FirstColumn;
                for MEANum = 1:length(TempData.(ExpermCond{Conditions}).(DIVs{DIV}).MEA)
                    MEAColumn{1} = TempData.(ExpermCond{Conditions}).(DIVs{DIV}).MEA(MEANum).ID;
                    ColumnData = TempData.(ExpermCond{Conditions}).(DIVs{DIV}).MEA(MEANum).Data.(ExportParams{ParamInt});
                    Column2Add = [MEAColumn num2cell(ColumnData)'];
                    Array      = [Array Column2Add'];
                end
                ExcelFileName = strcat(Save2, ...       %This is the file path
                     'NeuronalCulture_', DIVs{DIV},'_',ExpermCond{Conditions},'_DATA','.xlsx');
                xlswrite(ExcelFileName, Array, ExportParams{ParamInt})
            else
                for ii = 1:NumChan
                    FirstColumn{ii+1} = sprintf('Chan ID %d',ii);
                end
                FirstColumn{1}         = 'Channel';
                FirstColumn{NumChan+2} = 'Mean';
                FirstColumn = FirstColumn';
                Array = FirstColumn;
                for MEANum = 1:length(TempData.(ExpermCond{Conditions}).(DIVs{DIV}).MEA)
                    MEAColumn{1} = TempData.(ExpermCond{Conditions}).(DIVs{DIV}).MEA(MEANum).ID;
                    ColumnData = TempData.(ExpermCond{Conditions}).(DIVs{DIV}).MEA(MEANum).Data.(ExportParams{ParamInt});
                    Column2Add = [MEAColumn num2cell(ColumnData)'];
                    Array      = [Array Column2Add'];
                end
                ExcelFileName = strcat(Save2, ...       %This is the file path
                     'NeuronalCulture_', DIVs{DIV},'_',ExpermCond{Conditions},'_DATA','.xlsx');
                xlswrite(ExcelFileName, Array, ExportParams{ParamInt})
            end
            clear FirstColumn Array MEAColumn ColumnData Column2Add ExcelFileName
        end
    end
end
clear Conditions DIV CultureInd NCFileInd ParamInt MEANum
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MEAExporting = DataExcel(CultureFileNC)
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
%LAST MODIFIED: May 28, 2019
load(fullfile(CultureFileNC.folder,  CultureFileNC.name))
%Define the three key varibales that will determine the size of all of the
%matrices to be created for the variables of interest
NumChannels      = 60;
SpikeFreqThreshold  = 0.2*MEANCData.Time(1).RecordingTime;

NumSpikes         = zeros(NumChannels,1);
GlobalFreq        = zeros(NumChannels,1);
InstFreq          = zeros(NumChannels,1);
NumBurst          = zeros(NumChannels,1);
FreqBurst         = zeros(NumChannels,1);
AveBurstLength    = zeros(NumChannels,1);
AveSpikesPerBurst = zeros(NumChannels,1);
NumActiveChannels = 0;
for Chan = 1:NumChannels
    
    if MEANCData.PeakData(Chan).SpikeNum >= SpikeFreqThreshold
        NumSpikes(Chan)...
            = MEANCData.PeakData(Chan).SpikeNum;
        GlobalFreq(Chan)...
            = MEANCData.PeakData(Chan).GlobFreq;
        InstFreq(Chan) ...
            = MEANCData.PeakData(Chan).InstFreq; 
        NumBurst(Chan)...
            = MEANCData.BurstData(Chan).NumBurst;
        FreqBurst(Chan)...
            = MEANCData.BurstData(Chan).FrequencyofBurst;
        AveBurstLength(Chan)...
            = MEANCData.BurstData(Chan).AveBL;
        AveSpikesPerBurst(Chan)...
            = MEANCData.BurstData(Chan).AveSPB;
        NumActiveChannels = NumActiveChannels+1;
    else
        NumSpikes(Chan)        = NaN;
        GlobalFreq(Chan)       = NaN;
        InstFreq(Chan)         = NaN; 
        NumBurst(Chan)         = NaN;
        FreqBurst(Chan)        = NaN;
        AveBurstLength(Chan)   = NaN;
        AveSpikesPerBurst(Chan)= NaN;
    end
end
NumSpikes(NumChannels+1)         = nanmean(NumSpikes);
GlobalFreq(NumChannels+1)        = nanmean(GlobalFreq);
InstFreq(NumChannels+1)          = nanmean(InstFreq);
NumBurst(NumChannels+1)          = nanmean(NumBurst);
FreqBurst(NumChannels+1)         = nanmean(FreqBurst);
AveBurstLength(NumChannels+1)     = nanmean(AveBurstLength);
AveSpikesPerBurst(NumChannels+1) = nanmean(AveSpikesPerBurst);
ArrayData                 = NC_ArrayBurstDetector(MEANCData);
%Export Variables
MEAExporting.NumSpikes         = NumSpikes;
MEAExporting.GlobalFreq        = GlobalFreq;
MEAExporting.InstFreq          = InstFreq;
MEAExporting.NumBurst          = NumBurst;
MEAExporting.FreqBurst         = FreqBurst;
MEAExporting.AveBurstLength    = AveBurstLength;
MEAExporting.AveSpikesPerBurst = AveSpikesPerBurst;
MEAExporting.ArrayData         = ArrayData;
MEAExporting.NumActiveChannels =  NumActiveChannels;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DATAVIEWERNC(CultureFileNC, SaveDirectory, CultureID)
%This function will generate figures of the raw traces, spikes and
%rasterplots of the data.
%
%Inputs
%   CultureFileNC - Data contaning Neuronal Culture Data
%   SaveDirectory - Directory in which data will be saved in
%   CultureID     - Information on how to identify the culture being
%                   analysed
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN:  May 15, 2019
%LAST MODIFIED: May 28, 2019
Culture = strcat('MEA',CultureFileNC.MEA(5:end));
DIV     = CultureFileNC.DIV;
load(fullfile(CultureFileNC.folder,  CultureFileNC.name))
RAWNAME = strcat('RAW_',CultureFileNC.name);
load(fullfile(CultureFileNC.folder, RAWNAME));
SavingInfo.SaveDirectory = fullfile(SaveDirectory, CultureID, Culture(1:8),DIV,'/');
mkdir(SavingInfo.SaveDirectory);
 %close all other figures to ensure the save correctly
 close all
%  tic
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Plot the Raw Trace
for ChanNum = 1:60
    %obtain the relevant data to generate the plots
    figure
    hold on
    Fs = MEAData.RawData.SamplingRate;
    %obtain the time and voltage trace
    x = 1:length(MEAData.RawData.VoltageTrace);
    x = x / (Fs);
    VoltageTrace = MEAData.RawData.VoltageTrace(:,ChanNum);
    
    %Parameters of the filter
    Fstop = 75;
    Fpass = 100;
    Astop = 65;
    Apass = 0.5;
    d = designfilt('highpassiir','StopbandFrequency',Fstop ,...
      'PassbandFrequency',Fpass,'StopbandAttenuation',Astop, ...
      'PassbandRipple',Apass,'SampleRate',Fs,'DesignMethod','butter');

    %Filter the Data
    VoltageTrace = filtfilt(d, VoltageTrace);
    
    %calculate the root mean square(used as the threshold for spike
    %detection)
    RMS_Signal = rms(VoltageTrace);
    Mean       = mean(VoltageTrace);
    %Set threshold to be 5 times the Root Mean Square
    SetThreshold = 5*RMS_Signal;
    %generate a horizontal line
    Line = ones(1,length(VoltageTrace));
    %plot the data
    plot(x,VoltageTrace,'k')
    %plot the thresholds
%     plot(x,Line*SetThreshold,'r')
%     plot(x,Line*SetThreshold*-1,'r')
    %plot the mean
%     plot(x,Line*Mean, 'b')
    %Figure and axis titles
    T1 = sprintf('%s Trace of Chan %d', Culture, ChanNum);
    title(T1)
    xlabel('Time (s)')
    ylabel('Voltage (µV)')
    ylim([-60 75])
    hold off
end
%save the data
for ii = 1:60
numb = ii;
figure(ii)
NAME = sprintf('%s Channel ID %d', Culture, numb);
FigureName = strcat(SavingInfo.SaveDirectory,'RawData_', NAME);
saveas(gcf,FigureName,'tiffn')
end
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Spike plotting
for ChanNum = 1:60
    if MEANCData.PeakData(ChanNum).SpikeNum >= 1
        CTrace = [];
        figure
        hold on
        %obtain the trace of each spike
        for Spike =  1:MEANCData.PeakData(ChanNum).SpikeNum
  
                Peak = MEANCData.PeakData(ChanNum).SpikePeakLoc(Spike);
                leftedge = Peak - 0.001*Fs;
                rightedge = Peak + 0.003*Fs;
                %make sure that the trace before and after the spike can be
                %plotted
                if leftedge >= 1 
                    if rightedge <= length(MEAData.RawData.VoltageTrace(:,ChanNum))
                        Trace = MEAData.RawData.VoltageTrace(leftedge:rightedge,ChanNum);
                        time = rightedge-leftedge + 1;
                        plot((1:time)/Fs, Trace,'k')
                        CTrace(Spike,:) = Trace; 
                    else
                        P1 = sprintf('Not Plotted last spike of Chan %d',ChanNum);
                        disp(P1)
                    end
                else
                    P1 = sprintf('Not Plotted first spike of Chan %d',ChanNum);
                    disp(P1)
                end
            
        end
        %Plot the average spike
        plot((1:time)/Fs,mean(CTrace),'r','LineWidth', 6)
        %Figure and axis titles
        T1 = sprintf('Spike Traces of %s Chan %d', Culture, ChanNum);
        title(T1)
        xlabel('Time (s)')
        ylabel('Voltage (µV)')
    end
end
%Save the figures
for ii = 1:60
numb = ii;
figure(ii)
NAME = sprintf('%s Channel ID %d', Culture, numb);
FigureName = strcat(SavingInfo.SaveDirectory,'SpikeData_', NAME);
saveas(gcf,FigureName,'tiffn')
end
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot the raster plot
NC_RasterPlot(MEANCData)
FigureName = strcat(SavingInfo.SaveDirectory,'RasterPlot_', Culture,'_',DIV);
saveas(gcf,FigureName,'tiffn')
close all
disp('Finished');

% toc
end