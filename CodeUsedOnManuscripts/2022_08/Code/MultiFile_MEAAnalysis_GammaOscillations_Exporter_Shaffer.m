%%
%Get the data ready to be exported into excel or to be copy and  pasted to
%prism
%Select the data type to be analyzed
BandPowerData = true; PeriodogramData = false;

Data2Choose = 20;
%20 = ShafferCut

%50 = Assay Dev data
%-------------------------------------------------------------------------%
%Add the directory or directories of the data here
if Data2Choose <= 30 && Data2Choose >= 20
MainDirectory = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\ShafferCut';
elseif Data2Choose <= 60 && Data2Choose >= 50
MainDirectory = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\AssayDev';
end
%--------------------`-----------------------------------------------------%
%This is a check to make sure setting are set up correctly
if BandPowerData == 1 && PeriodogramData == 1
    disp('Selected both BandPowerData and PeriodogramData')
    disp('Please select one')
    return
elseif BandPowerData == 0 && PeriodogramData == 0
    disp('Did not select a data type')
    disp('Please select one')
    return 
%-------------------------------------------------------------------------%
elseif BandPowerData == 1
%set up the directory
MainDirectory = [MainDirectory '\BandPower'];
MainDirectory = [MainDirectory '\RAW_Power\'];%RAW_Power or Relative_Power
DataType = 'RawPower_'; %RawPower_ or RelativePower_
%----------------------------------------------------%
elseif PeriodogramData == 1
%Select the files to be loaded based on their conditions to be used for
%Periodogram Summary analysis
%set up the directory];
MainDirectory = [MainDirectory '\PowerSummaryData\'];%PowerSummaryData or AllKainatePowerSummaryData
DataType = 'PowSumD_'; %PowSumD_ or PowAllD_
%----------------------------------------------------%
end

SliceNames = {...
    '2022-01-16_Slice1','2022-01-17_Slice1'...%ShafferCut
    '2022-03-03_Slice1', '2022-03-03_Slice2',...
    '2021-11-14_Slice4',... %filler
    '2018-08-27_Slice2','2018-08-28_Slice1', ... %Bicuculline Slices
    '2018-09-12_Slice1', '2018-09-12_Slice2', '2018-09-14_Slice1',...
    '2018-09-14_Slice2', '2020-02-27_Slice1', '2020-02-27_Slice2',...
    '2018-09-18_Slice1','2018-09-18_Slice2',... %MK-801 Slices
    '2020-08-19_Slice2';... %DMSO
    };

%load the channel labels
load('G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\Slice_ChannelLabelInfo.mat');
%find the file and channel information
[FileNames,SlicesUsed,Channels] = DataSetChooser(Data2Choose,MainDirectory,DataType,SliceNames,SliceChannelInfo.Label);
TotalCA1 = 0; TotalCA3 = 0;
%load the file data
for FileIndex = 1:length(FileNames)
    AllData(FileIndex) = load(FileNames{FileIndex});
    TotalCA3 = TotalCA3+length(Channels(FileIndex).CA3Channels);
    TotalCA1 = TotalCA1+length(Channels(FileIndex).CA1Channels);
    Channels(FileIndex).SlicesUsed = SlicesUsed{FileIndex};
end

%Export area
%---------------------------------%
% Variables to specify
Time1 = 8/60;Time2 = 60+60;%in minutes
deltaT = 0.5;%in seconds
PlotFigsCA1 = false;PlotFigsCA3 = false;
%---------------------------------%
%determine the time inforation to be used to extract data
T1 = Time1*60; T2 = Time2*60;
DesiredTimes = T1:deltaT:T2;
[~,NumTimePoints] = size(DesiredTimes);

%preallocate
if strcmp(DataType, 'RawPower_') || strcmp(DataType, 'RelativePower_')
ExportCA1Gamma  = zeros(TotalCA1,NumTimePoints-1);
ExportCA3Gamma  = zeros(TotalCA3,NumTimePoints-1);
ExportCA1HGamma = zeros(TotalCA1,NumTimePoints-1);
ExportCA3HGamma = zeros(TotalCA3,NumTimePoints-1);
ExportCA1TGamma = zeros(TotalCA1,NumTimePoints-1);
ExportCA3TGamma = zeros(TotalCA3,NumTimePoints-1);
ExportCA1Times  = zeros(TotalCA1,NumTimePoints-1);
ExportCA3Times  = zeros(TotalCA3,NumTimePoints-1);
elseif strcmp(DataType, 'PowSumD_MEAData_')
ExportPower = struct; 
end
CA3ElectrodeCount= 0;CA1ElectrodeCount =0;
for FileIndex = 1:length(FileNames)
    if strcmp(DataType, 'RawPower_')
        %get the data for CA3
        if ~isempty(Channels(FileIndex).CA3Channels) == 1
            for ChannelI = 1:length(Channels(FileIndex).CA3Channels)
                ChannelIndex     = Channels(FileIndex).CA3Indexes(ChannelI);
                ChanGammaPower   = AllData(FileIndex).MEA_PowerData(ChannelIndex).PowerData.Raw.LowGamma; 
                ChanHGammaPower  = AllData(FileIndex).MEA_PowerData(ChannelIndex).PowerData.Raw.HighGamma; 
                ChanBetaPower    = AllData(FileIndex).MEA_PowerData(ChannelIndex).PowerData.Raw.Beta; 
                ChanThetaPower   = AllData(FileIndex).MEA_PowerData(ChannelIndex).PowerData.Raw.Theta; 
                ChanDeltaPower   = AllData(FileIndex).MEA_PowerData(ChannelIndex).PowerData.Raw.Delta; 
                ChanTotalPower   = AllData(FileIndex).MEA_PowerData(ChannelIndex).PowerData.Raw.Total;  
                Time             = AllData(FileIndex).MEA_PowerData(ChannelIndex).Time;
                if PlotFigsCA3 == true
                    figure
                    hold on
                    plot((Time./60),ChanGammaPower,'b','LineWidth', 1)
%                     plot((Time./60),ChanBetaPower,'r','LineWidth', 1)
%                     plot((Time./60),ChanThetaPower,'c','LineWidth', 1)
%                     plot((Time./60),ChanTotalPower,'k','LineWidth', 2)
                    %make the figure and axes titles
                    set(findall(gca, '-property', 'FontSize'),'FontSize',24);
                    xlabel('Time (min)');ylabel('Power in Gamma Band (�V^2)');
                    Title1 = sprintf('Power in CA3 Gamma \nfor electrode %d \nfrom %s',...
                        Channels(FileIndex).CA3Channels(ChannelI),SlicesUsed{FileIndex});
                    Title1 =splitlines(Title1);
                    title(Title1,'Interpreter', 'none')

                end
                CA3ElectrodeCount = CA3ElectrodeCount+1;
                DesiredIndexes = find(Time > T1 & Time < T2);
                ExportCA3Gamma(CA3ElectrodeCount,:) = ChanGammaPower(DesiredIndexes);
                ExportCA3HGamma(CA3ElectrodeCount,:) = ChanHGammaPower(DesiredIndexes);
                ExportCA3TGamma(CA3ElectrodeCount,:) = ChanGammaPower(DesiredIndexes)+ChanHGammaPower(DesiredIndexes);
                ExportCA3Times(CA3ElectrodeCount,:)  = Time(DesiredIndexes);
            end
        end
        %get the data for CA1
        if ~isempty(Channels(FileIndex).CA1Channels) == 1
            for ChannelI = 1:length(Channels(FileIndex).CA1Channels)
                ChannelIndex     = Channels(FileIndex).CA1Indexes(ChannelI);
                ChanGammaPower   = AllData(FileIndex).MEA_PowerData(ChannelIndex).PowerData.Raw.LowGamma; 
                ChanHGammaPower  = AllData(FileIndex).MEA_PowerData(ChannelIndex).PowerData.Raw.HighGamma; 
                ChanBetaPower    = AllData(FileIndex).MEA_PowerData(ChannelIndex).PowerData.Raw.Beta; 
                ChanThetaPower   = AllData(FileIndex).MEA_PowerData(ChannelIndex).PowerData.Raw.Theta; 
                ChanDeltaPower   = AllData(FileIndex).MEA_PowerData(ChannelIndex).PowerData.Raw.Delta; 
                ChanTotalPower   = AllData(FileIndex).MEA_PowerData(ChannelIndex).PowerData.Raw.Total;  
                Time             = AllData(FileIndex).MEA_PowerData(ChannelIndex).Time;
                if PlotFigsCA1 == true
                    figure
                    hold on
                    plot((Time./60),ChanGammaPower,'b','LineWidth', 1)
%                     plot((Time./60),ChanTotalPower,'k','LineWidth', 2)
                    %make the figure and axes titles
                    set(findall(gca, '-property', 'FontSize'),'FontSize',24);
                    xlabel('Time (min)');ylabel('Power in Gamma Band (�V^2)');
                    Title1 = sprintf('Power in CA1 Gamma \nfor electrode %d \nfrom %s',...
                        Channels(FileIndex).CA1Channels(ChannelI),SlicesUsed{FileIndex});
                    Title1 =splitlines(Title1);
                    title(Title1,'Interpreter', 'none')

                end
                CA1ElectrodeCount = CA1ElectrodeCount+1;
                DesiredIndexes = find(Time > T1 & Time < T2);
                ExportCA1Gamma(CA1ElectrodeCount,:) = ChanGammaPower(DesiredIndexes);
                ExportCA1HGamma(CA1ElectrodeCount,:) = ChanHGammaPower(DesiredIndexes);
                ExportCA1TGamma(CA1ElectrodeCount,:) = ChanGammaPower(DesiredIndexes)+ChanHGammaPower(DesiredIndexes);
                ExportCA1Times(CA1ElectrodeCount,:)  = Time(DesiredIndexes);
            end
        end
    elseif strcmp(DataType,'RelativePower_')
        %get the data for CA3
        if ~isempty(Channels(FileIndex).CA3Channels) == 1
            for ChannelI = 1:length(Channels(FileIndex).CA3Channels)
                ChannelIndex    = Channels(FileIndex).CA3Indexes(ChannelI);
                ChanGammaPower  = AllData(FileIndex).MEA_RelativeData(ChannelIndex).RelativePower.LowGamma; 
                ChanBetaPower   = AllData(FileIndex).MEA_RelativeData(ChannelIndex).RelativePower.Beta; 
                ChanThetaPower  = AllData(FileIndex).MEA_RelativeData(ChannelIndex).RelativePower.Theta; 
                ChanDeltaPower  = AllData(FileIndex).MEA_RelativeData(ChannelIndex).RelativePower.Delta; 
                ChanTotalPower  = AllData(FileIndex).MEA_RelativeData(ChannelIndex).RelativePower.Total;  
                Time            = AllData(FileIndex).MEA_RelativeData(ChannelIndex).Time;
                if PlotFigsCA3 == true
                    figure
                    hold on
                    plot((Time./60),ChanGammaPower,'b','LineWidth', 1)
                    plot((Time./60),ChanBetaPower,'r','LineWidth', 1)
                    plot((Time./60),ChanThetaPower,'c','LineWidth', 1)
                    plot((Time./60),ChanDeltaPower,'m','LineWidth', 1)
%                     plot((Time./60),ChanTotalPower,'k','LineWidth', 2)
                    %make the figure and axes titles
                    set(findall(gca, '-property', 'FontSize'),'FontSize',24);
                    xlabel('Time (min)');ylabel('Power in Gamma Band (�V^2)');
                    Title1 = sprintf('Power in CA3 Gamma \nfor electrode %d \nfrom %s',...
                        Channels(FileIndex).CA3Channels(ChannelI),SlicesUsed{FileIndex});
                    Title1 =splitlines(Title1);
                    title(Title1,'Interpreter', 'none')

                end
                CA3ElectrodeCount = CA3ElectrodeCount+1;
                DesiredIndexes = find(Time > T1 & Time < T2);
                ExportCA3Gamma(CA3ElectrodeCount,:) = ChanGammaPower(DesiredIndexes);
                ExportCA3Times(CA3ElectrodeCount,:)  = Time(DesiredIndexes);
            end
        end
        %get the data for CA1
        if ~isempty(Channels(FileIndex).CA1Channels) == 1
            for ChannelI = 1:length(Channels(FileIndex).CA1Channels)
                ChannelIndex    = Channels(FileIndex).CA1Indexes(ChannelI);
                ChanGammaPower  = AllData(FileIndex).MEA_RelativeData(ChannelIndex).RelativePower.LowGamma; 
                ChanBetaPower   = AllData(FileIndex).MEA_RelativeData(ChannelIndex).RelativePower.Beta; 
                ChanThetaPower  = AllData(FileIndex).MEA_RelativeData(ChannelIndex).RelativePower.Theta; 
                ChanDeltaPower  = AllData(FileIndex).MEA_RelativeData(ChannelIndex).RelativePower.Delta; 
                ChanTotalPower  = AllData(FileIndex).MEA_RelativeData(ChannelIndex).RelativePower.Total;  
                Time            = AllData(FileIndex).MEA_RelativeData(ChannelIndex).Time;
                if PlotFigsCA1 == true
                    figure
                    hold on
                    plot((Time./60),ChanGammaPower,'b','LineWidth', 1)
                    plot((Time./60),ChanTotalPower,'k','LineWidth', 2)
                    %make the figure and axes titles
                    set(findall(gca, '-property', 'FontSize'),'FontSize',24);
                    xlabel('Time (min)');ylabel('Power in Gamma Band (�V^2)');
                    Title1 = sprintf('Power in CA1 Gamma \nfor electrode %d \nfrom %s',...
                        Channels(FileIndex).CA1Channels(ChannelI),SlicesUsed{FileIndex});
                    Title1 =splitlines(Title1);
                    title(Title1,'Interpreter', 'none')

                end
                CA1ElectrodeCount = CA1ElectrodeCount+1;
                DesiredIndexes = find(Time > T1 & Time < T2);
                ExportCA1Gamma(CA1ElectrodeCount,:) = ChanGammaPower(DesiredIndexes);
                ExportCA1Times(CA1ElectrodeCount,:)  = Time(DesiredIndexes);
            end
        end
    elseif strcmp(DataType, 'PowSumD_')
        if ~isempty(Channels(FileIndex).CA3Channels) == 1
            for ChannelI = 1:length(Channels(FileIndex).CA3Channels)
                ChannelIndex    = Channels(FileIndex).CA3Indexes(ChannelI);
                ChanPowerStruct  = AllData(FileIndex).MEAPowerSum.Data(ChannelIndex); 
                CA3ElectrodeCount = CA3ElectrodeCount+1;
                
                ChanPowerStruct.Info.Label = Channels(FileIndex).CA3Channels(ChannelI);
                ChanPowerStruct.Info.Slice = SlicesUsed{FileIndex};
                %get the data
                ExportPower.CA3(CA3ElectrodeCount) = ChanPowerStruct;
            end
        end
        if ~isempty(Channels(FileIndex).CA1Channels) == 1
            for ChannelI = 1:length(Channels(FileIndex).CA1Channels)
                ChannelIndex    = Channels(FileIndex).CA1Indexes(ChannelI);
                ChanPowerStruct  = AllData(FileIndex).MEAPowerSum.Data(ChannelIndex); 
                CA1ElectrodeCount = CA1ElectrodeCount+1;
                
                ChanPowerStruct.Info.Label = Channels(FileIndex).CA1Channels(ChannelI);
                ChanPowerStruct.Info.Slice = SlicesUsed{FileIndex};

                ExportPower.CA1(CA1ElectrodeCount) = ChanPowerStruct;

            end
        end
    end
end
%clean up the workspace
clear ChannelI ChannelIndex DesiredIndexes DesiredTimes NumTimePoints T1 T2 Time1 Time2 TotalCA1 TotalCA3 
%%
%plot a figure with the area for CA3 and CA1 gamma band power across the
%recording using Time0 as the 0 time in the figure
Time0 = 60;
LimitX = [-10 60];LimitY = [0 20];
figure; hold on
MeanandSEMPlotter(ExportCA3Gamma,CA3ElectrodeCount,ExportCA3Times,'CA3',Time0,'k')
xlim(LimitX);ylim(LimitY);
figure; hold on
MeanandSEMPlotter(ExportCA1Gamma,CA1ElectrodeCount,ExportCA1Times,'CA1',Time0,'k')
xlim(LimitX);ylim(LimitY);
%%
%----------------------------%
%Rise Time Calculation Set up
%This will calculate a fall time for the selected data
Times2Analyze = [58, 120]; %in minutes
SmoothIt = true; fitWindow = 30;%in seconds
PlotCA1 = true; PlotCA3 = true;
%----------------------------%
Times2Analyze = Times2Analyze*60;%convert to seconds
[CA3RiseTime, CA3StepData] = RiseTimeCalc(ExportCA3Gamma,ExportCA3Times(1,:),Times2Analyze,CA3ElectrodeCount,SmoothIt,fitWindow,PlotCA3);
% [CA1RiseTime, CA1StepData] = RiseTimeCalc(ExportCA1Gamma,ExportCA1Times(1,:),Times2Analyze,CA1ElectrodeCount,SmoothIt,fitWindow,PlotCA1);
%%
%----------------------------%
%Fall Calculation Set up
%This will calculate a fall time for the selected data
Times2Analyze = [86, 88+10]; %in minutes
PlotCA1 = true; PlotCA3 = true;
%----------------------------%

Times2Analyze = Times2Analyze*60;%convert to seconds
CA3FallTimeData = FallTimeCalc(ExportCA3Gamma,ExportCA3Times(1,:),Times2Analyze,CA3ElectrodeCount,PlotCA3);
CA1FallTimeData = FallTimeCalc(ExportCA1Gamma,ExportCA1Times(1,:),Times2Analyze,CA1ElectrodeCount,PlotCA1);
%%
%plot a figure with the area for CA3 and CA1 for two groups together
%NOTE: To get this need to run the top section twice and save the data in
%the different seperate structures using the Groups here.
%-----------------------------%
%Plot set up
Time0 = 60;
Group1 = 'Wild Type';Col1 = [0,0,0];
Group3 = 'GluN2A+/-';Col3 = [0,0,1];
Group2 = 'Shaffer Cut';Col2 = [1,0,0];
Group4 = 'MK801 D7'; Col4 = [1,.5,.1];
Color = {Col1,Col2};%
LimitX = [-10 60];LimitY = [0 30];
V1 = G1;
V2 = G2;
%-----------------------------%
figure(1002);hold on
MeanandSEMPlotter(V1.ExportCA3Gamma,V1.CA3ElectrodeCount,V1.ExportCA3Times,'CA3',Time0,Color{1})
MeanandSEMPlotter(V2.ExportCA3Gamma,V2.CA3ElectrodeCount,V2.ExportCA3Times,'CA3',Time0,Color{2})
xlim(LimitX);ylim(LimitY);
hold off
figure(1003);hold on
MeanandSEMPlotter(V1.ExportCA1Gamma,V1.CA1ElectrodeCount,V1.ExportCA1Times,'CA1',Time0,Color{1})
MeanandSEMPlotter(V2.ExportCA1Gamma,V2.CA1ElectrodeCount,V2.ExportCA1Times,'CA1',Time0,Color{2})
xlim(LimitX);ylim(LimitY);
 hold off
%=========================================================================%
% The next two section are for seeting up this section
%%
G1.ExportCA1Gamma = ExportCA1HGamma;G1.ExportCA3Gamma = ExportCA3HGamma;
G1.CA1ElectrodeCount = CA1ElectrodeCount;G1.CA3ElectrodeCount = CA3ElectrodeCount;
G1.ExportCA3Times=ExportCA3Times;G1.ExportCA1Times=ExportCA1Times;
clearvars -except G1 G2 G3 G4
%%
G2.ExportCA1Gamma = ExportCA1HGamma;G2.ExportCA3Gamma = ExportCA3HGamma;
G2.CA1ElectrodeCount = CA1ElectrodeCount;G2.CA3ElectrodeCount = CA3ElectrodeCount;
G2.ExportCA3Times=ExportCA3Times;G2.ExportCA1Times=ExportCA1Times;
clearvars -except G1 G2 G3 G4
%%
G3.ExportCA1Gamma    = ExportCA1Gamma;   G3.ExportCA3Gamma = ExportCA3Gamma;
G3.CA1ElectrodeCount = CA1ElectrodeCount;G3.CA3ElectrodeCount = CA3ElectrodeCount;
G3.ExportCA3Times    = ExportCA3Times;   G3.ExportCA1Times=ExportCA1Times;
clearvars -except G1 G2 G3 G4
%%
G4.ExportCA1Gamma = ExportCA1Gamma;G4.ExportCA3Gamma = ExportCA3Gamma;
G4.CA1ElectrodeCount = CA1ElectrodeCount;G4.CA3ElectrodeCount = CA3ElectrodeCount;
G4.ExportCA3Times=ExportCA3Times;G4.ExportCA1Times=ExportCA1Times;
clearvars -except G1 G2 G3 G4
%=========================================================================%
%%
%----------------------------%
%Use this section for Periodogram data only
Region = {'CA1','CA3'}; SmoothIt=true;
PlotCA1 =false; PlotCA3 = false;
Stages2Use = 1:2; % 1 = Baseline; 2 = Kainate; 3 = Kainate+Drug
%----------------------------%

[CA1GammaPower,CA1NormGammaPower,CA1Qvalues,CA1PeakValue,CA1HalfBandwidth] = PowerParamsCalculator(ExportPower,Region{1},CA1ElectrodeCount,SmoothIt,PlotCA1,Stages2Use);
[CA3GammaPower,CA3NormGammaPower,CA3Qvalues,CA3PeakValue,CA3HalfBandwidth] = PowerParamsCalculator(ExportPower,Region{2},CA3ElectrodeCount,SmoothIt,PlotCA3,Stages2Use);

%%
function [FileNames,SlicesUsed,Channels] = DataSetChooser(Data2Choose,MainDirectory,DataType,SliceNames,Labels)
%this function will load the desired data set

%NOTE: These indexes need to be updated every time the data set is updated
ShafferCut = [1:4];
Glun2A01 = [];
Glun2A11 = [];
AssayDev = [6:16];
%20 = Glun 2A-/-

if Data2Choose == 20
    %homo data
    for ii = 1:length(ShafferCut)
    FileNames{ii}  = [MainDirectory DataType SliceNames{ShafferCut(ii)}];
    SlicesUsed{ii} = SliceNames{ShafferCut(ii)};
    [Channels(ii).CA1Channels, Channels(ii).CA3Channels, Channels(ii).CA1Indexes, Channels(ii).CA3Indexes] ...
        = CA1andCA3_Loader(SliceNames{ShafferCut(ii)},Labels);
    end
elseif Data2Choose == 21 
    % Het data
    for ii = 1:length(Glun2A01)
    FileNames{ii}  = [MainDirectory DataType SliceNames{Glun2A01(ii)}];
    SlicesUsed{ii} = SliceNames{Glun2A01(ii)};
    [Channels(ii).CA1Channels, Channels(ii).CA3Channels, Channels(ii).CA1Indexes, Channels(ii).CA3Indexes] ...
        = CA1andCA3_Loader(SliceNames{Glun2A01(ii)},Labels);
    end
elseif Data2Choose == 22
    %wt data
    for ii = 1:length(Glun2A11)
    FileNames{ii}  = [MainDirectory DataType SliceNames{Glun2A11(ii)}];
    SlicesUsed{ii} = SliceNames{Glun2A11(ii)};
    [Channels(ii).CA1Channels, Channels(ii).CA3Channels, Channels(ii).CA1Indexes, Channels(ii).CA3Indexes] ...
        = CA1andCA3_Loader(SliceNames{Glun2A11(ii)},Labels);
    end
elseif Data2Choose == 50
    %Assay Dev data
    for ii = 1:length(AssayDev)
    FileNames{ii}  = [MainDirectory DataType SliceNames{AssayDev(ii)}];
    SlicesUsed{ii} = SliceNames{AssayDev(ii)};
    [Channels(ii).CA1Channels, Channels(ii).CA3Channels, Channels(ii).CA1Indexes, Channels(ii).CA3Indexes] ...
        = CA1andCA3_Loader(SliceNames{AssayDev(ii)},Labels);
    end
end

end %function end
%-------------------------------------------------------------------------%
function [CA1Channels, CA3Channels, CA1Indexes, CA3Indexes] = CA1andCA3_Loader(SliceName,Labels)
%this function will output the labels and corresponding indexes for CA1 and
%CA3 electrodes for the slices . 
%-------------------------------------------------------------------------%
%Schaffer Cut
if isequal(SliceName,'2022-01-16_Slice1')
    CA1Channels =[34,25,15,24,23,22,11,33,32,31,43];
    CA3Channels =[61,62,71,72,81,82,83,101,92,93,104,84,74,75];   
elseif isequal(SliceName,'2022-01-17_Slice1')
    CA1Channels =[52,51,53,63,61,62,71,72,73,81,82,92,102,93,103,104,94,105,95,106,96,84,85];
    CA3Channels =[55,44,35,34,16,25,15,24,23,12,22,21];   
elseif isequal(SliceName,'2022-03-03_Slice1')
    CA1Channels =[56,55,46,45,44,36,35,34,26,16,25,15,24,23,12,22,11,21,32,43];
    CA3Channels =[72,81,82,83,101,92,102,103,104,94,105,95,106,96,84,85]; 
elseif isequal(SliceName,'2022-03-03_Slice2')
    CA1Channels =[34,26,16,25,15,24,23,12,22,11,21,33,32,31,43,42,41];
    CA3Channels =[71,72,73,81,82,83,103,105,106,96,84,85,86,74,75,76,66]; 
%-------------------------------------------------------------------------%
%AssayDevData
elseif isequal(SliceName,'2018-08-27_Slice2')
    CA1Channels =[85 84 96 106 95 105 94 104 103 93 102 101 83 82 73];%15 electrodes
    CA3Channels =[72 63 53 52 42 43 ]; %6 electrodes
elseif isequal(SliceName,'2018-08-28_Slice1')
    CA1Channels =[85 84 96 106 95 105 91 83 73 81];%10 electrodes
    CA3Channels =[71 62 61 52 42 43 21 11 12];%9 electrodes
elseif isequal(SliceName,'2018-09-12_Slice1')
    CA1Channels =[];
    CA3Channels =[76 75 74 86 85 84 104 61];%8 electrodes
elseif isequal(SliceName,'2018-09-12_Slice2')
    CA1Channels =[31 33 21 11 22 12 13 24 15 25 16 26 34 35]; %14 electrodes
    CA3Channels =[76 75 74 86 85 84 96 106 95 105 94 104 103 93 102 92 101 91 83 82 81 73 72 71 62 61 63 53 52];%29 electrodes
elseif isequal(SliceName,'2018-09-14_Slice1')
    CA1Channels =[15 16 26]; %3 electrodes
    CA3Channels =[86 96 81]; %4 electrodes removed 92
elseif isequal(SliceName,'2018-09-14_Slice2')
    CA1Channels =[75 74 86 85 84 96 106 95 105 94 104 101 91 83 81 92];%18 electrodes %removed 93 102
    CA3Channels =[64 72 73 71 62 61 63 53 52 43 45 46 16 26 35 54];%17 electrodes %removed 12
elseif isequal(SliceName,'2020-02-27_Slice1')
    CA1Channels =[62 61];%2 electrodes
    CA3Channels =[51 52 41 42 43 31 32 33 21 23 24 15 25 26 34 35 44]; %17 electrodes
elseif isequal(SliceName,'2020-02-27_Slice2')
    CA1Channels =[85 84 96 95 94 83 82 81 73 72 63];%12 electrodes
    CA3Channels =[];
elseif isequal(SliceName,'2018-09-18_Slice1')
    CA1Channels =[64 42 43 25 55 56 54];%7 electrodes
    CA3Channels =[74 86 85 84 96 106 95 105 94 104 103 102 101 91 83 82 81 73 72 71 62 61 63 53 51 52];%26 electrodes
elseif isequal(SliceName,'2018-09-18_Slice2')  
    CA1Channels =[43 23 13];%3 electrodes
    CA3Channels =[86 74 85 84 96 106 95 105 94 104 83 82 81 73 72 71 61 63 53 51 52];%21 electrodes

elseif isequal(SliceName,'2020-08-19_Slice2')
    CA1Channels =[85 84 96 106 95 105 94 104 103 93 92 72];%12 electrodes
    CA3Channels =[71 62 53 52 41 42 43 31 32 21 11 22 23 13 24 34 35 36 44 26];%20 electrodes
    
elseif isequal(SliceName, '2020-08-19_Slice1')
    CA1Channels =[43 33 12 13 15 25 16 26 34 35 36];%11 electrodes
    CA3Channels =[86 85 84 96 106 95 105 94 104 103 93 83 81 73 62 52];%16 electrodes
elseif isequal(SliceName,'2020-08-21_Slice1')
    CA1Channels =[];
    CA3Channels =[];
elseif isequal(SliceName,'2017-12-01_Slice2')
    CA1Channels =[53 63 61 71 72 73 83 91 102 93 104 105 106 96 86 74 75 66];
    CA3Channels =[34 15 13 12 11 33 32 31];
elseif isequal(SliceName,'2017-12-06_Slice3')
    CA1Channels =[72 71 62 61 43 33 21 22 12 24 25];
    CA3Channels =[104 103 102 83 82];
elseif isequal(SliceName,'2018-01-17_Slice2')
    CA1Channels =[];
    CA3Channels =[83 91 93 103 104 94 105 95 106 96 75];
elseif isequal(SliceName,'2018-01-11_Slice2')
    CA1Channels =[25 26];
    CA3Channels =[75 76 104 84 106 95 105 103 93 102 92 101 91];
elseif isequal(SliceName,'2018-03-02_Slice2')
    CA1Channels =[76 75 74 86 85 84 96 106 95 105 94 103 102 92 91 83 81 73 72];
    CA3Channels =[65 12 62 63 51 52 43 32 11 13 54 15 16];
elseif isequal(SliceName,'2017-12-08_Slice3')
    CA1Channels =[106 95 105 94 92 63];
    CA3Channels =[41 43 31 33 21 12 23 13 24 15 25 46 54];
elseif isequal(SliceName,'2017-12-14_Slice1')
    CA1Channels =[106 105 104 103 92 91 54];
    CA3Channels =[43 31 32 33 21 22 23 24 15 25 46];
elseif isequal(SliceName,'2017-12-14_Slice2')
    CA1Channels =[76 75 85 106 105 104 103 92 91 63];
    CA3Channels =[41 43 31 32 33 21 22 23 13 24 15 25];
elseif isequal(SliceName,'2018-02-22_Slice1')
    CA1Channels =[51 52 34 26 16 25 15 24 23 12 22 11 21 33 32 31 43 42 41];
    CA3Channels =[71 82 83 91 92 102 103 104 105];
elseif isequal(SliceName,'2018-02-22_Slice3')
    CA1Channels =[62 71 73 81 82 83 91 92 102 93 104 94 105 95 106 96 84 85 74 75];
    CA3Channels =[44 34 24 23 22 21 33 32 31 43 42 41 52 51];   
%-------------------------------------------------------------------------%
end
    
%preallocate the variables
CA1Indexes = zeros(length(CA1Channels),1);
CA3Indexes = zeros(length(CA3Channels),1);
for ChannelsIndex = 1:60 
    %set upt the search variables
    ChannelLabel = str2double(Labels{ChannelsIndex});
    ChaninCA1 = false;
    %check if the channel is in ca1
    for CA1_Chan = 1:length(CA1Channels)   
        if CA1Channels(CA1_Chan) == ChannelLabel
            CA1Indexes(CA1_Chan) = ChannelsIndex;
            ChaninCA1 = true;
        end
    end
    %Determine if the channel is not in CA1 to then check if it is in CA3
    if ChaninCA1 == false
        for CA3_Chan = 1:length(CA3Channels) 
            if CA3Channels(CA3_Chan) == ChannelLabel
                CA3Indexes(CA3_Chan) = ChannelsIndex;
            end
        end
    end
end

end%end of function
%-------------------------------------------------------------------------%
function MeanandSEMPlotter(ExportGammaData,ElectrodeCount,ExportTimes,RegionString,Time0,ColorChosen)
areaPower = [];
averagePower = mean(ExportGammaData);
semPower     = std(ExportGammaData)/sqrt(ElectrodeCount);
areaPower(1,:)    = averagePower+semPower; areaPower(2,:)    = averagePower-semPower; 
TimesAxis = ExportTimes(1,:)./60;
TimesAxis = TimesAxis - Time0;%resets the time to start at 0
% if strcmp(ColorChosen,'k')
%     patch([TimesAxis fliplr(TimesAxis)], [areaPower(1,:) fliplr(areaPower(2,:))], [0.9 0.9 0.9])
% elseif strcmp(ColorChosen,'r')
%     patch([TimesAxis fliplr(TimesAxis)], [areaPower(1,:) fliplr(areaPower(2,:))], [1 0.6 0.6])
% else
%     patch([TimesAxis fliplr(TimesAxis)], [areaPower(1,:) fliplr(areaPower(2,:))], [0.9 0.9 0.9])
% end
plot(TimesAxis,areaPower(1,:),'Color',ColorChosen,'LineWidth',0.1)
plot(TimesAxis,areaPower(2,:),'Color',ColorChosen,'LineWidth',0.1)
plot(TimesAxis,averagePower,'Color',ColorChosen,'LineWidth',2)
xlabel('Time (min)');ylabel('Power in Gamma Band (�V^2)');
T1 = sprintf('Power in %s Gamma', RegionString);
title(T1,'Interpreter', 'none')
end %end of function
%-------------------------------------------------------------------------%
function [RiseTime, StepData] = RiseTimeCalc(GammaDataMatrix,TimeVector,Times2Analyze,TotalElectrodes,SmoothIt,fitWindow,PlotIt)
%
dt = TimeVector(2)-TimeVector(1);
Window4fit = fitWindow/dt;
Window4SteadyState = floor(5*60/dt);%desired time window converted to seconds rounded down
DesiredIndexes = find(TimeVector > Times2Analyze(1) & TimeVector < Times2Analyze(2));
RiseTime = struct(); StepData = struct();
for ElectrodeNum = 1:TotalElectrodes
    if SmoothIt == 1
        Data2Analyze = smoothdata(GammaDataMatrix(ElectrodeNum,:),'movmedian',Window4fit);
        Data2Analyze = Data2Analyze(DesiredIndexes);
    else
        Data2Analyze = GammaDataMatrix(ElectrodeNum,DesiredIndexes);
    end
    [r,lt,ut,ll,ul] = risetime(Data2Analyze...
                    ,TimeVector(DesiredIndexes),'Tolerance', 7.5);
    RiseTime(ElectrodeNum).RiseTime = r;
    RiseTime(ElectrodeNum).LowerCrossing = lt;
    RiseTime(ElectrodeNum).UpperCrossing = ut;
    RiseTime(ElectrodeNum).LowerLevel = ll;
    RiseTime(ElectrodeNum).UpperLevel = ul;

    SteadyStateResponse = mean(Data2Analyze((end-Window4SteadyState):end));
    SID = stepinfo(Data2Analyze...
        ,TimeVector(DesiredIndexes)...
        ,SteadyStateResponse,'SettlingTimeThreshold',0.05);
    StepData(ElectrodeNum).RiseTime = SID.RiseTime;
    StepData(ElectrodeNum).SettlingTime = SID.SettlingTime;

    %plot the figure if desired
    if PlotIt == 1
        figure(500+ElectrodeNum)
        hold on
        plot((TimeVector(DesiredIndexes)./60),Data2Analyze,'r')
        if SmoothIt == 1
            plot((TimeVector(DesiredIndexes)./60),GammaDataMatrix(ElectrodeNum,DesiredIndexes),'k')
        end
        if ~isempty(r)
            if length(lt) > 1
                llv = ones(1,length(lt));llv=llv*ll;
                ulv = ones(1,length(lt));ulv=ulv*ul;
            else
                llv = ll; ulv = ul;
            end
            plot(lt./60,llv,'ro')
            plot(ut./60,ulv,'ro')
            T1 =sprintf('RiseTime = %d',r);
            yl = ylim;

            plot([SID.RiseTime SID.RiseTime]./60,yl,'-.r')
            plot([SID.SettlingTime SID.SettlingTime]./60,yl,'--m')
            title(T1)
            xlim(Times2Analyze./60)
            hold off
        else
            fprintf('No risetime detected for electrode id %d \n',ElectrodeNum)
        end
        clear llv ulv r lt ut ll ul
    end
end
end %end of function
%-------------------------------------------------------------------------%
function FallTimeData = FallTimeCalc(GammaDataMatrix,TimeVector,Times2Analyze,TotalElectrodes,PlotIt)

DesiredIndexes = find(TimeVector > Times2Analyze(1) & TimeVector < Times2Analyze(2));
for ElectrodeNum = 1:TotalElectrodes
[f,lt,ut,ll,ul] = falltime(GammaDataMatrix(ElectrodeNum,DesiredIndexes)...
                ,TimeVector(DesiredIndexes),'Tolerance', 5);
            FallTimeData(ElectrodeNum).FallTime = f;
            FallTimeData(ElectrodeNum).LowerCrossing = lt;
            FallTimeData(ElectrodeNum).UpperCrossing = ut;
            FallTimeData(ElectrodeNum).LowerLevel = ll;
            FallTimeData(ElectrodeNum).UpperLevel = ul;
            
            %plot the figure if desired
            if PlotIt == 1
                
                figure(1000+ElectrodeNum)
                hold on
                plot((TimeVector(DesiredIndexes)./60),GammaDataMatrix(ElectrodeNum,DesiredIndexes))
                if ~isempty(f)
                    if length(lt) > 1
                        llv = ones(1,length(lt));llv=llv*ll;
                        ulv = ones(1,length(lt));ulv=ulv*ul;
                    else
                        llv = ll; ulv = ul;
                    end
                    plot(lt./60,llv,'ro')
                    plot(ut./60,ulv,'ro')
                    T1 =sprintf('FallTime = %d',f);
                    title(T1)
                    xlim(Times2Analyze./60)
                    hold off
                else
                    title('No Fall Time Detected')
                end
                
            end
            
end
end %end of function
%-------------------------------------------------------------------------%
function PlotPeriodogram(MEAPowerSum)
plotInd = false; plotAll = true;SmoothIt = true; 
DesiredFreqWindow = 1;

for ElectrodeIndex = 50
    ElectrodeNumber = str2double(SliceChannelInfo.Label{ElectrodeIndex});
    if plotAll == true
        figure
        xlabel('Frequency(Hz)')
        ylabel('Power (dB)')
        T1 = sprintf('Channel %d',ElectrodeNumber);
        title(T1)
        hold on
    end
    for StageIndex = 1:3
        %get the data for the specific stage
        if StageIndex == 1    
            y=MEAPowerSum.Data(ElectrodeIndex).Baseline.Power;
            x=MEAPowerSum.Data(ElectrodeIndex).Baseline.FreqData;
        elseif StageIndex == 2 
            y=MEAPowerSum.Data(ElectrodeIndex).Kainate.Power;
            x=MEAPowerSum.Data(ElectrodeIndex).Kainate.FreqData;
        elseif StageIndex == 3 
            y=MEAPowerSum.Data(ElectrodeIndex).KainateDrug.Power;
            x=MEAPowerSum.Data(ElectrodeIndex).KainateDrug.FreqData;
        end
        deltax = x(2)-x(1);
        window = floor(DesiredFreqWindow/deltax);
        if plotInd == 1
            figure
            if SmoothIt == 0
                plot(x,10*log10(y))
            elseif SmoothIt == 1
                [ys, window]= smoothdata(10*log10(y),'movmean',window);    
                plot(x,ys)
                clear ys 
            end
        end
        if plotAll == 1
            if SmoothIt == 0
                plot(x,10*log10(y))
                title(T1)
            elseif SmoothIt == 1
                [ys, window] = smoothdata(10*log10(y),'movmean',window);       
                plot(x,ys)
                clear ys 
            end
        end
    end
end
end%end of function
%-------------------------------------------------------------------------%
function [GammaPower,NormGammaPower,Qvalue,PeakValue,HalfBandwidth] = PowerParamsCalculator(ExportPower,Region,NumElectrodes,SmoothIt,plotIt,Stages2Use)
%Preallocate the variables to be filled later on in the algorithms
NumStages = length(Stages2Use);
GammaPower = zeros(NumElectrodes,NumStages); NormGammaPower = zeros(NumElectrodes,NumStages);
GammaLF = [25,59];FreqBounds = [5 100];DesiredFreqWindow = 1;
Qvalue=zeros(NumElectrodes,1);PeakValue=zeros(NumElectrodes,1);HalfBandwidth=zeros(NumElectrodes,1);

%Do the Power params calculation for each electrode
for ElectrodeIndex = 1:NumElectrodes
    DataPoints = length(ExportPower.(Region)(ElectrodeIndex).Baseline.Power);
    PowerData = zeros(3,DataPoints);Frequencies = zeros(3,DataPoints);
    %run the calculations for each stage
    for StageIndex = Stages2Use
        if StageIndex == 1    
            PowerData(StageIndex,:)=ExportPower.(Region)(ElectrodeIndex).Baseline.Power;
            Frequencies(StageIndex,:)=ExportPower.(Region)(ElectrodeIndex).Baseline.FreqData;
            GammaPower(ElectrodeIndex,StageIndex) = AUC((PowerData(StageIndex,:))', Frequencies(StageIndex,:), GammaLF);
            NormGammaPower(ElectrodeIndex,StageIndex) = (GammaPower(ElectrodeIndex,StageIndex)-GammaPower(ElectrodeIndex,1))./GammaPower(ElectrodeIndex,1); 
        elseif StageIndex == 2 
            PowerData(StageIndex,:)=ExportPower.(Region)(ElectrodeIndex).Kainate.Power;
            Frequencies(StageIndex,:)=ExportPower.(Region)(ElectrodeIndex).Kainate.FreqData;
            GammaPower(ElectrodeIndex,StageIndex) = AUC((PowerData(StageIndex,:))', Frequencies(StageIndex,:), GammaLF);
            NormGammaPower(ElectrodeIndex,StageIndex) = (GammaPower(ElectrodeIndex,StageIndex)-GammaPower(ElectrodeIndex,1))./GammaPower(ElectrodeIndex,1);
            if SmoothIt == 0
                Power4Qval = PowerData(StageIndex,:)';
            elseif SmoothIt == 1
                deltax = Frequencies(StageIndex,2)-Frequencies(StageIndex,1);
                window = floor(DesiredFreqWindow/deltax);
%                 [Power4Qval2, ~]= smoothdata(PowerData(StageIndex,:)','movmean',window); 
%                 [Power4Qval1, ~]= smoothdata((PowerData(1,:)'),'movmean',window);
%                 Power4Qval = Power4Qval2./Power4Qval1;
                Power4Qval = ((PowerData(StageIndex,:))'-(PowerData(1,:)'))./(PowerData(1,:)');
                Power4Qval = smoothdata(Power4Qval,'movmean',window);
%                 if plotIt == 1
%                 figure(2000+ElectrodeIndex)
%                 plot(Frequencies(StageIndex,:),Power4Qval)
%                 T1 =sprintf('Periodogram of Chan %d from %s',ExportPower.(Region)(ElectrodeIndex).Info.Label,ExportPower.(Region)(ElectrodeIndex).Info.Slice);
%                 title(T1)
%                 end
            end
            [Qvalue(ElectrodeIndex), PeakValue(ElectrodeIndex), HalfBandwidth(ElectrodeIndex)] =...
                M_QvalueCalculator(Power4Qval,Frequencies(StageIndex,:),FreqBounds); 
        elseif StageIndex == 3 
            PowerData(StageIndex,:)=ExportPower.(Region)(ElectrodeIndex).KainateDrug.Power;
            Frequencies(StageIndex,:)=ExportPower.(Region)(ElectrodeIndex).KainateDrug.FreqData;
            GammaPower(ElectrodeIndex,StageIndex) = AUC((PowerData(StageIndex,:))', Frequencies(StageIndex,:), GammaLF);
            NormGammaPower(ElectrodeIndex,StageIndex) = (GammaPower(ElectrodeIndex,StageIndex)-GammaPower(ElectrodeIndex,2))./GammaPower(ElectrodeIndex,2);
        end
    end
    if plotIt==1
        figure; hold on
        for ii = Stages2Use   
            if SmoothIt == 0
                y = 10*log10(PowerData');
            elseif SmoothIt == 1
                deltax = Frequencies(ii,2)-Frequencies(ii,1);
                window = floor(DesiredFreqWindow/deltax);
                [y(ii,:), ~]= smoothdata(10*log10(PowerData(ii,:)'),'movmean',window); 
            end
        plot(Frequencies(ii,:),y(ii,:));
        end
        T1 =sprintf('Periodogram of Chan %d from %s',ExportPower.(Region)(ElectrodeIndex).Info.Label,ExportPower.(Region)(ElectrodeIndex).Info.Slice);
        title(T1)
        hold off
        clear y PowerData PowerData
    end
end
end %end of function
%-------------------------------------------------------------------------%
function [AUCPower] = AUC(Power, Frequency, FreqBoundaries) 
%This function will estimate the Peak power at the desired frequency range
%using the area under the curve. The area under the curve is estimated by
%using MATLAB's trapz function. 
%Inputs
%   Power     - Power calculated by the FFT
%   Frequency - Frequency 
%   Freq1     - Lower end of the desired frequency range
%   Freq2     - Higher end of the desired frequency range
%Output
%   AUCPower - Area under the curve of the designated frequency ranges

%Define the power and the frequency 
x = Frequency; 
%Determine the indexes for the data set that focuses on the desired
%frequency range
margingOfError = x(2)-x(1);
LeftI  = find(x >= (FreqBoundaries(1)- margingOfError) & x <= (FreqBoundaries(1)+ margingOfError)); LeftI  = LeftI(1);
RightI = find(x >= (FreqBoundaries(2)- margingOfError) & x <= (FreqBoundaries(2)+ margingOfError)); RightI = RightI(end); 
[~,NumTimePoints] = size(Power);
AUCPower=zeros(1,NumTimePoints);
for TimePoint = 1:NumTimePoints
    y = (Power(:,TimePoint));
    %Calculate the area under the curve
    AUCPower(TimePoint) = trapz(x(LeftI:RightI), y(LeftI:RightI));
end
end%end of function
%-------------------------------------------------------------------------%
function [QValue, PeakFreq, HalfBandwidth] = M_QvalueCalculator(Power,Frequency,FreqBoundaries)
%This function is used to calcuate a q value using the peak power and the
%half bandwidth. 
%the definition of this was taken from C.E. Lemercier et al. / Schizophrenia Research 188 (2017) 118�124
%Inputs
%   Power
%   Frequency
%   FreqUpper
%   FreqLower
%Output
%   QPower
%
%NOTE:
%Standard frequency ranges based on literature and multichannel systems 
%Delta : 0.5 - 4 Hz
%Theta : 4   - 8 Hz
%Alpha : 8   - 13 Hz
%Beta  : 13  - 30 Hz
%Gamma : 30  - 80 Hz
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN:
%LAST MODIFIED: January 3, 2018
%v1.0


%Define the power and the frequency 
x = Frequency; 
%Determine the indexes for the data set that focuses on the desired
%frequency range
margingOfError = x(2)-x(1);
LowerEnd  = find(x >= (FreqBoundaries(1)- margingOfError) & x <= (FreqBoundaries(1)+ margingOfError)); LowerEnd  = LowerEnd(1);
UpperEnd = find(x >= (FreqBoundaries(2)- margingOfError) & x <= (FreqBoundaries(2)+ margingOfError)); UpperEnd = UpperEnd(end); 
%Determine peak frequency of the Power region of interest
if UpperEnd > 60
PowerInterest     = Power(LowerEnd:UpperEnd); 
FrequencyInterest = Frequency(LowerEnd:UpperEnd);
LowerElim  = find(FrequencyInterest >= (59- margingOfError) & FrequencyInterest <= (59 + margingOfError)); LowerEnd  = LowerEnd(1);
UpperElim = find(FrequencyInterest >= (61- margingOfError) & FrequencyInterest <= (61+ margingOfError)); UpperEnd = UpperEnd(end); 
PowerInterest(LowerElim:UpperElim) = 0;
else
PowerInterest     = Power(LowerEnd:UpperEnd);
FrequencyInterest = Frequency(LowerEnd:UpperEnd);
end
[PeakFreqPower, PeakIndex] = max(PowerInterest);

%Determine half bandwidth
HalfPeakPower = PowerInterest(PeakIndex)/2;

%Find the two values that will be used to caluclate the bandwidth
 BandwidthPositions =[];
 IndexCont = 1;
 for ii = 1:length(PowerInterest)
     if PowerInterest(ii) > HalfPeakPower
         BandwidthPositions(IndexCont) = ii;
         IndexCont = IndexCont + 1;
     end
 end


%Determine the Half Bandwidth
HalfBandwithL  = FrequencyInterest(BandwidthPositions(1));
HalfBandwidthU = FrequencyInterest(BandwidthPositions(end));

HalfBandwidth = HalfBandwidthU - HalfBandwithL;
%Determine the Preak Frequency
TempFreq = Frequency(LowerEnd:UpperEnd);
PeakFreq = TempFreq(PeakIndex);
%Calculate the Q value
QValue = PeakFreq/HalfBandwidth;
end %end of function
