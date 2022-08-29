%Gamma Coherence Analysis Export and Figure Making Script
%This script is used to take the coh files and export them as tables or
%figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Start here
%This section is ran to establish import the coherence data in a desired
%format
%-------------------------------------------------------------------------%
%Define these variables
MainDriveFolder = 'G:\.shortcut-targets-by-id\0B4kHX-ERHmBENE1GUXJHYUQyX0U\Jean Rodriguez Diaz\';

folderdir = [MainDriveFolder 'Data Analysis\SliceData\CoherenceData\'];
DecimatedData = true;
UseJenkinsData = 0;
UseShafferCut = 0;%2 is for loadding it wiht the assay dev data

%Select the time window for the coherence data to be averaged on
TimeWindow2Use = 5; %in minutes
TimeWindow2Use4Rolling = [10,60];
Phases2Analyze = [1,2];         %Baseline Phase = 1; Kainate = 2; Drug = 3;
RollingTime = 0;
RollingWindow = 30; %in seconds
DesiredLayout = 1;                     %1 = SliceOrder; 2 = SliceIndexOrder
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%load the channel data
if UseJenkinsData == 1
    folderdir = [folderdir 'COH_Jenkins'];
    load([folderdir filesep 'Jenkins_Channels.mat'])
elseif UseShafferCut == 1
    folderdir = [folderdir 'COH_Shaffer'];
    load([folderdir filesep 'Shaffer_Channels.mat'])
elseif UseShafferCut == 2
    folderdir = [folderdir 'COH_AssayDev_Shaffer'];
    ShafferChannels = load([folderdir filesep 'Shaffer_Channels.mat']);
    AssayDevChannels = load([folderdir filesep 'AssayDev_Channels.mat']);
    Channels = [AssayDevChannels.Channels ShafferChannels.Channels];
    [~,NumAssayDevChan] = size(AssayDevChannels.Channels);
    [~,NumShafferChans] = size(ShafferChannels.Channels);
elseif UseJenkinsData == 0 && UseShafferCut == 0
    folderdir = [folderdir 'COH_AssayDev'];
    load([folderdir filesep 'AssayDev_Channels.mat'])
end
%-----------------------------------%
DesiredFreq = [25 59]; %this is the gamma band to be used
MEA_SliceIndexOrder = LoadLoadout(2);
for Phase2Analyze = Phases2Analyze
    %--------------------%
    if Phase2Analyze == 1
        DesiredTime(2) = [10];%in minutes
    elseif Phase2Analyze == 2
        DesiredTime(2) =[60];%in minutes    
    end
    if RollingTime == 0
        DesiredTime(1) = DesiredTime(2)-TimeWindow2Use; 
    elseif RollingTime == 1
        DesiredTime(1) = DesiredTime(2)-TimeWindow2Use4Rolling(Phase2Analyze); 
    end
    desiredTimeInSeconds = DesiredTime*60;%Change from minutes to seconds
    Conditions = {'Baseline','Kainate','Drug'};
    
    LayoutMEAInfo = LoadLoadout(DesiredLayout);
    [~,NumofFiles] = size(Channels);
    
    for FileIndex = 1:NumofFiles
    %make sure the data desired is used by using the correct file name
    MEAFileName = ['COH_MEAData_' Channels(FileIndex).SlicesUsed '_' Conditions{Phase2Analyze}];
    %load the coherence file
    load([folderdir filesep MEAFileName '.mat']);
    %-----------------------------------%
    desiredTimeInSeconds(1) = desiredTimeInSeconds(1)+t(1);%Correct for the time gap caused but using a time window for coherence calculation
    %t(1) should roughly equal to half the time window used for in the
    %coherence calculation
%     t(1);
    if RollingTime == 0
        coh2Export = zeros(60,60);
        time_ind  = find(t>desiredTimeInSeconds(1)&t<desiredTimeInSeconds(2));
        coh2Export = mean(coh(time_ind,:,:),1);
        coh2Export = squeeze(coh2Export);
        coh2Export = coh2Export' + coh2Export;
    elseif RollingTime == 1
        Times2Use = (desiredTimeInSeconds(1)+RollingWindow):RollingWindow:desiredTimeInSeconds(2);
        sizeofTimes = length(Times2Use);
        coh2Export = zeros(60,60,sizeofTimes);
        timeIndexCounter = 1;
        for rollingTime = Times2Use
            time_ind  = find(t>(rollingTime-RollingWindow)&t<rollingTime);
            coh2Export(:,:,timeIndexCounter) = mean(coh(time_ind,:,:),1);
            timeIndexCounter = timeIndexCounter+1;
        end
    end
    %incorporate the COH data into a structure that can be used to make figures
    CoherenceData{FileIndex}.(Conditions{Phase2Analyze}) = coh2Export;
    %Save some important setting used to calculate this dataset
    if RollingTime == 0
        TimeWindowDur2Exp = TimeWindow2Use;
    elseif RollingTime == 1
        TimeWindowDur2Exp = TimeWindow2Use4Rolling(Phase2Analyze);
    end
    CoherenceData{FileIndex}.Settings.(Conditions{Phase2Analyze}).TimeWindowDuration = TimeWindowDur2Exp;
    CoherenceData{FileIndex}.Settings.(Conditions{Phase2Analyze}).TimeWindowUsed     = desiredTimeInSeconds;
    CoherenceData{FileIndex}.Settings.(Conditions{Phase2Analyze}).Rolling            = RollingTime;
    
    if RollingTime == 1
        CoherenceData{FileIndex}.Settings.(Conditions{Phase2Analyze}).RollingTime = RollingWindow;
        CoherenceData{FileIndex}.Settings.(Conditions{Phase2Analyze}).Timevector  = Times2Use;  
    end
    clear coh coh2Export Times2Use sizeofTimes time_ind timeIndexCounter TimeWindowDur2Exp
    end
end
%clean workspace

clear DesiredTime desiredTimeInSeconds FileIndex MEAFileName NumofFiles Phase2Analyze RollingWindow t time_ind 
%%
%-------------------------------------------------------------------------%
%
%THIS SEGMENT IS NOT FOR ROLLING TIME
%
%-------------------------------------------------------------------------%
%This segment of the scirpt will organize the data to be in a table format.
%-------------------------------------------------------------------------%
[~,NumFiles] = size(Channels);
%------------------------%
%Set these variables
Files2Use = 1:NumFiles; %Use NumFiles if table should be made out of all files
%------------------------%
IndexCounter = 0;
for FileIndex = Files2Use
    %make sure the settings are correct and rolling time wasn't used
    if CoherenceData{FileIndex}.Settings.Baseline.Rolling == 1||CoherenceData{FileIndex}.Settings.Kainate.Rolling == 1
        disp('WARNING: Rolling Time Window was used.')
        disp('Please use the Rolling Time Window analysis segment further down or change data to not have a rolling time window')
        return
    end
    %Generate the vectors for the variables of interest and organize them
    %in a table
    %Go through CA3 - CA3 coherences
    if ~isempty(Channels(FileIndex).CA3Channels)
        [~,NumCA3] = size(Channels(FileIndex).CA3Channels);
        for ChanCA3Index1 = 1:NumCA3
            for ChanCA3Index2 = (ChanCA3Index1+1):NumCA3
                IndexCounter = IndexCounter + 1;
                Elec1 = Channels(FileIndex).CA3Indexes(ChanCA3Index1);
                Elec2 = Channels(FileIndex).CA3Indexes(ChanCA3Index2);

                CoherenceVector(IndexCounter) = CoherenceData{FileIndex}.Kainate(Elec1,Elec2);
                CoherenceBVector(IndexCounter) = CoherenceData{FileIndex}.Baseline(Elec1,Elec2);
                Elec1Vector{IndexCounter} = sprintf('S%d_C%d',FileIndex,(ChanCA3Index1));
                Elec2Vector{IndexCounter}= sprintf('S%d_C%d',FileIndex,(ChanCA3Index2));
                ElectrodeNum1Vector{IndexCounter} = Channels(FileIndex).CA3Channels(ChanCA3Index1);
                ElectrodeNum2Vector{IndexCounter} = Channels(FileIndex).CA3Channels(ChanCA3Index2);
                SliceVector(IndexCounter)= FileIndex;
                RegionVector{IndexCounter}= 'CA3';
                if UseJenkinsData == 1
                    if FileIndex >5
                        ModelVector{IndexCounter}= 'MT';
                    else
                        ModelVector{IndexCounter}= 'WT';
                    end
                elseif UseShafferCut == 2
                    if FileIndex > NumAssayDevChan
                        ModelVector{IndexCounter}= 'Cut';
                    else
                        ModelVector{IndexCounter}= 'NoCut';
                    end
                end
            end
        end
    end
    %Go through CA1 - CA1 coherences
    if ~isempty(Channels(FileIndex).CA1Channels)
        [~,NumCA1] = size(Channels(FileIndex).CA1Channels);
        for ChanCA1Index1 = 1:NumCA1
            for ChanCA1Index2 = (ChanCA1Index1+1):NumCA1
                IndexCounter = IndexCounter + 1;
                Elec1 = Channels(FileIndex).CA1Indexes(ChanCA1Index1);
                Elec2 = Channels(FileIndex).CA1Indexes(ChanCA1Index2);
                CoherenceVector(IndexCounter) = CoherenceData{FileIndex}.Kainate(Elec1,Elec2);
                CoherenceBVector(IndexCounter) = CoherenceData{FileIndex}.Baseline(Elec1,Elec2);
                Elec1Vector{IndexCounter} = sprintf('S%d_A%d',FileIndex,(ChanCA1Index1));
                Elec2Vector{IndexCounter}= sprintf('S%d_A%d',FileIndex,(ChanCA1Index2));
                ElectrodeNum1Vector{IndexCounter} = Channels(FileIndex).CA1Channels(ChanCA1Index1);
                ElectrodeNum2Vector{IndexCounter} = Channels(FileIndex).CA1Channels(ChanCA1Index2);
                SliceVector(IndexCounter)= FileIndex;
                RegionVector{IndexCounter}= 'CA1';
                if UseJenkinsData == 1
                    if FileIndex >5
                        ModelVector{IndexCounter}= 'MT';
                    else
                        ModelVector{IndexCounter}= 'WT';
                    end
                elseif UseShafferCut == 2
                    if FileIndex > NumAssayDevChan
                        ModelVector{IndexCounter}= 'Cut';
                    else
                        ModelVector{IndexCounter}= 'NoCut';
                    end
                end
            end
        end
    end
    %Go through CA3 - CA1 coherences
    if ~isempty(Channels(FileIndex).CA3Channels)&&~isempty(Channels(FileIndex).CA1Channels)
        for ChanCA3Index = 1:NumCA3
            Elec1 = Channels(FileIndex).CA3Indexes(ChanCA3Index);
            for ChanCA1Index = 1:NumCA1
                IndexCounter = IndexCounter + 1;
                Elec2 = Channels(FileIndex).CA1Indexes(ChanCA1Index);
                CoherenceVector(IndexCounter) = CoherenceData{FileIndex}.Kainate(Elec1,Elec2);
                CoherenceBVector(IndexCounter) = CoherenceData{FileIndex}.Baseline(Elec1,Elec2);
                Elec1Vector{IndexCounter} = sprintf('S%d_C%d',FileIndex,(ChanCA3Index));
                Elec2Vector{IndexCounter}= sprintf('S%d_A%d',FileIndex,(ChanCA1Index));
                ElectrodeNum1Vector{IndexCounter} = Channels(FileIndex).CA3Channels(ChanCA3Index);
                ElectrodeNum2Vector{IndexCounter} = Channels(FileIndex).CA1Channels(ChanCA1Index);
                SliceVector(IndexCounter)= FileIndex;
                RegionVector{IndexCounter}= 'CA1_3';
                if UseJenkinsData == 1
                    if FileIndex >5
                        ModelVector{IndexCounter}= 'MT';
                    else
                        ModelVector{IndexCounter}= 'WT';
                    end
                elseif UseShafferCut == 2
                    if FileIndex > NumAssayDevChan
                        ModelVector{IndexCounter}= 'Cut';
                    else
                        ModelVector{IndexCounter}= 'NoCut';
                    end
                end              
            end
        end
    end
    clear ChanCA3Index1 ChanCA3Index2 ChanCA1index1 ChanCA1index2 NumCA1 NumCA3
end
%Set the Variable Names Desired

Regions = RegionVector';
Slice   = SliceVector';
Elec1   = Elec1Vector';
Elec2   = Elec2Vector';
ElectrodeNum1 = ElectrodeNum1Vector';
ElectrodeNum2 = ElectrodeNum2Vector';
Coherence = CoherenceVector./CoherenceBVector;Coherence=Coherence';
CoherenceKainate = CoherenceVector';
CoherenceBaseline = CoherenceBVector';
DeltaCoherence = CoherenceKainate - CoherenceBaseline;

if UseJenkinsData == 1 || UseShafferCut == 2
     Model = ModelVector';          
end
if UseJenkinsData == 1 || UseShafferCut == 2
    CoherenceTable = table(Model,Regions,Slice,Elec1,Elec2,Coherence,CoherenceKainate,CoherenceBaseline,DeltaCoherence,ElectrodeNum1,ElectrodeNum2);
else
    CoherenceTable = table(Regions,Slice,Elec1,Elec2,Coherence,CoherenceKainate,CoherenceBaseline,DeltaCoherence,ElectrodeNum1,ElectrodeNum2);
end

%calculate the distance of each electrode pair
for DataIndex = 1:height(CoherenceTable)
Electrode1 = CoherenceTable.ElectrodeNum1(DataIndex);
Electrode2 = CoherenceTable.ElectrodeNum2(DataIndex);
ElectrodeDistanceStandard = 100;

FirstNum(1) = floor(Electrode1{1}/10);
FirstNum(2) = floor(Electrode2{1}/10);
SecondNum(1) = Electrode1{1} - (FirstNum(1)*10);
SecondNum(2) = Electrode2{1} - (FirstNum(2)*10);
HorDist  = abs(diff(FirstNum));
VertDist = abs(diff(SecondNum));

Distance = sqrt(((HorDist)^2)+((VertDist)^2));
CoherenceTable.ElecPairDistance(DataIndex)  = Distance;

end

[MModel1, MModel2,D1] = CoherenceMixedModelF(CoherenceTable,UseJenkinsData,UseShafferCut);

%clean up the workspace
clear IndexCounter FileIndex ChanCA3Index1 ChanCA3Index2 ChanCA1index1 ChanCA1index2
%%
%-------------------------------------------------------------------------%
%this section will make the coherence heatmap figures
%-------------------------------------------------------------------------%
%Define these variables
Files2Use = [1:4];
Phase2Analyze = [2];
Conditions = {'Baseline','Kainate','Drug'};
SaveFig = true; %Set true or false SPECIAL CASE 2 = Save figure with the Channel Labels
DesiredFormat = 'svg'; %Set desired format for the automatic figure saving
%-------------------------------------------------------------------------%

if UseJenkinsData == 1
    Fig2SaveDir = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\CoherenceData\Figures\Jenkins\';
elseif UseShafferCut == 1
    Fig2SaveDir = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\CoherenceData\Figures\Shaffer\';
elseif UseJenkinsData == 0 && UseShafferCut == 0
	Fig2SaveDir = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\CoherenceData\Figures\AssayDev\';
end
if ~exist(Fig2SaveDir, 'dir')
   mkdir(Fig2SaveDir)
end
    
for FileIndex = Files2Use
    if CoherenceData{FileIndex}.Settings.Baseline.Rolling == 1||CoherenceData{FileIndex}.Settings.Kainate.Rolling == 1
        disp('WARNING: Rolling Time Window was used.')
        disp('Please use the Rolling Time Window analysis segment further down or change data to not have a rolling time window')
        return
    end
    
    clear ChannelsWithData ChanNumWithData coh NewCOH T2 CA1Electrodes CA3Electrodes
    for Phases2Analyze = Phase2Analyze
        %-----------------------------------------------------------------%
        %set up or get the necessary variables
        coh = CoherenceData{FileIndex}.(Conditions{Phases2Analyze});
        CA1Electrodes = Channels(FileIndex).CA1Indexes;
        CA3Electrodes = Channels(FileIndex).CA3Indexes;
        ChannelsWithData = [CA1Electrodes' CA3Electrodes'];
        ChanNumWithData = [Channels(FileIndex).CA1Channels Channels(FileIndex).CA3Channels];
        NewCOH = coh([ChannelsWithData],[ChannelsWithData]);
        %%----------------------------------%
        %----Make the coherence heatmap----%
        f1 = figure;
        hm = heatmap(NewCOH);
        %load the heatmap color blind friendly color scheme
        inferno=inferno();colormap inferno;
        caxis([0 1]) %set the c axis limits
        T2 = [Conditions{Phases2Analyze}];
        title(T2);
        %-------------%
        %If desired this section will save the figures automatically in
        %desired format
        if SaveFig == 1
            ToSaveFig = [Fig2SaveDir 'COHFig_' Channels(FileIndex).SlicesUsed '_' T2 '.' DesiredFormat];
            saveas(f1,ToSaveFig, DesiredFormat)
        elseif SaveFig == 2
            CustomXLabels = string(ChanNumWithData);
            hm.XDisplayLabels = CustomXLabels;
            hm.YDisplayLabels = CustomXLabels;
            ToSaveFig = [Fig2SaveDir 'FigWithChanLabel' filesep 'COHFig_' Channels(FileIndex).SlicesUsed '_' T2 'ChanLabel.' DesiredFormat];
            saveas(f1,ToSaveFig, DesiredFormat)
        end
    end
end
%clean up the workspace
clear ChannelsWithData ChanNumWithData coh NewCOH T2 CA1Electrodes CA3Electrodes DesiredFormat 
%%
%this will make the spatial coherence map with the location of the electrodes
%-------------------------------------------------------------------------%
%Define these variables
Files2Use = [6];
Phase2Analyze = [1:2];
OnlyShowRegions = true;
WeightMultiplier = 3; %this will change the thickness of the lines to facilitate visualization
WeightsThreshold = .45;
GammaVisScale = 1; %change this number to facilitate the visualization of power in the different regions
Conditions = {'Baseline','Kainate','Drug'};
Thresholds = [3,40]; %Use htese to determine the minimum and maximum size of the marker.
SaveFig = true; %Set true or false 
DesiredFormat = 'emf'; %Set desired format for the automatic figure saving
CA1Color = [1 .5 0]; CA3Color = [0 .2 1]; %set the color in rgb
%-------------------------------------------------------------------------%
MainDirectory = [MainDriveFolder 'Data Analysis\SliceData\AssayDev'];
Fig2SaveDir = [MainDirectory '\CoherenceData\Figures\AssayDev\'];
if UseJenkinsData == 1
    MainDirectory = [ MainDriveFolder 'Data Analysis\SliceData\Jenkins'];
    Fig2SaveDir = [ MainDriveFolder 'Data Analysis\SliceData\CoherenceData\Figures\Jenkins\'];
elseif UseShafferCut == 1
    MainDirectory = [MainDriveFolder 'Data Analysis\SliceData\ShafferCut'];
    Fig2SaveDir = [MainDriveFolder 'Data Analysis\SliceData\CoherenceData\Figures\Shaffer\'];
end
if ~exist(Fig2SaveDir, 'dir')
   mkdir(Fig2SaveDir)
end
MainDirectory = [MainDirectory '\PowerSummaryData\'];
DataType = 'PowSumD_';

LowGammaFreq = [25,59];
%----------------------------------------------------%
SlicePower = struct; 
%Get the power in the gamma band for CA1 and CA3 regions from each slice.
for FileIndex = Files2Use
    if CoherenceData{FileIndex}.Settings.Baseline.Rolling == 1||CoherenceData{FileIndex}.Settings.Kainate.Rolling == 1
        disp('WARNING: Rolling Time Window was used.')
        disp('Please use the Rolling Time Window analysis segment further down or change data to not have a rolling time window')
        return
    end
    SliceFileName = [MainDirectory DataType Channels(FileIndex).SlicesUsed '.mat'];
    AllData = load(SliceFileName);
    %Get the data for all both Baseline and Kainate phase
    for Phases2Analyze = Phase2Analyze
        %Get the CA1 Data
        if ~isempty(Channels(FileIndex).CA1Channels) == 1
            for ChannelI = 1:length(Channels(FileIndex).CA1Channels)
                ChannelIndex     = Channels(FileIndex).CA1Indexes(ChannelI);
                %get the data
                AllPower = AllData.MEAPowerSum.Data(ChannelIndex).(Conditions{Phases2Analyze}).Power;
                AllFreq  = AllData.MEAPowerSum.Data(ChannelIndex).(Conditions{Phases2Analyze}).FreqData;
                [GammaPower] = AUC(AllPower, AllFreq, LowGammaFreq); 
                SlicePower(FileIndex).CA1(ChannelI,Phases2Analyze) = GammaPower;
            end
        end
        %Get the CA3 Data
        if ~isempty(Channels(FileIndex).CA3Channels) == 1
            for ChannelI = 1:length(Channels(FileIndex).CA3Channels)
                ChannelIndex     = Channels(FileIndex).CA3Indexes(ChannelI);              
                %get the data
                AllPower = AllData.MEAPowerSum.Data(ChannelIndex).(Conditions{Phases2Analyze}).Power;
                AllFreq  = AllData.MEAPowerSum.Data(ChannelIndex).(Conditions{Phases2Analyze}).FreqData;
                [GammaPower] = AUC(AllPower, AllFreq, LowGammaFreq); 
                SlicePower(FileIndex).CA3(ChannelI,Phases2Analyze) = GammaPower;
            end
        end
    end
    clear AllData
end
%clean up the workspace
clear FileIndex Phases2Analyze ChannelI ChannelIndex GammaPower AllPower AllFreq
%-------------------------------------------------------------------------%
%set up or get the necessary variables for the figure
ChannelLabelLocation = [MainDriveFolder '\Data Analysis\SliceData\Slice_ChannelLabelInfo.mat'];
load(ChannelLabelLocation)
ChanLabels = SliceChannelInfo.Label;
for FileIndex = Files2Use
    for Phases2Analyze = Phase2Analyze
        %get the data for the figure
        coh = CoherenceData{FileIndex}.(Conditions{Phases2Analyze});
        %get the channel information
        CA3Channels = Channels(FileIndex).CA3Indexes; CA1Channels = Channels(FileIndex).CA1Indexes;
        %prep the figure
        figure
        FirstElec = []; SecondElec = []; weights = [];
        if OnlyShowRegions == 0
            for FirElecIndex = 1:59
                for SecElecIndex = FirElecIndex+1:60
                    FirstElec = [FirstElec FirElecIndex];
                    SecondElec = [SecondElec SecElecIndex];
                    %Append the coherence between elec1 and elec2 to the weights vector
                    weights = [weights coh(FirElecIndex,SecElecIndex)];
                end
            end
        elseif OnlyShowRegions == 1
            RegionElecs = [CA3Channels' CA1Channels'];
            for FirElecIndex = 1:59
                for SecElecIndex = FirElecIndex+1:60
                    FirstElec = [FirstElec FirElecIndex];
                    SecondElec = [SecondElec SecElecIndex];
                    %Append the coherence between elec1 and elec2 to the weights vector
                    if ismember(FirElecIndex,RegionElecs) && ismember(SecElecIndex,RegionElecs)
                        weights = [weights coh(FirElecIndex,SecElecIndex)];
                    else 
                        weights = [weights (1/10000000)]; 
                    end
                end
            end
        end
        
        %plot the graph with the coherence strenghts 
        %DOM HELP THIS DOESN"T WORK
        G = graph(FirstElec,SecondElec,weights,ChanLabels);
        P = plot(G,'k');
        %--------------------------------------------------------%
        %adjust the width of the lines representing the coherence
        %first get rid of all the lines that represent little to no coherence
        indx = find(weights<WeightsThreshold);
        P.LineWidth = (weights*WeightMultiplier); %this thickness of the line
        dontplot =  cellstr(repmat('-',length(weights),1))';

        for ii = 1:length(indx)
            dontplot{indx(ii)} = 'none';
        end
        P.LineStyle = dontplot;

        x = zeros(length(P.XData),1);
        y = x;
        MEA_SliceIndexOrder2 = flipud(MEA_SliceIndexOrder);
        %this section will change the location of the nodes to their actual
        %physical one
        for ii = 1:length(P.XData)
            [xx,yy] = find(MEA_SliceIndexOrder2==ii);
            if ii == 15
                [xx,yy] = find(MEA_SliceIndexOrder2==0);
            end
            x(ii) = yy;
            y(ii) = xx;
        P.XData = x;
        P.YData = y;
        end
        hold on;
        
        %Make the markers for CA1 and CA3 visible 
        for ChannelInd = 1:length(CA3Channels)
            Val = SlicePower(FileIndex).CA3(ChannelInd,Phases2Analyze)*GammaVisScale;
            %this will ensure that the markers are of a atleast threshold 1
            %but no larger than threshold 2
            if Val > Thresholds(2)
                Val = Thresholds(2);
            elseif Val < Thresholds(1)
                Val = Thresholds(1);
            end
            u= plot(P.XData(CA3Channels(ChannelInd)),P.YData(CA3Channels(ChannelInd)),'o','MarkerEdgeColor',CA3Color,'MarkerFaceColor',CA3Color,'MarkerSize',Val);
        end
        for ChannelInd = 1:length(CA1Channels)
            Val = SlicePower(FileIndex).CA1(ChannelInd,Phases2Analyze)*GammaVisScale;
            %this will ensure that the markers are of a atleast threshold 1
            %but no larger than threshold 2
            if Val > Thresholds(2)
                Val = Thresholds(2);
            elseif Val < Thresholds(1)
                Val = Thresholds(1);
            end
            u= plot(P.XData(CA1Channels(ChannelInd)),P.YData(CA1Channels(ChannelInd)),'o','MarkerEdgeColor',CA1Color,'MarkerFaceColor',CA1Color,'MarkerSize',Val);
        end
%         for ChannelInd = 1:length(CA3Channels)
%             Val = floor(SlicePower(FileIndex).CA3(ChannelInd,Phases2Analyze)*256/5);
%             if Val > 256
%                 Val = 256;
%             end
%             u= plot(P.XData(CA3Channels(ChannelInd)),P.YData(CA3Channels(ChannelInd)),'o','MarkerEdgeColor',CA3Color,'MarkerFaceColor',inferno(Val,:),'MarkerSize',MarkerSize,'LineWidth',2);
%         end
%         for ChannelInd = 1:length(CA1Channels)
%             Val = floor(SlicePower(FileIndex).CA1(ChannelInd,Phases2Analyze)*256/5);
%             if Val > 256
%                 Val = 256;
%             end
%             u= plot(P.XData(CA1Channels(ChannelInd)),P.YData(CA1Channels(ChannelInd)),'o','MarkerEdgeColor',CA1Color,'MarkerFaceColor',inferno(Val,:),'MarkerSize',MarkerSize,'LineWidth',2);
%         end
        hold off;
        T1 = [Conditions{Phases2Analyze}];
        title(T1)
        if SaveFig == 1
            GammaVisScript = string(GammaVisScale); 
            ToSaveFig = strcat(Fig2SaveDir, 'COHMAPFig_', Channels(FileIndex).SlicesUsed, '_', T1, 'MAX50_VF1','.', DesiredFormat);
            saveas(P,ToSaveFig, DesiredFormat)
        end
        
        %clear workspace for next figure
        clear coh AllFreq AllPower DesiredFreq dontplot TimeWindow2Use x y xx yy ii FirElecIndex SecElecIndex indx u G weights FirstElec SecondElec T1 ToSaveFig CustomXLabels
    end
end
%Clean the workspace 
clear FileIndex Files2Use Phases2Analyze Phase2Analyze factor4Plot GammaVisScale SaveFig DesiredFormat CA3Channels CA1Channels MEA_SliceIndexOrder2 CA1Color CA3Color
%%
%-------------------------------------------------------------------------%
%this will format the data into tables that can be exported easily into
%graphing software such as GraphPad Prism or Excel
%-------------------------------------------------------------------------%
if CoherenceData{1}.Settings.Baseline.Rolling == 1||CoherenceData{1}.Settings.Kainate.Rolling == 1
    disp('WARNING: Rolling Time Window was used.')
    disp('Please use the Rolling Time Window analysis segment further down or change data to not have a rolling time window')
    return
end
ResultsBMean = varfun(@mean,CoherenceTable,'InputVariables','CoherenceBaseline','GroupingVariables',{'Regions', 'Slice'});
ResultsBSTD = varfun(@std,CoherenceTable,'InputVariables','CoherenceBaseline','GroupingVariables',{'Regions', 'Slice'});
ResultsBMeanSTD = ResultsBMean;ResultsBMeanSTD.std_Coherence = ResultsBSTD.std_CoherenceBaseline;

ResultsKMean = varfun(@mean,CoherenceTable,'InputVariables','CoherenceKainate','GroupingVariables',{'Regions', 'Slice'});
ResultsKSTD = varfun(@std,CoherenceTable,'InputVariables','CoherenceKainate','GroupingVariables',{'Regions', 'Slice'});
ResultsKMeanSTD = ResultsKMean;ResultsKMeanSTD.std_Coherence = ResultsKSTD.std_CoherenceKainate;

ResultsMean = varfun(@mean,CoherenceTable,'InputVariables','Coherence','GroupingVariables',{'Regions', 'Slice'});
ResultsSTD = varfun(@std,CoherenceTable,'InputVariables','Coherence','GroupingVariables',{'Regions', 'Slice'});
ResultsMeanSTD = ResultsMean;ResultsMeanSTD.std_Coherence = ResultsSTD.std_Coherence;

LargestCount = max(ResultsMean.GroupCount);
[NumSlices,~] =size(unique(CoherenceTable.Slice));[NumRegions,~] = size(unique(CoherenceTable.Regions));
ResultsIndvCoherence = zeros(LargestCount,(NumSlices*NumRegions));ResultsIndvCoherence=ResultsIndvCoherence/0;
ResultsIndvCoherenceBaseline = ResultsIndvCoherence; ResultsIndvCoherenceKainate = ResultsIndvCoherence;
RegionsUsed = unique(CoherenceTable.Regions); 
SlicesNumbered = unique(CoherenceTable.Slice);
CoherenceTable.Regions = categorical(CoherenceTable.Regions);
for RegionIndx = 1:NumRegions
    for SliceIndex = 1:NumSlices
        Slice1DataTable = CoherenceTable(CoherenceTable.Slice == SlicesNumbered(SliceIndex)& CoherenceTable.Regions == RegionsUsed(RegionIndx), :);
        ColumnIndex = (NumSlices*(RegionIndx-1))+SliceIndex;
        CoherenceVector = Slice1DataTable.Coherence;
        CoherenceBVector = Slice1DataTable.CoherenceBaseline;
        CoherenceKVector = Slice1DataTable.CoherenceKainate;
        [NumDatas,~] = size(CoherenceVector);
        ResultsIndvCoherence(1:NumDatas,ColumnIndex)         = CoherenceVector;
        ResultsIndvCoherenceBaseline(1:NumDatas,ColumnIndex) = CoherenceBVector;
        ResultsIndvCoherenceKainate(1:NumDatas,ColumnIndex)  = CoherenceKVector;
    end
end
T1 = sprintf('\n %s ',RegionsUsed); T1 = ['Region order is :' T1];
disp(T1)
%clean up workspace
clear VariableofInterest NumSlices CoherenceVector CoherenceBVector CoherenceKVector NumDatas ColumnIndex T1 SlicesNumbered RegionsUsed
%%
%-------------------------------------------------------------------------%
%This segment is for plotting the cohernce with respect to electrode
%distance
%-------------------------------------------------------------------------%
ByRegion = 1;
if CoherenceData{1}.Settings.Baseline.Rolling == 1||CoherenceData{1}.Settings.Kainate.Rolling == 1
    disp('WARNING: Rolling Time Window was used.')
    disp('Please use the Rolling Time Window analysis segment further down or change data to not have a rolling time window')
    return
end

if ByRegion == 0
    if UseJenkinsData == 1
        T1 = CoherenceTable(strcmp(CoherenceTable.Model, 'WT'),:);
        T2 = CoherenceTable(strcmp(CoherenceTable.Model, 'MT'),:);
        figure;hold on
        plot((T1.ElecPairDistance)*100,T1.Coherence,'ok')
        plot((T2.ElecPairDistance)*100,T2.Coherence,'or')

%         fit1 = polyfit((T1.ElecPairDistance)*100,T1.Coherence,2);
        
        exponFit = fit((T1.ElecPairDistance)*100,T1.Coherence,'exp1');
%         xIndv = CA1CT.ElecPairDistance*100; yIndv = CA1CT.Coherence;
        exponFit
        plot(exponFit,'k')
        
        exponFit = fit((T2.ElecPairDistance)*100,T2.Coherence,'exp1');
%         xIndv = CA1CT.ElecPairDistance*100; yIndv = CA1CT.Coherence;
        exponFit
        plot(exponFit,'r')
        
%         x = 50:50:1100; yest = polyval(fit1,x);
%         plot(x,yest,'k')
%         fit2 = polyfit((T2.ElecPairDistance)*100,T2.Coherence,2);
%         yest2 = polyval(fit2,x);
%         
%         plot(x,yest2,'r')
        title('Ank Model')
    elseif UseShafferCut == 2
        T1 = CoherenceTable(strcmp(CoherenceTable.Model, 'NoCut'),:);
        T2 = CoherenceTable(strcmp(CoherenceTable.Model, 'Cut'),:);
        figure;hold on
        plot((T1.ElecPairDistance)*100,T1.Coherence,'ok')
        plot((T2.ElecPairDistance)*100,T2.Coherence,'om')

        fit1 = polyfit((T1.ElecPairDistance)*100,T1.Coherence,2);
        x = 50:50:1100; yest = polyval(fit1,x);
        plot(x,yest,'k')
        fit2 = polyfit((T2.ElecPairDistance)*100,T2.Coherence,2);
        yest2 = polyval(fit2,x);
        plot(x,yest2,'m')
        title('Schaffer Lesion Model')
    else
        figure;hold on
        plot((CoherenceTable.ElecPairDistance)*100,CoherenceTable.Coherence,'ok')
        fit1 = polyfit((CoherenceTable.ElecPairDistance)*100,CoherenceTable.Coherence,2);
        x = 50:50:1100; yest = polyval(fit1,x);
        plot(x,yest,'k')
    end
elseif ByRegion == 1
    if UseJenkinsData == 1
       Models = {'WT','MT'};Colors = {'k','r'};
        ColorMods = {[0,0,0];[0,0,0];[0,0,0]};
        figure;hold on
        for modelIndex = 1:2
        %get the data for the model
        TableData = CoherenceTable(strcmp(CoherenceTable.Model, Models{modelIndex}),:);
        %for each model get the region data
        
        [CA1CT, CA3CT,CA13CT] = DataTableGetter(TableData);
        
        x = 50:50:1100; 
        %plot each dataset
        subplot(1,3,1);title('within CA1');hold on;
        plot((CA1CT.ElecPairDistance)*100,CA1CT.Coherence,'o', 'Color',   (Colors{modelIndex}))
%         fitCA1 = polyfit((CA1CT.ElecPairDistance)*100,CA1CT.Coherence,2); yestCA1 = polyval(fitCA1,x);
        exponFit = fit((CA1CT.ElecPairDistance)*100,CA1CT.Coherence,'exp1');
        xIndv = CA1CT.ElecPairDistance*100; yIndv = CA1CT.Coherence;
        exponFit
        plot(exponFit,Colors{modelIndex})
        xlabel('Electrode Distance(µm)'); ylabel('Coherence');ylim([0.5 3]);xlim([0 1100]);hold off
        
        subplot(1,3,2);title('within CA3');hold on;
        plot((CA3CT.ElecPairDistance)*100,CA3CT.Coherence,'s', 'Color', (Colors{modelIndex}))
%         fitCA3 = polyfit((CA3CT.ElecPairDistance)*100,CA3CT.Coherence,2);yestCA3 = polyval(fitCA3,x);
        
        exponFit = fit((CA3CT.ElecPairDistance)*100,CA3CT.Coherence,'exp1');
        xIndv = CA3CT.ElecPairDistance*100; yIndv = CA3CT.Coherence;
        plot(exponFit,Colors{modelIndex})
        xlabel('Electrode Distance(µm)'); ylabel('Coherence');ylim([0.5 3]);xlim([0 1100]);hold off
        exponFit
        subplot(1,3,3);title('between CA1 and CA3');hold on;
        plot((CA13CT.ElecPairDistance)*100,CA13CT.Coherence,'s', 'Color', (Colors{modelIndex}))
        fitCA13 = polyfit((CA13CT.ElecPairDistance)*100,CA13CT.Coherence,2);
        yestCA13 = polyval(fitCA13,x);
        exponFit = fit((CA13CT.ElecPairDistance)*100,CA13CT.Coherence,'exp1');
        xIndv = CA13CT.ElecPairDistance*100; yIndv = CA13CT.Coherence;
        plot(exponFit,Colors{modelIndex})
        xlabel('Electrode Distance(µm)'); ylabel('Coherence');ylim([0.5 3]);xlim([0 1100]);hold off
        exponFit
        end
        sgtitle('Ank Model')
    elseif UseShafferCut == 2
        Models = {'NoCut','Cut'};Colors = {'k','m'};
        ColorMods = {[0,0,0];[0,0,0];[0,0,0]};
        figure;hold on
        for modelIndex = 1:2
        %get the data for the model
        TableData = CoherenceTable(strcmp(CoherenceTable.Model, Models{modelIndex}),:);
        %for each model get the region data
        
        [CA1CT, CA3CT,CA13CT] = DataTableGetter(TableData);
        
        x = 50:50:1100; 
        %plot each dataset
        subplot(1,3,1);title('within CA1');hold on;
        plot((CA1CT.ElecPairDistance)*100,CA1CT.Coherence,'o', 'Color',   (Colors{modelIndex}))
%         fitCA1 = polyfit((CA1CT.ElecPairDistance)*100,CA1CT.Coherence,2); yestCA1 = polyval(fitCA1,x);
        exponFit = fit((CA1CT.ElecPairDistance)*100,CA1CT.Coherence,'exp1');
        xIndv = CA1CT.ElecPairDistance*100; yIndv = CA1CT.Coherence;
        exponFit
        
        plot(exponFit,Colors{modelIndex})
        xlabel('Electrode Distance(µm)'); ylabel('Coherence');ylim([0.5 3]);xlim([0 1100]);hold off
        
        subplot(1,3,2);title('within CA3');hold on;
        plot((CA3CT.ElecPairDistance)*100,CA3CT.Coherence,'s', 'Color', (Colors{modelIndex}))
%         fitCA3 = polyfit((CA3CT.ElecPairDistance)*100,CA3CT.Coherence,2);yestCA3 = polyval(fitCA3,x);
        
        exponFit = fit((CA3CT.ElecPairDistance)*100,CA3CT.Coherence,'exp1');
        xIndv = CA3CT.ElecPairDistance*100; yIndv = CA3CT.Coherence;
        plot(exponFit,Colors{modelIndex})
        xlabel('Electrode Distance(µm)'); ylabel('Coherence');ylim([0.5 3]);xlim([0 1100]);hold off
        
        exponFit
        subplot(1,3,3);title('between CA1 and CA3');hold on;
        plot((CA13CT.ElecPairDistance)*100,CA13CT.Coherence,'s', 'Color', (Colors{modelIndex}))
        fitCA13 = polyfit((CA13CT.ElecPairDistance)*100,CA13CT.Coherence,2);
        yestCA13 = polyval(fitCA13,x);
        exponFit = fit((CA13CT.ElecPairDistance)*100,CA13CT.Coherence,'exp1');
        xIndv = CA13CT.ElecPairDistance*100; yIndv = CA13CT.Coherence;
        plot(exponFit,Colors{modelIndex})
        exponFit
        xlabel('Electrode Distance(µm)'); ylabel('Coherence');ylim([0.5 3]);xlim([0 1100]);hold off
        end
        sgtitle('Schaffer Lesion Model')
        %clean up
        clear fitCA1 fitCA3 fitCA13 yestCA1 yestCA3 yestCA13 CA1CT CA13CT CA3CT
    else
        
        Colors = [0,0,0];
        ColorMods = {[1,0,0];[0,0,1];[1,0,1]};
        figure;hold on
        %get the data for the model
        TableData = CoherenceTable;
        %for each model get the region data
        
        [CA1CT, CA3CT,CA13CT] = DataTableGetter(TableData);
        x = 50:50:1100; 
        %plot each dataset
        subplot(1,3,1);title('within CA1');hold on;
        plot((CA1CT.ElecPairDistance)*100,CA1CT.Coherence,'o', 'Color',   (Colors + ColorMods{1}))
        fitCA1 = polyfit((CA1CT.ElecPairDistance)*100,CA1CT.Coherence,2);
        yestCA1 = polyval(fitCA1,x);
        plot(x,yestCA1,'-','Color',  (Colors  + ColorMods{1}))
        xlabel('Electrode Distance(µm)'); ylabel('Coherence');ylim([0.5 3]);xlim([0 1100]);hold off
        
        subplot(1,3,2);title('within CA3');hold on;
        plot((CA3CT.ElecPairDistance)*100,CA3CT.Coherence,'s', 'Color',   (Colors  + ColorMods{2}))
        fitCA3 = polyfit((CA3CT.ElecPairDistance)*100,CA3CT.Coherence,2);
        yestCA3 = polyval(fitCA3,x);
        plot(x,yestCA3,':','Color',  (Colors  + ColorMods{2}))
        xlabel('Electrode Distance(µm)'); ylabel('Coherence');ylim([0.5 3]);xlim([0 1100]);hold off
        
        subplot(1,3,3);title('between CA1 and CA3');hold on;
        plot((CA13CT.ElecPairDistance)*100,CA13CT.Coherence,'d', 'Color', (Colors  + ColorMods{3}))
        fitCA13 = polyfit((CA13CT.ElecPairDistance)*100,CA13CT.Coherence,2);
        yestCA13 = polyval(fitCA13,x);
        plot(x,yestCA13,'--','Color',(Colors  + ColorMods{3}))
        xlabel('Electrode Distance(µm)'); ylabel('Coherence');ylim([0.5 3]);xlim([0 1100]);hold off
    end
end

clear x fit1 fit2 yest yest2 T1 T1 
%%
%-------------------------------------------------------------------------%
%
%THIS SEGMENT FOR ROLLING TIME
%
%-------------------------------------------------------------------------%



%this segment will make coherece figures showing the baseline then the
%progression trhought the recording. 

%Define these variables
[~,NumFiles] = size(Channels);
Files2Use = 1:NumFiles;
Phase2Analyze = 1:2;
Conditions = {'Baseline','Kainate','Drug'};
DesiredFormat = 'png'; %Set desired format for the automatic figure saving
%-------------------------------------------------------------------------%
for FileIndex = Files2Use
    if CoherenceData{FileIndex}.Settings.Baseline.Rolling == 0||CoherenceData{FileIndex}.Settings.Kainate.Rolling == 0
        disp('WARNING: Rolling Time Window was not used.')
        disp('Please use the Non Rolling Time Window analysis segment further up or change data to have a rolling time window')
        return
    end
    
    clear ChannelsWithData ChanNumWithData coh NewCOH T2 CA1Electrodes CA3Electrodes    
    CA1Electrodes = Channels(FileIndex).CA1Indexes; NumCA1 = length(CA1Electrodes);
    CA3Electrodes = Channels(FileIndex).CA3Indexes; NumCA3 = length(CA3Electrodes);
    totalTimePoints = 0;
    TimeVectorAll = [];
    for Phases2Analyze = Phase2Analyze
        [~,~,timePointer] = size(CoherenceData{FileIndex}.(Conditions{Phases2Analyze}));
        Timer(Phases2Analyze) = timePointer;
        if Phases2Analyze == 1
            TimePairs(Phases2Analyze,1) = 1;
        else
            TimePairs(Phases2Analyze,1) = totalTimePoints+1;
        end
        totalTimePoints = timePointer + totalTimePoints;
        TimePairs(Phases2Analyze,2) = totalTimePoints;
        coh = CoherenceData{FileIndex}.(Conditions{Phases2Analyze});
        coher = permute(coh,[2,1,3]);
        coherence(Phases2Analyze).Phase = coh + coher;
        TimeWindowUsed = CoherenceData{FileIndex}.Settings.(Conditions{Phases2Analyze}).TimeWindowUsed;
        Time2Add = TimeWindowUsed(1):CoherenceData{FileIndex}.Settings.(Conditions{Phases2Analyze}).RollingTime:TimeWindowUsed(2);
        TimeVectorAll = [TimeVectorAll Time2Add];
    end
    %-----------------------------------------------------------------%
    vectorPreallc = zeros(1,totalTimePoints);
    %-----------start of CA3 section---------------%
    if ~isempty(Channels(FileIndex).CA3Channels)
    CoherenceMatrixCA3 = vectorPreallc;
    IndexCounter = 0;
    for ChanCA3Index1 = 1:NumCA3
        for ChanCA3Index2 = (ChanCA3Index1+1):NumCA3
            IndexCounter = IndexCounter + 1;
            Elec1 = CA3Electrodes(ChanCA3Index1);
            Elec2 = CA3Electrodes(ChanCA3Index2);
            Elec1Vector{IndexCounter} = sprintf('S%d_C%d',FileIndex,(ChanCA3Index1));
            Elec2Vector{IndexCounter}= sprintf('S%d_C%d',FileIndex,(ChanCA3Index2));
            ElectrodeNum1Vector{IndexCounter} = Channels(FileIndex).CA3Channels(ChanCA3Index1);
            ElectrodeNum2Vector{IndexCounter} = Channels(FileIndex).CA3Channels(ChanCA3Index2);
            for PhaseIndex = Phase2Analyze
                CoherenceMatrixCA3(IndexCounter,TimePairs(PhaseIndex,1):TimePairs(PhaseIndex,2)) = coherence(PhaseIndex).Phase(Elec1,Elec2,:);
            end
            
        end
    end
    ExportingCoherence(FileIndex).CA3ElecVector1  = Elec1Vector;
    ExportingCoherence(FileIndex).CA3ElecVector2  = Elec2Vector;
    ExportingCoherence(FileIndex).CA3ElectrodeNum1Vector  = ElectrodeNum1Vector;
    ExportingCoherence(FileIndex).CA3ElectrodeNum2Vector  = ElectrodeNum2Vector;
    ExportingCoherence(FileIndex).CA3Coherence  = CoherenceMatrixCA3;
    end
    clear ElecVector1 ElecVector2 ElectrodeNum1Vector ElectrodeNum2Vector
    %-----------end of CA3 section---------------%
    %-----------start of CA1 section---------------%
    if ~isempty(Channels(FileIndex).CA1Channels)
    CoherenceMatrixCA1 = vectorPreallc;
    IndexCounter = 0;
    for ChanCA1Index1 = 1:NumCA1
        for ChanCA1Index2 = (ChanCA1Index1+1):NumCA1
            IndexCounter = IndexCounter + 1;
            Elec1 = CA1Electrodes(ChanCA1Index1);
            Elec2 = CA1Electrodes(ChanCA1Index2);
            for PhaseIndex = Phase2Analyze
                CoherenceMatrixCA1(IndexCounter,TimePairs(PhaseIndex,1):TimePairs(PhaseIndex,2)) = coherence(PhaseIndex).Phase(Elec1,Elec2,:);
            end
            Elec1Vector{IndexCounter} = sprintf('S%d_A%d',FileIndex,(ChanCA1Index1));
            Elec2Vector{IndexCounter}= sprintf('S%d_A%d',FileIndex,(ChanCA1Index2));
            ElectrodeNum1Vector{IndexCounter} = Channels(FileIndex).CA1Channels(ChanCA1Index1);
            ElectrodeNum2Vector{IndexCounter} = Channels(FileIndex).CA1Channels(ChanCA1Index2);
        end
    end 
    ExportingCoherence(FileIndex).CA1ElecVector1  = Elec1Vector;
    ExportingCoherence(FileIndex).CA1ElecVector2  = Elec2Vector;
    ExportingCoherence(FileIndex).CA1ElectrodeNum1Vector  = ElectrodeNum1Vector;
    ExportingCoherence(FileIndex).CA1ElectrodeNum2Vector  = ElectrodeNum2Vector;
    ExportingCoherence(FileIndex).CA1Coherence  = CoherenceMatrixCA1;
    end
    
    clear ElecVector1 ElecVector2 ElectrodeNum1Vector ElectrodeNum2Vector
    
    %-----------end of CA1 section---------------%
    %-----------start of CA1 - CA3 section---------------%
    if ~isempty(Channels(FileIndex).CA3Channels)&&~isempty(Channels(FileIndex).CA1Channels)
    CoherenceMatrixCA13 = vectorPreallc;
    IndexCounter = 0;
    for ChanCA3Index = 1:NumCA3
        for ChanCA1Index = 1:NumCA1
            IndexCounter = IndexCounter + 1;
            Elec1 = CA3Electrodes(ChanCA3Index);
            Elec2 = CA1Electrodes(ChanCA1Index);
            for PhaseIndex = Phase2Analyze
                CoherenceMatrixCA13(IndexCounter,TimePairs(PhaseIndex,1):TimePairs(PhaseIndex,2)) = coherence(PhaseIndex).Phase(Elec1,Elec2,:);
            end
            Elec1Vector{IndexCounter} = sprintf('S%d_C%d',FileIndex,(ChanCA3Index));
            Elec2Vector{IndexCounter}= sprintf('S%d_A%d',FileIndex,(ChanCA1Index));
            ElectrodeNum1Vector{IndexCounter} = Channels(FileIndex).CA3Channels(ChanCA3Index);
            ElectrodeNum2Vector{IndexCounter} = Channels(FileIndex).CA1Channels(ChanCA1Index);
        end
    end
    ExportingCoherence(FileIndex).CA13ElecVector1  = Elec1Vector;
    ExportingCoherence(FileIndex).CA13ElecVector2  = Elec2Vector;
    ExportingCoherence(FileIndex).CA13ElectrodeNum1Vector  = ElectrodeNum1Vector;
    ExportingCoherence(FileIndex).CA13ElectrodeNum2Vector  = ElectrodeNum2Vector;
    ExportingCoherence(FileIndex).CA13Coherence = CoherenceMatrixCA13;
    end
    clear ElecVector1 ElecVector2 ElectrodeNum1Vector ElectrodeNum2Vector IndexCounter Elec1 Elec2
    
    %-----------end of CA1 - CA3 section---------------%
    %Exporting data to larger Structure
    if UseJenkinsData == 1
        if FileIndex >5
            Model= 'MT';
        else
            Model= 'WT';
        end
    elseif UseShafferCut == 2
        if FileIndex > NumAssayDevChan
            Model= 'Cut';
        else
            Model= 'NoCut';
        end
    end
    ExportingCoherence(FileIndex).Model = Model;
    ExportingCoherence(FileIndex).Time = TimeVectorAll;
    %This segment will make the graph desired
    clear CoherenceMatrixCA1 CoherenceMatrixCA3 CoherenceMatrixCA13 Model
end
%clean up the workspace
clear ChannelsWithData ChanNumWithData coh NewCOH T2 CA1Electrodes CA3Electrodes DesiredFormat 
%%
%-------------------------------------------------------------------------%
%this section will make the coherence heatmap figures
%-------------------------------------------------------------------------%
%Define these variables
Files2Use = [1:4];
Phase2Analyze = 1:2;
Conditions = {'Baseline','Kainate','Drug'};
DesiredFormat = 'gif'; %Set desired format for the automatic figure saving
%-------------------------------------------------------------------------%

if UseJenkinsData == 1
    Fig2SaveDir = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\CoherenceData\CoherenceTimeSeries\Jenkins\';
elseif UseShafferCut == 1
    Fig2SaveDir = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\CoherenceData\CoherenceTimeSeries\Shaffer\';
elseif UseJenkinsData == 0 && UseShafferCut == 0
	Fig2SaveDir = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\CoherenceData\CoherenceTimeSeries\AssayDev\';
end
if ~exist(Fig2SaveDir, 'dir')
   mkdir(Fig2SaveDir)
end
    
for FileIndex = Files2Use
    if CoherenceData{FileIndex}.Settings.Baseline.Rolling == 0||CoherenceData{FileIndex}.Settings.Kainate.Rolling == 0
        disp('WARNING: Rolling Time Window was not used.')
        disp('Please use the Non Rolling Time Window analysis segment further up or change data to have a rolling time window')
        return
    end
    
    clear ChannelsWithData ChanNumWithData coh NewCOH T2 CA1Electrodes CA3Electrodes    
    for Phases2Analyze = Phase2Analyze
        %-----------------------------------------------------------------%
        %set up or get the necessary variables
        coh = CoherenceData{FileIndex}.(Conditions{Phases2Analyze});
        CA1Electrodes = Channels(FileIndex).CA1Indexes;
        CA3Electrodes = Channels(FileIndex).CA3Indexes;
        ChannelsWithData = [CA1Electrodes' CA3Electrodes'];
        ChanNumWithData = [Channels(FileIndex).CA1Channels Channels(FileIndex).CA3Channels];
        [~,~,FramesinPhase] = size(CoherenceData{FileIndex}.(Conditions{Phases2Analyze}));
        %-----------------------------------------------------------------%
        for FrameIndexLoop = 1:FramesinPhase
            %get the coherence of desired electrodes
            NewCOH = coh(ChannelsWithData,ChannelsWithData,FrameIndexLoop);
            NewCOH = NewCOH + NewCOH';
            %----Make the coherence heatmap----%
            f1 = figure(1000);
            hm = heatmap(NewCOH);
            %load the heatmap color blind friendly color scheme
            inferno=inferno();colormap inferno;
            caxis([0 1]) %set the c axis limits
            T2 = [Conditions{Phases2Analyze}];
            title(T2);
            drawnow 
            % Capture the plot as an image 
            frame = getframe(f1); 
            im = frame2im(frame); 
            [imind,cm] = rgb2ind(im,256); 
            %If desired this section will save the figures automatically in
            %desired format
            ToSaveFig = [Fig2SaveDir 'HighResCOHMovie_' Channels(FileIndex).SlicesUsed '.' DesiredFormat];
            if strcmp (DesiredFormat,'gif')
                if FrameIndex == 1
                    imwrite(imind,cm,ToSaveFig,DesiredFormat, 'Loopcount',inf); 
                else
                    imwrite(imind,cm,ToSaveFig,DesiredFormat,'WriteMode','append'); 
                end
            elseif strcmp (DesiredFormat,'avi')
                if FrameIndex == 1
                    v = VideoWriter(ToSaveFig);
                    open(v);
                    writeVideo(v,frame);
                else 
                    writeVideo(v,frame);
                end
            end
            FrameIndex = FrameIndex + 1;
        end
    end
end
%clean up the workspace
disp('Finished Making Movies')
clear ChannelsWithData ChanNumWithData coh NewCOH T2 CA1Electrodes CA3Electrodes DesiredFormat 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %functions%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LayoutMEAInfo = LoadLoadout(DesiredLayout)
%This function will either load the slice order using the actual numbers of
%the electrode layout or use a matrix with the matrix that has the Indexes
%in order

MEA_SliceOrder = ...
    [11,21,31,41,51,61,71,81,91,101,...
    12,22,32,42,52,62,72,82,92,102,...
    13,23,33,43,53,63,73,83,93,103,...
    0, 24,34,44,54,64,74,84,94,104,...
    15,25,35,45,55,65,75,85,95,105,...
    16,26,36,46,56,66,76,86,96,106];

MEA_SliceIndexOrder = ...
    [20,21,24,27,29,32,34,37,40,41;...
    18,19,23,26,28,33,35,38,42,43;...
    16,17,22,25,30,31,36,39,44,45;...
    0,14,9,6,1,60,55,52,47,46;...
    13,12,8,5,3,58,56,53,49,48;...
    11,10,7,4,2,59,57,54,51,50];

if DesiredLayout == 1
    LayoutMEAInfo = MEA_SliceOrder;
elseif DesiredLayout == 2
    LayoutMEAInfo = MEA_SliceIndexOrder;
end

end %end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Model1, Model2,D1] = CoherenceMixedModelF(DataTable,UseJenkinsData,UseShafferCut)
%This function runs the linear mixed model using the given data

%convert categorical data into categorical vaariable type
if UseJenkinsData == 1
    %This seection will turn the model category into a specified number
    %that way the control is the 0 in the model
    %revise this section
    DataTable.Model = categorical(DataTable.Model);
    DataTable.Model = renamecats(DataTable.Model,'WT','0');
    DataTable.Model = renamecats(DataTable.Model,'MT','1');
    DataTable.Model = str2num(char(DataTable.Model));
    DataTable.Model = categorical(DataTable.Model);
elseif UseShafferCut == 2
    DataTable.Model = categorical(DataTable.Model);
    DataTable.Model = renamecats(DataTable.Model,'NoCut','0');
    DataTable.Model = renamecats(DataTable.Model,'Cut','1');
    DataTable.Model = str2num(char(DataTable.Model));
    DataTable.Model = categorical(DataTable.Model);
end
DataTable.Regions   = categorical(DataTable.Regions);
DataTable.Regions   = reordercats(DataTable.Regions,{'CA3','CA1','CA1_3'});
DataTable.Slice     = nominal(DataTable.Slice);
DataTable.Elec1     = categorical(DataTable.Elec1);
DataTable.Elec2     = categorical(DataTable.Elec2);
% DataTable.Coherence = table2array(DataTable.Coherence);
if UseJenkinsData == 1 || UseShafferCut == 2
    Model1 = fitlme(DataTable,'Coherence~Model*Regions+(1|Slice)');
    Model2 = fitlme(DataTable,'Coherence~Model*Regions+(1|Slice)+(1|Slice:Elec1)');
else
    Model1 = fitlme(DataTable,'Coherence~Regions+(1|Slice)');
    Model2 = fitlme(DataTable,'Coherence~Regions+(1|Slice)+(1|Slice:Elec1)');
end
D1 = DataTable;
end% end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CA1Data, CA3Data,CA13Data] = DataTableGetter(CohTable)

%get the data for each specific region and create new subtables with it
CA1Data  = CohTable(strcmp(CohTable.Regions, 'CA1'),:);
CA3Data  = CohTable(strcmp(CohTable.Regions, 'CA3'),:);
CA13Data = CohTable(strcmp(CohTable.Regions, 'CA1_3'),:);
end% end of function