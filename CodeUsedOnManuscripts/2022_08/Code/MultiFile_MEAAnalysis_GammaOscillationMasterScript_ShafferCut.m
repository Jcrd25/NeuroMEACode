%MEA Analysis Master Script
%Import Data from HDF5 File
%Select HDF5 File and create a struct that contains the files to be
%analyzed
%You may select mulitple files but then you need to adjust the following
%script segments to analyze the desired channels at a time.

%run this first
NumFilesToStudy = 2;

for ii = 1:NumFilesToStudy
    FileG(ii).HDF5DataFile = uipickfiles;  
    T1 = sprintf('Selected File Set %d',ii);
    disp(T1);
end
%%
%%-----------------------------------------------------------------------%%
%Run this after running the top section to generate the raw filtered data
%files
for FG = 1:NumFilesToStudy
    
HDF5DataFile = FileG(FG).HDF5DataFile;
    
SavePath = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\ShafferCut\FilteredRawData\';
if ~exist(SavePath, 'dir')
       mkdir(SavePath)
end
Fpass = 100;
Fstop = 120;
Apass = 0.5;
Astop = 65;
Fs = 20000;

d = designfilt('lowpassiir', ...
  'PassbandFrequency',Fpass,'StopbandFrequency',Fstop, ...
  'PassbandRipple',Apass,'StopbandAttenuation',Astop, ...
  'DesignMethod','butter','SampleRate',Fs);
clear Fpass Fstop Apass Astop

[RecordingLength, Settings] = M_HDF5LengthRecordingCalculator(HDF5DataFile);
Traces.Data = struct; Traces.Fs = struct;
for ChannelNumber = 1:60
    %Preallocate and set up variables
    ChannelTrace = zeros(RecordingLength,1); Fs = 20000; DesiredFs = 1000;
    MinTime = 1; MaxTime = 0;
    %Load the complete trace
    for FileIndex = 1:length(HDF5DataFile)
        [ChannelVoltageTrace] = HDF5_ChannelDataExtractor(HDF5DataFile{FileIndex},ChannelNumber);
        [NumTime,~]=size(ChannelVoltageTrace);
        MaxTime = MaxTime+NumTime;
        ChannelTrace(MinTime:MaxTime) = ChannelVoltageTrace;
        MinTime = MaxTime+1;
    end
    %filter and downsample the data
    ChannelVoltageTrace = filtfilt(d,ChannelVoltageTrace);
    [ChannelTrace,NewFs] = M_Downsampler(ChannelTrace, Fs, DesiredFs);
    Traces(ChannelNumber).Data = ChannelTrace;
    Traces(ChannelNumber).Fs   = NewFs;
end
MEAData.RawData  = Traces;
MEAData.Settings = Settings;
MEAData.Filter = d;
clear Traces

%Save the MEAData Struct
[NAME] = M_MATLABsaver(Settings{1}.FileName);
NAME = char(NAME);


Name = strcat(SavePath, NAME);
Name = char(Name);
save(Name, 'MEAData','-v7.3'); 
T1 = ['Saved' ' ' NAME];
disp(T1)
end
%-------------------------------------------------------------------------%
%%
%Calculate the spectograms and periodograms
SavedPath   = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\ShafferCut\FilteredRawData\';
SavePathDir = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\Jenkins\MultiTaperSpectrumData\';
listing = dir(SavedPath);
counter = 0;
tic
for ii = 1:length(listing)
    if listing(ii).bytes > 1000
        counter = counter+1;
        MEAFileNames{counter} = [listing(ii).name];
    end
end

for MEAFileIndx = 1:counter
    FileFullName = [SavedPath MEAFileNames{MEAFileIndx}];
    load(FileFullName)
    %Define the parameters to be used
    %From Chronux Toolbox: 
    %params: structure with fields tapers, pad, Fs, fpass, err, trialave
    params.tapers   =[5 9];
    params.pad      = 0;
    params.Fs       = MEAData.RawData(1).Fs;
    params.fpass    = [0 100];%frequency to analyze
    params.err      = 0; %
    params.trialave = 0;
    movingwin       = [10 .5];
for ChannelID = 1:60
    
[S,t,f]    = mtspecgramc(MEAData.RawData(ChannelID).Data,movingwin,params);
MEAGamma(ChannelID).Frequencies = f;
MEAGamma(ChannelID).Time        = t;
MEAGamma(ChannelID).PowerSpec   = S;
%save the parameters used to do the calculations with the data
MEAGamma(ChannelID).params           = params; 
MEAGamma(ChannelID).params.movingwin = movingwin;
clear S t f 
MEAData.RawData(ChannelID).Data =[];
end
%Save the MEAData Struct
Name = strcat(SavePathDir,'MTS_' ,MEAFileNames{MEAFileIndx});
Name = char(Name);
save(Name, 'MEAGamma','-v7.3'); 
T1 = ['Saved' ' ' MEAFileNames{MEAFileIndx} 'in_' ];
disp(T1);toc

end
disp('Finished Multitaper analysis in_')
toc
%%
%Calculate the spectograms and periodograms adn save the data in channel
%form
%---------------------------Run this second-------------------------------%
SavedPath   = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\ShafferCut\FilteredRawData\';
Dir2Save = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\ShafferCut\MultiTaperSpectrumChannelData\';
listing = dir(SavedPath);
counter = 0;
tic
for ii = 1:length(listing)
    if listing(ii).bytes > 1000
        counter = counter+1;
        MEAFileNames{counter} = [listing(ii).name];
    end
end

for MEAFileIndx = 1:counter
    FileFullName = [SavedPath MEAFileNames{MEAFileIndx}];
    load(FileFullName)
    %Define the parameters to be used
    %From Chronux Toolbox: 
    %params: structure with fields tapers, pad, Fs, fpass, err, trialave
    params.tapers   =[5 9];
    params.pad      = 0;
    params.Fs       = MEAData.RawData(1).Fs;
    params.fpass    = [0 100];%frequency to analyze
    params.err      = 0; %
    params.trialave = 0;
    movingwin       = [15 .5];
    
    MEAGamma = struct;
    
for ChannelID = 1:60
    
[S,t,f]    = mtspecgramc(MEAData.RawData(ChannelID).Data,movingwin,params);

MEAGamma.Frequencies = f;
MEAGamma.Time        = t;
MEAGamma.PowerSpec   = S;

MEAGamma.Params      = params; MEAGamma.Params.movingwin = movingwin;
clear S t f
%save the data file for each channel
%make a series of folders to for the channel
ChanName = sprintf('MEAGamma_Chan%d',ChannelID);
SliceName = MEAFileNames{MEAFileIndx}(9:(end-4));
SubFile = [Dir2Save SliceName '\' ];
if ~exist(SubFile, 'dir')
    mkdir(SubFile)
end
Name = [SubFile ChanName '.mat'];
save(Name, 'MEAGamma','-v7.3'); 
MEAData.RawData(ChannelID).Data =[];
end
clear MEAGamma Name SubFile SliceName ChanName
T1 = ['Saved' ' ' MEAFileNames{MEAFileIndx} ' in_' ];
disp(T1);toc

end
disp('Finished MultitaperSpectrum analysis in_')
toc

%%
%get the data using wavelets for the spectograms
%Calculate the spectograms and periodograms
SavedPath   = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\FilteredRawData\';
Dir2Save = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\WaveletsSpectrumData\';
listing = dir(SavedPath);
counter = 0;
tic
for ii = 1:length(listing)
    if listing(ii).bytes > 1000
        counter = counter+1;
        MEAFileNames{counter} = [listing(ii).name];
    end
end

for MEAFileIndx = 1:counter
    FileFullName = [SavedPath MEAFileNames{MEAFileIndx}];
    load(FileFullName)
    %Define the parameters to be used
    freqVec = [1:1:100];  width = 7; 
for ChannelID = 1:60
    
Fs = MEAData.RawData(ChannelID).Fs;
[TFR,timeVec,freqVec] = traces2TFR(MEAData.RawData(ChannelID).Data,freqVec,Fs,width);
MEAGamma(ChannelID).Frequencies = freqVec;
MEAGamma(ChannelID).Time        = timeVec;
MEAGamma(ChannelID).PowerSpec   = TFR;

MEAGamma(ChannelID).params.width = width;
clear S t f 
MEAData.RawData(ChannelID).Data =[];
end
%Save the MEAData Struct
Name = strcat(Dir2Save,'WaveW7_' ,MEAFileNames{MEAFileIndx});
Name = char(Name);
save(Name, 'MEAGamma','-v7.3'); 
clear MEAGamma Name
T1 = ['Saved' ' ' MEAFileNames{MEAFileIndx} 'in_' ];
disp(T1);toc

end
disp('Finished Wavelets analysis in_')
toc

%%
%-------------------------------------------------------------------------%
%---------------------------Run this third--------------------------------%
ChanDataPathDir = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\ShafferCut\MultiTaperSpectrumChannelData\';
Save2Directory = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\ShafferCut\';
listing = dir(ChanDataPathDir); listing = listing(3:end-1);
UsedB = false;
for SliceIndex = 1:length(listing)
    MEA_PowerData = struct; MEA_RelativeData = struct;
    for chanID = 1:60
        ChannelName = sprintf('MEAGamma_Chan%d', chanID);
        ChannelFile = [ChanDataPathDir listing(SliceIndex).name '\' ChannelName];
        load(ChannelFile)
        No60Hz = 1;
        if UsedB == 1
            Powers = 10*log10(MEAGamma.PowerSpec');
        else
            Powers = MEAGamma.PowerSpec';
        end
        [PowerData, Relative] = PowerPCalc(Powers,MEAGamma.Frequencies, No60Hz);
        MEA_PowerData(chanID).PowerData = PowerData;
        MEA_RelativeData(chanID).RelativePower = Relative;
        
        MEA_PowerData(chanID).Time    = MEAGamma.Time;
        MEA_RelativeData(chanID).Time = MEAGamma.Time;
        clear ChannelName ChannelFile
    end
    clear MEAData
    %save the data
    if UsedB == 1
        Save2Dir = [Save2Directory 'dBBandPower\'];
        NamePD = strcat(Save2Dir,'RAW_dBPower\');
        NameRD = strcat(Save2Dir,'Relative_dBPower\');
        NamePD = char(NamePD); NameRD = char(NameRD);
    else
        Save2Dir = [Save2Directory 'BandPower\'];
        NamePD = strcat(Save2Dir,'RAW_Power\');
        NameRD = strcat(Save2Dir,'Relative_Power\');
        NamePD = char(NamePD); NameRD = char(NameRD);
    end
    
    
    if ~exist(NamePD, 'dir')
       mkdir(NamePD)
    end
    if ~exist(NameRD, 'dir')
        mkdir(NameRD)
    end
    NamePD = [NamePD 'RawPower_' listing(SliceIndex).name];
    NameRD = [NameRD 'RelativePower_' listing(SliceIndex).name];
    save(NamePD, 'MEA_PowerData','-v7.3'); 
    save(NameRD, 'MEA_RelativeData','-v7.3'); 
    T1 = ['Saved Power Band Data for ' listing(SliceIndex).name];
    disp(T1);
    clear MEA_PowerData MEA_RelativeData NamePD NameRD
end
disp('Finished Band Analysis')
%%
%-------------------------------------------------------------------------%
%---------------------------Run this Fourth-------------------------------%
%This section is to Prepare the data to calculate the average Peak
%Frequency and the Q values
SavedPath   = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\ShafferCut\FilteredRawData\';
Save2DirSum = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\ShafferCut\PowerSummaryData\';
Save2DirAll = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\ShafferCut\AllKainatePowerSummaryData\';

KainateEnd = 90; KainateStart = 60; %set in minutes
DrugDuration = 10;
Phases2Analyze = 1:2;
TimeWindowAnalysis = 5;
DoAll = false;
%-------------------------%

Times4All = [KainateStart, KainateEnd;KainateEnd, (KainateEnd+DrugDuration)];
Time4Focused = [(KainateStart-TimeWindowAnalysis),KainateStart;(KainateEnd-TimeWindowAnalysis), KainateEnd;...
(KainateEnd+DrugDuration-TimeWindowAnalysis),(KainateEnd+DrugDuration)];
Times4All = Times4All(Phases2Analyze);
Time4Focused = Time4Focused(Phases2Analyze,:);
if ~exist(Save2DirSum, 'dir')
    mkdir(Save2DirSum)
end
if DoAll == 1
    if ~exist(Save2DirAll, 'dir')
        mkdir(Save2DirAll)
    end
end
listing = dir(SavedPath);
counter = 0; NumSegs = size(Time4Focused);NumSegs = NumSegs(1);
Num4AllSegs = size(Times4All);Num4AllSegs = Num4AllSegs(1);
Times4All = Times4All*60;Time4Focused=Time4Focused*60;%convert 2 seconds
%Define the parameters to be used
%From Chronux Toolbox: 
%params: structure with fields tapers, pad, Fs, fpass, err, trialave
params.tapers   =[5 9];
params.pad      = 0;
params.fpass    = [0 100];%frequency to analyze
params.err      = 0; %
params.trialave = 0;

tic
%load the data
for ii = 1:length(listing)
    if listing(ii).bytes > 1000
        counter = counter+1;
        MEAFileNames{counter} = [listing(ii).name];
    end
end

for FileNum = 1:counter
    File2LoadName = [SavedPath MEAFileNames{FileNum}];
    load(File2LoadName)
    Channel = struct;
    for ElectIndex = 1:60
        %calcualte the power spectrum for each data segment of interest
        for DataSegment = 1:NumSegs
            params.Fs      = MEAData.RawData(1).Fs;
            NumDatPoints   = length(MEAData.RawData(ElectIndex).Data);
            TimeVect       = 1/params.Fs:1/params.Fs:NumDatPoints/params.Fs;
            DesiredIndexes = find(TimeVect > Time4Focused(DataSegment,1) & TimeVect < Time4Focused(DataSegment,2));
            clear TimeVect NumDatPoints;
            data2Analyze   = MEAData.RawData(ElectIndex).Data(DesiredIndexes);
            [PowData,FreqData] = mtspectrumc(data2Analyze,params);
            %save the data in the appropriate segment
            if DataSegment == 1 
                Channel(ElectIndex).Baseline.Power      = PowData;
                Channel(ElectIndex).Baseline.FreqData   = FreqData;
%                 Channel(ChannelIndex).Baseline.PowerError = PowerError;
            elseif DataSegment == 2
                Channel(ElectIndex).Kainate.Power      = PowData;
                Channel(ElectIndex).Kainate.FreqData   = FreqData;
%                 Channel(ChannelIndex).Kainate.PowerError = PowerError;
            elseif DataSegment == 3
                Channel(ElectIndex).KainateDrug.Power      = PowData;
                Channel(ElectIndex).KainateDrug.FreqData   = FreqData;
%                 Channel(ChannelIndex).KainateDrug.PowerError = PowerError;
            end
            clear data2Analyze DesiredIndexes Power FreqData PowerError     
        end
    end
    MEAPowerSum.Data = Channel;
    clear Channel;
    SliceName = MEAFileNames{FileNum}(9:(end));
    Name4File = [Save2DirSum 'PowSumD_' SliceName];
    save(Name4File, 'MEAPowerSum','-v7.3'); 
    clear Name4File
    clear MEAPowerSum
    T1 = ['Saved Power Band Windowed sumary Data for ' MEAFileNames{FileNum}];
    disp(T1)
    if DoAll == 1
    %calcualte the average for the entire experimental phase
    for ElectIndex = 1:60
        %calcualte the power spectrum for each period of the experiment
        for DataSegment = 1:Num4AllSegs
            NumDatPoints   = length(MEAData.RawData(ElectIndex).Data);
            params.Fs      = MEAData.RawData(1).Fs;
            TimeVect       = 1/params.Fs:1/params.Fs:NumDatPoints/params.Fs;
            DesiredIndexes = find(TimeVect > Times4All(DataSegment,1) & TimeVect < Times4All(DataSegment,2));
            data2Analyze   = MEAData.RawData(ElectIndex).Data(DesiredIndexes);
            [PowData,FreqData]=mtspectrumc(data2Analyze,params);
            if DataSegment == 1 
                Channel(ElectIndex).Kainate.Power      = PowData;
                Channel(ElectIndex).Kainate.FreqData   = FreqData;
%                 Channel(ChannelIndex).Kainate.PowerError = PowerError;
            elseif DataSegment == 2
                Channel(ElectIndex).KainateDrug.Power      = PowData;
                Channel(ElectIndex).KainateDrug.FreqData   = FreqData;
%                 Channel(ChannelIndex).KainateDrug.PowerError = PowerError;
            end
            clear data2Analyze DesiredIndexes Power FreqData PowerError TimeVect NumDatPoints    
        end
    end
    MEA_AllPower.Data = Channel;
    clear Channel;
    SliceName = MEAFileNames{FileNum}(9:(end));
    Name4File = [Save2DirAll 'PowAllD_' SliceName];
    save(Name4File, 'MEA_AllPower','-v7.3'); 
    T1 = ['Saved Power Band sumary Data for ' SliceName];
    clear MEA_AllPower
    disp(T1);
    end
end
%%
%-------------------------------------------------------------------------%
%---------------------------Run this Fifth--------------------------------%
%-----------------------%
%Define these parameters
SavedPath   = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\ShafferCut\FilteredRawData\';
% ChanDataPathDir = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\MK801_SzModel\MultiTaperSpectrumChannelData\';
Save2Directory = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\ShafferCut\Periodograms\';
%Define these variables
Phase2Analyze = [1:2];
Conditions = {'Baseline','Kainate','Bicuculline'};
ConditionsTimes = [55,60;85,90;93,98];
SaveFig = 1; %Set true or false SPECIAL CASE 2 = Wont save figure but will display it with a minor pause
DesiredFormat = 'png'; %Set desired format for the automatic figure saving
%-----------------------%
listing = dir(SavedPath); listing = listing(3:end-1);
UsedB = true;
%Define the parameters to be used
%From Chronux Toolbox: 
%params: structure with fields tapers, pad, Fs, fpass, err, trialave
params.tapers   =[5 9];
params.pad      = 0;
params.fpass    = [0 100];%frequency to analyze
params.err      = 0; %
params.trialave = 0;
movingwin       = [15 .5];

for SliceIndex = 1:length(listing)
% for SliceIndex = 5
    %load the data
    File2Load = [SavedPath listing(SliceIndex).name];
    load(File2Load)
    params.Fs       = MEAData.RawData(1).Fs;
    SliceChannelInfo = MEAData.Settings{1,1}.InfoChannel.Label;
    for ChannelID = 1:60
        ChannelName = sprintf('MEAGamma_Chan%d', ChannelID);
        TimesOfInterest = ConditionsTimes*60;%change it to seconds
        f1 = figure;
        hold on
        %plot the periodogram for each stage
        for StageID = Phase2Analyze
            Fs = MEAData.RawData(1).Fs;
            TimeL = TimesOfInterest(StageID, 1)*Fs;TimeR = TimesOfInterest(StageID,2)*Fs;
            [S,f]=mtspectrumc(MEAData.RawData(ChannelID).Data(TimeL:TimeR),params);
            SS = smoothdata(S);
            plot(f,10*log10(SS))
            %add the axis limits and labels
            ylabel('Power (dB)')
            xlabel('Frequency (Hz)')
            xlim([0 100])
            %Title of the figure
            ChanNumber = str2double(SliceChannelInfo{ChannelID,1});
            T1 = sprintf('Channel %d',ChanNumber);
            title(T1)
            % ylim([1 100])
        end
        %Add legend
        legend(Conditions{Phase2Analyze})
        hold off
        if SaveFig == 1
            %Save the periodogram for this dataset
            ChanName = sprintf('Chan%d',ChannelID);
            SliceName = listing(SliceIndex).name(1:end-4);
            SubFile = [Save2Directory SliceName '\' ];
            if ~exist(SubFile, 'dir')
                mkdir(SubFile)
            end
            Sub2File = [Save2Directory SliceName '\MatFiles\' ];
            if ~exist(Sub2File, 'dir')
                mkdir(Sub2File)
            end
            Name = [SubFile ChanName 'NoBic.' DesiredFormat];
            saveas(f1,Name,DesiredFormat)
            Name2 = [Sub2File ChanName '.fig'];
            saveas(f1,Name2,'fig')
        elseif SaveFig ==2
            pause(2)
        end
    %clear workspace
    clear Name Name2 ChanName ChannelID ChannelName ChanNumber f f1 Fs S T1 SliceName SubFile TimeL TimeR TimesOfInterest
    close all
    end
    T1 = ['Saved Periodograms for ' listing(SliceIndex).name];
    disp(T1);
    clear MEAData SliceChannelInfo ChannelID
end
disp('Finished saving All Periodograms')
%%
%Run this segment to generate the spectograms utilizing the wavelets

%Wavelet section for the generation of spectograms
SavedPath   = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\AssayDev\FilteredRawData\';
% ChanDataPathDir = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\MK801_SzModel\MultiTaperSpectrumChannelData\';
Save2Directory = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\AssayDev\Spectograms\';
%Define these variables
plotIt =0;MovieIt =1;
Time0 = 60;
TimeWindowFromTime0 = [-5 29]; %in minutes
TimeWindow4Movie = 1;% this is the time step for each frame
WindowSize= 10; %this is the window the data will be displayed
SaveFig = 1; %Set true or false SPECIAL CASE 2 = Wont save figure but will display it with a minor pause
DesiredFormat = 'png'; %Set desired format for the automatic figure saving
%-----------------------%
listing = dir(SavedPath); listing = listing(3:end-1);
%50 is a good channel;
ChansInterest = [22];
for SliceIndex = 10
    %load the data
    File2Load = [SavedPath listing(SliceIndex).name];
    load(File2Load)
    params.Fs       = MEAData.RawData(1).Fs;
    SliceChannelInfo = MEAData.Settings{1,1}.InfoChannel.Label;
    for ChannelID = ChansInterest
        ChannelName = sprintf('MEAGamma_Chan%d', ChannelID);
        %plot the periodogram for each stage
        Fs = MEAData.RawData(1).Fs;
        TimeL = (Time0 + TimeWindowFromTime0(1)-1)*Fs*60;
        TimeR = (Time0 + TimeWindowFromTime0(2))*Fs*60;
        Trace = MEAData.RawData(ChannelID).Data;
        TraceTrim = Trace(TimeL:TimeR);
        freqVec = [1:0.5:100];  width = 25;
        [TFR,timeVec,freqVec] = traces2TFR(Trace,freqVec,Fs,width);
        if plotIt == 1
            
            f1 = figure;
            hold on
            % imagesc(timeVec, freqVec, TFR);
            imagesc(((timeVec./60)-Time0), freqVec, 10*log10(TFR));
            ylabel('Frequency (Hz)')
            xlabel('Time (min)')
            c = colorbar;
            % xlim([0 120])
            ylim([1 100])
            set(gca,'YDir','normal')
            caxis([-5 15]);
            xlim([TimeWindowFromTime0])
            ylabel(c, 'Power (dB)')
            colormap inferno
            ChanNumber = str2double(SliceChannelInfo{ChannelID});
            T1 = sprintf('Channel %d',ChanNumber);
            title(T1)
            T1 = sprintf('Channel%d',ChanNumber);
            ExpFileName = [Save2Directory T1];
            if SaveFig == 1
        %Save the periodogram for this dataset
                ChanName = sprintf('Chan%d',ChannelID);
                SliceName = listing(SliceIndex).name(1:end-4);
                SubFile = [Save2Directory SliceName '\' ];
                if ~exist(SubFile, 'dir')
                    mkdir(SubFile)
                end
                Sub2File = [Save2Directory SliceName '\MatFiles\' ];
                if ~exist(Sub2File, 'dir')
                    mkdir(Sub2File)
                end
                Name = [SubFile 'Specto_' ChanName '.' DesiredFormat];
                saveas(f1,Name,DesiredFormat)
                Name2 = [Sub2File 'Spectogram_' ChanName '.fig'];
                %saveas(f1,Name2,'fig')
            elseif SaveFig ==2
                pause(2)
            end
            
            clear Trace TFR timeVec freqVec f
            close all
       elseif MovieIt == 1 
            f1 = figure;
            hold on
            % imagesc(timeVec, freqVec, TFR);
            imagesc(((timeVec)-(Time0*60)), freqVec, 10*log10(TFR));
            ylabel('Frequency (Hz)')
            xlabel('Time (s)')
            c = colorbar;

            ylim([1 100])
            set(gca,'YDir','normal')
            caxis([-5 15]);
            ylabel(c, 'Power (dB)')
            colormap inferno
            ChanNumberLabel = str2double(SliceChannelInfo{ChannelID});
             
            T1 = sprintf('Channel %d',ChanNumberLabel);
            title(T1)
            ChanNumber = sprintf('%d',ChannelID);
            
            T1 = ['Channel' ChanNumber];
            ExpFileName = [Save2Directory T1];
            ChanName = sprintf('Chan%d',ChannelID);
            SliceName = listing(SliceIndex).name(1:end-4);
            SubFile = [Save2Directory SliceName '\' ];
            TimeBins = (TimeWindowFromTime0(1)*60):TimeWindow4Movie:(TimeWindowFromTime0(2)*60); %in seconds
            NumFrames = length(TimeBins);
            Dirname = [SubFile 'SpectoMovie\'];
            if ~exist(Dirname, 'dir')
                   mkdir(Dirname)
            end
            Windower = floor(WindowSize/TimeWindow4Movie);
            Name = [Dirname 'SpectogramMovie_Chan' ChanNumber '.gif'];
            for framIndex = Windower:NumFrames
                TimeL = TimeBins(framIndex-Windower+1);
                TimeR = TimeBins(framIndex);
                xlim([TimeL TimeR])

                  % Capture the plot as an image 
                  frame = getframe(f1); 
                  im = frame2im(frame); 
                  [imind,cm] = rgb2ind(im,256); 
                  % Write to the GIF File 
                  if framIndex == Windower 
                      imwrite(imind,cm,Name,'gif','Loopcount',inf); 
                  else 
                      imwrite(imind,cm,Name,'gif','DelayTime', 0.1, 'WriteMode','append'); 
                  end 
            end
          close all
          clear TFR timeVec freVec Fs TimeL TimeR TraceTrim Trace 
        end
    end
end
%%
%-------------------------------------------------------------------------%
function [ChannelVoltageTrace] = HDF5_ChannelDataExtractor(HDF5DataFile,ChannelNumber)
%This function is meant to be used to import the MEA data from HDF5 files
%exported through MultiChannel DataManager from multichannel systems. It
%will create structure in MatLab. 
%The AnalogStream is important since it might be a filtered data set or the 
%raw trace.
%This updated version is meant to reduce the size of the traces analyzed.
%It is used to only look at downsampled data instead of the whole data
%file.
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN: April 12, 2017
%LAST MODIFIED: August 17, 2020
%v1.0

%Preallocate for ChannelData
Temp = struct;
%Import the data into each struct
%Data are in int 32, to analyze them much change the to double by using
%double(data)
%The ChannelData is no longer being stored in the structure in order to
%save memory.
%define the objects as RawData and Settings


Temp.ChannelData        = h5read(HDF5DataFile, '/Data/Recording_0/AnalogStream/Stream_0/ChannelData');
MEAData.Settings.InfoChannel    = h5read(HDF5DataFile, '/Data/Recording_0/AnalogStream/Stream_0/InfoChannel');
MEAData.Settings.DataType       = h5readatt(HDF5DataFile, '/Data/Recording_0/AnalogStream/Stream_0','Label');
MEAData.Settings.Date           = h5readatt(HDF5DataFile, '/Data','Date');
MEAData.Settings.MeaName        = h5readatt(HDF5DataFile, '/Data','MeaName');
%Get the file name
Temporary_Variable              = h5info(HDF5DataFile);
MEAData.Settings.FileName       = Temporary_Variable.Filename;
%Make sure the conversion factor is the same for all channels
for iii = 1:size(MEAData.Settings.InfoChannel.ConversionFactor, 1)
    if MEAData.Settings.InfoChannel.ConversionFactor(1) ~= MEAData.Settings.InfoChannel.ConversionFactor(iii)
        fprintf('WARNING: The conversion factors are not the same for all channels please review InfoChannel \n')
    end
end
%Convert data into a Voltage trace
%the conversion factor is ~59605. Included the change into µV. This info was obtained from the InfoChannel ConversionFactor.
VoltageTrace = (double(Temp.ChannelData)*double(MEAData.Settings.InfoChannel.ConversionFactor(1)))/(10^6); 

ChannelVoltageTrace = VoltageTrace(:,ChannelNumber);
end
%-------------------------------------------------------------------------%
function [RecordingLength, Settings] = M_HDF5LengthRecordingCalculator(HDF5DataFile)
RecordingLength = 0;

for ii = 1:length(HDF5DataFile)
    InfoStruct = h5info(HDF5DataFile{ii}, '/Data/Recording_0/AnalogStream/Stream_0/ChannelData');
    if ii ==1
        Settings{1}.InfoChannel = h5read(HDF5DataFile{ii}, '/Data/Recording_0/AnalogStream/Stream_0/InfoChannel');
        Settings{1}.DataType       = h5readatt(HDF5DataFile{ii}, '/Data/Recording_0/AnalogStream/Stream_0','Label');
        Settings{1}.Date           = h5readatt(HDF5DataFile{ii}, '/Data','Date');
        Settings{1}.MeaName        = h5readatt(HDF5DataFile{ii}, '/Data','MeaName');
    end
    [FileLength] = InfoStruct.Dataspace.Size(1);
    RecordingLength = RecordingLength + FileLength;
    %Get the file name
    Temporary_Variable = h5info(HDF5DataFile{ii});
    Settings{ii}.FileName  = Temporary_Variable.Filename;
end
end
%-------------------------------------------------------------------------%
function [DownSampledData, UpdatedSamplingRate] = M_Downsampler(VoltageTrace, SamplingRate, DesiredSamplingRate)
%Use this function to downsample the data recroded from a MEA. It will the
%recorded frequency to a desired frequency. Also known as decimating.
%
%Inputs
%   MEAlData            - The voltage trace of the desired channel in µV
%   SamplingRate        - Sampling rate used in the data adquisition
%   DesiredSamplingRate - Desired end sampling rate
%Output
%   DownSampledData     - The Voltage trace downsampled to the desired
%                         rate
%   UpdatedSamplingRate - The new sampling rate at which the downsampled
%                         data was downsampled to
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN: May 13, 2017
%LAST MODIFIED: May 13, 2017
%v2.0


DownSampledValue = floor(SamplingRate/DesiredSamplingRate); 
% DownSampledData = zeros(
%Determine the amount of channels
SizeVariable = size(VoltageTrace);
NumChannels = SizeVariable(2);
TDownSampledData = [];

%Downsample the data to the desired frequency for all Channels
for i = 1:NumChannels
    Trace =  VoltageTrace(:,i);             
    DownSampledChannelData = Trace(1:DownSampledValue:end,:);    
    %Add the ith channel data to the downsampled data set.
        %OPTIMIZATION:It works through concatanation, this is not the fastest 
        %nor RAM saving way to do it
    TDownSampledData = [TDownSampledData DownSampledChannelData];
end

%Set the new sampling rate
DownSampledData     = TDownSampledData;
UpdatedSamplingRate = SamplingRate/DownSampledValue;    
end
%-------------------------------------------------------------------------%
function [NAME] = M_MATLABsaver(FileName)
%This code is meant to be used in order to determine the name of the file
%you want to save
%MEAData

SplitFilePathUsed = strsplit(FileName,'\');

%Convert the last cell array that contains the File Name
CharSFPU = char(SplitFilePathUsed(end));

SplitFileName = strsplit(CharSFPU,'_');

%Obtain the date and the Slice number from the file name
Date     = char(SplitFileName(1));
SliceNum = char(SplitFileName(3));

NAME = strcat('MEAData_', Date,'_',SliceNum);
end
%-------------------------------------------------------------------------%
function [PowerData, Relative] = PowerPCalc(Power,Frequency,No60Hz)
%This function is used in order to calculate the power in the frequencies

%Inputs
%   Power at all frequencies
%   Frequency
%Output
%   PowerData = structure with Raw and relative powers
%
%NOTE:
%Standard frequency ranges based on literature and multichannel systems 
%Delta : 0.5 - 4 Hz
%Theta : 4   - 8 Hz
%Alpha : 8   - 13 Hz
%Beta  : 13  - 30 Hz
%Gamma : 30  - 80 Hz

%Define the gamma frequency range

DeltaF = [2,4];%adjusted for the low freq mechanical noise
ThetaF = [5,8];
AlphaF = [9,12];
BetaF  = [13,24];
if No60Hz == 1
    GammaLF = [25,59];
elseif No60Hz == 0
    GammaLF = [25,60];
end
GammaHF = [61,100];

%Calculate the Area under the curve as an estimate of the total power in
%the gamma band

[Delta]      = AUC(Power, Frequency, DeltaF);
[Theta]      = AUC(Power, Frequency, ThetaF);
[Alpha]      = AUC(Power, Frequency, AlphaF);
[Beta]       = AUC(Power, Frequency,BetaF);
[HighGamma]  = AUC(Power, Frequency,GammaHF);
[LowGamma]   = AUC(Power, Frequency,GammaLF);
TotalPower   = Delta+Theta+Alpha+Beta+LowGamma+HighGamma;

%Export
PowerData.Raw.Delta     = Delta;
PowerData.Raw.Theta     = Theta;
PowerData.Raw.Alpha     = Alpha;
PowerData.Raw.Beta      = Beta; 
PowerData.Raw.HighGamma = HighGamma;
PowerData.Raw.LowGamma  = LowGamma; 
PowerData.Raw.Total     = TotalPower;

Relative.HighGamma = HighGamma./TotalPower;
Relative.LowGamma  = LowGamma./TotalPower; 
Relative.Delta     = Delta./TotalPower;
Relative.Theta     = Theta./TotalPower;
Relative.Alpha     = Alpha./TotalPower;
Relative.Beta      = Beta./TotalPower; 
Relative.Total     = TotalPower./TotalPower;
end
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
margingOfError = 0.06103516;
LeftI  = find(x >= (FreqBoundaries(1)- margingOfError) & x <= (FreqBoundaries(1)+ margingOfError));
RightI = find(x >= (FreqBoundaries(2)- margingOfError) & x <= (FreqBoundaries(2)+ margingOfError));
[~,NumTimePoints] = size(Power);
AUCPower=zeros(1,NumTimePoints);
for TimePoint = 1:NumTimePoints
    y = (Power(:,TimePoint));
    %Calculate the area under the curve
    AUCPower(TimePoint) = trapz(x(LeftI:RightI), y(LeftI:RightI));
end
end
%-------------------------------------------------------------------------%
%this is the code that will run the wavelets for the spectograms
function [TFR,timeVec,freqVec] = traces2TFR(S,freqVec,Fs,width)
% function [TFR,timeVec,freqVec] = traces2TFR(S,freqVec,Fs,width);
%
% Calculates the average of a time-frequency energy representation of
% multiple trials using a Morlet wavelet method.                            
%
% Input
% -----
% S    : signals = time x Trials      
% freqVec    : frequencies over which to calculate TF energy        
% Fs   : sampling frequency
% width: number of cycles in wavelet (> 5 advisable)  
%
% Output
% ------
% t    : time
% f    : frequency
% B    : phase-locking factor = frequency x time
%     
%------------------------------------------------------------------------
% Ole Jensen, Brain Resarch Unit, Low Temperature Laboratory,
% Helsinki University of Technology, 02015 HUT, Finland,
% Report bugs to ojensen@neuro.hut.fi
%------------------------------------------------------------------------?
%    Copyright (C) 2000 by Ole Jensen 
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You can find a copy of the GNU General Public License
%    along with this package (4DToolbox); if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


S = S';
LNR50 = 1; 

timeVec = (1:size(S,2))/Fs;  

B = zeros(length(freqVec),size(S,2)); 

for SampleIndex=1:size(S,1)          
    fprintf(1,'%d ',SampleIndex); 
    for FreqIndex=1:length(freqVec)
        if LNR50
            B(FreqIndex,:) = energyvec(freqVec(FreqIndex),lnr50(detrend(S(SampleIndex,:)),Fs),Fs,width) + B(FreqIndex,:);
        else
            B(FreqIndex,:) = energyvec(freqVec(FreqIndex),detrend(S(SampleIndex,:)),Fs,width) + B(FreqIndex,:);
        end
    end
end
TFR = B/size(S,1);     
end

function y = energyvec(f,s,Fs,width)
% function y = energyvec(f,s,Fs,width)
%
% Return a vector containing the energy as a
% function of time for frequency f. The energy
% is calculated using Morlet's wavelets. 
% s : signal
% Fs: sampling frequency
% width : width of Morlet wavelet (>= 5 suggested).
%
% 

dt = 1/Fs;
sf = f/width;
st = 1/(2*pi*sf);

t=-3.5*st:dt:3.5*st;
m = morlet(f,t,width);
y = conv(s,m);
y = (2*abs(y)/Fs).^2;
y = y(ceil(length(m)/2):length(y)-floor(length(m)/2));
end

function y = morlet(f,t,width)
% function y = morlet(f,t,width)
% 
% Morlet's wavelet for frequency f and time t. 
% The wavelet will be normalized so the total energy is 1.
% width defines the ``width'' of the wavelet. 
% A value >= 5 is suggested.
%
% Ref: Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997)
%
% See also: PHASEGRAM, PHASEVEC, WAVEGRAM, ENERGY 
%
% Ole Jensen, August 1998 

sf = f/width;
st = 1/(2*pi*sf);
A = 1/(st*sqrt(2*pi));

y = A*exp(-t.^2/(2*st^2)).*exp(1i*2*pi*f.*t);
end
%-------------------------------------------------------------------------%
function Sc = lnr50(S,Fs)
% function Sc = lnr50Hz(S,Fs)
%
% Line noise reduction (50 Hz) The amplitude and phase of the line noise is 
% estimated. A sinusoide with these characteristics is then subtracted from
% the signal. In many cases this approach is better than a notch filter, since 
% it does not reduce remove brain responses in the 50 Hz band. Works excellent
% on EFs with strong line noise.
%
% S   : time x trials    
% Fs  : sampling frequency
% Sc  : S 'cleaned' from 50 Hz
%------------------------------------------------------------------------
% Ole Jensen and Cristina Simoes, Brain Resarch Unit, Low Temperature 
% Laboratory, Helsinki University of Technology, 02015 HUT, Finland,
% Report bugs to ojensen@neuro.hut.fi
%------------------------------------------------------------------------

%    Copyright (C) 2000 by Ole Jensen 
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You can find a copy of the GNU General Public License
%    along with this package (4DToolbox); if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


fNoise=60.0;
tv = (0:length(S)-1)/Fs;
S = S';

transP = 0;
if size(S,2) == 1
    S = S';
    transP = 1;
end

for k=1:size(S,1)
   Sft = ft(S(k,:),fNoise,Fs);
   Sc(k,:) = S(k,:) - abs(Sft)*cos(2*pi*fNoise*tv - angle(Sft));
end

Sc = Sc';

if transP 
    Sc = Sc';
end
end
function S=ft(s, f, Fs)

tv = (0:length(s)-1)/Fs;
tmp = exp(i*2*pi*f*tv);
S= 2*sum(s.*exp(i*2*pi*f*tv))/length(s);

end