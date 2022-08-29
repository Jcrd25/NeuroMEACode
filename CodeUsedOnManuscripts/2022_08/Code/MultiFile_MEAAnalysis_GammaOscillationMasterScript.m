%MEA Analysis Master Script
%Import Data from HDF5 File
%Select HDF5 File and create a struct that contains the files to be
%analyzed
%You may select mulitple files but then you need to adjust the following
%script segments to analyze the desired channels at a time.

NumFilesToStudy = 3;

for ii = 1:NumFilesToStudy
    FileG(ii).HDF5DataFile = uipickfiles;  
    T1 = sprintf('Selected File Set %d',ii);
    disp(T1);
end
%%
for FG = 1:NumFilesToStudy
    
HDF5DataFile = FileG(FG).HDF5DataFile;
    
SavePath = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\ShafferCut\FilteredRawData';
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

%%
%Calculate the spectograms and periodograms
SavedPath   = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\FilteredRawData\Jenkins_AnkG_Data\';
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
SavedPath   = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\FilteredRawData\Jenkins_AnkG_Data\';
Dir2Save = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\Jenkins\MultiTaperSpectrumChannelData\';
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
T1 = ['Saved' ' ' MEAFileNames{MEAFileIndx} 'in_' ];
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
ChanDataPathDir = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\ShafferCut\MultiTaperSpectrumChannelData\';
Save2Dir = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\ShafferCut\BandPower\';
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
    NamePD = strcat(Save2Dir,'RAW_Power\');
    NameRD = strcat(Save2Dir,'Relative_Power\');
    NamePD = char(NamePD); NameRD = char(NameRD);
    
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
%This section is to Prepare the data to calculate the average Peak
%Frequency and the Q values
SavedPath   = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\FilteredRawData\2Analyze\';
Save2DirSum = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\PowerSummaryData\';
Save2DirAll = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\SliceData\AllKainatePowerSummaryData\';

KainateEnd = 120; KainateStart = 60; %set in minutes
DrugDuration = 16;
TimeWindowAnalysis = 5;
%-------------------------%

Times4All = [KainateStart, KainateEnd;KainateEnd, (KainateEnd+DrugDuration)];
Time4Focused = [(KainateStart-TimeWindowAnalysis),KainateStart;(KainateEnd-TimeWindowAnalysis), KainateEnd;...
(KainateEnd+DrugDuration-TimeWindowAnalysis),(KainateEnd+DrugDuration)];
if ~exist(Save2DirSum, 'dir')
    mkdir(Save2DirSum)
end
if ~exist(Save2DirAll, 'dir')
    mkdir(Save2DirAll)
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
    Name4File = [Save2DirSum 'PowSumD_' MEAFileNames{FileNum}];
    save(Name4File, 'MEAPowerSum','-v7.3'); 
    clear Name4File
    clear MEAPowerSum
    
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
    Name4File = [Save2DirAll 'PowAllD_' MEAFileNames{FileNum}];
    save(Name4File, 'MEA_AllPower','-v7.3'); 
    T1 = ['Saved Power Band sumary Data for ' MEAFileNames{FileNum}];
    clear MEA_AllPower
    disp(T1);
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
ThetaHippoF = [4,12];
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
[ThetaHippo]      = AUC(Power, Frequency, ThetaHippoF);

%Export
PowerData.Raw.Delta     = Delta;
PowerData.Raw.Theta     = Theta;
PowerData.Raw.Alpha     = Alpha;
PowerData.Raw.Beta      = Beta; 
PowerData.Raw.HighGamma = HighGamma;
PowerData.Raw.LowGamma  = LowGamma; 
PowerData.Raw.Total     = TotalPower;

PowerData.Raw.ThetaHippo  = ThetaHippo; 

Relative.HighGamma = HighGamma./TotalPower;
Relative.LowGamma  = LowGamma./TotalPower; 
Relative.Delta     = Delta./TotalPower;
Relative.Theta     = Theta./TotalPower;
Relative.Alpha     = Alpha./TotalPower;
Relative.Beta      = Beta./TotalPower; 
Relative.Total     = TotalPower./TotalPower;
Relative.ThetaHippo     = ThetaHippo./TotalPower;
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


