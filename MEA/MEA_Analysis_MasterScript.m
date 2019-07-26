%MEA Analysis Master Script
%Import Data from HDF5 File
%Select HDF5 File and create a struct that contains the files to be
%analyzed
%You may select mulitple files but then you need to adjust the following
%script segments to analyze the desired channels at a time.
  
HDF5DataFile = uipickfiles;    
MEAData = HDF5_MEA_data_structuremaker(HDF5DataFile);
%%
%Downsample the data from all the files
for n = 1:length(MEAData) 
     
    tic
    [MEAData(n).RawData.DownSampledData, MEAData(n).RawData.DownSamplingRate] = MEA_downsampler(MEAData(n),20000, 20000);
    toc
end
%%
%This part if for plotting all 60 channels
%Define n as the file index of the file you want to look at
for n = 1
    
    MEADataPlotter(MEAData(n))
    
end
%%
%this functions gives out the periodograms for the specified channel number
%REMEMBER the channel number and ID are not necesarily the same
%There are three settings, see MEAPowerSpectrumCalc help for more
%information

%normally I use Settings 1 or 2 to look at the spectograms of the channels 
%in order to determine which channels I want to later look at the 
%spectogram using setting 3
Setting = 4;
Plot = 'S';
%Select the file index of the file you want to look at
% for n = 1:length(MEAData) 
    
    for Chan = 1:60
        %select the channel index you wish to look at (The channel ID can be
        %checked in InfoChannel in MEAData) 

         figure
         hold on
        for n = [1:5]
            MEAPowerSpectrumCalc(MEAData(n), Chan, Setting,Plot);
    %         ylim([0 1])
    %          caxis([0 20]);

        end
        hold off
%         P1 = sprintf('Downsampled at %d', rate);
%         title(P1)

    end

%%
MEA_FFT_SpectCont(MEAData, SamplingRate,InterFileInterval)

%%
%Spike Detection
for n = 1:5
    ChannelNumber = 50;
    MEA_RasterPlotter(MEAData(n), ChannelNumber, 'Y')
end

%%
tic
NumFiles            = length(MEAData);  
NumChannels         = 60;
MEA_CountSpike(MEAData, NumFiles, NumChannels)
toccl
   


%%
for ii = 1:length(MEAData)
    MEAData(ii).AnalyzedData.Gpresence = zeros(60, 1);
end
%%
for Channel= 1:60
    for ii = 15:length(MEAData)
        
        if MEAData(ii).AnalyzedData.GammaPower(Channel) < MEAData(1).AnalyzedData.GammaPower(Channel)
            MEAData(ii).AnalyzedData.GPresence(Channel) = 1;
%             fprintf('Channel found to have gamma %d', Channel)
        else
            MEAData(ii).AnalyzedData.GPresence(Channel) = 0 ;
        end
    end
end


%%
%Calculate the Power at the gamma band (30- 80 Hz)
for n = 1:length(MEAData)
    for Chan = 1:60
        [MEAData(n).AnalyzedData.GammaPower(Chan),...
            MEAData(n).AnalyzedData.GammaPowerHigh(Chan),...
            MEAData(n).AnalyzedData.GammaPowerLow(Chan)]...
            = MEA_GammaPCalc(MEAData(n), Chan);
    end
MEAData(n).AnalyzedData.GammaPower     = MEAData(n).AnalyzedData.GammaPower';
MEAData(n).AnalyzedData.GammaPowerHigh = MEAData(n).AnalyzedData.GammaPowerHigh';
MEAData(n).AnalyzedData.GammaPowerLow  = MEAData(n).AnalyzedData.GammaPowerLow';
end
%%
BasFile = 1;
ChemInd = 3;
Cond1   = 5;

%Positive Channels
PosChan = struct;
n3=1;
for Channel = [12,31,39]
%    figure
%    hold on
%    n=1;
   
%     for ii=[1,10,15]
%             y(n) = MEAData(ii).AnalyzedData.GammaPower(Channel);
%             n=n+1;
%     end
%         plot(1:3, y)
%     
%     hold off
%     ChannelDouble = str2double(MEAData(1).Settings.InfoChannel.Label(Channel));
%     P1 = sprintf('Channel %d', ChannelDouble);
%     title(P1)
    
    n2 = 1;
    
    PosChan.TGP(n3).GP1 = MEAData(BasFile).AnalyzedData.GammaPower(Channel)/MEAData(1).AnalyzedData.GammaPower(Channel);
    PosChan.TGP(n3).GP2 = MEAData(ChemInd).AnalyzedData.GammaPower(Channel)/MEAData(1).AnalyzedData.GammaPower(Channel);
    PosChan.TGP(n3).GP3 = MEAData(Cond1).AnalyzedData.GammaPower(Channel)/MEAData(1).AnalyzedData.GammaPower(Channel);
    
    PosChan.HGP(n3).GP1 = MEAData(BasFile).AnalyzedData.GammaPowerHigh(Channel)/MEAData(1).AnalyzedData.GammaPowerHigh(Channel);
    PosChan.HGP(n3).GP2 = MEAData(ChemInd).AnalyzedData.GammaPowerHigh(Channel)/MEAData(1).AnalyzedData.GammaPowerHigh(Channel);
    PosChan.HGP(n3).GP3 = MEAData(Cond1).AnalyzedData.GammaPowerHigh(Channel)/MEAData(1).AnalyzedData.GammaPowerHigh(Channel);
    
    PosChan.LGP(n3).GP1 = MEAData(BasFile).AnalyzedData.GammaPowerLow(Channel)/MEAData(1).AnalyzedData.GammaPowerLow(Channel);
    PosChan.LGP(n3).GP2 = MEAData(ChemInd).AnalyzedData.GammaPowerLow(Channel)/MEAData(1).AnalyzedData.GammaPowerLow(Channel);
    PosChan.LGP(n3).GP3 = MEAData(Cond1).AnalyzedData.GammaPowerLow(Channel)/MEAData(1).AnalyzedData.GammaPowerLow(Channel);
        
    ChannelDouble = str2double(MEAData(1).Settings.InfoChannel.Label(Channel));
    P1 = sprintf('Channel %d', ChannelDouble);
    PosChan.Label(n3).ChanID = P1;
    
    n3= n3+1;
end

%Negative Channels
NegChan = struct;
n3=1;

for Channel = [4,17,18,19,20,21,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,42,43,44,48,49,50,51,52,53,54,56,57,58,59,60]   
    
    NegChan.TGP(n3).GP1 = MEAData(BasFile).AnalyzedData.GammaPower(Channel)/MEAData(BasFile).AnalyzedData.GammaPower(Channel);
    NegChan.TGP(n3).GP2 = MEAData(ChemInd).AnalyzedData.GammaPower(Channel)/MEAData(BasFile).AnalyzedData.GammaPower(Channel);
    NegChan.TGP(n3).GP3 = MEAData(Cond1  ).AnalyzedData.GammaPower(Channel)/MEAData(BasFile).AnalyzedData.GammaPower(Channel);
    
    NegChan.HGP(n3).GP1 = MEAData(BasFile).AnalyzedData.GammaPowerHigh(Channel)/MEAData(BasFile).AnalyzedData.GammaPowerHigh(Channel);
    NegChan.HGP(n3).GP2 = MEAData(ChemInd).AnalyzedData.GammaPowerHigh(Channel)/MEAData(BasFile).AnalyzedData.GammaPowerHigh(Channel);
    NegChan.HGP(n3).GP3 = MEAData(Cond1  ).AnalyzedData.GammaPowerHigh(Channel)/MEAData(BasFile).AnalyzedData.GammaPowerHigh(Channel);
    
    NegChan.LGP(n3).GP1 = MEAData(BasFile).AnalyzedData.GammaPowerLow(Channel)/MEAData(BasFile).AnalyzedData.GammaPowerLow(Channel);
    NegChan.LGP(n3).GP2 = MEAData(ChemInd).AnalyzedData.GammaPowerLow(Channel)/MEAData(BasFile).AnalyzedData.GammaPowerLow(Channel);
    NegChan.LGP(n3).GP3 = MEAData(Cond1  ).AnalyzedData.GammaPowerLow(Channel)/MEAData(BasFile).AnalyzedData.GammaPowerLow(Channel);
    
        
    ChannelDouble = str2double(MEAData(1).Settings.InfoChannel.Label(Channel));
    P1 = sprintf('Channel %d', ChannelDouble);
    NegChan.Label(n3).ChanID = P1;
    
    n3= n3+1;
end
% Organize the data in a manner that is easy to save and export to other
% data plotting tools such as GraphPad Prism
GP = cell(60,length(MEAData));
for chan = 1:60
    for ii= 1:length(MEAData)
        GP{chan,ii} = MEAData(ii).AnalyzedData.GammaPower(chan);
    end
end
SummaryAllData.GammaPower = GP;

%Save data
MData.MEAData = MEAData;
MData.Summary.GammaRange.AllData     = SummaryAllData;
MData.Summary.GammaRange.PosChannels = PosChan;
MData.Summary.GammaRange.NegChannels = NegChan;
%%
clearvars -except MEAData MData