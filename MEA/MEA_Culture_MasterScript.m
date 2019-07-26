%MEA Analysis Master Script
%Import Data from HDF5 File
%Select HDF5 File and create a struct that contains the files to be
%analyzed
%You may select mulitple files but then you need to adjust the following
%script segments to analyze the desired channels at a time.

HDF5DataFile = uipickfiles;
MEAData = MEA_culturedata_structuremaker(HDF5DataFile);


%%
%Downsample the data from all the files

for n = 1:length(MEAData)
    tic
    [MEAData(n).DownSampledData, MEAData(n).UpdatedSamplingRate] = MEA_downsampler(MEAData(n),20000, 10000);
    toc
end
%%
%This part if for plotting all 60 channels
%Define n as the file index of the file you want to look at
for n = [1:12]
   
    
    MEADataPlotter(MEAData(n))
    
end

%%
%Spike Detection
for n = 1:12;
    ChannelNumber = 1:60;
    MEA_RasterPlotter(MEAData(n), ChannelNumber, 'N')
end

%%
%This function might not be needed in this script


%this functions gives out the periodograms for the specified channel number
%REMEMBER the channel number and ID are not necesarily the same
%There are three settings, see MEAPowerSpectrumCalc help for more
%information

%normally I use Settings 1 or 2 to look at the spectograms of the channels in order to determine which channels I want to later look at the spectogram using setting 3
Setting = 2;
%Select the file index of the file you want to look at
for n = 1:12
    %select the channel index you wish to look at (The channel ID can be
    %checked in InfoChannel in MEAData) 
    for Chan = 50
        
        MEAPowerSpectrumCalc(MEAData(n), Chan, Setting)
    
    end
end