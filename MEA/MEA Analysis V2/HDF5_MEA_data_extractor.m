function MEAData = HDF5_MEA_data_extractor(HDF5DataFile)
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
%LAST MODIFIED: May 9, 2018
%v1.0

%Create struct
MEAData = struct;
%Preallocate for ChannelData
Temp = struct;

    %Import the data into each struct
    %Data are in int 32, to analyze them much change the to double by using
    %double(data)
    %The ChannelData is no longer being stored in the structure in order to
    %save memory.
    %define the objects as RawData and Settings
    
    
    Temp.ChannelData                = h5read(HDF5DataFile, '/Data/Recording_0/AnalogStream/Stream_0/ChannelData');
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
    %the conversion factor is ~59605. Included the change into �V. This info was obtained from the InfoChannel ConversionFactor.
    VoltageTrace = (double(Temp.ChannelData)*double(MEAData.Settings.InfoChannel.ConversionFactor(1)))/(10^6); 
    [MEAData.RawData.VoltageTrace, MEAData.RawData.SamplingRate] = M_downsampler(VoltageTrace, 20000, 5000);
    



