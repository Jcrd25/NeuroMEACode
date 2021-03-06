function MEAData = MEA_culturedata_structuremaker(HDF5DataFile)
%This function is meant to be used to import the MEA data from HDF5 files
%exported through MultiChannel DataManager from multichannel systems. It
%will create structure in MatLab. 
%The AnalogStream is important since it might be a filtered data set or the 
%raw trace. 
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN: April 12, 2017
%LAST MODIFIED: July 11, 2017
%v1.0

files_imported = 0;

%Create struct
MEAData = struct;
%Preallocate for ChannelData
Temp = struct;
for ii = 1:length(HDF5DataFile)
    %Import the data into each struct
    %Data are in int 32, to analyze them much change the to double by using
    %double(data)
    %The ChannelData is no longer being stored in the structure in order to
    %save memory.
    
    Temp(ii).ChannelData                = h5read(HDF5DataFile{ii}, '/Data/Recording_0/AnalogStream/Stream_0/ChannelData');
    MEAData(ii).Settings.InfoChannel    = h5read(HDF5DataFile{ii}, '/Data/Recording_0/AnalogStream/Stream_0/InfoChannel');
    MEAData(ii).Settings.DataType       = h5readatt(HDF5DataFile{ii}, '/Data/Recording_0/AnalogStream/Stream_0','Label');
    MEAData(ii).Settings.Date           = h5readatt(HDF5DataFile{ii}, '/Data','Date');
    MEAData(ii).Settings.MeaName        = h5readatt(HDF5DataFile{ii}, '/Data','MeaName');
    %Get the file name
    Temporary_Variable                  = h5info(HDF5DataFile{ii});
    MEAData(ii).FileName       = Temporary_Variable.Filename;
    %Make sure the conversion factor is the same for all channles
    for iii = 1:size(MEAData(ii).Settings.InfoChannel.ConversionFactor, 1)
        if MEAData(ii).Settings.InfoChannel.ConversionFactor(1) ~= MEAData(ii).Settings.InfoChannel.ConversionFactor(iii)
            fprintf('WARNING: The conversion factors are not the same for all channels please review InfoChannel \n')
        end
    end
    %Convert data into a Voltage trace
    %the conversion factor is ~59605. Included the change into �V. This info was obtained from the InfoChannel ConversionFactor.
    MEAData(ii).RawData.VoltageTrace = (double(Temp(ii).ChannelData)*double(MEAData(ii).Settings.InfoChannel.ConversionFactor(1)))/(10^6); 
    files_imported = files_imported + 1;
end

fprintf('Imported %d files\n', files_imported)
