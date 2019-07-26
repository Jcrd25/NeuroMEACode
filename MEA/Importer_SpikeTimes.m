function [SpikeTimes] = Importer_SpikeTimes(ImportedVariableName1, ImportedVariableName2, ImportedVariableName3)
%Import Spike Times from Text File
%Data is in seconds




ChannelNumber = [12,13,14,15,16,17,21,22,23,24,25,26,27,28,31,32,33,34,35,36,37,38,41,42,43,44,45,46,47,48,51,52,53,54,55,56,57,58,61,62,63,64,65,66,67,68,71,72,73,74,75,76,77,78,82,83,84,85,86,87];

%Get the Data for all of the Channels
%Determine where the strings 't' are located. normally they will be found
%2 lines before the time data is displayed. 
Indexofts = find(strncmp(ImportedVariableName1,'t',1));
    
        
for ii = 1:length(ChannelNumber)
    %Use the 't' as a marker of where the data should be.
    
        
       
                


