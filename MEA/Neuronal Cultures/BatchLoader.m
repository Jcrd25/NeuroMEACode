function MatlabFiles = BatchLoader(SavedDataDirectory)
%This is set for the section of the plotting the different conditions of
%interest so you can get a figure based on the data


%obtain the name of the different folders in the directory
listing = dir(SavedDataDirectory);

for ii = 1:length(listing)
    %get the name of each local directory
    %determine if the name of the folders inside the specified directory
    if listing(ii).isdir == 1
        if contains(listing(ii).name,'NC') == 1
            NamedDirectory = fullfile(SavedDataDirectory,listing(ii).name);
            MatlabFiles.(listing(ii).name) = dir(fullfile(NamedDirectory, '*.mat'));
        end
    end 
end
