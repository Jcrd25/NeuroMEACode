%%This script is meant to be used to generate plots for recordings done in
%%multiple neuronal cultures in MEAs

%Define all important variables for analysis

%Select the conditions that will be used in this case
Condition.DrugTreatment{1} = 'Vehicle';
Condition.DrugTreatment{2} = 'MK801';

Condition.DIV = [10 ,17, 24];

%Which neuronal cultures got each treatment

%August Culture
VehicleCultures = {'MEA-31991','MEA-31999', 'MEA-32000', 'MEA-32985', ...
                   'MEA-33075', 'MEA-33076'};
TreatedCultures = {'MEA-31990', 'MEA-33759', 'MEA-32001'};
% %October Culture
% VehicleCultures = {'MEA-31990','MEA-32001', 'MEA-32000', 'MEA-32985', ...
%                    'MEA-33075', 'MEA-33076'};
% TreatedCultures = { 'MEA-33759','MEA-31999','MEA-31991'};
% 
% %November Culture
% VehicleCultures = {'MEA-32985','MEA-31999'};
% TreatedCultures = {'MEA-31990', 'MEA-31991'};

SavedPathway = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\NeuronalCulture\NC_2018_08_31\NC Data\';
SavingInfo.SaveDirectory = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\NeuronalCulture\Excel Files\NC_2018_08_31\';
SavingInfo.VehicleName = VehicleCultures;
SavingInfo.TreatedName = TreatedCultures;
listing = dir(SavedPathway);
%%
CultureName = [];
Culture = [];
%Make a structure containing all of the vehicle treated cultures
for CultureIndex = 1:length(VehicleCultures)
    Counter = 1;
    %Get the complete file name for all Vehicle Treated Cultures
    for FileNum = 1:length(listing)
        String = sprintf(listing(FileNum).name);
        String2 = sprintf(VehicleCultures{CultureIndex});
        %If it's the MEA of interest, save the Data file path
        if isempty(strfind(String,String2)) == 0
            CultureName(CultureIndex).FileName{Counter} =  ...
                strcat(SavedPathway,'\',listing(FileNum).name);
            Counter = Counter + 1;
        end
    end
    
end

%plot the MEA mean data
for CultureID = 1:length(CultureName)
    for DIV_ID = 1:length(CultureName(CultureID).FileName)
        load(CultureName(CultureID).FileName{DIV_ID});
        Culture(CultureID).DIV(DIV_ID).BurstData = MEANCData.BurstData;
        Culture(CultureID).DIV(DIV_ID).PeakData  = MEANCData.PeakData;
        Culture(CultureID).Time                  = MEANCData.Time(1).RecordingTime;
    end
end
SavingInfo.Treatment = 'Vehicle';
% NC_CulturePlotter2(Culture,Condition.DIV,SavingInfo);
NC_BurstPlotter(Culture,Condition.DIV,SavingInfo);   
% [ISIExported,Val2Kev] = NC_ISIHistogram(Culture,Condition.DIV,SavingInfo, VehicleCultures);          
 
 CultureName = [];
 Culture = [];
 %Make a structure containing all of the Drug treated cultures
for CultureIndex = 1:length(TreatedCultures)
    Counter = 1;
    %Get the complete file name for all Drug Treated Cultures
    for FileNum = 1:length(listing)
        String = sprintf(listing(FileNum).name);
        String2 = sprintf(TreatedCultures{CultureIndex});
        %If it's the MEA of interest, save the Data file path
        if isempty(strfind(String,String2)) == 0
            CultureName(CultureIndex).FileName{Counter} =  ...
                strcat(SavedPathway,'\',listing(FileNum).name);
            Counter = Counter + 1;
        end
        
    end
    
end

%plot the MEA mean data
for CultureID = 1:length(CultureName)
    for DIV_ID = 1:length(CultureName(CultureID).FileName)
        load(CultureName(CultureID).FileName{DIV_ID});
        Culture(CultureID).DIV(DIV_ID).BurstData = MEANCData.BurstData;
        Culture(CultureID).DIV(DIV_ID).PeakData  = MEANCData.PeakData;
        Culture(CultureID).Time                  = MEANCData.Time(1).RecordingTime;
    end
end
SavingInfo.Treatment = 'MK801';
% NC_CulturePlotter2(Culture,Condition.DIV,SavingInfo);
NC_BurstPlotter(Culture,Condition.DIV,SavingInfo);   
% [ISIExported,Val2Kev] = NC_ISIHistogram(Culture,Condition.DIV,SavingInfo, TreatedCultures);


