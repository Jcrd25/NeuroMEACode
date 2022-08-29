%% preprocessing step
%-------------------------------------%
%define the parameters for the coherence data analysis
%data = struct2cell(ForDom2.Kainate);
addpath(genpath('/home/jcrd/CustomMatlabCode/'));
folder = '/home/jcrd/DevData/';
folder2save = '/scratch/kevjon_root/kevjon99/jcrd/Results/DevData/';
folder2SaveCOH = '/home/jcrd/Dev_COH/';
MEAFileName2use = {...
 'MEAData_2021-07-27_Slice1', 'MEAData_2021-07-27_Slice2',...
 'MEAData_2021-07-29_Slice1', 'MEAData_2021-07-29_Slice2',...
 'MEAData_2021-07-31_Slice1', 'MEAData_2021-07-31_Slice2',...
 'MEAData_2021-08-04_Slice1', 'MEAData_2021-08-04_Slice2',...
 'MEAData_2021-08-06_Slice1', 'MEAData_2021-08-06_Slice2',...
 'MEAData_2021-08-07_Slice1', 'MEAData_2021-08-11_Slice1',...
 'MEAData_2021-08-11_Slice2', 'MEAData_2021-08-13_Slice1',...
 'MEAData_2021-08-13_Slice2', 'MEAData_2021-08-14_Slice3',...
 'MEAData_2021-08-17_Slice1', 'MEAData_2021-08-17_Slice2',...
 'MEAData_2021-08-19_Slice1', 'MEAData_2021-08-19_Slice2',...
 'MEAData_2021-08-20_Slice1', 'MEAData_2021-08-20_Slice2'};


%fileID2use = 9;
%TimeWindow2Use = 10; %in minutes
%Phase2Analyze = [1]; 
%Baseline Phase = 1; Kainate = 2; Drug = 3;
% fileID2use     = str2num(getenv('fileID2use')); %use this for single batch submission
fileID2use     = str2num(getenv('SLURM_ARRAY_TASK_ID')); %use this for array batch submissions
TimeWindow2Use = str2num(getenv('TimeWindow2Use'));
Phase2Analyze  = str2num(getenv('Phase2Analyze'));
DesiredWorker  = str2num(getenv('DesiredWorker'));
SetWorker = true; 
%-------------------------------------%
%define the parameters for the spectral analysis in this section
params.tapers = [5 9];
params.Fs = 1000;
params.fpass = [0 100];
movingwin = [10 0.05];
GammaBandLims = [25 59];
%-------------------------------------%
if Phase2Analyze == 1
    folder2save  = [folder2save  'Baseline' filesep];
    FileType = 'Baseline';
elseif Phase2Analyze == 2
    folder2save  = [folder2save  'Kainate' filesep];
    FileType = 'Kainate';
elseif Phase2Analyze == 3
    folder2save  = [folder2save  'Drug' filesep];
    FileType = 'Drug';
end
%-------------------------------------------------------------------------%
fnum = fileID2use;
clear MEAData
%get the data's location
shortname = MEAFileName2use{fnum};
tic
%load the data
load([folder shortname '.mat']);
disp('load time:');toc
tic
%determine save location
mkdir(folder2save, shortname);
file_folder = [folder2save shortname];

if ~exist(folder2SaveCOH, 'dir')
    mkdir(folder2SaveCOH)
end


%decimating save location
mkdir(file_folder, 'decimated');
file_folder = [file_folder filesep 'decimated'];
     
data = struct2cell(MEAData.RawData);
Fs = data{2,1,1};
hr = 3600 * Fs ;
TimeFactor4Kainate = 1.4667; %at what time did the Kainate portion end (ie 120min of recording then the factor is 2; if its 28 min the factor is 1.4667 as it is calculated by [(DurationofBaseline+DurationofKainate)/DurationofBaseline])
TimeWindowinSec = 60*Fs*TimeWindow2Use; %define the desidered time window
% get the indexes for the desired time window
if Phase2Analyze == 1
    sec_hr = (hr - TimeWindowinSec):hr; 
elseif Phase2Analyze == 2
    sec_hr = ((hr*TimeFactor4Kainate) - TimeWindowinSec):(hr*TimeFactor4Kainate); 
elseif Phase2Analyze == 3
%     sec_hr = (hr - TimeWindowinSec):hr;
end

dataInd2analyize = sec_hr;
%select the time segment of data to be analyzed
data = [data{1,1,:}];
data = data(dataInd2analyize,:); 
% end data selection
%------------------------------------------------------------------------%
%decimation step
n = 10; % average every n values
s1 = size(data, 1);      % Find the next smaller multiple of n
m  = s1 - mod(s1, n);
y  = reshape(data(1:m,:), n, m/n, []);     % Reshape x to a [n, m/n] matrix
Avg = sum(y, 1) / n; %mean accross column

data = squeeze(Avg);
params.Fs = Fs/n; %New decimated time
clear m y Avg

num_elec = size(data,2);
disp('serial rest')
toc
%num_elec = 5;
%determine the size of the data to preallocate 
[Cmats,~,~,~,~,t,f2]=cohgramc(data(:,1),data(:,2),movingwin,params);
CmatsSizes(1) = size(Cmats,1);CmatsSizes(2) = size(Cmats,2);
AllTimer = tic;
[~,TimefileSize] = size(t);
%preallocate the coh file
coh = zeros(TimefileSize,60,60);
clear Cmats 
%=========================================================================%
%=========================================================================%
%Parallel sections
%We get the number of tasks from the job and use that for the
%the number of workers
%-------------------------------------------------------------------------%
if isempty(getenv('SLURM_NTASKS'))
    nWorkers = 1;
elseif SetWorker == 1
   nWorkers = DesiredWorker;
else
   nWorkers = str2double(getenv('SLURM_NTASKS'))*str2double(getenv('SLURM_CPUS_PER_TASK'));
end
%Set up the Matlab cluster object
theCluster = parcluster('local');
%Create the pool of workers
thePool = parpool(theCluster, nWorkers);
%That worked, right?  If not, exit
if isempty(thePool)
    error('pctexample:backslashbench:poolClosed', ...
         ['This example requires a parallel pool. ' ...
          'Manually start a pool using the parpool command or set ' ...
          'your parallel preferences to automatically start a pool.']);
    exit
end
%Calculate the coherence between different electrodes
parfor electNum1 = 1:num_elec-1
%     clear Cmats PHImats Smats S12 S1 tnot f2not
    %preallocate the variables
    CmatsLopp = zeros(CmatsSizes(1),CmatsSizes(2),num_elec-1); 
    PHImats = CmatsLopp;
    Smats = zeros(CmatsSizes(1),CmatsSizes(2));
    %Get the coherence for the chosen electrode1 with the rest
    for electNum2 = electNum1+1:num_elec
    %calculate the coherence        
    [CmatsLopp(:,:,electNum2) ,PHImats(:,:,electNum2),S12,S1,Smats,tnot,f2not]=...
        cohgramc(data(:,electNum1),data(:,electNum2),movingwin,params);
    end
    %save the data in the a file that will store the important parameters
    %for the coherence analysis for the chosen electrode
%     save([file_folder filesep 'hippo' num2str(electNum1) '.mat'], '-v7.3', 'CmatsLopp', 'Smats', 'PHImats', 'f2', 't');
%     ParForSaver(CmatsLopp,Smats,PHImats,f2,t,file_folder, electNum1)
    disp(['hippo' num2str(electNum1)]);toc(AllTimer)
    %calculate coh
    gamma_ind = find(f2not>GammaBandLims(1)&f2not<GammaBandLims(2));
    coh(:,electNum1,:) = mean(CmatsLopp(:,gamma_ind,:),2);
end

SaveName = [folder2SaveCOH 'COH_' shortname '_' FileType];
save(SaveName, '-v7.3', 'coh', 't');
T2 = spritnf('Finished Saving coherence data for %s which is ArrayIndex %d',shortname,fileID2use);
disp(T2); 
disp(shortname)
disp('total time:')
toc(AllTimer)
delete(thePool)
exit
%-----------------------------------%
function ParForSaver(CmatsLopp,Smats,PHImats,f2,t,file_folder, electNum1)
save([file_folder filesep 'hippo' num2str(electNum1) '.mat'], '-v7.3', 'CmatsLopp', 'Smats', 'PHImats', 'f2', 't');
end
