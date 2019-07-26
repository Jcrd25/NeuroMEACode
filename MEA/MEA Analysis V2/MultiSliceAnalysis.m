%This script is meant to be used to generate figures containing all
%generated from files used using one condition.
%MK801 injection
FileName{1} = 'H:\Documents\Analyzed_MATLAB\MK801 inj\MEAData_2018-04-19_Slice1.mat';
FileName{2} = 'H:\Documents\Analyzed_MATLAB\MK801 inj\MEAData_2018-04-19_Slice2.mat';
FileName{3} = 'H:\Documents\Analyzed_MATLAB\MK801 inj\MEAData_2018-04-10_Slice1.mat';
%%
%PBS injection
FileName{1} = 'H:\Documents\Analyzed_MATLAB\MK801 inj\MEAData_2018-04-11_Slice2.mat';
% FileName{2} = 'H:\Documents\Analyzed_MATLAB\MK801 inj\MEAData_2018-04-03_Slice1.mat';
FileName{2} = 'H:\Documents\Analyzed_MATLAB\MK801 inj\MEAData_2018-05-03_Slice2.mat';
%%
%Jenkins WT
FileName{1} = 'H:\Documents\Analyzed_MATLAB\Jenkins\MEAData_2017-12-01_Slice2.mat';
FileName{2} = 'H:\Documents\Analyzed_MATLAB\Jenkins\MEAData_2017-12-06_Slice3.mat';
FileName{3} = 'H:\Documents\Analyzed_MATLAB\Jenkins\MEAData_2018-03-02_Slice2.mat';
FileName{4} = 'H:\Documents\Analyzed_MATLAB\Jenkins\MEAData_2018-01-17_Slice2.mat';
FileName{5} = 'H:\Documents\Analyzed_MATLAB\Jenkins\MEAData_2018-01-11_Slice2.mat';
%%
%Jenkins Mutant
FileName{1} = 'H:\Documents\Analyzed_MATLAB\Jenkins\MEAData_2017-12-08_Slice3.mat';
FileName{2} = 'H:\Documents\Analyzed_MATLAB\Jenkins\MEAData_2017-12-14_Slice1.mat';
FileName{3} = 'H:\Documents\Analyzed_MATLAB\Jenkins\MEAData_2017-12-14_Slice2.mat';
FileName{4} = 'H:\Documents\Analyzed_MATLAB\Jenkins\MEAData_2018-02-22_Slice1.mat';
FileName{5} = 'H:\Documents\Analyzed_MATLAB\Jenkins\MEAData_2018-02-22_Slice3.mat';


%%
%Bicuculline block
FileName{1} = 'H:\Documents\Analyzed_MATLAB\Bicuculline Recrodings\MEAData_2018-08-27_Slice2.mat';
% FileName{2} = 'H:\Documents\Analyzed_MATLAB\Bicuculline Recrodings\MEAData_2018-08-28_Slice1.mat';
FileName{2} = 'H:\Documents\Analyzed_MATLAB\Bicuculline Recrodings\MEAData_2018-09-12_Slice1.mat';
FileName{3} = 'H:\Documents\Analyzed_MATLAB\Bicuculline Recrodings\MEAData_2018-09-12_Slice2.mat';
FileName{4} = 'H:\Documents\Analyzed_MATLAB\Bicuculline Recrodings\MEAData_2018-09-14_Slice1.mat';
FileName{5} = 'H:\Documents\Analyzed_MATLAB\Bicuculline Recrodings\MEAData_2018-09-14_Slice2.mat';
%%
%MK801 Block

FileName{1} = 'H:\Documents\Analyzed_MATLAB\MK801 Bath application\MEAData_2018-09-18_Slice1.mat';
FileName{2} = 'H:\Documents\Analyzed_MATLAB\MK801 Bath application\MEAData_2018-09-18_Slice2.mat';

%%
%Caculate the average power spectrum for each brain region during the last
%kainate application for each slice
LastBaseline = 30;
LastKainate  = 60;
% LastBaseline = 1;
% LastKainate  = 31;

FirstFreq    = 1;
SecondFreq   = 100;

for File = 1:length(FileName)
%  for File = 1:5
    Files(File) = MEA_BatchDataExtractor(FileName{File},LastBaseline, FirstFreq, SecondFreq);
         
 end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NORMALIZED DATA SEGMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Plot the compounded figures for CA1
%Define the number of files that are analyzed
NumofFiles = 38;
FirstFile = LastBaseline;
%Begin by setting up the time (x values in minutes) 
x = (0:NumofFiles-1)*2;
yy = zeros(1, NumofFiles);
y=[];
for File = 1:length(FileName)
    if isfield(Files(File).ChannelID, 'CA1Channels') == 1
    %empty y and indy
    indy=[];
    %Get the number of active channels
        NumChan = length(Files(File).ChannelID.CA1Channels);
        %Get all of the data from each time point
        
        for ii = 1:NumofFiles
            FN = ii + FirstFile -1 ;
            %define the first timepoint (baseline) as 1
                if ii ==1
                    y(File,ii) = 1;
                    for n = 1:NumChan
                        indy(ii,n) = 1;
                    end
                else
                    y(File,ii) = Files(File).Time(FN).CA1MeanLG;

                    %Get the Low gamma for all the individual channels
                    for n = 1:NumChan
                        indy(ii,n) = Files(File).Time(FN).CA1LG(n);
                    end
                end

        end
        yy = yy+y(File,:);
        figure
        hold on
        for n = 1:NumChan
            plot(x,indy(:,n),'Color', [0 0 1], 'LineWidth', 1)
        end
        plot(x,y(File,:),'k','LineWidth', 6)
        set(findall(gca, '-property', 'FontSize'),'FontSize',24);
        xlabel('Time (min)')
        ylabel('Change in power from baseline')
        ExtractedName = M_AnalyzedNamer(FileName{File});
        P1 = sprintf('Gamma Power of CA1 in %s', ExtractedName);
        title(P1,'Interpreter', 'none')
    end
end
%plot the average
figure 
yy = yy*(1/length(FileName));
plot(x,yy,'k','LineWidth', 6)
set(findall(gca, '-property', 'FontSize'),'FontSize',24);
xlabel('Time (min)')
ylabel('Change in power from baseline')
P1 = sprintf('Average Gamma Power of CA1'); 
title(P1)
%%
%Plot the compounded figures for CA3
%Define the number of files that are analyzed
NumofFiles = 45;
FirstFile = LastBaseline;
%Begin by setting up the time (x values in minutes) 
x = (0:NumofFiles-1)*2;
yy = zeros(1, NumofFiles);
y=[];
for File = 1:length(FileName)
    if isfield(Files(File).ChannelID, 'CA3Channels') == 1
    %empty y and indy
    indy=[];
    %Get the number of active channels
        NumChan = length(Files(File).ChannelID.CA3Channels);
        %Get all of the data from each time point
        
        for ii = 1:NumofFiles
            FN = ii + FirstFile -1 ;
            %define the first timepoint (baseline) as 1
                if ii ==1
                    y(File,ii) = 1;
                    for n = 1:NumChan
                        indy(ii,n) = 1;
                    end
                else
                    y(File,ii) = Files(File).Time(FN).CA3MeanLG;

                    %Get the Low gamma for all the individual channels
                    for n = 1:NumChan
                        indy(ii,n) = Files(File).Time(FN).CA3LG(n);
                    end
                end

        end
        yy = yy+y(File,:);
        figure
        hold on
        for n = 1:NumChan
            plot(x,indy(:,n),'Color', [0 0 1], 'LineWidth', 1)
        end
        plot(x,y(File,:),'k','LineWidth', 6)
        set(findall(gca, '-property', 'FontSize'),'FontSize',24);
        xlabel('Time (min)')
        ylabel('Change in power from baseline')
        ExtractedName = M_AnalyzedNamer(FileName{File});
        P1 = sprintf('Gamma Power of CA3 in %s', ExtractedName);
        title(P1,'Interpreter', 'none')
    end
end
%plot the average
figure 
yy = yy*(1/length(FileName));
plot(x,yy,'k','LineWidth', 6)
set(findall(gca, '-property', 'FontSize'),'FontSize',24);
xlabel('Time (min)')
ylabel('Change in power from baseline')
P1 = sprintf('Average Gamma Power of CA3'); 
title(P1)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SONOGRAM SEGMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CA1 Figure Sonogram
FileToPlot = 1;
NumofFiles = 38;
FirstFile = LastBaseline;
SumZ = zeros(100,NumofFiles);

for File = 1:length(FileName)

   
    if isfield(Files(File).ChannelID, 'CA1Channels') == 1
    figure
    x = (0:NumofFiles-1)*2;
    y = 1:100;
    z=[];
    for ii = 1:NumofFiles
        FN = ii + FirstFile -1 ;
        for iii = 1:100
            if ii ==1
                z(iii,ii) = 1;
            else
                z(iii,ii) = Files(File).Time(FN).MeanCA1(iii);
            end
        end
    end
    imagesc(x,y,z)
    ylabel('Frequency (Hz)')
    xlabel('Time (min)')
    c = colorbar;
    ylabel(c, 'Fold Change')
    colormap jet
    ylim([0 100])
    caxis([0 7]);
    ham = gca;
    set(ham, 'YDir', 'normal')
    set(findall(gca, '-property', 'FontSize'),'FontSize',24);
    ExtractedName = M_AnalyzedNamer(FileName{File});
    P1 = sprintf('Power Spectrum of CA1 in %s', ExtractedName); 
    title(P1,'Interpreter', 'none')
    SumZ = SumZ+z;
    end
end

%%

%CA1 figure Average
figure
x = (0:NumofFiles-1)*2;
y = 1:100;
z = SumZ*(1/length(FileName));
imagesc(x,y,z)
ylabel('Frequency (Hz)')
xlabel('Time (min)')
c = colorbar;
ylabel(c, 'Fold Change')
colormap jet
ylim([0 100])
caxis([0 7]);
ham = gca;
set(ham, 'YDir', 'normal')
set(findall(gca, '-property', 'FontSize'),'FontSize',24);
P1 = sprintf('Average Sonogram of CA1'); 
title(P1,'Interpreter', 'none')
%%
%CA3 Plot

FileToPlot = 38;
SumZ = zeros(100,NumofFiles);
for File = 1:length(FileName)
    if isfield(Files(File).ChannelID, 'CA3Channels') == 1
    figure
    x = (0:NumofFiles-1)*2;
    y = 1:100;
    z=[];
    for ii = 1:NumofFiles
        FN = ii + FirstFile - 1;
        for iii = 1:100
            if ii ==1
                z(iii,ii) = 1;
            else
                z(iii,ii) = Files(File).Time(FN).MeanCA3(iii);
            end
        end
    end
    imagesc(x,y,z)
    ylabel('Frequency (Hz)')
    xlabel('Time (min)')
    c = colorbar;
    ylabel(c, 'Fold Change')
    colormap jet
    ylim([0 100])
    caxis([0 7.2]);
    ham = gca;
    set(ham, 'YDir', 'normal')
    set(findall(gca, '-property', 'FontSize'),'FontSize',24);
    ExtractedName = M_AnalyzedNamer(FileName{File});
    P1 = sprintf('Power Spectrum of CA3 in %s', ExtractedName); 
    title(P1,'Interpreter', 'none')
    SumZ = SumZ+z;
    end
end
%%
figure
x = (0:NumofFiles-1)*2;
y = 1:100;
z = SumZ*(1/length(FileName));
imagesc(x,y,z)
ylabel('Frequency (Hz)')
xlabel('Time (min)')
c = colorbar;
ylabel(c, 'Fold Change')
colormap jet
ylim([0 100])
caxis([0 3]);
ham = gca;
set(ham, 'YDir', 'normal')
set(findall(gca, '-property', 'FontSize'),'FontSize',24);
P1 = sprintf('Average Sonogram of CA3'); 
title(P1,'Interpreter', 'none')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RAW POWER SEGMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%CA1 Figure Sonogram
FileToPlot = 1;
NumofFiles = 68;
SumZ = zeros(100,NumofFiles);
NumZ = 0;
for File = 1
    if isfield(Files(File).ChannelID, 'CA1Channels') == 1
        for ChannelI = 1:length(Files(File).ChannelID.CA1Channels)
            Channel = Files(File).ChannelID.CA1Channels(ChannelI); 

            figure
            x = (0:NumofFiles-1)*2;
            y = 1:100;
           close all z=[];
            for ii = 1:NumofFiles
                for iii = 1:100
                    z(iii,ii) = Files(File).Time(ii).RawPowerCA1(ChannelI,iii);
                end
            end
            imagesc(x,y,10*log10(z))
            ylabel('Frequency (Hz)')
            xlabel('Time (min)')
            c = colorbar;
            ylabel(c, 'Fold Change')
            colormap jet
            ylim([0 100])
            caxis([0 10]);
            ham = gca;
            set(ham, 'YDir', 'normal')
            set(findall(gca, '-property', 'FontSize'),'FontSize',24);
            ExtractedName = M_AnalyzedNamer(FileName{File});
            P1 = sprintf('Power Spectrum of CA1 in Channel %d in file %s', Channel, ExtractedName); 
            title(P1,'Interpreter', 'none')
            NumZ = NumZ +1;
            SumZ = SumZ+z;
        end
    end
end
%%
figure
x = (0:NumofFiles-1)*2;
y = 1:100;
z = SumZ*(1/length(FileName));
imagesc(x,y,10*log10(z))
ylabel('Frequency (Hz)')
xlabel('Time (min)')
c = colorbar;
ylabel(c, 'Power (dB/Hz)')
colormap jet
ylim([0 100])
%caxis([0 20]);
ham = gca;
set(ham, 'YDir', 'normal')
set(findall(gca, '-property', 'FontSize'),'FontSize',24);
P1 = sprintf('Sonogram of Average CA1'); 
title(P1,'Interpreter', 'none')
%%
%CA3 Figure Sonogram
FileToPlot = 1;
NumofFiles = 74;
SumZ = zeros(100,NumofFiles);
NumZ = 0;
for File = 1:length(Files)
    if isfield(Files(File).ChannelID, 'CA3Channels') == 1
        for ChannelI = 1:length(Files(File).ChannelID.CA3Channels)
            Channel = Files(File).ChannelID.CA3Channels(ChannelI); 
        figure

        x = (0:NumofFiles-1)*2;
        y = 1:100;
        z=[];
        for ii = 1:NumofFiles
            for iii = 1:100

                    z(iii,ii) = Files(File).Time(ii).RawPowerCA3(ChannelI,iii);

            end
        end
        imagesc(x,y,10*log10(z))
        ylabel('Frequency (Hz)')
        xlabel('Time (min)')
        c = colorbar;
        ylabel(c, 'Fold Change')
        colormap jet
        ylim([0 100])
    %     caxis([0 10]);
        ham = gca;
        set(ham, 'YDir', 'normal')
        set(findall(gca, '-property', 'FontSize'),'FontSize',24);
        ExtractedName = M_AnalyzedNamer(FileName{File});
        P1 = sprintf('Power Spectrum of CA3 in Channel %d in file %s', Channel, ExtractedName); 
        title(P1,'Interpreter', 'none')
        NumZ = NumZ +1;
        SumZ = SumZ+z;
        end
    end
end

figure
x = (0:NumofFiles-1)*2;
y = 1:100;
z = SumZ*(1/length(FileName));
imagesc(x,y,10*log10(z))
ylabel('Frequency (Hz)')
xlabel('Time (min)')
c = colorbar;
ylabel(c, 'Fold Change')
colormap jet
ylim([0 100])
%caxis([0 20]);
ham = gca;
set(ham, 'YDir', 'normal')
set(findall(gca, '-property', 'FontSize'),'FontSize',24);
P1 = sprintf('Sonogram of Average CA3'); 
title(P1,'Interpreter', 'none')

%%
%Get the peak data information fof the last kainate for both CA1 and CA3
%CA1 peak data
LastCondition = LastKainate;
%Generate or reset the variables before analysis
AllQValue = []; AllPeakFreq = []; AllHBandwidth = [];
AverageQ  = []; AveragePF   = []; AverageHB     = [];
Counter = 1;
for File = 1:length(Files)
    %Reset the variables for the single brain section analysis
    SumQ  = 0;
    SumPF = 0;
    SumHB = 0;
    NumChannels = 0;
    %Make sure the current brain section has channels of interest
    if isfield(Files(File).ChannelID, 'CA1Channels') == 1
        for ChannelI = 1:length(Files(File).ChannelID.CA1Channels)
            %Extract the data of interest 
            PowerChanInterest = Files(File).Time(LastCondition).PowerAllCA1(ChannelI,1:100);
            FreqChanInterest  = Files(File).Time(LastCondition).Frequencies(1:100);
            %Calculate the desired Qvalue, Peak frequency and the
            %HalfBandwidth
            [QValue, PeakFreq,HalfBandwidth] = M_QvalueCalculator(PowerChanInterest, FreqChanInterest, 100, 15);
            AllQValue(Counter,1)   = QValue;
            AllPeakFreq(Counter,1) = PeakFreq;
            AllHBandwidth(Counter,1) = HalfBandwidth;
            Counter = Counter + 1;
            %Add the calculated values to determine the average of the
            %current brain section
            SumQ  = SumQ + QValue; 
            SumPF = SumPF + PeakFreq; 
            SumHB = SumHB + HalfBandwidth;
            NumChannels = NumChannels + 1;
        end
    end
    %Calculate the average of the current brain section data
    AverageQ(File,1)  = SumQ/NumChannels;
    AveragePF(File,1) = SumPF/NumChannels;
    AverageHB(File,1) = SumHB/NumChannels;
end
%%
% CA3 Analyzer
LastCondition = LastKainate;
%Generate or reset the variables before analysis
AllQValue = []; AllPeakFreq = []; AllHBandwidth = [];
AverageQ  = []; AveragePF   = []; AverageHB     = [];
Counter = 1;
for File = 1:length(Files)
    %Reset the variables for the single brain section analysis
    SumQ  = 0;
    SumPF = 0;
    SumHB = 0;
    NumChannels = 0;
    %Make sure the current brain section has channels of interest
    if isfield(Files(File).ChannelID, 'CA3Channels') == 1
        for ChannelI = 1:length(Files(File).ChannelID.CA3Channels)
            %Extract the data of interest 
            PowerChanInterest = Files(File).Time(LastCondition).PowerAllCA3(ChannelI,1:100);
            FreqChanInterest  = Files(File).Time(LastCondition).Frequencies(1:100);
            %Calculate the desired Qvalue, Peak frequency and the
            %HalfBandwidth
            [QValue, PeakFreq,HalfBandwidth] = M_QvalueCalculator(PowerChanInterest, FreqChanInterest, 100, 15);
            AllQValue(Counter,1)   = QValue;
            AllPeakFreq(Counter,1) = PeakFreq;
            AllHBandwidth(Counter,1) = HalfBandwidth;
            Counter = Counter + 1;
            %Add the calculated values to determine the average of the
            %current brain section
            SumQ  = SumQ + QValue; 
            SumPF = SumPF + PeakFreq; 
            SumHB = SumHB + HalfBandwidth;
            NumChannels = NumChannels + 1;
        end
    end
    %Calculate the average of the current brain section data
    AverageQ(File,1)  = SumQ/NumChannels;
    AveragePF(File,1) = SumPF/NumChannels;
    AverageHB(File,1) = SumHB/NumChannels;
end