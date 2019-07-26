
%Neuronal Culture Master Script for spike analysis

%Identify the desired files to be analyzed.
 uiopen;    
 
 %%
Culture = 'MEA33075b';
DIV     = 14;
DIVSave = sprintf('DIV%d',DIV);
SaveDirectory = 'G:\My Drive\Jean Rodriguez Diaz\Data Analysis\NeuronalCulture\Shigeki Collab\Figures\';
SavingInfo.SaveDirectory = strcat(SaveDirectory, Culture,'\',DIVSave,'\');
mkdir(SavingInfo.SaveDirectory);
 
 tic
for ChanNum = 1:60
    %obtain the relevant data to generate the plots
    figure
    hold on
    Fs = MEAData.RawData.SamplingRate;
    %obtain the time and voltage trace
    x = 1:length(MEAData.RawData.VoltageTrace);
    x = x / (Fs);
    VoltageTrace = MEAData.RawData.VoltageTrace(:,ChanNum);
    
    %Parameters of the filter
    Fstop = 75;
    Fpass = 100;
    Astop = 65;
    Apass = 0.5;
    d = designfilt('highpassiir','StopbandFrequency',Fstop ,...
      'PassbandFrequency',Fpass,'StopbandAttenuation',Astop, ...
      'PassbandRipple',Apass,'SampleRate',Fs,'DesignMethod','butter');

    %Filter the Data
    VoltageTrace = filtfilt(d, VoltageTrace);
    
    %calculate the root mean square(used as the threshold for spike
    %detection)
    RMS_Signal = rms(VoltageTrace);
    Mean       = mean(VoltageTrace);
    %Set threshold to be 5 times the Root Mean Square
    SetThreshold = 5*RMS_Signal;
    %generate a horizontal line
    Line = ones(1,length(VoltageTrace));
    %plot the data
    plot(x,VoltageTrace)
    %plot the thresholds
    plot(x,Line*SetThreshold,'r')
    plot(x,Line*SetThreshold*-1,'r')
    %plot the mean
    plot(x,Line*Mean, 'b')
    %Figure and axis titles
    T1 = sprintf('%s Trace of Chan %d', Culture, ChanNum);
    title(T1)
    xlabel('Time (s)')
    ylabel('Voltage (µV)')
    hold off
end

for ii = 1:60
numb = ii;
figure(ii)
NAME = sprintf('%s Channel ID %d', Culture, numb);
FigureName = strcat(SavingInfo.SaveDirectory,'RawData_', NAME);
saveas(gcf,FigureName,'tiffn')
end
close all
%Spike plotting
for ChanNum = 1:60
    if MEANCData.PeakData(ChanNum).SpikeNum >= 1
        CTrace = [];
        figure
        hold on
        %obtain the trace of each spike
        for Spike =  1:MEANCData.PeakData(ChanNum).SpikeNum
  
                Peak = MEANCData.PeakData(ChanNum).SpikePeakLoc(Spike);
                leftedge = Peak - 0.001*Fs;
                rightedge = Peak + 0.003*Fs;
                %make sure that the trace before and after the spike can be
                %plotted
                if leftedge >= 1 
                    if rightedge <= length(MEAData.RawData.VoltageTrace(:,ChanNum))
                        Trace = MEAData.RawData.VoltageTrace(leftedge:rightedge,ChanNum);
                        time = rightedge-leftedge + 1;
                        plot((1:time)/Fs, Trace,'k')
                        CTrace(Spike,:) = Trace; 
                    else
                        P1 = sprintf('Not Plotted last spike of Chan %d',ChanNum);
                        disp(P1)
                    end
                else
                    P1 = sprintf('Not Plotted first spike of Chan %d',ChanNum);
                    disp(P1)
                end
            
        end
        %Plot the average spike
        plot((1:time)/Fs,mean(CTrace),'r','LineWidth', 6)
        %Figure and axis titles
        T1 = sprintf('Spike Traces of %s Chan %d', Culture, ChanNum);
        title(T1)
        xlabel('Time (s)')
        ylabel('Voltage (µV)')
    end
end
%Save the figures
for ii = 1:60
numb = ii;
figure(ii)
NAME = sprintf('%s Channel ID %d', Culture, numb);
FigureName = strcat(SavingInfo.SaveDirectory,'SpikeData_', NAME);
saveas(gcf,FigureName,'tiffn')
end
close all

NC_RasterPlot(MEANCData)
DIVName = sprintf('DIV%d', DIV);
FigureName = strcat(SavingInfo.SaveDirectory,'RasterPlot_', Culture,'_',DIVName);
saveas(gcf,FigureName,'tiffn')
close all
disp('Finished');

toc 