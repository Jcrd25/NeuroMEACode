function PlotPeakDetection(Time, Voltage, LocOfPeaks)
%This funtction is meant to plot the detected peaks marked on the data to
%visually confirm that all desired peaks were detected and that no
%undesired peaks were accidentally tagged. 



figure
hold on
plot(Time,Voltage)
plot(Time(LocOfPeaks), Voltage(LocOfPeaks), 'rv', 'MarkerFaceColor', 'r')
xlabel('Time (ms)')
ylabel('Voltage (mV)')
title('Detected Peaks')
