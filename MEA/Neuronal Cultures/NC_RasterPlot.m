function NC_RasterPlot(MEANCData)
%Make rasterplot of the data.
Time = 300;
Fs = 20000;
figure
hold on
for chan = 1:60
    x = (MEANCData.PeakData(chan).SpikePeakLoc)/Fs;
    y = ones(size(x))*chan;
    plot(x, y, 'k.', 'MarkerSize', 6, 'LineStyle', 'none')
end 
ylim([0.5 60.5])
xlim([0 Time])
ylabel('Channel ID')
xlabel('Time (s)')
title('Raster Plot')
hold off







