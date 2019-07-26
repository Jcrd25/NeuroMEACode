function MEA_RasterPlotter(MEAData, ChannelsToPlot, ShowDetectedPeaks)
%This funciton is meant to calculate determine the spike times of sevaral
%channles and plot a RasterPlot.
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN: June7, 2017
%LAST MODIFIED: June 7, 2017
%v1.0



NumberChannels = size(ChannelsToPlot, 2); 
%Preallocate for speed and memory
DataRasterPlot = struct;


%Get the Spike data for each channel and save it to the struct
for ii = 1:NumberChannels
   Temp_PeakData = MEA_Spike_Detector(MEAData, ChannelsToPlot(ii), ShowDetectedPeaks);
   %Define hte parameters in the struct
   DataRasterPlot.PeakData(ii)    = Temp_PeakData;
      
end

 
%Plot Raster Plot of all channels
   figure
   hold on
   for ii = 1:length(ChannelsToPlot)
       ChannelNumberXAxis = (ones( 1, length(DataRasterPlot.PeakData(ii).SpikePeakTime)))*ii;
       scatter(DataRasterPlot.PeakData(ii).SpikePeakTime(:, :), ChannelNumberXAxis, '*k')
   end
   xlim([0 , DataRasterPlot.PeakData(1).Time(end)])
   xlabel('Time (s)')
   set(gca,'YMinorTick','on')
  hold off