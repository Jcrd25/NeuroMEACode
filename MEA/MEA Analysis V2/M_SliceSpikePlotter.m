function [indData] = M_SliceSpikePlotter(SpikeData, PlotFields)


rows    = length(SpikeData.Time(1).(PlotFields));
columns = length(SpikeData.Time);
indData = zeros(rows, columns);
for NumberOfChan = 1:rows
  for Time = 1:columns
    indData(NumberOfChan,Time) = SpikeData.Time(Time).(PlotFields)(NumberOfChan);
  end
end