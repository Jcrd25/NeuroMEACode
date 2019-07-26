%Cumulativ eProbability Code to be used with the Neuronal Culture Plotting
%Script
SpikeFreqThreshold  = 0.2*Culture(1).Time;
for DayID = 1:length(Condition.DIV)
    ISIs = [];
    for CultureID = 1:length(Culture)
        for Chan = 1:60
            if Culture(CultureID).DIV(DayID).PeakData(Chan).SpikeNum >= SpikeFreqThreshold
            ISIs = [ISIs Culture(CultureID).DIV(DayID).PeakData(Chan).ISI'];
            end
        end
    end
    Treatment(1).Prep(2).Day(DayID).ISIs = ISIs;
end

ISIExcel = [];
Day = 3;
for ii = 1:3
    if ii == 2
        ISIExcel = [ISIExcel Treatment.Prep(ii).Day(Day).ISIs];
    else
      ISIExcel = [ISIExcel Treatment.Prep(ii).Day(Day).ISIs];  
    end
end
[Exporty,Exportx] = ecdf(ISIExcel);