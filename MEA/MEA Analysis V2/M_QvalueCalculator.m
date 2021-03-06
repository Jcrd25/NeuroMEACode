function [QValue, PeakFreq, HalfBandwidth] = M_QvalueCalculator(Power,Frequency, FreqUpper, FreqLower)
%This function is used to calcuate a q value using the peak power and the
%half bandwidth. 
%the definition of this was taken from C.E. Lemercier et al. / Schizophrenia Research 188 (2017) 118�124
%Inputs
%   Power
%   Frequency
%   FreqUpper
%   FreqLower
%Output
%   QPower
%
%NOTE:
%Standard frequency ranges based on literature and multichannel systems 
%Delta : 0.5 - 4 Hz
%Theta : 4   - 8 Hz
%Alpha : 8   - 13 Hz
%Beta  : 13  - 30 Hz
%Gamma : 30  - 80 Hz
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN:
%LAST MODIFIED: January 3, 2018
%v1.0


%Find the peaks in the frequency range of interest

%Determine the indexed of the frequencies desired
LowerEnd  = find(Frequency == FreqLower);
UpperEnd  = find(Frequency == FreqUpper);
%Determine peak frequency of the Power region of interest
PowerInterest = Power(LowerEnd:UpperEnd);
[PeakFreqPower, PeakIndex] = max(PowerInterest);

%Determine half bandwidth

HalfPeakPower = PowerInterest(PeakIndex)/2;

%Find the two values that will be used to caluclate the bandwidth

 BandwidthPositions =[];
 IndexCont = 1;
 for ii = 1:length(PowerInterest)
     if PowerInterest(ii) > HalfPeakPower
         BandwidthPositions(IndexCont) = ii;
         IndexCont = IndexCont + 1;
     end
 end


%Determine the Half Bandwidth
HalfBandwithL = BandwidthPositions(1);
HalfBandwidthU = BandwidthPositions(end);

HalfBandwidth = HalfBandwidthU - HalfBandwithL;
%Determine the Preak Frequency
TempFreq = Frequency(LowerEnd:UpperEnd);
PeakFreq = TempFreq(PeakIndex);
%Calculate the Q value
QValue = PeakFreq/HalfBandwidth;
