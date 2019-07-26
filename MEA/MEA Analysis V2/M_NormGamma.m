function NormData = M_NormGamma(Baseline, Condition)
%This function is used to normalize the data to the percent increase
%compared to baseline
%
%Inputs
%   Baseline - Baseline recodring that is used to normalize the data
%   Condition - This is the desired part that iis going to be normalized
%
%Outputs
%   NormData - The Data for 0 - FreqMax Hz normalized to the baseline
%              recording
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN: January 23, 2018
%LAST MODIFIED: January 23, 2018
%v1.0

FreqMax = length(Baseline);
NormData = zeros(1, FreqMax);

for ii = 1:FreqMax
    NormData(ii) = (Condition(ii) - Baseline(ii))/Baseline(ii);
end