function [AUCPower] = M_AUC(Power, Frequency, Freq1, Freq2) 
%This function will estimate the Peak power at the desired frequency range
%using the area under the curve. The area under the curve is estimated by
%using MATLAB's trapz function. 
%Inputs
%   Power     - Power calculated by the FFT
%   Frequency - Frequency 
%   Freq1     - Lower end of the desired frequency range
%   Freq2     - Higher end of the desired frequency range
%Output
%   AUCPower - Area under the curve of the designated frequency ranges
%
%NOTE:
%Standard frequency ranges based on literature and multichannel systems 
%Delta : 0.5 - 4 Hz
%Theta : 4   - 8 Hz
%Alpha : 8   - 13 Hz
%Beta  : 13  - 30 Hz
%Gamma : 30  - 100 Hz
%
%AUTHOR: Jean Carlos Rodriguez Diaz 
%EMAIL:  jcrd@umich.edu
%WRITTEN: September 15, 2017
%LAST MODIFIED: Septermber 15, 2017
%v1.0

%Define the power and the frequency 
x = Frequency; 
y = (Power);
%Determine the indexes for the data set that focuses on the desired
%frequency range
LeftI  = find(x == Freq1);
RightI = find(x == Freq2);
%Calculate the area under the curve
AUCPower = trapz(x(LeftI:RightI), y(LeftI:RightI));
