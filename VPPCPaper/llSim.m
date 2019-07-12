clear;clc;close all;

load('fcFit');
aux = 1 - .014;
sc = 1 - 0.002;
dcdc = 1 - 0.03;

averagePower = 20.7;
motorPower = 77;
motorDuty = averagePower / motorPower;
%NO LL-------------------------
simLen = 1e6;

p = zeros(simLen, 1);
idx = 1:length(p);
%p(1:round(simLen*motorDuty)) = motorPower; %no LL

%p = idx / simLen * 13; %13 watts ripple %SUPERCAP LL
p = p - mean(p) + averagePower; %set mean power

v = polyval(vFitC, p, [], vFitMu);
i = p ./ v;

eOut = sum(p);
eIn = eOut + sum(polyval(lossFitCoeffs, p));
eff = eOut / eIn;

figure;
plot(i);