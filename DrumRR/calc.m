clear;clc;close all;

noLoad = 11.4;
radAxle = 0.185;
radApp = 0.46;
g = 9.8;
armWeight = 2.11;
motorSpeed = 2500;%rpm
gearRatio = 14/120;
drumRad = 0.1524;
linearVelo = motorSpeed * gearRatio ./ 60 .* 2*pi * drumRad;

scale = [5.0, 8.1, 11.4, 11.4];%40 psi
power = [19.8, 23.7, 28.1, 27.0];

scale = [scale, 13.3, 2.2, 4.1, 7.0, 8.3, 10.9, 13.1, 13.1, 13.1, 13.1];%70 psi
power = [power, 23.6, 14.9, 16.5, 18.5, 19.4, 21.3, 23.4, 22.6, 23.0, 24.4];

power = power - noLoad;

fNormal = (armWeight + scale) * (radApp / radAxle) * g;
fDrag = power ./ linearVelo;
plot(fNormal, fDrag, '.'); hold on;

crrModel = @(crr, fn) crr .* fn;
coeffs = lsqcurvefit(crrModel, [0.01], fNormal, fDrag, [], [], optimset('Display','off'))

fSweep = linspace(0, 500, 1000);
plot(fSweep, crrModel(coeffs, fSweep));