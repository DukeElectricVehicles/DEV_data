% Patrick Grady
% analyze fuel cell testing data

clear; close all; clc;
%% import data

leakRate = 0.010; %mg per sec

energyDensityH2 = 119.93e3; %joules per g

data = importdata('base_IV2.txt');
%data = data(1:4400, :);

flow = data(:,1);%milligrams per second
power = data(:,2);
eff = data(:,3);
health = data(:,4);
voltage = data(:,5);
current = data(:,6);
time = data(:,8);
leak = data(:,7);
total = data(:,9);
temp = data(:,13);
pres = data(:,12);

power = smooth(power, 21);
flow = smooth(flow, 21);

h2energy = total* energyDensityH2;

h2power = flow .* energyDensityH2 / 1000;
leakPower = leakRate * energyDensityH2 / 1000;

voltageEff = voltage ./ 20 ./ 1.229;
totalEff = power ./ h2power;

flowByCurrent = current / 1.60217e-19 / 2 / 6.02214e23 * 2.01588 * 1e3 * 20; 
powerByCurrent = flowByCurrent * energyDensityH2 / 1000;


%plot(current, voltageEff); hold on;
%plot(power, power ./ (powerByCurrent + leakPower)); hold on;
%plot(power, power ./ (powerByCurrent));
cutPlotIdx = (power > 18) | (power < 16);
cutFitIdx = ((power > 18) | (power < 16)) & (power > 5);

fitCoeffs = polyfit(power(cutFitIdx), totalEff(cutFitIdx), 3);
powerSweep = linspace(10, 100, 1000);


scatter(power(cutPlotIdx), totalEff(cutPlotIdx), 20, 'filled'); hold on;
scatter(powerSweep, polyval(fitCoeffs, powerSweep), 10, 'filled'); hold on;


%plot(power, flow);
legend('Measured efficiency', 'Fit')

ylim([0.45 0.65]);
grid on;
xlabel("Output power, watts")
ylabel("Efficiency");
title("Horizon H-100 Stack-Only Efficiency")

%figure();
%plot(flowByCurrent + leakRate); hold on;
%plot(flow);