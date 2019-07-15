% Patrick Grady
% analyze fuel cell testing data

clear; close all; clc;
%% import data

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

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
%leakPower = leakRate * energyDensityH2 / 1000;

%voltageEff = voltage ./ 20 ./ 1.229;
totalEff = power ./ h2power;
lossPower = h2power - power;

%flowByCurrent = current / 1.60217e-19 / 2 / 6.02214e23 * 2.01588 * 1e3 * 20; 
%powerByCurrent = flowByCurrent * energyDensityH2 / 1000;


%plot(current, voltageEff); hold on;
%plot(power, power ./ (powerByCurrent + leakPower)); hold on;
%plot(power, power ./ (powerByCurrent));
cutPlotIdx = (power > 18) | (power < 16);
cutFitIdx = ((power > 18) | (power < 16)) & (power > 5) & (power < 40);

%effFitCoeffs = polyfit(power(cutFitIdx), totalEff(cutFitIdx), 3);
lossFitCoeffs = polyfit(power(cutFitIdx), lossPower(cutFitIdx), 2);
powerOutSweep = linspace(00, 100, 100);

vcutFitIdx = ((power > 18) | (power < 16)) & (power > 2) & (power < 100);
[vFitC, vFitS, vFitMu] = polyfit(power(vcutFitIdx), voltage(vcutFitIdx), 5);

effLossModelPin = powerOutSweep + polyval(lossFitCoeffs, powerOutSweep);
effLossModel = (powerOutSweep) ./ effLossModelPin;

yyaxis left
scatter(power(cutPlotIdx), totalEff(cutPlotIdx), 6, 'filled'); hold on;
%scatter(powerOutSweep, polyval(effFitCoeffs, powerOutSweep), 10, 'filled'); hold on;
plot(powerOutSweep, effLossModel); hold on;
ylim([0.45 0.75]);
ylabel("Efficiency");

yyaxis right
IVcutPlotIdx = ((power > 18) | (power < 16)) & (power > 1);
scatter(power(IVcutPlotIdx), voltage(IVcutPlotIdx), 6, 'filled'); hold on;
plot(powerOutSweep, polyval(vFitC, powerOutSweep, [], vFitMu));
ylim([8 19]);
ylabel("Voltage (V)");

%plot(power, flow);
legend('Measured', 'Fit')


grid on;
xlabel("Output power (W)")

title("Horizon H-100 Stack-Only Efficiency")

% figure;
% plot(power, lossPower); hold on;
% plot(powerOutSweep, polyval(lossFitCoeffs, powerOutSweep));

%figure();
%plot(flowByCurrent + leakRate); hold on;
%plot(flow);

save('fcFit','vFitC','vFitMu', 'lossFitCoeffs');