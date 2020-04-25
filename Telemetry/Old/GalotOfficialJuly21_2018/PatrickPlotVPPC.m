clear; clc; close all;

FILENAME = 'officialDay2Run3.TXT';

ACCEL_WINDOW = 100;

c = 190;

data = importdata(FILENAME);

data = data(4976:27315, :); %official day 2 run 3

vBMS = data(:, 1);
iBMS = data(:, 2);
energy = data(:, 5);
vFC = data(:, 7);
iFC = data(:, 8);
eFC = data(:, 9);
tempFC = data(:, 13);
pressFC = data(:, 14);
flow = data(:, 15);
totalFlow = data(:, 16);
instantEff = data(:, 17);

power = data(:, 3);
velo = data(:, 4);
dist = data(:, 6);
elapsed = data(:, 10) ./ 1000;
lat = data(:, 11);
lon = data(:, 12);


totalFlow = totalFlow - totalFlow(1);
eFC = eFC - eFC(1);
energy = energy - energy(1);
dist = dist - dist(1);
elapsed = elapsed - elapsed(1);
h2Energy = totalFlow .* 1000 .* 119.93;
instantEff = movmean(instantEff, 21);

capEnergy = 0.5 .* c .* (vFC.^2 - vFC(1).^2);
totalFCEff = (eFC) ./ h2Energy;
totalBMSEffComp = (energy + capEnergy) ./ h2Energy;


%velo = smooth(velo, 21);
ke = 0.5 * (26 + 40) .* (velo .^2);
miles = dist ./ 1609;
totalEnergy = energy - ke;

kwh = (energy) ./ (1000 * 3600);
kwhTE = totalEnergy ./ (1000 * 3600);

mpkwh = miles ./ kwh;
mpkwhTE = miles ./ kwhTE;

litersGasEq = (h2Energy - (capEnergy ./ 0.6)) ./ 1000 ./ (42.9 * 0.74262);

kmpl = dist ./ litersGasEq;

figure(1);
yyaxis left;
plot(elapsed, velo);
ylim([0 10]);
grid on;
ylabel('Velocity in m/s');
xlabel('Time in minutes');
yyaxis right;
plot(elapsed, power);
ylabel('Power in watts');
ylim([0 200]);

figure;
plot(mpkwh); hold on;
plot(mpkwhTE);
legend('no TE', 'TEC');
ylim([550 850]);

figure; clf;
plot(dist ./ elapsed); hold on;
plot(velo);
ylim([0 10]);

% figure(4);
% scatter3(lat, lon, 1:length(lat));
% xlim([35.315, 35.33])
% ylim([-78.513 -78.508])

figure;
plot(elapsed, instantEff, '.'); hold on;
plot(elapsed, (totalFCEff), '.'); 
%plot(elapsed, tempFC ./ 100);
plot(elapsed, totalBMSEffComp, '.');
ylim([0.5 0.65]); grid on;
legend('Instantaneous FC Efficiency', 'Integrated FC Efficiency', 'Total System Efficiency');
title('Fuel Cell Efficiency during the World Record Attempt');
xlabel('Time in seconds');
ylabel('Efficiency');

%figure;
%plot(eFC); hold on;
%plot(energy);
%plot((energy + capEnergy));

%figure;
%plot(smooth(pressFC, 100), smooth(instantEff, 21), 'o');

figure;
pFC = vFC .* iFC;
pFC = smooth(pFC, 11);
instantEff = smooth(instantEff, 11);
scatter(pFC, instantEff, 1); grid on; hold on;

figure;
[foo, startLap] = min(abs(dist - 2700));
[foo, endLap] = min(abs(dist - 4700));
lapTime = elapsed(startLap:endLap) - elapsed(startLap);
lapPower = power(startLap:endLap);
lapPower = smooth(lapPower, 5);

plot(lapTime, lapPower, 'LineWidth',1); grid on; hold on;
xlabel('Time, s');
ylabel('Power, W');
title('Motor power');

iSC = iBMS - iFC;
rSC = 10e-3;
pSC = iSC .^2 * rSC;
eSC = cumsum(pSC .* gradient(elapsed));
figure;
plot(elapsed, eSC ./ elapsed);

close all;
figure;
sub1 = subplot(2, 1, 1);
plot(elapsed, velo * 3.6, 'LineWidth',2.0); grid on;
ylim([0 30]);
ylabel('Velocity (kph)');
title('Velocity and Power during World Record Attempt');

sub2 = subplot(2, 1, 2);
plot(elapsed, power, 'Color', [1 0.3 0.0], 'LineWidth',2.0);
ylabel('Power (W)');
xlabel('Time (s)');
linkaxes([sub1, sub2], 'x'); grid on;