clear; close all; clc;


FILENAME = 'cleaned.txt';
%FILENAME = '../GalotOfficialJuly21_2018/officialDay2Run3.txt';
data = importdata(FILENAME);

data = data(46700:59900, :);
%data = data(6200:26760, :);

voltage = data(:, 1);
current = data(:, 2);
power = data(:, 3);
velo = data(:, 4) .* 8 ./ 9;
energy = data(:, 5);
dist = data(:, 6);
elapsed = data(:, 10) ./ 1000;

plot(velo); hold on;
plot(smooth(current, 41));

integratedDist = zeros(size(dist));

for i = 2 : length(dist)
    dt = elapsed(i) - elapsed(i - 1);
    integratedDist(i) = integratedDist(i-1) + velo(i) * dt;
end

figure;
plot(dist); hold on;
plot(integratedDist);

averageVelo = mean(velo)
rms = sqrt(mean(velo .^2))

rms / averageVelo

mean(velo .^2) / mean(velo).^2