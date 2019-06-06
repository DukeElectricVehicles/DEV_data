%%  Gerry Chen
clear;

load 2019-06-05_spindown
FILENAME = 'motorCC2A_2019-06-05.TXT';

data = importdata(FILENAME);
mass = 20+67;

% endInd = find(data(:,9) < 2000);
% endInd(endInd<1000) = [];
% data = data(1:endInd,:);
data = data(2500:10500,:);

voltage = data(:, 1);
current = data(:, 2);
power = data(:, 3);
velo = data(:, 4);
energy = data(:, 5);
dist = data(:, 6);
time = data(:, 9) ./ 1000;
accel = zeros(size(velo));

accel = gradient(velo)./gradient(time);

motorAccel = accel - f(coeffs, velo);

mPower = motorAccel*mass;
ePower = voltage.*current;
mPower = smooth(mPower,50);
ePower = smooth(ePower,50);

figure(1);clf;
plot(time, voltage); ylabel('Voltage');
yyaxis right
plot(time, current); ylabel('Current');

figure(2);clf;
plot(time, mPower); hold on;
plot(time, ePower);
legend('mPower','ePower');
ylim([0,30]);