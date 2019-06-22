%% Patrick Grady
%  Nov. 2017

clear; %clc; close all;

% filenames = sprintfc('../installationLap.TXT',0);
% filenames = sprintfc('../installationLap_cut.mat',0);
filenames = sprintfc('../flyinglaps1_cut.mat',0);

data = zeros(1,12);
for i = 1:length(filenames)
%     datanew = importdata(filenames{i});
    datanew = load(filenames{i});
    lapInds = datanew.lapInds;
    datanew = datanew.data;
    % data = data(16830:end, :);
    %datanew = [datanew(:, 1:6) datanew(:, 6:end)];%insert extra column because I messed up format
    datanew(:,10) = datanew(:,10) - datanew(1,10) + data(end,10) + 10000; % add time
    datanew(:,6) = datanew(:,6) - datanew(1,6) + data(end,6) + 1000; % add time
    data = [data; datanew];
end

data = data(2:end,:);

%%
%data = importdata('WRRun.TXT');
%data = data(6100:26460, :);

%data = importdata('thirdRuns.TXT');
%data = data(2680:end, :);

ACCEL_WINDOW = 100;

%CAR MODEL---------------------------
crr = 0.0015;
mass = (50 + 21);

densityAir = 1.225;
frontalArea = 0.353;
cd = 0.135;
rollingForce = crr * mass * 9.8;
airForce = @(v) 0.5 * cd * densityAir * frontalArea * v.^2;

accelModel = @(v) -(rollingForce + airForce(v)) / mass;

%END CAR MODEL-----------------------


voltage = data(:, 1);
current = data(:, 2);
power = data(:, 3);
velo = data(:, 4);
energy = data(:, 5);
dist = data(:, 6);
dpsCurrent = data(:, 9);
elapsed = data(:, 10) ./ 1000;
lat = data(:, 11);
lon = data(:, 12);
accel = zeros(size(velo));
deltaTE = zeros(size(velo));

dist = dist - dist(1);
energy = energy - energy(1);
elapsed = elapsed - elapsed(1);

wheelDia = .475;
kv = 26.6;
%kv = 25.5;
omega = velo ./ (wheelDia / 2);
omega = omega * 60 / (2*pi);
emf = omega / kv;

altRTK = lookupElev(lat, lon);
altRTK = smooth(altRTK, 11);
[x, y, z] = geodetic2ned(lat, lon, zeros(size(lat)), 35.322, -78.5111, 0, referenceEllipsoid('GRS80','m'));
x = -x;%idk why but it works
x(abs(x) > 3000) = 0;
y(abs(y) > 3000) = 0;

velo = smooth(velo, 11);
ke = 0.5 * mass .* (velo .^2);
pe = mass * 9.8 * altRTK;

te = ke + pe;
te = te - te(1);
mipkwhTEC = (dist ./ 1609) ./ ((energy-te) ./ 3.6e6);
mipkwh = (dist ./ 1609) ./ ((energy) ./ 3.6e6);


%windowPoints = PatrickWindow(velo, power, elapsed);
% windowPoints = [];

dv = gradient(velo);
dt = gradient(elapsed);
de = gradient(te);
accel = dv ./ dt;
deltaTE = de ./ dt;
accel = smooth(accel, ACCEL_WINDOW);
deltaTE = smooth(deltaTE, ACCEL_WINDOW);
% for i = ACCEL_WINDOW + 1: length(velo) - ACCEL_WINDOW
%    dv = velo(i + ACCEL_WINDOW) - velo(i - ACCEL_WINDOW);
%    dt = elapsed(i + ACCEL_WINDOW) - elapsed(i - ACCEL_WINDOW);
%    de = te(i + ACCEL_WINDOW) - te(i - ACCEL_WINDOW);
%    
%    accel(i) = dv / dt;
%    deltaTE(i) = de / dt;%power in watts
% end

accelComp = deltaTE ./ (velo * mass);

%% Plot start/stop lines--------------------------------------------
% figure(1); clf;
% plot(elapsed, velo); hold on;
% for i = 1:size(windowPoints, 1) 
%     start = elapsed(windowPoints(i, 1));
%     stop = elapsed(windowPoints(i, 2));
%     
%     line([start, start], [0, 10], 'Color', 'black', 'LineWidth', 3);
%     line([stop, stop], [0, 10], 'Color', 'red', 'LineWidth', 3);
% end
% title('Start/Stop Lines');xlabel('time (s)');ylabel('velocity (m/s)');
% ylim([0 10]);
% %accelComp = accelComp .* mass .* velo;
% %% Plot decelerations------------------------------------------------
% figure(2); clf;
% for window = 1 : size(windowPoints, 1)
%    start = windowPoints(window, 1);
%    stop = windowPoints(window, 2);
%    
%    %plot(velo(start:stop), accel(start: stop), 'o','MarkerSize',2); hold on;
%    yyaxis left
%    plot(velo(start:stop), accelComp(start: stop), 'o','MarkerSize',2,'DisplayName',sprintf('window %d',window)); hold on;
%    ylim([-0.2 0.1]); grid on;
%    yyaxis right
%    plot(velo(start:stop), smooth(deltaTE(start:stop), 50), 'o','MarkerSize',2,'DisplayName',sprintf('window %d',window));
%    ylim([-25,0]); grid on;
% end
% 
% % Plot fit
% yyaxis left
% veloSweep = linspace(0, 10, 1000);
% decelModel = accelModel(veloSweep);
% plot(veloSweep, decelModel, 'DisplayName','target');
% title('Acceleration vs Velocity');
% xlabel('Velocity (m/s)');ylabel('Acceleration(m/s^2)');
% legend show
% %legend('Uphill','Uphill TEC', 'Downhill', 'Downhill TEC', 'Drag model');
% 
% %% Plot decel vs map
% 
% % figure(3); clf;
% % for window = 1 : size(windowPoints, 1)
% %    
% %    start = windowPoints(window, 1);
% %    stop = windowPoints(window, 2);
% %    
% %     wLat = lat(start : stop);
% %     wLon = lon(start : stop);
% %     wAccel = accelComp(start : stop);
% %     wVelo = velo(start : stop);
% % 
% %     %scatter(wLat, wLon, 1, wAccel); %hold on;
% %     
% %     scatter3(wLat, wLon, wAccel, 3, 'filled'); hold on;
% %     %scatter3(wLat, wLon, accelModel(wVelo), 2, 'filled'); hold on;
% %     
% %     xlim([35.315 35.33]);
% %     ylim([-78.5135 -78.51]);
% %     zlim([-0.07 0]);
% % end

figure(1);clf;
scatter(lat,lon); hold on;
plot(lat(1),lon(1),'r*');
plot(lat(end),lon(end),'r*');

figure(4); clf;
p1 = plot(dist, mipkwhTEC);
ylim([0 1000]); ylabel('Score');
hold on; grid on;
plot(dist, mipkwh);
legend('Total Energy Compensation','raw');
fprintf('mipkwhTEC: %.2f\n',mipkwhTEC(end));

teC = energy - te; %total energy consumed
tpW = 20;
totalPower = gradient(teC)./gradient(elapsed);
totalPower = smooth(totalPower, tpW*2+1);

figure(5); clf;
a = subplot(3, 1, 1);
plot(elapsed, totalPower); hold on;
plot(elapsed, smooth(power, 51)); grid on;
ylabel('Power'); legend('Total Power','Motor power');
ylim([-50 100]);
b = subplot(3, 1, 2);
plot(elapsed, velo); grid on;
ylabel('Velocity');
c = subplot(3,1,3);
plot(elapsed, mipkwhTEC);
ylabel('Score (total energy compensated)');
linkaxes([a,b,c], 'x');
xlabel('time');

figure(6);
title('Power, velocity*10, motor power');
scatter3(x, y, totalPower); hold on;
scatter3(x, y, velo * 10);
scatter3(x, y, power);
zlim([-50 100]);

figure(7);
plot(current, voltage - emf, '.'); hold on; grid on;
xlabel('Current in A');
ylabel('BMS voltage - emf');
title('Motor Series Resistance');
%plot(emf);

%figure;
%plot(dist ./ elapsed);