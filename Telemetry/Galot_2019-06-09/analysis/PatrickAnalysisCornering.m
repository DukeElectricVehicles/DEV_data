%% Patrick Grady
%  Nov. 2017

clear; clc; close all;

% filenames = sprintfc('../spindowns%d.TXT',4);
% filenames = sprintfc('../racesim1.TXT',0);
filenames = sprintfc('../cornering1.TXT',0);

data = zeros(1,12);
for i = 1:length(filenames)
    datanew = importdata(filenames{i});
    % data = data(16830:end, :);
    datanew = [datanew(:, 1:6) datanew(:, 6:end)];%insert extra column because I messed up format
    datanew(:,10) = datanew(:,10) - datanew(1,10) + data(end,10) + 10000; % add time
    datanew(:,6) = datanew(:,6) - datanew(1,6) + data(end,6) + 1000; % add time
    data = [data; datanew];
end

%data = data(1693:13410,:);

%% racesim1
% data = data(8600:end,:);

% beginInd = find(data(:,4) > 6.7);
% beginInd = beginInd(1);
% endInd = 14440;
% data = data(beginInd:endInd,:);

%%
%data = importdata('WRRun.TXT');
%data = data(6100:26460, :);

%data = importdata('thirdRuns.TXT');
%data = data(2680:end, :);

DIST_WINDOW = 200;

%CAR MODEL---------------------------
%crr = 0.0015;
crr = 0.0015;
mass = (67 + 21);
ca = 130;
cornerRadius = DIST_WINDOW / pi / 2;

densityAir = 1.225;
frontalArea = 0.353;
cd = 0.135;
rollingForce = crr * mass * 9.8;
airForce = @(v) 0.5 * cd * densityAir * frontalArea * v.^2;
alpha = @(v) (mass .* v.^2 ./ cornerRadius) ./ ca; %tire slip angle, degrees
corneringDragForce = @(v) ca .* alpha(v) .^ 2 .* pi ./ 180; %traction force needed


accelModel = @(v) -(rollingForce + airForce(v) + corneringDragForce(v)) / mass;

%END CAR MODEL-----------------------


voltage = data(:, 1);
current = data(:, 2);
power = data(:, 3);
velo = data(:, 4);
energy = data(:, 5);
dist = data(:, 6);
elapsed = data(:, 10) ./ 1000;
lat = data(:, 11);
lon = data(:, 12);
accel = zeros(size(velo));

dist = dist - dist(1);
energy = energy - energy(1);
elapsed = elapsed - elapsed(1);


velo = smooth(velo, 11);
mipkwh = (dist ./ 1609) ./ (energy ./ 3.6e6);

windowPoints = PatrickWindow(velo, power, elapsed);
% windowPoints = [];

%dv = gradient(velo);
%dt = gradient(elapsed);
%de = gradient(te);
%accel = dv ./ dt;
%deltaTE = de ./ dt;
%accel = smooth(accel, ACCEL_WINDOW);
%deltaTE = smooth(deltaTE, ACCEL_WINDOW);
for i = 100:length(velo-100)
    curDist = dist(i);
    [resid, indFirst] = min(abs(dist + DIST_WINDOW / 2 - curDist));
    [resid, indLast] = min(abs(dist - DIST_WINDOW / 2 - curDist));
    
    
    dv = velo(indLast) - velo(indFirst);
    dt = elapsed(indLast) - elapsed(indFirst);
   
    accel(i) = dv / dt;
end

%% Plot start/stop lines--------------------------------------------
figure(1); clf;
plot(elapsed, velo); hold on;
for i = 1:size(windowPoints, 1) 
    start = elapsed(windowPoints(i, 1));
    stop = elapsed(windowPoints(i, 2));
    
    line([start, start], [0, 10], 'Color', 'black', 'LineWidth', 3);
    line([stop, stop], [0, 10], 'Color', 'red', 'LineWidth', 3);
end
title('Start/Stop Lines');xlabel('time (s)');ylabel('velocity (m/s)');
ylim([0 10]);
%accelComp = accelComp .* mass .* velo;
%% Plot decelerations------------------------------------------------
figure(2); clf;
for window = 1 : size(windowPoints, 1)
   start = windowPoints(window, 1);
   stop = windowPoints(window, 2);
   
   %plot(velo(start:stop), accel(start: stop), 'o','MarkerSize',2); hold on;
   yyaxis left
   plot(velo(start:stop), accel(start: stop), 'o','MarkerSize',2,'DisplayName',sprintf('window %d',window)); hold on;
   ylim([-0.1 0.0]); grid on;
   yyaxis right
   ylim([-25,0]); grid on;
end

% Plot fit
yyaxis left
veloSweep = linspace(0, 10, 1000);
decelModel = accelModel(veloSweep);
plot(veloSweep, decelModel, 'DisplayName','target');
title('Acceleration vs Velocity');
xlabel('Velocity (m/s)');ylabel('Acceleration(m/s^2)');
legend show
%legend('Uphill','Uphill TEC', 'Downhill', 'Downhill TEC', 'Drag model');

%% Plot decel vs map

figure(3); clf;
for window = 1 : size(windowPoints, 1)
   
   start = windowPoints(window, 1);
   stop = windowPoints(window, 2);
   
    wLat = lat(start : stop);
    wLon = lon(start : stop);
    wAccel = accel(start : stop);
    wVelo = velo(start : stop);

    %scatter(wLat, wLon, 1, wAccel); %hold on;
    
    scatter3(wLat, wLon, wAccel, 3, 'filled'); hold on;
    %scatter3(wLat, wLon, accelModel(wVelo), 2, 'filled'); hold on;
    
    xlim([35.315 35.33]);
    ylim([-78.5135 -78.51]);
    zlim([-0.07 0]);
end

figure;
plot(dist, accel);
