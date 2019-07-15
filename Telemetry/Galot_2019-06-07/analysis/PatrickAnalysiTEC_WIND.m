%% Patrick Grady
%  Nov. 2017

clear; clc; close all;

filenames = sprintfc('../spindowns%d.TXT',4);
%filenames = sprintfc('../racesim1.TXT',0);
%filenames = sprintfc('../cornering1.TXT',0);

data = zeros(1,12);
for i = 1:length(filenames)
    datanew = importdata(filenames{i});
    % data = data(16830:end, :);
    datanew = [datanew(:, 1:6) datanew(:, 6:end)];%insert extra column because I messed up format
    datanew(:,10) = datanew(:,10) - datanew(1,10) + data(end,10) + 10000; % add time
    datanew(:,6) = datanew(:,6) - datanew(1,6) + data(end,6) + 1000; % add time
    data = [data; datanew];
end
    
%% racesim1
% data = data(8600:end,:);
% 
% beginInd = find(data(:,4) > 6.7);
% beginInd = beginInd(1);
% endInd = 14640;
% data = data(beginInd:endInd,:);

%%
% data = importdata('WRRun.TXT');
% data = data(6100:26460, :);

%data = importdata('thirdRuns.TXT');
%data = data(2680:end, :);

ACCEL_WINDOW = 50;

%CAR MODEL---------------------------
mass = (67 + 21);

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
deltaTE = zeros(size(velo));

dist = dist - dist(1);
energy = energy - energy(1);
elapsed = elapsed - elapsed(1);

altRTK = lookupElev(lat, lon);
%altRTK = smooth(altRTK, 11);
[x, y, z] = geodetic2ned(lat, lon, zeros(size(lat)), 35.322, -78.5111, 0, referenceEllipsoid('GRS80','m'));
x = -x;
x(abs(x) > 3000) = 0;
y(abs(y) > 3000) = 0;

%velo = smooth(velo, 11);
ke = 0.5 * mass .* (velo .^2);
pe = mass * 9.8 * altRTK;

te = ke + pe;
mipkwh = (dist ./ 1609) ./ (energy ./ 3.6e6);

windowPoints = PatrickWindow(velo, power, elapsed);
direction = zeros(size(windowPoints(:, 1)));
% windowPoints = [];

dv = gradient(velo);
dt = gradient(elapsed);
de = gradient(te);
accel = dv ./ dt;
deltaTE = de ./ dt;
accel = smooth(accel, ACCEL_WINDOW);
deltaTE = smooth(deltaTE, ACCEL_WINDOW);
accelComp = deltaTE ./ (velo * mass);

%% Plot start/stop lines--------------------------------------------
figure(1); clf;
plot(velo); hold on;
for i = 1:size(windowPoints, 1) 
    start = (windowPoints(i, 1));
    stop = (windowPoints(i, 2));
    direction(i) = x(start) > x(stop);
    
    line([start, start], [0, 10], 'Color', 'black', 'LineWidth', 3);
    line([stop, stop], [0, 10], 'Color', 'red', 'LineWidth', 3);
end
title('Start/Stop Lines');xlabel('ticks');ylabel('velocity (m/s)');
ylim([0 10]); grid on;

direction(direction == 0) = -1;

%% Plot decelerations------------------------------------------------
figure(2); clf;
for window = 1 : size(windowPoints, 1)
   start = windowPoints(window, 1);
   stop = windowPoints(window, 2);
   
   %plot(velo(start:stop), accel(start: stop), 'o','MarkerSize',2); hold on;
   
%    if direction(window) > 0
%        yyaxis left
%    else
%        yyaxis right
%    end
   plot(velo(start:stop), accelComp(start: stop), 'o','MarkerSize',1,'DisplayName',sprintf('window %d',window)); hold on;
   %ylim([-0.08 0.00]); grid on;
   %yyaxis right
   %plot(velo(start:stop), smooth(deltaTE(start:stop), 50), 'o','MarkerSize',2,'DisplayName',sprintf('window %d',window));
   %ylim([-25,0]); grid on;
end

% Plot fit
yyaxis left
veloSweep = linspace(0, 10, 1000);
decelModel = LossModel(veloSweep, 0, 0, mass, 0.0015, 0.033);
plot(veloSweep, decelModel, 'DisplayName','target');
title('Acceleration vs Velocity');
xlabel('Velocity (m/s)');ylabel('Acceleration(m/s^2)');
legend show
%legend('Uphill','Uphill TEC', 'Downhill', 'Downhill TEC', 'Drag model');

%SOLVE
x0 = [-4, 0.002, 0.033];
lossSpindown = @(coeffs) lossWind(coeffs, windowPoints, direction, velo, accel, mass);
dragCoeffs = fminsearch(lossSpindown, x0)
%dragCoeffs = x0;

decelModel = LossModel(veloSweep, dragCoeffs(1), 1, mass, dragCoeffs(2), dragCoeffs(3));
plot(veloSweep, decelModel, 'DisplayName','upwind');
decelModel = LossModel(veloSweep, dragCoeffs(1), -1, mass, dragCoeffs(2), dragCoeffs(3));
plot(veloSweep, decelModel, 'DisplayName','downwind');

%% Plot decel vs map

figure(3); clf;
for window = 1 : size(windowPoints, 1)
   
   start = windowPoints(window, 1);
   stop = windowPoints(window, 2);
   
    wLat = x(start : stop);
    wLon = y(start : stop);
    wAccel = accelComp(start : stop);
    wVelo = velo(start : stop);

    %scatter(wLat, wLon, 1, wAccel); %hold on;
    
    scatter3(wLat, wLon, wAccel, 3, 'filled'); hold on;
    %scatter3(wLat, wLon, accelModel(wVelo), 2, 'filled'); hold on;
    
    %xlim([35.315 35.33]);
    %ylim([-78.5135 -78.51]);
    zlim([-0.07 0]);
    xlabel('X m'); ylabel('Y m');zlabel('Accel');
end

function loss = lossWind(coeffs, windowPoints, direction, velo, accel, mass)
    loss = 0;
    for window = 1 : size(windowPoints, 1)
        start = windowPoints(window, 1);
        stop = windowPoints(window, 2);
        veloCut = velo(start:stop);
        accelCut = accel(start: stop);
        
        accelPred = LossModel(veloCut, coeffs(1), direction(window), mass, coeffs(2), coeffs(3));
        resid = (accelPred - accelCut);
        
        %resid(abs(resid) > 0.03) = 0;
        
        loss = loss + sum(resid.^2);
    end
end

