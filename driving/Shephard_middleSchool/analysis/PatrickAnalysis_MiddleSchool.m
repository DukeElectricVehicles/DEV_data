%% Patrick Grady
%  Nov. 2017

clear; clc; %close all;

FILENAME1 = 'spindowns_good_2019-06-05.TXT';
FILENAME2 = 'spindowns_good2_2019-06-05.TXT';
FILENAMEs = sprintfc('maxwell/%02d.TXT',1:3);
% FILENAME = '../Creepin-at-a-Middle-School/01.TXT';
% FILENAME2= '../Creepin-at-a-Middle-School/02.TXT';
% FILENAME3= '../Creepin-at-a-Middle-School/03.TXT';
ACCEL_WINDOW = 200;

linecolors = {'b','k','r'};

data = importdata(FILENAME1);
% data = data(18000:20500 , :);
data(:,end+1) = 1;
data2 = importdata(FILENAME2);
data2(:,end+1) = 2;
data = [data; data2];
for filename = FILENAMEs
    data2 = importdata(filename{1});
    data2(:,9) = [];
    data2(:,12:13) = [];
    data2(:,9) = data2(:,9) + data(end,9) + 1500*1000;
    data2(:,end+1) = 3;
    data = [data; data2];
end

voltage = data(:, 1);
current = data(:, 2);
power = data(:, 3);
velo = data(:, 4);
energy = data(:, 5);
dist = data(:, 6);
elapsed = data(:, 9) ./ 1000;
lat = data(:,10);
lon = data(:,11);
accel = zeros(size(velo));
dataID = data(:,end);

velo = smooth(velo, 21);
ke = 0.5 * (60.1 + 21.1 + 0.8) .* (velo .^2);

windowPoints = PatrickWindow(velo, power, elapsed);

% for i = ACCEL_WINDOW + 1: length(velo) - ACCEL_WINDOW
%    dv = velo(i + ACCEL_WINDOW) - velo(i - ACCEL_WINDOW);
%    dt = elapsed(i + ACCEL_WINDOW) - elapsed(i - ACCEL_WINDOW);
%    
%    accel(i) = dv / dt;
%    data(i, 11) = dv / dt;
% end
accel = gradient(velo)./gradient(elapsed);
assert(~any(isnan(accel)), 'some nan''s');
accel(isnan(accel)) = 0;
accel = smooth(accel, ACCEL_WINDOW);
data(:,7) = accel;

%% Plot start/stop lines--------------------------------------------
figure(1); clf;
plot(elapsed, velo); hold on;
for i = 1:size(windowPoints) 
    start = elapsed(windowPoints(i, 1));
    stop = elapsed(windowPoints(i, 2));
    
    line([start, start], [0, 10], 'Color', 'black', 'LineWidth', 3);
    line([stop, stop], [0, 10], 'Color', 'red', 'LineWidth', 3);
end
title('Start/Stop Lines');xlabel('time (s)');ylabel('velocity (m/s)');
ylim([-0.15 12]);
%% Plot decelerations------------------------------------------------
figure(2); clf;
veloAll=0;
accelAll=0;
windowTrim = 200;
for window = 1 : size(windowPoints)
   start = windowPoints(window, 1) + windowTrim;
   stop = windowPoints(window, 2) - windowTrim;
   
   plot(velo(start:stop), accel(start: stop), ...
       'Color',linecolors{dataID(start)}, ...
       'DisplayName',['Data ',num2str(window)]); hold on;
   ylim([-0.15 0]);
   
   veloAll = [veloAll;velo(start:stop)];
   accelAll = [accelAll;accel(start: stop)];
end

% Plot fit
%coeffs = [-0.00035 0 -0.025];
veloSweep = linspace(0, 10, 1000);
%decelFit = polyval(coeffs, veloSweep);
coeffs = polyfit(veloAll,accelAll,2);

f = @(a, x) a(1) + a(2).*x.^2;
fSSR = @(a, xm, ym) sum((ym - f(a,xm)).^2);
coeffs = fminsearch(@(a) fSSR(a,veloAll,accelAll), [0,0])
estCRR = -coeffs(1)/9.81
% estCD

St = sum(( accelAll - mean(accelAll) ).^2);
Sr = sum(( accelAll - f(coeffs,veloAll) ).^2);
r2 = (St - Sr) / St

decelFit = f(coeffs,veloSweep);%polyval(coeffs,veloSweep);
plot(veloSweep, decelFit,'DisplayName','fit');
title('Acceleration vs Velocity');
xlabel('Velocity (m/s)');ylabel('Acceleration(m/s^2)');
legend show

%% Plot decel vs dist
figure(3);clf;

for window = 1 : length(windowPoints)
   start = windowPoints(window, 1);
   stop = windowPoints(window, 2);
   
   if window == 1
       runData = data(start:stop, :);
       runData(:, 6) = runData(:, 6) - runData(1, 6) + 25;
   else
       runData = data(fliplr(start:stop), :);
       runData(:, 6) = -runData(:, 6) + runData(1, 6);
   end
   
   dist = runData(:, 6);
   accel = runData(:, 7);
   velo = runData(:, 4);
   lat = runData(:,10);
   lon = runData(:,11);
   
   subtractedAccel = accel - f(coeffs, velo);
   
   scatter(lat, lon); hold on;
%    yyaxis left;
%    plot(dist, subtractedAccel,'DisplayName',sprintf('window %d',window)); hold on;
%    %plot(dist, 0);
%    yyaxis right;
%    plot(dist, velo,'DisplayName',sprintf('window %d',window));
end
% xlim(mean(lat)+.1*[-1,1]);
% ylim(mean(lon)+.1*[-1,1]);
grid on
% yyaxis left; ylabel('Acceleration (m/s^2)'); %ylim([-10,10]);
% yyaxis right; ylabel('Velocity (m/s)'); ylim([0,10]);