%%  Gerry Chen
%   Jun 22

clear;

filename = '../flyingLaps1.TXT';

%% import
data = importdata(filename);

voltage = data(:, 1);
current = data(:, 2);
power = data(:, 3);
velo = data(:, 4);
energy = data(:, 5);
distTraveled = data(:, 6);
dpsCurrent = data(:, 9);
elapsed = data(:, 10) ./ 1000;
lat = data(:, 11); lat(lat==0) = nan; lat(isnan(lat)) = mean(lat,'omitnan');
lon = data(:, 12); lon(lon==0) = nan; lon(isnan(lon)) = mean(lon,'omitnan');
timeRaw = elapsed;
timeCum = timeRaw; 
indsInc = find(diff(timeCum)<0);
for i = indsInc'
    timeCum(i:end) = timeCum(i:end)+timeCum(i-1) + 100;
end

%% plot
figure(1);clf;
subplot(4,1,1);
plot(elapsed); ylabel('BMS time');
subplot(4,1,2);
plot(velo);    ylabel('velo');
subplot(4,1,3);
plot(power);   ylabel('power');
subplot(4,1,4);
plot(lat-mean(lat)); hold on;
plot(lon-mean(lon)); legend('lat','lon');
xlabel('index');

figure(2);clf;
scatter(lon,lat,3,timeCum);colorbar
title('Position');

%% get user click - from https://www.mathworks.com/matlabcentral/answers/263209-how-to-get-data-cursor-values-on-user-input
fprintf('Use the cursor to select the starting point then click spacebar when complete\n');
fprintf('Recommended to use the Power plot in Figure 1\n');
figure(1);
set(gcf,'CurrentCharacter',char(1));
h=datacursormode;
set(h,'DisplayStyle','datatip','SnapToData','off');
drawnow();
waitfor(gcf,'CurrentCharacter',char(32));
s = getCursorInfo(h);
indStart = fix(s.Position(1));

%% find points
dist = @(indStart,ind) sqrt((lat(ind)-lat(indStart)).^2 + (lon(ind)-lon(indStart)).^2);
figure(3);clf;
plot(indStart:length(lat),dist(indStart,indStart:length(lat))); hold on;

startLoc = [lat(indStart),lon(indStart)];
lapInds = indStart;
minDist = 1e99;
samePlace = true;
for i = indStart:length(timeCum)
    if (dist(indStart,i) < 2e-4);
        if (samePlace)
            if (dist(indStart,i) < minDist)
                minDist = dist(indStart,i);
                lapInds(end) = i;
            end
            continue;
        else
            lapInds = [lapInds,i];
            samePlace = true;
            minDist = dist(indStart,i);
        end
    else
        samePlace = false;
    end
%     if ((lat(i)-lat(indStart)^2 - (lon(i)-lon(lapInds(end)))^2 < .00001)
%     if (samePlace)
end
plot(lapInds,dist(indStart,lapInds),'r*');
drawnow();
fprintf('how many laps?\n');
numLaps = input('');
lapInds = lapInds(1:numLaps+1);
indEnd = lapInds(end);

%% find places with same velocity
perturbInds = fix(mean(diff(lapInds))/3);
perturbInds = -perturbInds:perturbInds;
veloSpline = spline(1:length(velo),velo);
[dun, duninds] = unique(distTraveled);
antiDistSpline = spline(dun, duninds);
indOffset = fminsearch(@(i) abs(ppval(veloSpline,antiDistSpline(i)+indStart)-ppval(veloSpline,antiDistSpline(i)+indEnd)), 0)
indOffset = fix(indOffset);
indOffset = antiDistSpline(i)
assert(abs(elapsed(indOffset+indStart)-elapsed(indStart))<30,'couldn''t find location with same velocity very nearby');
figure(2); hold on;
plot(lon(lapInds),lat(lapInds),'b*');
plot(lon(lapInds+indOffset),lat(lapInds+indOffset),'r*');
legend('path','selected pos','matching velocity pos');
velo(lapInds+indOffset)

%%
data = data(indStart+indOffset:indEnd+indOffset,:);
lapInds = lapInds-indStart+1;
newfilename = [filename(1:end-4),'_cut'];
save(newfilename,'data','lapInds');