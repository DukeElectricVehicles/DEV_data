%%  Gerry Chen
%   Jun 22

clear;

filename = '../installationLap.TXT';

%% import
data = importdata(filename);

voltage = data(:, 1);
current = data(:, 2);
power = data(:, 3);
velo = data(:, 4);
energy = data(:, 5);
dist = data(:, 6);
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
subplot(3,1,1);
plot(elapsed); ylabel('BMS time');
subplot(3,1,2);
plot(velo);    ylabel('velo');
subplot(3,1,3);
plot(power);   ylabel('power');
xlabel('index');

figure(2);clf;
subplot(2,1,1);
plot(lat-mean(lat)); hold on;
plot(lon-mean(lon)); legend('lat','lon');
subplot(2,1,2);
scatter(lat,lon,3,timeCum);colorbar
title('Position');

%% get user click - from https://www.mathworks.com/matlabcentral/answers/263209-how-to-get-data-cursor-values-on-user-input
fprintf('Use the cursor to select the starting point then click spacebar when complete\n');
fprintf('Recommended to use the Power plot in Figure 1\n');
figure(1);
set(gcf,'CurrentCharacter',char(1));
h=datacursormode;
set(h,'DisplayStyle','datatip','SnapToData','off');
waitfor(gcf,'CurrentCharacter',char(32));
s = getCursorInfo(h);
indStart = s.Position(1);

%% find points
dist = @(indStart,ind) sqrt((lat(ind)-lat(indStart)).^2 + (lon(ind)-lon(indStart)).^2)
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

%%
data = data(indStart:lapInds(1));
newfilename = [filename(1:end-4),'new'];
save(filename,'data');