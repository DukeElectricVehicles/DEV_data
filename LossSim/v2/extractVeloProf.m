%%  Gerry Chen
%   Extract velocity profile

% basic stuff
voltage = data(:, 1);
current = data(:, 2);
power = data(:, 3);
velo = data(:, 4);
energy = data(:, 5);
distTraveled = data(:, 6);
SCvoltage = data(:, 7);
dpsCurrent = data(:, 9);
elapsed = data(:, 10) ./ 1000;
% lat/lon
lat = data(:, 11); lat(lat==0) = nan; lat = deleteoutliers(lat,.1,1); %lat(isnan(lat)) = mean(lat,'omitnan');
lon = data(:, 12); lon(lon==0) = nan; lon = deleteoutliers(lon,.1,1); %lon(isnan(lon)) = mean(lon,'omitnan');
[x, y, z] = geodetic2ned(lat, lon, zeros(size(lat)), 35.322, -78.5111, 0, referenceEllipsoid('GRS80','m'));
assert(abs(mean(x,'omitnan') - mean(track.x)) < 0.01*range(x) , 'coordinate frames don''t match');
assert(abs(mean(y,'omitnan') - mean(track.y)) < 0.01*range(y) , 'coordinate frames don''t match');
x = spline(1:length(x),x,1:length(x));
y = spline(1:length(y),y,1:length(y));
% FC stuff
vFC = data(:, 7);
iFC = data(:, 8);
eFC = data(:, 9);
tempFC = data(:, 13);
pressFC = data(:, 14);
flow = data(:, 15);
totalFlow = data(:, 16);
instantEff = data(:, 17);
h2Energy = totalFlow .* 1000 .* 119.93;
h2Power = smooth(gradient(h2Energy)./gradient(elapsed),50);
% SCenergy = .5*191*SCvoltage.^2;
SCenergy = cumtrapz(elapsed,vFC.*iFC - power);
FCloss = smooth((gradient(h2Energy) - gradient(eFC)) ./ gradient(elapsed), 1000);

%% find points
dist = @(ind) sqrt((x(ind)-track.x(1)).^2 + (y(ind)-track.y(1)).^2);
figure(3);clf;
plot(dist(1:length(x))); hold on;

startLoc = [track.x(1),track.y(1)];
lapIndsOld = lapInds;
lapInds = [];
minDist = 1e99;
samePlace = true;
for i = 1:lapIndsOld(end)
    if (dist(i) < 5)
        if (samePlace)
            if (dist(i) < minDist)
                minDist = dist(i);
                lapInds(end) = i;
            end
            continue;
        else
            lapInds = [lapInds,i];
            samePlace = true;
            minDist = dist(i);
        end
    else
        samePlace = false;
    end
%     if ((lat(i)-lat(indStart)^2 - (lon(i)-lon(lapInds(end)))^2 < .00001)
%     if (samePlace)
end
plot(lapInds,dist(lapInds),'r*');
drawnow();
fprintf('how many laps?\n');
% numLaps = input('');
numLaps = length(lapInds)-1;fprintf('%d\n',numLaps);
lapInds = lapInds(1:numLaps+1);

%%
ds = sqrt(gradient(x).^2 + gradient(y).^2);
s = cumtrapz(ds);

xnew = zeros(length(track.x),numLaps);
ynew = zeros(length(track.y),numLaps);
snew = zeros(length(track.s),numLaps);
vnew = zeros(length(track.x),numLaps);
Inew = zeros(length(track.x),numLaps);
FClossnew = zeros(length(track.x),numLaps);
H2powernew = zeros(length(track.x),numLaps);
H2energynew = zeros(length(track.x),numLaps);
SCenergynew = zeros(length(track.x),numLaps);
[sun,iun] = unique(s);
for i = 1:numLaps
    snew(:,i) = track.s * ((s(lapInds(i+1))-s(lapInds(i))) / track.totalS) + s(lapInds(i));
    xnew(:,i) = spline(sun,x(iun), snew(:,i));
    ynew(:,i) = spline(sun,y(iun), snew(:,i));
    vnew(:,i) = spline(sun,velo(iun), snew(:,i));
    Inew(:,i) = spline(sun,current(iun), snew(:,i));
    FClossnew(:,i) = spline(sun,FCloss(iun), snew(:,i));
    H2powernew(:,i) = spline(sun,h2Power(iun), snew(:,i));
    H2energynew(:,i) = spline(sun,h2Energy(iun), snew(:,i));
    SCenergynew(:,i) = spline(sun,SCenergy(iun), snew(:,i));
end

figure(4);clf;
for i = 1:numLaps
    scatter(xnew(:,i),ynew(:,i),fix(Inew(:,i))+1,vnew(:,i)); hold on;
end

inds = 2:4;
v = mean(vnew(:,inds),2);
I = mean(Inew(:,inds),2);
FCloss = mean(FClossnew(:,inds),2);
H2power = mean(H2powernew(:,inds),2);
for i = 1:numLaps
    H2energygrad(:,i) = gradient(H2energynew(:,i));
    SCenergygrad(:,i) = gradient(SCenergynew(:,i));
end
H2energy = cumtrapz(mean(H2energygrad(:,inds),2));
SCenergy = cumtrapz(mean(SCenergygrad(:,inds),2));
% ind = 1;
% v = vnew(:,ind);
% I = Inew(:,ind);
% FCloss = FClossnew(:,ind);
% H2power = H2powernew(:,ind);
% H2energy = H2energynew(:,ind);
% SCenergy = SCenergynew(:,ind);