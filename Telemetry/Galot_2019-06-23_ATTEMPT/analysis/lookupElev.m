function [elev,enofix] = lookupElev(lat,lon)

    map = importdata('galotMap.csv');
    corrections = load('elevCorr.mat');
    
    mapLat = map(:, 1);
    mapLon = map(:, 2);
    mapElev = map(:, 3);
    
    mapElev = mapElev - mean(mapElev);
    
    mapLat = smooth(mapLat,10);
    mapLon = smooth(mapLon,10);
    mapElev = smooth(mapElev,10);
    
    E = scatteredInterpolant([mapLon, mapLat], mapElev,'nearest','nearest');
    elev = E(lon,lat);
    enofix = elev;
    
    Efix = scatteredInterpolant([corrections.lon, corrections.lat], corrections.elevCorr,'nearest','nearest');
    elev = elev + Efix(lon,lat);
    
%     figure(3);clf;
%     plot(elev); hold on;
%     plot(E(lon,lat));
%     plot(Efix(lon,lat));
    
%     figure(1);clf;
%     plot3(lon,lat,elev,'k.'); hold on;
%     plot3(mapLon,mapLat,mapElev,'r-');
% 
%     figure(2);clf;
%     [lonp,latp] = meshgrid(linspace(min(lon),max(lon)),linspace(min(lat),max(lat)));
%     surfc(lonp,latp,E(lonp,latp));
%     hold on;
%     plot3(lon,lat,elev,'k-');
%     plot3(mapLon,mapLat,mapElev,'r-');

%     elev = interp2(mapLat,mapLon,mapElev, lat,lon);
%     elev = zeros(size(lat));
%     
%     for i = 1:length(lat)
%         bestDist = 100;
%         bestElev = 0;
%         
%         for j = 1:length(mapLat)
%             dist = [lat(i) - mapLat(j); lon(i) - mapLon(j)];
%             dist = norm(dist);
%             
%             if dist < bestDist
%                bestDist = dist;
%                bestElev = mapElev(j);
%             end
%         end
%         
%         elev(i) = bestElev;
%     end
end

