function elev = lookupElev(lat,lon)

    map = importdata('galotMap.csv');
    
    mapLat = map(:, 1);
    mapLon = map(:, 2);
    mapElev = map(:, 3);
    
    mapElev = mapElev - mean(mapElev);
    
    elev = zeros(size(lat));
    
    for i = 1:length(lat)
        bestDist = 100;
        bestElev = 0;
        
        for j = 1:length(mapLat)
            dist = [lat(i) - mapLat(j); lon(i) - mapLon(j)];
            dist = norm(dist);
            
            if dist < bestDist
               bestDist = dist;
               bestElev = mapElev(j);
            end
        end
        
        elev(i) = bestElev;
    end
end

