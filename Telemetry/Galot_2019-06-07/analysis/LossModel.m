function accel = LossModel(v, wind, dir, mass, crr, cdA)

    densityAir = 1.225;
    rollingForce = crr * mass * 9.8;
    airForce = @(v) 0.5 * cdA * densityAir * v.^2;

    lossPoly_aeroAndBearing = [-1.06527e-08	-6.50352e-07	-1.02305e-04	1.12781e-03];
    lossPoly_eddy = [-1.11897e-08	-5.30369e-06	-6.07395e-03	3.03073e-02] - lossPoly_aeroAndBearing;
    lossPoly_aeroAndBearing(end) = 0;
    lossPoly_eddy(end) = 0;
    PlossMag_W = @(v) -polyval(lossPoly_eddy, v / (.475/2) * 60/(2*pi));
    %PlossMech_W = @(v) -polyval(lossPoly_aeroAndBearing, v / (.475/2) * 60/(2*pi));
    forceMag = (PlossMag_W(v) ./ v);
    forceMag(isnan(forceMag)) = 0;
    
    accel = -(rollingForce + airForce(v + wind * dir) + forceMag) ./ mass;
end

