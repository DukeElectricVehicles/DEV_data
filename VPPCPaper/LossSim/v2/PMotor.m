% Yukai Qian
% Duke Electric Vehicles

function [PTOT, POUT, PLOSS] = PMotor(F, V, U)
% PMotor  Total power, output power, power loss and efficiency of motor.
%
%   PTOT = PMotor(F, V)
%   PTOT = PMotor(F, V, U)
%   [PTOT, POUT, PLOSS, ETA] = PMotor(___)
%
%   F       (N)     1-by-N vector of tractive forces.
%   V       (m/s)   1-by-N vector of horizontal velocities.
%   U       (V)     Motor voltage. 12 by default.
%   PTOT    (W)     1-by-N vector of total powers.
%   POUT    (W)     1-by-N vector of output powers.
%   PLOSS   (W)     1-by-N vector of power losses.
%
%   Assuming PLOSS equals resistance loss.

if nargin == 2
    U = 12;
end

r = 0.25; % Ohm, motor resistance

% Output power
POUT = F.*V;

% magnetic losses
lossPoly_aeroAndBearing = [-1.06527e-08	-6.50352e-07	-1.02305e-04	1.12781e-03];
lossPoly_eddy = [-1.11897e-08	-5.30369e-06	-6.07395e-03	3.03073e-02] - lossPoly_aeroAndBearing;
Pmag = max(abs(polyval(lossPoly_eddy, V/(0.475*pi)*60)),1);

% Current
i = (abs(POUT)+Pmag)/U;

% Power loss, total power, efficiency
PLOSS = r * i.^2 + Pmag;
PTOT = (POUT + PLOSS) .* (F>0);
