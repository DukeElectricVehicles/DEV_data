%%  Patrick Grady
%   From patrick's pie chart script

% @param v - nominal race speed in m/s
% @param kappa - turn radius in meters

% all loss terms are in watts
%   *edit* total loss is in force

massCar = 26; %mass of car in kg
massDriver = 45; %mass of driver in kg
g = 9.807; %acceleration of gravity
d_wheel = 0.475; %diameter of the wheel in m

massTotal = massCar + massDriver;

totalLosses = {};
totalLossesLabels = {};

%air drag------------------------------------------------
airCD = 0.13; %coefficient of drag
airFrontal = 0.353; %frontal area in m^2
airCdA = 0.0438;
airDensity = 1.225; %density of air at 25C and standard pressure, kg/m^3

airForce = @(v,vw) 0.5 * airCdA * airDensity * (v+vw).^2;
airPower = @(v) airForce * v;

totalLosses{end+1} = @(args) airForce(args.v,args.vw);
totalLossesLabels{end+1} = 'External air drag';

%rolling resistance--------------------------------------
rrCoeff = 0.0015; %coefficient of rolling resistance
rrForce = @(v) massTotal * g * rrCoeff; %drag force of rolling resistance in newtons
rrPower = @(v) rrForce * v; %power loss of rolling resistance in watts

totalLosses{end+1} = @(args) rrForce(args.v);
totalLossesLabels{end+1} = 'Tire rolling resistance';

%cornering losses----------------------------------------
ca = 120; %tire cornering stiffness, newtons per degree
trackLength = 1947.1; %track length in meters. Galot raceway in Benson, NC

alpha = @(v, kappa) (massTotal .* v.^2 ./ kappa) ./ ca; %tire slip angle, degrees
corneringDragForce = @(v, kappa) ca .* alpha(v, kappa) .^ 2 .* pi ./ 180; %traction force needed

corneringPower = @(v, kappa) corneringDragForce(v, kappa) .* v;

totalLosses{end+1} = @(args) corneringDragForce(args.v, args.kappa);
totalLossesLabels{end+1} = 'Tire cornering losses';

%wheel air drag------------------------------------------
wheelDia = 0.475; %diameter of the wheel in m
wheelOmega = @(v) v / (wheelDia / 2);
wheelCdA = 1.1e-3;
wheelAirLoss = @(v) 3 * 0.5 * airDensity * v.^3 * wheelCdA; %three wheels in the car

totalLosses{end+1} = @(args) wheelAirLoss(args.v) ./ args.v;
totalLossesLabels{end+1} = 'Internal wheel air drag';
%bearing drag--------------------------------------------
Tb = 4.9e-3; %bearing frictional moment, Nm
bearingLoss = @(v) 3 * v * Tb / (wheelDia/2); %three wheels in the car

totalLosses{end+1} = @(args) bearingLoss(args.v) ./ args.v;
totalLossesLabels{end+1} = 'Wheel bearing drag';
%freewheel drag------------------------------------------


%motor losses--------------------------------------------
motorGearRatio = 120/14;
motorRPM = @(v,motorCurrent) wheelOmega(v) * motorGearRatio * 60/(2*pi) .*(motorCurrent>.5);
% motorCurrent = 5;
% motorVoltage = 16;
motorWindingResistance = 0.186*1.34;
motorWindingInductance = 0.160e-3 * 1.2;
motorKv = 189;
motorKt = 1./(motorKv*2*pi/60);
load nonElectricalLosses
PvsRPM = PvsERPM .* 2.^(3:-1:0);

motorTorque = @(motorCurrent) (motorCurrent-.25)*motorKt;

motorResistanceLoss = @(motorCurrent) motorWindingResistance * motorCurrent.^2;
motorNoLoadLoss = @(v,motorCurrent) polyval(PvsRPM,motorRPM(v,motorCurrent));
motorControllerLoss = 0.6;
motorFakeLoss = @(v,motorCurrent)  6e-3 * motorTorque(motorCurrent).^2.*motorRPM(v,motorCurrent);
motorTotalLoss = @(v,motorCurrent) ...
                        motorResistanceLoss(motorCurrent) + ...
                        motorNoLoadLoss(v,motorCurrent) + ...
                        motorControllerLoss + ...
                        motorFakeLoss(v,motorCurrent);

totalLosses{end+1} = @(args) motorTotalLoss(args.v,args.motorCurrent)./args.v;
totalLossesLabels{end+1} = 'Motor losses';

% motorEff = @(v) 1 - (motorTotalLoss(v)) / sum(totalLosses);
%chain losses--------------------------------------------

% chainPower_0 = 3.0;
% chainRPM_0 = 2600;
% 
% chainLoss = @(v,motorCurrent) chainPower_0 * motorRPM(v,motorCurrent) / chainRPM_0;
% totalLosses{end+1} = @(args) chainLoss(args.v,args.motorCurrent) ./ args.v;
% totalLossesLabels{end+1} = 'Chain losses';

% figure;
% 
% H = pie(totalLosses);
% T = H(strcmpi(get(H,'Type'),'text'));
% P = cell2mat(get(T,'Position'));
% set(T,{'Position'},num2cell(P*0.6,2));
% P(1, 1) = -1.7;%hack to make text better
% text(P(:,1),P(:,2),totalLossesLabels(:));

%fuel cell-----------------------------------------------
% h2Eff = 0.584;
% electricalPower = 23; %sum(totalLosses);
% fuelCellLoss = @(v) electricalPower / h2Eff - electricalPower;

totalLosses{end+1} = @(args) args.FCloss ./ args.v;
totalLossesLabels{end+1} = 'Fuel cell losses';

% %plotting and scorekeeping------------------------------
% figure;
% H = pie(totalLosses);
% T = H(strcmpi(get(H,'Type'),'text'));
% P = cell2mat(get(T,'Position'));
% set(T,{'Position'},num2cell(P*0.6,2));
% text(P(:,1),P(:,2),totalLossesLabels(:));
% 
% 
% totalPower = sum(totalLosses);
% electricScoreMetric = v ./ (electricalPower ./ 3600); %electric score in km per kWh
% electricScoreEnglish = electricScoreMetric ./ 1.609; %electric score in miles per kWh
% joulesPerLiterGas = 42.9e6 * 0.7646; %using constants from 2018 Eco-Marathon rules
% scoreMetric = joulesPerLiterGas .* v ./ (totalPower .* 1000); %hydrogen score in km per liter of gas
% scoreEnglish = scoreMetric ./ 1.609 .* 3.78541; % miles per gallon
% 
% %note, actual world record hydrogen score was 14,573 MPG, 6196 km/L
% 
% fprintf('Predicted electric score: %.1f mi/kWh, Hydrogen score: %.1f km/L\n', electricScoreEnglish, scoreMetric);

%%
