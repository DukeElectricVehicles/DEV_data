clear; clc; % close all;

ROT_INERTIA = 0.8489;
HALLTEETH = 48;
FLYTEETH = 54;
MOTORTEETH = 72;

isRegen = 0;

load ../spindown/spindown_yesRotor_jun14_before % PARASITIC LOSSES
% data = importdata('./12V_t5A_-2A_FF_H_0.txt');
data = importdata('./12V5A_DCDC_FF_5khz_0.txt');

if (isRegen)
    data = data(data(:,2)<-.1,:); % current < -.1
else
    data = data(data(:,2)>.1,:); % current > .1
end
badInds = find(data(:,4)>1000); % sometimes noise causes bad readings at start
badInds = [1; badInds];
startInd = badInds(end)+50;
data = data(startInd:end-5,:);

time = data(:,6)/1000;
throttle = data(:, 5);
voltage = data(:, 1);
current = data(:, 2);
rpm = data(:, 4);

for i = 1:length(rpm) - 2%fix glitches in rpm readout
   if (rpm(i) > 0) && (rpm(i+2) > 0) && (rpm(i+1) == 0)
       rpm(i+1) = rpm(i);
   end
end
% glitches = find(abs(diff(rpm))>5);
% for glitch = glitches'
%     rpm(glitch+1) = rpm(glitch);%mean(rpm([glitch,glitch+2]));
% end
rpm = 1./smooth(1./rpm, HALLTEETH);
% voltage = smooth(voltage, 54*72/60);
% current = smooth(current, 54*72/60);

% rpm = smooth(time,rpm, 1001, 'sgolay', 5);

omega = rpm * 2 * pi / 60;
rpm_motor = rpm * FLYTEETH/MOTORTEETH;
velo_mph = rpm_motor*(18.7/12/5280*pi)*60;

ePower = voltage .* current;
ePower = smooth(ePower, fix(HALLTEETH * MOTORTEETH / FLYTEETH));
% ePower = smooth(time,ePower, 5);

accel = gradient(omega)./gradient(time);
accel = smooth(accel, 48);

% accel = smooth(rpm,accel, 20, 'sgolay');
accelOld = accel;
accelInterp = fit(omega, accel, 'smoothingspline', 'SmoothingParam', 0.99);
accel = accelInterp(omega);

accelComp = accel - polyval(PARASITIC_LOSSES_ACC_OF_FLYWHEEL_RPS, omega);


torque = ROT_INERTIA .* accelComp;
mPower = torque .* omega;
% mPower = smooth(time,mPower, 150, 'sgolay');

if (isRegen)
    eff = ePower ./ mPower;
else
    eff = mPower ./ ePower;
end
eff = smooth(time,eff, 104, 'sgolay');

%%
figure(1);clf;
subplot(3,1,1);
plot(time, rpm_motor);%, time, data(:,4)*54/72);
xlabel('Time (s)'); ylabel('Motor RPM');
subplot(3,1,2);
plot(rpm_motor, voltage, 'DisplayName','Voltage');
ylabel('Voltage (V)'); yyaxis right
plot(rpm_motor, current, 'DisplayName','Current');
ylabel('Current (A)'); xlabel('RPM of motor');
legend show
subplot(3,1,3);
plot(rpm_motor,ePower, 'DisplayName','Electrical Power');
hold on;
plot(rpm_motor,mPower, 'DisplayName','Mechanical Power');
if (isRegen)
    ylim([-100,0])
else
    ylim([0,100])
end
ylabel('Power (W)'); xlabel('RPM of motor');
legend show;
grid on;

figure(2);clf;
plot(rpm_motor, accelComp); hold on;
plot(rpm_motor, accelOld  - polyval(PARASITIC_LOSSES_ACC_OF_FLYWHEEL_RPS, omega));
xlabel('RPM'); ylabel('Acceleration');
if (isRegen)
    ylim([-10,0])
else
    ylim([0,10])
end
grid on;

figure(3);clf;
plot(rpm_motor, eff); hold on;
ylim([0.6, 1]); xlim([0,max(velo_mph)]);
ylabel('efficiency');
xlabel('speed (mph)');
grid on;
title('efficiency vs speed');

figure(4);clf;
plot(ePower, eff); hold on;
ylim([0.6, 1]);
if (isRegen)
    xlim([-100,0]);
else 
    xlim([0,100]);
end
ylabel('efficiency');
xlabel('Electrical power (W)');
grid on;
title('efficiency vs power');