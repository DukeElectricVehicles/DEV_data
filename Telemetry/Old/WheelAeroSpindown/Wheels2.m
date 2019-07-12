clear; clc; close all;

TICKS_REV = 8;
ACCEL_WINDOW = 4;

filesStruct = dir('*.txt');

targetVelo = 6.7;
wheelRad = 0.475 / 2;

radSec = targetVelo / wheelRad;

for i = 1:numel(filesStruct)
    filename = filesStruct(i).name;
    
    data = importdata(filename);
    
    time = data ./ 1000;
    ticks = 1:length(time);
    radians = ticks  ./ TICKS_REV .* 2 * pi;
    
    velo = zeros(size(ticks));
    accel = zeros(size(ticks));
    
    for i = 1:length(ticks) - TICKS_REV
        i2 = i + TICKS_REV;
        velo(i) = (radians(i2) - radians(i)) / (time(i2) - time(i)); 
    end
    
    velo = smooth(velo, 21);
    
    for i = 1:length(ticks) - ACCEL_WINDOW
        i2 = i + ACCEL_WINDOW;
        accel(i) = (velo(i2) - velo(i)) / (time(i2) - time(i)); 
    end
    
    accel = smooth(accel, 21);
    
    I = 0.021;
    
    torque = accel * I;
    power = torque .* velo;
    
    linearVelo = velo * 0.475/2;
    
    valid = accel < -0.002 & velo > 5 & velo < 60;%throw out extrema of data

    dragModel = @(k, v) k(1) + k(2) * v.^2;    
    coeffs = lsqcurvefit(dragModel, [-0.02, -0.001], velo(valid), torque(valid), [], [], optimset('Display','off'));
    xspaced = linspace(0, 60, 1000);
    
    
    figure(1);
    plot(linearVelo, accel, '.', 'DisplayName', filename); hold on;

    figure(2);
    plot(linearVelo, power, '.', 'DisplayName', filename); hold on;
    
    figure(3);
    s = sprintf('k_c=%0.2e k_q=%0.2e', coeffs(1), coeffs(2));
    plot(velo(valid), torque(valid), '.', 'DisplayName', filename); hold on;
    plot(xspaced, dragModel(coeffs, xspaced), '-', 'DisplayName', s); hold on;
end

figure(1);
xlabel('rad/sec');
ylabel('rad/sec^2');
legend(gca,'show');
ylim([-2 0]);
xlim([0 10]);

figure(2);
xlabel('m/sec');
ylabel('watts')
legend(gca,'show');
ylim([-4 0]);
xlim([0 10]);

figure(3);
xlabel('Wheel angular velocity in rad/sec');
ylabel('Drag torque in N/m');
legend(gca,'show');
ylim([-0.1 0]);