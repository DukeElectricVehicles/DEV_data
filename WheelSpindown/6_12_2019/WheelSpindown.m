clear; clc; close all;

ACCEL_WINDOW = 32;
wheelRad = 0.475 / 2; %NOTE THIS IS FOR MICHELIN TIRE

I = 0.570 * 0.227^2; %Michelin
I = I + 0.020 * 0.23^2; %Sealant
I = I + 0.380 * 0.20^2; %DX32 rim
I = I + 0.170 * 0.13^2; %Spokes

filesStruct = dir('*s*.txt');

for i = 1:numel(filesStruct)
    filename = filesStruct(i).name;
    data = importdata(filename);
    meterspersec = data(:, 4);
    elapsed = data(:, 10) ./ 1000;
    
    velo = meterspersec / wheelRad;%to rad/sec
    %velo = smooth(velo, 21);
    dv = gradient(velo);
    dt = gradient(elapsed);
    accel = smooth(dv ./ dt, ACCEL_WINDOW);
    
    valid = accel < -0.002 & velo > 2 & velo < 50;%throw out extrema of data
    
    torque = accel * I;
    power = torque .* velo;
    linearVelo = velo * wheelRad;
    
    dragModel = @(k, v) k(1) + k(2) * v.^2;
    
    coeffs = lsqcurvefit(dragModel, [-0.02, -0.001], velo(valid), torque(valid), [], [], optimset('Display','off'));
    xspaced = linspace(0, 60, 1000);
    
    figure(1);
    plot(velo(valid), accel(valid), '.', 'DisplayName', filename); hold on;

    figure(2);
    plot(linearVelo, power, '.', 'DisplayName', filename); hold on;
    %plot(linearVelo, dragModel([coeffs(1), 0], velo) .* velo, 'DisplayName', 'k_c drag');
    %plot(linearVelo, dragModel([0, coeffs(2)], velo) .* velo, 'DisplayName', 'k_q drag');
    
    figure(3);
    
    s = sprintf('k_c=%0.2e k_q=%0.2e', coeffs(1), coeffs(2));
    plot(velo, torque, '.', 'DisplayName', filename); hold on;
    plot(xspaced, dragModel(coeffs, xspaced), '.', 'DisplayName', s); hold on;
end

figure(1);
xlabel('rad/sec');
ylabel('rad/sec^2');
legend(gca,'show');

figure(2);
xlabel('Car velocity in m/s');
ylabel('Power dissipation per wheel, watts')
legend(gca,'show');
ylim([-2 0]);
xlim([0 10]);
grid on;

figure(3);
xlabel('Wheel angular velocity in rad/sec');
ylabel('Drag torque in Nm');
legend(gca,'show');
xlim([0 35]);
ylim([-0.1 0]); grid on;