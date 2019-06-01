clear; clc; %close all;
figure(1);clf;figure(2);clf;figure(3);clf;figure(4);clf;figure(5);clf;

ACCEL_WINDOW = 1;
ROT_INERTIA = 0.8489;% + 0.00745;

load ../spindown/spindown_noRotor_may30_after

filesStruct = dir('*.txt');

legendShow = 'on';
for i = 1:numel(filesStruct)
    filename = filesStruct(i).name;
    clear linecolor
    if (contains(filename,'7,25mm'))
        linecolor = 'k';
    elseif (contains(filename,'7,50mm'))
        linecolor = 'b';
    elseif (contains(filename,'7,75mm'))
        linecolor = 'c';
    elseif (contains(filename,'8,00mm_0advance'))
        linecolor = 'g';
    elseif (contains(filename,'8,25mm'))
        linecolor = 'm';
    elseif (contains(filename,'8,50mm'))
        linecolor = 'r';
%     elseif (contains(filename,'_8advance'))
%         linecolor = 'r';
%     elseif (contains(filename,'_Cadvance'))
%         linecolor = 'g';
    end
    if (~exist('linecolor','var'))
        continue
    end
    filePath = strcat(filesStruct(i).folder, '/', filename);
    
    data = importdata(filePath);
    data = data(data(:,2)>.1,:); % current > .1
    badInds = find(data(:,4)>1000); % sometimes noise causes bad readings at start
    badInds = [1;badInds];
    startInd = badInds(end)+50;
    data = data(startInd:end-5,:);
    
    throttle = data(:, 5);
    time = data(:, 6) ./ 1000;
    voltage = data(:, 1);
    current = data(:, 2);
    rpm_fly = data(:, 4);
    
    for i = 1:length(rpm_fly) - 2%fix glitches in rpm readout
       if (rpm_fly(i) > 0) && (rpm_fly(i+2) > 0) && (rpm_fly(i+1) == 0)
           rpm_fly(i+1) = rpm_fly(i);
       end
    end
%     glitches = find(abs(diff(rpm_fly))>5);
%     for glitch = glitches'
%         rpm_fly(glitch+1) = rpm_fly(glitch);
%     end
    rpm_fly = 1./smooth(1./rpm_fly, 54);

%     rpm_fly = smooth(time, rpm_fly, 1001, 'sgolay', 5);
    rpm_motor = rpm_fly * 60/72;
    omega_fly = rpm_fly * 2 * pi / 60;

    ePower = voltage .* current;
    ePower = smooth(ePower, 54*72/60);
%     ePower = smooth(ePower, 500, 'sgolay');

    accel = gradient(omega_fly)./gradient(time);
    smooth(accel, 54);

%     accel = smooth(time, accel, 401, 'sgolay');
    accelInterp = fit(omega_fly, accel, 'smoothingspline', 'SmoothingParam', 0.99);
    accel = accelInterp(omega_fly);

    accelComp = accel - polyval(PARASITIC_LOSSES_ACC_OF_FLYWHEEL_RPS, omega_fly);

    torque = ROT_INERTIA .* accelComp;
    mPower = torque .* omega_fly;
    eff = mPower ./ ePower;
%     eff = smooth(eff, 101, 'sgolay');
    
    filename = strrep(filename,'_',' ');
    filename = strrep(filename,',','.');
    
    figure(1);
%     yyaxis left
    plot(rpm_motor, eff, [linecolor,'-'], 'DisplayName', filename); hold on;
%     yyaxis right
%     plot(rpmMotor, mPower, '-');
    
    figure(2);
    scatter3(rpm_motor, mPower, eff, [linecolor,'.'], 'DisplayName', filename); hold on;
    
    figure(3);
    ax1 = subplot(2,1,1);
    plot(rpm_motor, voltage, [linecolor,'.'], 'DisplayName', filename); hold on;
   	ax2 = subplot(2,1,2);
    plot(rpm_motor, current, [linecolor,'.'], 'DisplayName', filename); hold on;
    linkaxes([ax1,ax2],'x');
    
    figure(4);
    plot(rpm_motor, accel, linecolor, 'DisplayName',filename); hold on;
    
    currentPlot = smooth(current, 54*72/60);
    figure(5);
    yyaxis left
    plot(currentPlot, rpm_motor, [linecolor,'-'], 'DisplayName','RPM','HandleVisibility',legendShow); hold on;
    yyaxis right
    plot(currentPlot, torque/9.81*100, [linecolor,'^'], 'DisplayName','Torque (kgf.cm)','HandleVisibility',legendShow); hold on;
    plot(currentPlot, eff*100, [linecolor,'.'], 'DisplayName','Efficiency (%)','HandleVisibility',legendShow);
    legendShow = 'off';
end

figure(1);
legend(gca,'show','Location','South');
% yyaxis left
xlabel('RPM'); ylabel('efficiency'); title('efficiency vs speed (DEV Controller)');
grid on;
ylim([0.6, 1]);
% yyaxis right
% ylabel('Power');

figure(2);
legend(gca,'show');
xlabel('RPM'); ylabel('Power'); zlabel('Efficiency'); title('Efficiency Map (DEV Controller)');
grid on;
zlim([0.6, 1]);
ylim([0, 100]);
xlim([0, 300]);

figure(3);
subplot(2,1,1);
legend(gca,'show');
ylabel('Voltage'); title('Voltage and Current vs Speed (DEV Controller)');
ylim([0,20]);
subplot(2,1,2);
legend show
xlabel('RPM'); ylabel('Current');
ylim([0,20]);
grid on;

figure(4);
legend show
xlabel('RPM'); ylabel('Acceleration');
ylim([0,10]);
grid on;

%%
figure(5);
xlabel('Current'); title('Mitsuba datasheet graph (DEV Controller)'); grid on
xlim([0,18]);
legend show
yyaxis left
ylabel('Speed (RPM)'); ylim([0,500]);
yyaxis right
ylabel('Torque (kgf.cm) and Efficiency (%)'); ylim([0,100]);
rectangle('Position',[12.5,15,5,20],'FaceColor','w');
text(13, 30, '0 advance','Color','k','FontSize',16);
text(13, 25, '50 advance','Color','b','FontSize',16);
text(13, 20, '100 advance','Color','r','FontSize',16);