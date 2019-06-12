clear; clc; %close all;
figure(1);clf;figure(2);clf;figure(3);clf;figure(4);clf;figure(5);clf;

ACCEL_WINDOW = 1;
ROT_INERTIA = 0.8489;% + 0.00745;

load ../spindown/spindown_noRotor_may30_after

filesStruct = dir('*.txt');

allPlotColors = {'k','b','c','g','y','r','m','k','b','c','g'};
filenameFormat = '8,00mm_(?<advance>-?\d+)advance_\d\.txt';

ismemberstruct = @(A, B) arrayfun( @(x) isequal( B, x ), A );
allParameters = [];
for i = 1:numel(filesStruct)
    filename = filesStruct(i).name;
    stuff = regexp(filename,filenameFormat,'names');
%     stuff = sscanf(filename,'8,00mm_%dadvance_%d.txt');
    if (length(stuff)~=1)
        continue;
    end
    if (str2num(stuff.advance)~=0)
        continue;
    end
    if (~any(ismemberstruct(allParameters,stuff)))
        allParameters = [allParameters, stuff];
    end
end

allRs = [];
allKv = [];

legendShow = 'on';
for i = 1:numel(filesStruct)
    filename = filesStruct(i).name;
    stuff = regexp(filename,filenameFormat,'names');
    try
        if (any(ismemberstruct(allParameters,stuff)))
            parameter = stuff;
        else
            continue
        end
    catch error
        continue
    end
    linecolor = allPlotColors{find(ismemberstruct(allParameters,stuff))};
    filePath = strcat(filesStruct(i).folder, '/', filesStruct(i).name);
    
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
    eff = smooth(eff, 501, 'sgolay');
    
    model = load('../MotorLossModel.mat');
    Padj = current.*voltage - ...
        model.PlossMag_W(rpm_motor) - ...
        model.PlossMech_W(rpm_motor) - ...
        model.PlossContr_W(voltage, 1, rpm_motor, 6e3);
    currentAdj = Padj ./ voltage;
    [CoRV,int,~,~,stats] = regress(currentAdj,[rpm_motor,voltage]);
    R2 = stats(1);
    assert (R2 > 0.95, sprintf('Kv calculation bad due to excessive nonlinearity: R2=%.4f\n',R2));
%     assert (abs(CoRV(3)) < 0, 'y intercept problem');
    Rs = 1/CoRV(2);
    Ke = -CoRV(1) * Rs;
    Kv = 1/Ke;
    
    fprintf('%s:\n',filename);
    fprintf('\tI = '); fprintf('(%.3f)RPM + (%.3f)V', CoRV);
    fprintf('(R2 = %.5f)\n', R2);
    fprintf('\tRs = %.6fohms\n\tKv = %.6f\n',Rs,Kv);
    
    allRs = [allRs; Rs];
    allKv = [allKv; Kv];
    
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
    plot(currentPlot(1:100:end), torque(1:100:end)/9.81*100, [linecolor,'^'], 'DisplayName','Torque (kgf.cm)','HandleVisibility',legendShow,'MarkerSize',5); hold on;
    plot(currentPlot, eff*100, [linecolor,'.'], 'DisplayName','Efficiency (%)','HandleVisibility',legendShow,'MarkerSize',1);
    legendShow = 'off';
end
fprintf('Kv = %.4f +/- %.4f\n', mean(allKv), std(allKv));
fprintf('Rs = %.4f +/- %.4f\n', mean(allRs), std(allRs));

%% plot model
model = load('../MotorLossModel.mat');
rpmVals = linspace(0,350,1000);
figure(1); plot(rpmVals,model.eff(12,1,rpmVals,6e3),'DisplayName','Motor Model');
figure(2); plot(rpmVals,model.Ptot_W(12,1,rpmVals,6e3),'DisplayName','Motor Model');
figure(3); plot(12.*ones(size(rpmVals)),model.Ptot_W(12,1,rpmVals,6e3)/12,'DisplayName','Motor Model');
figure(4); plot(rpmVals,model.torque_Nm(12,1,rpmVals,6e3)/ROT_INERTIA,'DisplayName','Motor Model');
figure(5);
    Ivals = model.Ptot_W(12,1,rpmVals,6e3)/12;
    yyaxis left;
    plot(Ivals,rpmVals,'DisplayName','Motor Model');
    yyaxis right;
    plot(Ivals, 100*model.eff(12,1,rpmVals,6e3),'DisplayName','Motor Model'); hold on;
    plot(Ivals, model.torque_Nm(12, 1, rpmVals, 6e3)/9.81*100,'DisplayName','Motor Model');

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

figure(5);
xlabel('Current'); title('Mitsuba datasheet graph (DEV Controller)'); grid on
xlim([0,18]);
legend show
yyaxis left
ylabel('Speed (RPM)'); ylim([0,500]);
yyaxis right
ylabel('Torque (kgf.cm) and Efficiency (%)'); ylim([0,100]);
rectangle('Position',[12.5,5,5,4*(1+length(allParameters))],'FaceColor','w');
for i = 1:length(allParameters)
    text(13, 5+4*(length(allParameters)-i+1), sprintf('%s advance',allParameters(i).advance),'Color',allPlotColors{i},'FontSize',12,'FontName','FixedWidth');
end