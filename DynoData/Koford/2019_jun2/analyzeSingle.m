function [Rs, Kv, mysteryLosses, R2] = analyzeSingle(filePath, linecolor, toPlot, stuff)

    global PARASITIC_LOSSES_ACC_OF_FLYWHEEL_RPS PARASITIC_LOSSES_POWER_OF_FLYWHEEL_RPM
    global ROT_INERTIA
    persistent legendShow
    if (isempty(legendShow))
        legendShow = 'on';
    end
    
    filename = split(filePath, '/');
    filename = filename{end};
    
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
    rpm_fly(isnan(rpm_fly)) = 0;
    rpm_fly = 1./smooth(1./rpm_fly, 54);
    rpm_fly = smooth(rpm_fly, 150, 'sgolay');

%     rpm_fly = smooth(time, rpm_fly, 1001, 'sgolay', 5);
    rpm_motor = rpm_fly * 60/13;
    omega_fly = rpm_fly * 2 * pi / 60;

    ePower = voltage .* current;
    ePower = smooth(ePower, 54*14/60);
%     ePower = smooth(ePower, 500, 'sgolay');

    accel = gradient(omega_fly)./smooth(gradient(time),100,'sgolay');
    smooth(accel, 54);

%     accel = smooth(time, accel, 401, 'sgolay');
%     accel(isnan(accel)) = 0;
%     accelInterp = fit(omega_fly, accel, 'smoothingspline', 'SmoothingParam', 0.99);
%     accel = accelInterp(omega_fly);

    accelComp = accel - polyval(PARASITIC_LOSSES_ACC_OF_FLYWHEEL_RPS, omega_fly);

    torque = ROT_INERTIA .* accelComp;
    mPower = torque .* omega_fly;
    eff = mPower ./ ePower;
    eff = smooth(eff, 501, 'sgolay');
    
    currCutoff = max(find(current>18));
%     currOfRPM = polyfit(rpm_motor(currCutoff:end), current(currCutoff:end), 1);
%     Rs = mean(voltage(currCutoff:end))/currOfRPM(2);
%     Ke = -currOfRPM(1)*Rs;
%     Kv = 1/Ke;
%     R2 = corrcoef(current(currCutoff:end), polyval(currOfRPM, rpm_motor(currCutoff:end)));
%     R2 = R2(1,2)^2;
    [CoRV,int,~,~,stats] = regress(current,[rpm_motor,voltage]);
    R2lin = stats(1);
%     assert (R2 > 0.99, 'Kv calculation bad due to excessive nonlinearity\n');
%     assert (abs(CoRV(3)) < 0, 'y intercept problem');
    Rs = 1/CoRV(2);
    Ke = -CoRV(1) * Rs;
    Kv = 1/Ke;
    
    toFit = current<18 & current>0;
    y = ePower(toFit) - current(toFit).^2*Rs - mPower(toFit);
    mysteryLosses = polyfit(rpm_motor(toFit), y, 3);
%     y2 = y./rpm_motor(toFit); y2(isnan(y2)) = 0;
%     mysteryLosses = polyfit(rpm_motor(toFit), y2, 2); mysteryLosses(end+1) = 0;
    y = smooth(y,100);
    R2 = 1 - ...
        sum((polyval(mysteryLosses,rpm_motor(toFit)) - y).^2) / sum((y-mean(y)).^2);
    
    fprintf('%s:\n',filename);
    fprintf('\tI = '); fprintf('(%.3f)RPM + (%.3f)V\t', CoRV);
    fprintf('(R2 = %.5f)\n', R2lin);
    fprintf('\tRs = %.6fohms\n\tKv = %.6f\n',Rs,Kv);
    fprintf('\tR2 = %.6f\n',R2);
    
    filename = strrep(filename,'_',' ');
    filename = strrep(filename,',','.');
    
    if (toPlot)
        figure(1);
    %     yyaxis left
        rpmCut = rpm_motor>1150;
        if (filename(end-4)=='0')
            visi = 'on';
        else
            visi = 'off';
        end
        plot(rpm_motor(currCutoff:end), smooth(eff(currCutoff:end),800), ['-'], 'DisplayName', sprintf('$%d^\circ{}$ advance\n',fix(str2num(stuff.advance)*.15)),'Color',linecolor,'HandleVisibility',visi); hold on;
    %     yyaxis right
    %     plot(rpmMotor, mPower, '-');

        figure(2);
        plot3(rpm_motor, mPower, eff, ['.'], 'Color', linecolor, 'DisplayName', filename); hold on;

        figure(3);
        ax1 = subplot(2,1,1);
        plot(rpm_motor, voltage, ['.'], 'Color', linecolor, 'DisplayName', filename); hold on;
        ax2 = subplot(2,1,2);
        plot(rpm_motor, current, ['.'], 'Color', linecolor, 'DisplayName', filename); hold on;
        plot(rpm_motor, [rpm_motor,voltage]*CoRV, [':'], 'Color', linecolor, 'DisplayName', filename);
        linkaxes([ax1,ax2],'x');
    
        figure(4);
        plot(rpm_motor, accel, 'Color', linecolor, 'DisplayName',filename); hold on;

        currentPlot = smooth(current, 54*13/60);
        figure(5);
        yyaxis left
        currentPlot2 = smooth(currentPlot,300);
        plot(currentPlot2(1:50:end), rpm_motor(1:50:end), ['-'], 'Color', linecolor, 'DisplayName','RPM','HandleVisibility',legendShow); hold on;
        yyaxis right
        cToPlot = linspace(0,18,50);% +rand(1)*18/50;
        t1 = smooth(torque,50);
        torqueToPlot = spline(currentPlot(1:100:end),13/60*t1(1:100:end),cToPlot);
        torqueToPlot = deleteoutliers(torqueToPlot,[],1);
        plot(cToPlot, torqueToPlot*100, ['^'], 'Color', linecolor, 'DisplayName','Torque (N.cm)','HandleVisibility',legendShow,'MarkerSize',5); hold on;
%         plot(currentPlot(1:200:end), torque(1:200:end)*100, ['^'], 'Color', linecolor, 'DisplayName','Torque (N.cm)','HandleVisibility',legendShow,'MarkerSize',5); hold on;
        plot(currentPlot(1:10:end), eff(1:10:end)*100, ['.'], 'Color', linecolor, 'DisplayName','Efficiency (\%)','HandleVisibility',legendShow,'MarkerSize',1);
        
        if (filename(end-4)=='0')
            visi = 'on';
        else
            visi = 'off';
        end
        figure(6);
        plot(currentPlot, smooth(ePower - mPower,500), '-', 'Color',linecolor,'DisplayName', filename(3:5),'HandleVisibility',visi); hold on;
        
%         eff2 = 1 - (currentPlot(toFit).^2*Rs + polyval(mysteryLosses,rpm_motor(toFit)))./ePower(toFit);
%         plot(currentPlot(toFit), eff2*100, ['<'], 'Color', linecolor, 'DisplayName','Efficiency (%)','HandleVisibility',legendShow,'MarkerSize',10);
        legendShow = 'off';
    end
end