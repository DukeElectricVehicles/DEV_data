%% motor model standalone
clear;
load ../noLoad/nonElectricalLosses
load ../spindown/chain+mag_sketchy
% mystery = mean(allMys(allR2>.35,:),1);
modelKv = 189;
modelLs = 0.160e-3 * 1.75;
modelRs = 0.186 * 1.34;
% modelRs = 0.186*1.5;
modelKt = 1./(modelKv*2*pi/60);
paramSets = table();
    
%     V = str2num(paramSets.voltage(i,:));
    
    [torqueVals, RPMvals] = meshgrid(linspace(0,1,300),linspace(000,4000,300));
    IVals = torqueVals / modelKt + 0.25;
    modelRsTmp = sqrt(modelRs.^2 + (modelLs*RPMvals*2 / 60 * 2 * pi).^2);
    V = RPMvals / modelKv + IVals.*modelRsTmp;
    modelRPM = RPMvals;
    modelTorque = torqueVals;
    
%     modelRsTmp = sqrt(modelRs.^2 + (modelLs*modelRPM*2 / 60 * 2 * pi).^2);
%     modelRPM = (V - IVals.*modelRsTmp) * modelKv;
%     modelRsTmp = sqrt(modelRs.^2 + (modelLs*modelRPM*2 / 60 * 2 * pi).^2);
%     modelRPM = (V - IVals.*modelRsTmp) * modelKv;
    
%     modelTorque = (IVals-0.25) * modelKt;
%     PvsERPM(1) = -PvsERPM(1);
    modelLosses(:,:,1) = IVals.^2.*modelRs;            % I2R
    modelLosses(:,:,2) = polyval(PvsERPM,modelRPM*2);  % nonelectrical
    modelLosses(:,:,3) = 0.3*ones(size(modelRPM));     % controller
    modelLosses(:,:,4) = 6e-3*modelTorque.^2.*modelRPM;% chain
    modelLosses(:,:,5) = 0.3*ones(size(modelRPM));     % chain
%     modelLosses = [IVals.^2.*modelRsTmp, polyval(PvsERPM,modelRPM*2), -polyval(CHAIN_LOSSES_POWER_OF_FLYWHEEL_RPM,modelRPM*14/72/2)];
%     modelLosses = [IVals.^2*paramSets.Rs(i), polyval(paramSets.Mys(i,:),modelRPM)];
    IVals = IVals + .03; % controller
    modelEff = max(0,1 - sum(modelLosses,3) ./ (IVals .* V));
    
    % impose limits
    modelEff(torqueVals.*RPMvals/60*2*pi > 349) = nan;
    modelEff(torqueVals > 1.01686) = nan;
    modelEff(RPMvals > 4200) = nan;
    
    figure(1);clf;
    colormap jet
%     pcolor(modelTorque',modelRPM',modelEff'); hold on;
%     shading interp;
    contourf(modelRPM',modelTorque',modelEff',200,'LineStyle','none'); hold on;
    [C,h] = contour(modelRPM',modelTorque',modelEff',[.1:.1:.7,.75,.8,.85:.01:.9],'LineColor','k');
    clabel(C,h,'LabelSpacing',200);
    ylabel('Torque (Nm)');
    xlabel('Speed (RPM)');
    legend('Efficiency');
    title('Motor Efficiency Map');
    
    figure(2);clf;
    colormap jet
    [C,h] = contour(modelTorque',modelEff',modelRPM',5);
    clabel(C,h);
    xlabel('Torque (Nm)');
    ylabel('Efficiency (\%)');
    legend('Speed (RPM)');
    
    figure(3);clf;
    colormap jet
    [C,h] = contour(modelRPM',modelEff',modelTorque',5);
    clabel(C,h);
    xlabel('Speed (Nm)');
    ylabel('Efficiency (\%)');
    legend('Torque (Nm)');
    
%     plot(modelTorque',modelEff');
    
%     legend(sprintfc('%d RPM',RPMvals(:,1)))