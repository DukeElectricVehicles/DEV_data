%% process Model

if (exist('p1','var'))
    delete(p0);
    delete(p1);
    delete(p2);
    delete(p3);
    delete(p4);
    delete(p5);
    delete(p6);
end
load ../noLoad/nonElectricalLosses
load ../spindown/chain+mag_sketchy
% mystery = mean(allMys(allR2>.35,:),1);
modelKv = mean(allKv) * 1.01 / 1.01;
modelZs = mean(allRs)*1.1;
modelLs = 0.160e-3 * 1.2;
% modelRs = 0.186 * 1.45;
% modelLs = 0;
modelLs = 0.160e-3 * 1.75;
modelRs = mean(allRs) * 0.9;
modelRs = 0.186 * 1.34;
% modelRs = 0.186*1.5;
modelKt = 1./(modelKv*2*pi/60);
paramSets = table();
for i = 1:length(allParameters)
    paramSets.Kv(i) = mean(allKv(allChar==i));
    paramSets.Rs(i) = mean(allRs(allChar==i)) * 1.15;
    paramSets.Mys(i,:) = mean(allMys(allChar==i,:),1);
    paramSets.Kt(i) = 1./(paramSets.Kv(i)*2*pi/60);
    paramSets.VofI(i,:) = mean(allVofI(allChar==i,:),1);
    
%     fields = fieldnames(allParameters(i));
%     vals = cellfun(@(f) getfield(allParameters(i),f),fields,'UniformOutput',false);
    for field = fieldnames(allParameters(i))'
        field = field{1};
        paramSets.(field)(i,:) = getfield(allParameters(i),field);
    end
    
%     V = str2num(paramSets.voltage(i,:));
    IVals = linspace(0,18,10*99)';
    V = polyval(paramSets.VofI(i,:),IVals);
%     modelRPM = (V - IVals*paramSets.Rs(i)) * paramSets.Kv(i);
%     modelTorque = IVals * paramSets.Kt(i) * 1;
%     modelLosses = [IVals.^2*paramSets.Rs(i), polyval(PvsERPM,modelRPM*2)]; 
    modelRPM = (V - IVals*modelZs) * modelKv;   
    modelRsTmp = sqrt(modelRs.^2 + (modelLs*modelRPM*2 / 60 * 2 * pi).^2);
    modelRPM = (V - IVals.*modelRsTmp) * modelKv;
    modelRsTmp = sqrt(modelRs.^2 + (modelLs*modelRPM*2 / 60 * 2 * pi).^2);
    modelRPM = (V - IVals.*modelRsTmp) * modelKv;
    
    modelTorque = (IVals) * modelKt - polyval(PvsERPM,modelRPM*2)./(modelRPM/60*2*pi);
    modelLosses = [IVals.^2.*modelRs, ...           % I2R
                   polyval(PvsERPM,modelRPM*2), ... % nonelectrical
                   0.3*ones(size(modelRPM)), ...    % controller
                   6e-3*modelTorque.^2.*modelRPM, ...% chain
                   0.3*ones(size(modelRPM)) ...     % chain
                   ];               
%     modelLosses = [IVals.^2.*modelRsTmp, polyval(PvsERPM,modelRPM*2), -polyval(CHAIN_LOSSES_POWER_OF_FLYWHEEL_RPM,modelRPM*14/72/2)];
%     modelLosses = [IVals.^2*paramSets.Rs(i), polyval(paramSets.Mys(i,:),modelRPM)];
    IVals = IVals + .03; % controller
    modelEff = 1 - sum(modelLosses,2) ./ (IVals .* V);
    modelEff(:,2) = 1 - sum(modelLosses(:,1:3),2) ./ (IVals.*V);
    
%     linecolor = allPlotColors{mod(i-1,length(allPlotColors))+1};
    linecolor = 'r';
    
    nominalVoltage = str2num(paramSets.voltage(i,:));
    
    figure(1);
    p0(i,1)=plot(modelRPM, modelEff(:,1), '-', 'Color',linecolor,'DisplayName', ['model - ',num2str(nominalVoltage),'V']); hold on;
    
    figure(3);
    subplot(2,1,1);
    p5(i,1)=plot(modelRPM, V,'-','Color',linecolor','DisplayName', ['model - ',num2str(nominalVoltage),'V']);
    subplot(2,1,2);
    p6(i,1)=plot(modelRPM, IVals,'-','Color',linecolor','DisplayName', ['model - ',num2str(nominalVoltage),'V']);
    
    sel = fix((linspace(0,1,40).^2 * (size(modelEff,1)-1) + 1));
    figure(5);
    yyaxis left
    p1(i,1)=plot(IVals(sel), modelRPM(sel), [linecolor,'-'],'HandleVisibility','off'); hold on;
    yyaxis right
    p2(i,1)=plot(IVals(1:2:end), modelTorque(1:2:end)*100, [linecolor,'-'], 'MarkerSize',5,'HandleVisibility','off'); hold on;
    p3(i,1)=plot(IVals(sel), modelEff(sel,1)*100, [linecolor,'o'], 'LineWidth',1.5,'HandleVisibility','off');
    legendShow = 'off';
    
    figure(6);
    p4(i,1)=plot(IVals, sum(modelLosses,2), '-', 'Color',linecolor,'DisplayName', ['model - ',num2str(nominalVoltage),'V']); hold on;
    
    
    linecolor = 'g';
    figure(1);
    p0(i,2)=plot(modelRPM, modelEff(:,2), '-', 'Color',linecolor,'DisplayName', ['model - ',num2str(nominalVoltage),'V']); hold on;
    
    figure(3);
    subplot(2,1,1);
    p5(i,2)=plot(modelRPM, V,'-','Color',linecolor','DisplayName', ['model - ',num2str(nominalVoltage),'V']);
    subplot(2,1,2);
    p6(i,2)=plot(modelRPM, IVals,'-','Color',linecolor','DisplayName', ['model - ',num2str(nominalVoltage),'V']);
    
    figure(5);
    yyaxis left
    p1(i,2)=plot(IVals, modelRPM, [linecolor,'--'],'HandleVisibility','off'); hold on;
    yyaxis right
    p2(i,2)=plot(IVals(1:2:end), modelTorque(1:2:end)*100, [linecolor,'--'], 'MarkerSize',5,'HandleVisibility','off'); hold on;
    p3(i,2)=plot(IVals(sel), modelEff(sel,2)*100, [linecolor,'o'], 'LineWidth',1.5,'HandleVisibility','off');
    legendShow = 'off';
    
    figure(6);
    p4(i,2)=plot(IVals, sum(modelLosses,2), '-', 'Color',linecolor,'DisplayName', ['model - ',num2str(nominalVoltage),'V']); hold on;
end
paramSets.Kv(length(allParameters)+1) = modelKv;
paramSets.Rs(length(allParameters)+1) = modelRs;
paramSets.Kt(length(allParameters)+1) = modelKt;
for field = fieldnames(allParameters(i))'
    n = 'model';
    paramSets.(field{1})(i+1,:) = n(1:length(paramSets.(field{1})(1,:)));
end