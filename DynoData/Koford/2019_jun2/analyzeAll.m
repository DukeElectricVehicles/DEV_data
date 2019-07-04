clear; clc; %close all;
clear analyzeSingle
figure(1);clf;figure(2);clf;figure(3);clf;figure(4);clf;figure(5);clf;figure(6);clf;

global PARASITIC_LOSSES_ACC_OF_FLYWHEEL_RPS PARASITIC_LOSSES_POWER_OF_FLYWHEEL_RPM
global ROT_INERTIA

ACCEL_WINDOW = 1;
ROT_INERTIA = 0.8489;% + 0.00745;

load ../spindown/spindown_noChain_jun2_before

filesStruct = dir('*.txt');

% filenameFormat = 'PS(?<voltage>\d+)V_D(?<duty>[01].\d+)_\d\.txt';
% filenameFormat = 'PS(?<voltage>\d+)V(?<comm>.*?)_D(?<duty>-?[01].\d+)(?<mode>sync)?_\d\.txt';
% filenameFormat = 'PS(?<voltage>\d+)V_D(?<duty>[01].\d+)(?<mode>sync)?_a(?<advance>-?\d+)_\d\.txt';
filenameFormat = 'PS(?<current>\d+)A(?<comm>.*?)_D(?<duty>[01].\d+)(?<mode>sync)?_\d\.txt';

ismemberstruct = @(A, B) arrayfun( @(x) isequal( B, x ), A );
allParameters = [];
% fileInds = zeros(length(filesStruct),1);
for i = 1:numel(filesStruct)
    filename = replace(filesStruct(i).name,',','.');
    stuff = regexp(filename,filenameFormat,'names');
    if (length(stuff)~=1)
        continue;
    end
    
    if (str2num(stuff.duty)~=1)
        continue;
    end
    if (length(stuff.comm)~=0)
        continue;
    end
%     if (str2num(stuff.duty)>=0)
%         continue;
%     end
%     if (abs(str2num(stuff.voltage)*str2num(stuff.duty) - 12) > .1 )
%         continue
%     end
%     if (contains(stuff.mode,'sync'))
%         continue;
%     end
%     if (~contains(stuff.mode,'sync') && (str2num(stuff.voltage)~=12) && (str2num(stuff.voltage)~=9))
%         continue
%     end
    if (~any(ismemberstruct(allParameters,stuff)))
        allParameters = [allParameters, stuff];
%         stuff.comm = '';
%         stuff.mode = '';
%         stuff.duty = sprintf('%.2f',-str2num(stuff.duty));
%         allParameters = [allParameters, stuff];
    end
end
allPlotColors = parula(length(allParameters)+2);

%%
allRs = [];
allKv = [];
allMys = [];
allR2 = [];
allChar = [];
for i = 1:numel(filesStruct)
    filename = replace(filesStruct(i).name,',','.');
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
    linecolor = allPlotColors(find(ismemberstruct(allParameters,stuff)),:);
    filePath = strcat(filesStruct(i).folder, '/', filesStruct(i).name);
    
    [Rs, Kv, mys, R2] = analyzeSingle(filePath, linecolor, true);
    allRs = [allRs; Rs];
    allKv = [allKv; Kv];
    allMys = [allMys; mys];
    allR2 = [allR2; R2];
    allChar = [allChar; find(ismemberstruct(allParameters,stuff))];
end

%%
motorModel
fprintf('Kv = %.4f +/- %.4f\n', mean(allKv), std(allKv));
fprintf('Rs = %.4f +/- %.4f\n', mean(allRs), std(allRs));

%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

figure(1);
legend(gca,'show','Location','SouthEast');
% yyaxis left
xlabel('RPM'); ylabel('Efficiency (\%)'); title('Efficiency vs Speed');
grid on;
ylim([0.75, 0.90]);
xlim([750,3000]);
% xlim([1000,4000]);
% yyaxis right
% ylabel('Power');

figure(2);
legend(gca,'show');
xlabel('RPM'); ylabel('Power'); zlabel('Efficiency'); title('Efficiency Map (DEV Controller)');
grid on;
zlim([0.6, 1]);
ylim([0, 300]);
xlim([0, 3500]);

figure(3);
subplot(2,1,1);
legend(gca,'show');
ylabel('Voltage'); title('Voltage and Current vs Speed (DEV Controller)');
ylim([0,20]);
subplot(2,1,2);
legend show
xlabel('RPM'); ylabel('Current');
ylim([0,20]);
xlim('auto');
grid on;

figure(4);
legend show
xlabel('RPM'); ylabel('Acceleration');
ylim([0,10]);
grid on;


figure(5);
xlabel('Current'); title('Speed, Torque, and Efficiency vs Current'); grid on
xlim([0,18]);
legend show
yyaxis left
ylabel('Speed (RPM)'); ylim([0,5000]);
yyaxis right
ylabel('Torque (N.cm) and Efficiency (\%)'); ylim([0,100]);
width = 4;
rectangle('Position',[12,5,width,2*(2+length(allParameters))],'FaceColor','w');
for i = 1:length(allParameters)
    fields = fieldnames(allParameters(i));
    vals = cellfun(@(f) getfield(allParameters(i),f),fields,'UniformOutput',false);
    fields = cellfun(@(f) f(1:end), fields,'UniformOutput',false);
    allP = {fields{[1]};vals{[1]}};
%     allP = {'advance',allParameters(i).advance};
    text(12+width/2, 5+2*(length(allParameters)-i+2), ...
        sprintf('%s=%sV', allP{:}),...
        'Color',allPlotColors(i,:),'FontSize',12,'FontName','FixedWidth','HorizontalAlignment','center');
end
text(12+width/2, 5+2*(1), 'motor model', ...
        'Color','r','FontSize',12,'FontName','FixedWidth','HorizontalAlignment','center');
    
figure(6);
xlabel('Current (A)'); ylabel('Power Loss (W)'); title('Power Loss vs Current'); grid on;
xlim([0,18]);ylim([0,100]);
legend('Location','NorthWest');