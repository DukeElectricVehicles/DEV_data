clear; clc; %close all;
figure(1);clf;figure(2);clf;figure(3);clf;figure(4);clf;figure(5);clf;

global PARASITIC_LOSSES_ACC_OF_FLYWHEEL_RPS PARASITIC_LOSSES_POWER_OF_FLYWHEEL_RPM
global ROT_INERTIA

ACCEL_WINDOW = 1;
ROT_INERTIA = 0.8489;% + 0.00745;

load ../spindown/spindown_noRotor_jun1_after

filesStruct = dir('*.txt');

allPlotColors = {'k','b','c','g','y','r','m','k','b','c','g'};
% filenameFormat = 'PS(?<voltage>\d+)V_D(?<duty>[01].\d+)_\d\.txt';
% filenameFormat = 'PS(?<voltage>\d+)V_D(?<duty>[01].\d+)(?<mode>sync)?_\d\.txt';

ismemberstruct = @(A, B) arrayfun( @(x) isequal( B, x ), A );
allParameters = [];
for i = 1:numel(filesStruct)
    filename = replace(filesStruct(i).name,',','.');
    if (~contains(filename, '8mm') && ~contains(filename,'PS12V_D1.00'))
        continue
    end
%     stuff = regexp(filename,filenameFormat,'names');
    stuff = struct('name',filename(1:end-6));
    if (length(stuff)~=1)
        continue;
    end
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
%         stuff.mode = '';
%         allParameters = [allParameters, stuff];
    end
end

for i = 1:numel(filesStruct)
    filename = replace(filesStruct(i).name,',','.');
%     stuff = regexp(filename,filenameFormat,'names');
    stuff = struct('name',filename(1:end-6));
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
    
    analyzeSingle(filePath, linecolor, true)
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
rectangle('Position',[12.5,5,5,4*(1+length(allParameters))],'FaceColor','w');
for i = 1:length(allParameters)
    fields = fieldnames(allParameters(i));
    vals = cellfun(@(f) getfield(allParameters(i),f),fields,'UniformOutput',false);
    allP = {fields{:};vals{:}};
    text(13, 5+4*(length(allParameters)-i+1), ...
        sprintf('%s=%s\t', allP{:}),...
        'Color',allPlotColors{i},'FontSize',12,'FontName','FixedWidth');
end