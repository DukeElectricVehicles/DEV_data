%%  Gerry Chen
%   PieChart2
%       Upgrade to the "pie chart" to generate itemized loss as a function
%       of position on the track

clear;

lossEquations

%% track stuff
load('Galot3.mat'); % stolen from Yukai's script
shiftAmt = 320; %320; %240;
track.dyaw = circshift(track.dyaw,shiftAmt);
track.yaw = circshift(track.yaw,shiftAmt);
track.x = circshift(track.x,shiftAmt);
track.y = circshift(track.y,shiftAmt);
track.z = circshift(track.z,shiftAmt);
track.s = cumtrapz(sqrt(gradient(track.x).^2+gradient(track.y).^2));
dhds = gradient(track.z) ./ gradient(track.s);
kappa = gradient(track.s) ./ abs(track.dyaw);
kappa = min(kappa,500);
% figure(1);clf; plot(kappa)

%% driver profile stuff
load('2018_officialDay2Run3_cut.mat');
extractVeloProf
assert(exist('v','var')==1,'didn''t set velocity profile');
assert(exist('I','var')==1,'didn''t set motor current profile');

%% calculate model
args = struct('v',v,...                     % velocity (m/s)
              'kappa',smooth(kappa,3),...   % radius of curvature (m)
              'motorCurrent',smooth(I,1),...% motor current (A)
              'FCloss',FCloss,...           % fuel cell loss (W)
              'vw',0*sin(track.yaw+pi/2)... % wind (m/s)
              );
totalLossesVals = zeros(length(track.s),size(totalLosses,2));
for i = 1:size(totalLossesVals,2)
    totalLossesVals(:,i) = totalLosses{i}(args);
end

%% TEC
totalEnergy = [-smooth(H2energy,15), ...
                smooth(SCenergy,1), ...
                smooth(track.z*massTotal*9.81,7), ...
                smooth(.5*massTotal*v.^2,9)];
for i = 1:size(totalEnergy,2)
    totalPower(:,i) = -gradient(totalEnergy(:,i))./(gradient(track.s)./v);
    totalPower(:,i) = smooth(totalPower(:,i),5);
end
totalPower(:,1) = circshift(totalPower(:,1),-8); % sketchy time offset
totalPower(:,2) = circshift(totalPower(:,2),3)+1.1;

%% plotting
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
order = [8,2,4,5,1,6,7,3];
assert(all(sort(order)==1:size(totalLossesVals,2)),'invalid order!');

figure(1);clf;
colormap jet
% dt = track.s ./ track.v
plot(cumtrapz(gradient(track.s)./v), smooth(sum(totalPower,2),9),'k-',...
    'LineWidth',1,'DisplayName','World record run data');
hold on;
area(cumtrapz(gradient(track.s)./v), totalLossesVals(:,order) .* v,'FaceColor','flat','FaceAlpha',0.75);
xlabel('Time (s)'); ylabel('Power (W)'); title('Power Loss Breakdown - Model and Experimental Data'); grid on;
legend('World record run data',totalLossesLabels{order})
xlim([0,trapz(gradient(track.s)./v)]); ylim([0,80]);
% annotation('textarrow',[.45,.57],[.75,.56],'String','4 turns per lap')
% annotation('arrow',[.455,.6],[.76,.68]);
% annotation('arrow',[.33,.245],[.75,.61]);
% annotation('arrow',[.325,.21],[.76,.62]);
figure(4);clf;
plot(totalPower)

figure(2);clf;
colormap jet
area(track.s, totalLossesVals(:,order),'FaceColor','flat','FaceAlpha',0.75);
xlabel('Distance along lap (m)'); ylabel('Force (N)'); title('Power Consumption Breakdown'); grid on;
legend(totalLossesLabels{order})
xlim([0,track.s(end)]);

%% statistics
totalPowerData = sum(totalPower,2);
totalPowerModel = sum(totalLossesVals,2).*v;
time = cumtrapz(gradient(track.s)./v);
EData = trapz(time,totalPowerData);
EModel = trapz(time,totalPowerModel);
error = abs(EModel - EData) ./ EData;
RMSerror = sqrt(mean((totalPowerData-totalPowerModel).^2));
fprintf('%% Error = %.1f%%\n',error*100);
fprintf('RMSerror = %.3fW\n',RMSerror);
fprintf('         = %.1f%%\n',RMSerror/mean(totalPowerData)*100);
38.3;
EModel / time(end)
EData / time(end)