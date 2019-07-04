%% no load losses analysis

%% quiescent controller power
clear;
data = importdata('quiescentPower_KofordController.txt');
voltage = data(:, 1);
current = data(:, 2);
% p = polyfit(voltage,current,2); % current is very constant wrt voltage
% figure(1);clf;
% plot(voltage,current); hold on;
% vVals = linspace(0,30);
% plot(vVals,polyval(p,vVals));
Vq = mean(voltage);
Iq = mean(current);
Pq = mean(voltage.*current);

%% data
data = importdata('PSrampdown_KofordController.txt');
voltage = data(:, 1);
current = data(:, 2) - Iq;
power = data(:,3) - Pq;
rpm = data(:, 4);
time = data(:, 6) ./ 1000;
energy = data(:,9);
power = gradient(energy) ./ gradient(time);

for i = 1:length(rpm) - 2%fix glitches in rpm readout
   if (rpm(i) > 0) && (rpm(i+2) > 0) && (rpm(i+1) == 0)
       rpm(i+1) = rpm(i);
   end
   if (rpm(i+1) > 210)
       rpm(i+1) = rpm(i);
   end
end

rpm = smooth(rpm, 21);
rpm = rpm*48; % because forgot to change sprocket ticks to 1

current = power ./ voltage;
power = smooth(power,51) - current.^2*0.24924;

PvsERPM = polyfit(rpm,power./rpm,2); PvsERPM(end+1) = 0;
% PvsERPM = polyfit(rpm,power,3);
R2 = 1 - sum((power - polyval(PvsERPM,rpm)).^2) ./ sum((power-mean(power)).^2)
ERPMvals = linspace(0,max(rpm)*1.1);

% save('nonElectricalLosses','PvsERPM');

%% plot
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

figure(1);clf;
subplot(3,1,1);
plot(time,voltage);
subplot(3,1,2);
plot(time,current);
subplot(3,1,3);
plot(time,power);

figure(2);clf;
plot(rpm,power,'k.'); hold on;
plot(ERPMvals,polyval(PvsERPM,ERPMvals),'k-');
grid on;
title('Non-electrical Losses');
ylabel('Power (W)'); xlabel('Speed (eRPM)');
legend('data','regression','Location','SouthEast');
text(50*48,5,sprintf('$P = av^3 + bv^2 + cv$\n~~~$a=$%+.3e\n~~~$b=$%+.3e\n~~~$c=$%+.3e\n$R^2=%.4f$',...
    PvsERPM(1:end-1),R2),'FontSize',10)

% print -dpng noLoadLosses_Koford