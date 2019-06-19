clear;

%% non-electrical losses
lossPoly_aeroAndBearing = [-1.06527e-08	-6.50352e-07	-1.02305e-04	1.12781e-03];
lossPoly_eddy = [-1.11897e-08	-5.30369e-06	-6.07395e-03	3.03073e-02] - lossPoly_aeroAndBearing;
lossPoly_aeroAndBearing(end) = 0;
lossPoly_eddy(end) = 0;

figure(1);clf;
rpmVals = linspace(0,350,1000);
plot(rpmVals, polyval(lossPoly_aeroAndBearing, rpmVals),...
    'DisplayName','Aero and Bearing Losses'); hold on;
plot(rpmVals, polyval(lossPoly_eddy, rpmVals), ...
    'DisplayName','Magnetic Losses');
plot(rpmVals, polyval(lossPoly_aeroAndBearing+lossPoly_eddy,rpmVals),...
    'DisplayName','Total non-electrical Losses');
legend show;
xlabel('motor mRPM'); ylabel('Power (W)'); title('Non-electrical Losses');

%% I2R losses
NUMPOLES = 16;
% Kv = 26.16; % from 314RPM
Kv = 25.3; % reasonable efficiency, also from Jun1 data
Rs = 0.25;
% Kv = 25.35; % these were extracted from 0advance data
% Rs = 0.245;
% Kv = 26.16; % use these for speed control
% Rs = 0.37;
L = 128.65e-6;
Vf = 0.3;
Kt_Nm_A = 1/(Kv/60*2*pi);
bemf_V =    @(rpm)          rpm/Kv;
Ibus_A =    @(Vbus,D,rpm)   (Vbus.*D - bemf_V(rpm))/Rs;
i_App =     @(Vbus,D,rpm,fs_Hz) (Vbus-bemf_V(rpm)).*(D./fs_Hz)/L;
Pelec_W =   @(Vbus,D,rpm)   Ibus_A(Vbus,D,rpm).*Vbus.*D;
body_W =    @(I, D)     (1-D) * I * Vf;
I2R_W =     @(Vbus,D,rpm)   -Ibus_A(Vbus,D,rpm).^2 * Rs .* sign(Ibus_A(Vbus,D,rpm));
I2Rs_W =    @(Vbus,D,rpm,fs_Hz) -(Ibus_A(Vbus,D,rpm).^2 + i_App(Vbus,D,rpm,fs_Hz).^2/12) * ...
                                    Rs .* sign(Ibus_A(Vbus,D,rpm));

%% all losses
PlossContr_W = @(Vbus,D,rpm,fs_Hz) 0.3.*ones(size(rpm));
PlossI2R_W = @(Vbus,D,rpm,fs_Hz) abs(I2Rs_W(Vbus,D,rpm,fs_Hz));
PlossMag_W = @(rpm) -polyval(lossPoly_eddy, rpm);
PlossMech_W = @(rpm) -polyval(lossPoly_aeroAndBearing, rpm);
PlossTot_W = @(Vbus,D,rpm,fs_Hz)    PlossContr_W(Vbus,D,rpm,fs_Hz) + ...
                                    abs(PlossI2R_W(Vbus,D,rpm,fs_Hz)) + ...
                                    0*abs(body_W(Ibus_A(Vbus,D,rpm), D)) + ...
                                    PlossMag_W(rpm) + ...
                                    PlossMech_W(rpm);
Ptot_W = @(Vbus, D, rpm, fs_Hz) Pelec_W(Vbus,D,rpm) + PlossContr_W(Vbus,D,rpm,fs_Hz);
eff = @(Vbus, D, rpm, fs_Hz) 1 - abs(PlossTot_W(Vbus,D,rpm,fs_Hz) ./ Ptot_W(Vbus,D,rpm,fs_Hz));

%% motor mechanical output power
Prot_W = @(Vbus, D, rpm, fs_Hz) Ptot_W(Vbus, D, rpm, fs_Hz)-PlossTot_W(Vbus, D, rpm, fs_Hz);
torque_Nm = @(Vbus, D, rpm, fs_Hz) Prot_W(Vbus, D, rpm, fs_Hz) ./ (rpm*2*pi/60);
torque_Nm = @(Vbus, D, rpm, fs_Hz) Kt_Nm_A * Ibus_A(Vbus,D,rpm);

save('MotorLossModel','PlossContr_W','PlossI2R_W','PlossMag_W','PlossMech_W','PlossTot_W','Ptot_W','eff','Prot_W','torque_Nm');

%% plot
clf(5);clf(6);
Ds = [1,.5];
Vbuss = 12 ./ Ds;
fs_Hz = 6e3;

for i = 1:length(Ds)
    D = Ds(i);
    Vbus = Vbuss(i);
    rpmVals = linspace(0,350,1000);
    figure(6);
    p1 = subplot(3,1,1);
    plot(rpmVals, PlossContr_W(Vbus,D,rpmVals,fs_Hz), 'DisplayName','Controller losses'); hold on;
    plot(rpmVals, PlossI2R_W(Vbus,D,rpmVals,fs_Hz), 'DisplayName','Electrical losses');
    plot(rpmVals, PlossMag_W(rpmVals), 'DisplayName','Magnetic losses');
    plot(rpmVals, PlossMech_W(rpmVals), 'DisplayName','Mechanical losses');
    plot(rpmVals, PlossTot_W(Vbus,D,rpmVals,fs_Hz), 'DisplayName','Total Losses');
    plot(rpmVals, Ptot_W(Vbus,D,rpmVals,fs_Hz), 'DisplayName','Total Power');
    legend show; grid on;
    ylabel('Power (W)'); xlabel('motor rpm');
    p2 = subplot(3,1,2);
    Ptot_W_Vals = abs(Ptot_W(Vbus,D,rpmVals,fs_Hz));
    plot(rpmVals, PlossContr_W(Vbus,D,rpmVals,fs_Hz)./Ptot_W_Vals, 'DisplayName','Controller losses'); hold on;
    plot(rpmVals, PlossI2R_W(Vbus,D,rpmVals,fs_Hz)./Ptot_W_Vals, 'DisplayName','Electrical losses');
    plot(rpmVals, PlossMag_W(rpmVals)./Ptot_W_Vals, 'DisplayName','Magnetic losses');
    plot(rpmVals, PlossMech_W(rpmVals)./Ptot_W_Vals, 'DisplayName','Mechanical losses');
    plot(rpmVals, PlossTot_W(Vbus,D,rpmVals,fs_Hz)./Ptot_W_Vals, 'DisplayName','Total Losses');
    plot(rpmVals, abs(Ptot_W(Vbus,D,rpmVals,fs_Hz))./Ptot_W_Vals, 'DisplayName','Total Power');
    legend('Location','SouthWest'); grid on;
    ylabel('proportion of total power'); xlabel('motor rpm');
    ylim([0,1]);
    p3 = subplot(3,1,3);
    plot(rpmVals, eff(Vbus,D,rpmVals,1e99),'DisplayName','no current ripple'); hold on;
    plot(rpmVals, eff(Vbus,D,rpmVals,fs_Hz),'DisplayName','current ripple');
    legend show; grid on;
    ylim([0.6,1]);

    figure(5); % mitsuba plot
    IVals = Ptot_W(Vbus, D, rpmVals, fs_Hz)/(Vbus*D);
    yyaxis left
    plot(IVals, rpmVals, 'DisplayName','RPM (motor)');
    ylim([0,500]);
    yyaxis right
    plot(IVals, 100*eff(Vbus,D,rpmVals,fs_Hz),'DisplayName','efficiency (%)'); hold on;
    plot(IVals, torque_Nm(Vbus, D, rpmVals, fs_Hz)/9.81*100,'DisplayName','torque (kg.cm)');
    legend show; grid on;
    ylim([0,100]); xlim([0,18]);
end

%% plot styles
linestyles = {'-','-.','--',':'};
a = figure(6);
linkaxes(a.Children(2:2:end),'x');
for ax = a.Children(2:2:end)'
    numPlotEl = length(ax.Children);
    colors = hsv(numPlotEl/i);
    for plotElInd = 1 : numPlotEl
        color = colors(mod(plotElInd,numPlotEl/i)+1,:);
        ax.Children(plotElInd).LineStyle = linestyles{fix(plotElInd/(numPlotEl/i))+1};
        ax.Children(plotElInd).Color = colors(mod(plotElInd,numPlotEl/i)+1,:);
        fprintf('%s: %d %d %d\n',ax.Children(plotElInd).DisplayName,ax.Children(plotElInd).Color);
    end
end
subplot(3,1,1);
legend show;
subplot(3,1,2);
legend show;
subplot(3,1,3);
legend show;