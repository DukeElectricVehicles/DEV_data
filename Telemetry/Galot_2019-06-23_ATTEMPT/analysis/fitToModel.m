%% model fit

getData

t = elapsed;
v = velo;
int_ = cumtrapz(t, v>0.02);
intV = cumtrapz(t, v);
intV2 = cumtrapz(t, v.^2);
int__ = cumtrapz(t, v>0.02);
intX = cumtrapz(t,dist);
intintV2 = 
intFind = @(coeffs,window) coeffs(1) .* int_(window(1):window(2)) + ...
                           coeffs(2) .* intV(window(1):window(2)) + ...
                           coeffs(3) .* intV2(window(1):window(2));
intF0 = @(coeffs,window) coeffs * [int_(window(1));...
                                   intV(window(1));...
                                   intV2(window(1));...
                                   1];
intFdef = @(coeffs,window) intFind(coeffs,window) - intF0(coeffs,window) + mass*v(window(1));

% dhdt = gradient(altRTK)./gradient(elapsed);
% grade = max(min(dhdt ./ velo,0.01),-0.01);
grade = gradient(altRTK) ./ gradient(dist);
grade(v < 1) = 0;
grade = max(min(grade, 0.01),-0.01);
grade(isnan(grade)) = 0;
grade(isinf(grade)) = 0;
dh = gradient(altRTK);

% effMotor = 1 - (current * 0.35) ./ 12; % IR / V

lossPoly_aeroAndBearing = [-1.06527e-08	-6.50352e-07	-1.02305e-04	1.12781e-03];
lossPoly_eddy = [-1.11897e-08	-5.30369e-06	-6.07395e-03	3.03073e-02] - lossPoly_aeroAndBearing;
wheelDia = 0.475;
omega = velo ./ (wheelDia / 2);
rpm = omega * 60 / (2*pi);
magLosses = polyval(lossPoly_eddy, rpm);

Fmotor = min((power - current.^2*.25 + magLosses) ./ velo,100);
Fmotor(isnan(Fmotor)) = 0;
Fmotor(isinf(Fmotor)) = 0;

%%
figure(1);clf;
coeffs = [0.0015*mass*9.81, 0, 0.0292, 0]
mv = @(coeffs,window) cumtrapz(t(window(1):window(2)), Fmotor(window(1):window(2))) - ...
                      intFdef(coeffs, window) - ...
                      mass*9.81*dh(window(1):window(2));
mx = @(coeffs,window) cumtrapz(t(window(1):window(2)), Fmotor(window(1):window(2))) - ...
                      intFdef(coeffs, window) - ...
                      mass*9.81*dh(window(1):window(2));
plot(mv(coeffs,[1,length(v)])); hold on;
plot(mass*v);

SSE = @(coeffs,window) sum((mv(coeffs,window) - mass*v(window(1):window(2))).^2);
coeffs = fminsearch(@(coeffs) SSE(coeffs,[1,length(v)]), coeffs)

plot(0*mv(coeffs,[1,length(v)]));

legend('model','mv','model2');

%%
spindownWindows = [];
spinningDown = false;
startInd = 1;
for i = 1:length(v)
    if (spinningDown)
        if (power(i) > 1 || v(i)<0.05)
            toAdd = [startInd,i];
            if (i-startInd > 100)
                spindownWindows = [spindownWindows; toAdd];
            end
            spinningDown = false;
        end
    else
        if (power(i) < 1 && v(i)>0.02)
            spinningDown = true;
            startInd = i;
        end
    end
end

allcoeffs = [];
allR2 = [];
for window = spindownWindows'
    [coeffs,sse,exitflag] = fmincon(@(coeffs) SSE(coeffs,window'), [0.0015*mass*9.81, 0, 0.0292, 0],...
        [],[],[],[],[0,0,1e-3,-1e9],[5,5,1,1e9]);
    mvaoeu = mass*v(window(1):window(2));
    SSR = sum((mvaoeu - mean(mvaoeu)).^2);
    coeffs
    R2 = 1-sse/SSR
    if (R2>0.95)
        allcoeffs = [allcoeffs;coeffs];
        allR2 = [allR2;R2];
        plot(window(1):window(2),mv(coeffs,window'),'k-');
    end
    fprintf('sse: %e\n',sse);
%     break;
end
allcoeffs,allR2
plot(power)
grid on;
ylim([-1,1]*max(mass*v)*2);

%%
figure(2);clf;
bar((allcoeffs ./ mean(allcoeffs,1))');

a = deleteoutliers(allcoeffs(:,1));
b = deleteoutliers(allcoeffs(:,2));
c = deleteoutliers(allcoeffs(:,3));

crr = mean(a) / mass / 9.81
CdA = mean(c) / 0.5 / 1.225
crrActual = 1.5e-3
CdAactual = 0.353*0.135