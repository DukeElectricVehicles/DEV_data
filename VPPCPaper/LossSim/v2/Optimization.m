% Yukai Qian
% Duke Electric Vehicles

clear;
% close all; format shortg
clearAllMemoizedCaches

tic

global track m totalS

% load 1906091056

%% Import data

track = load('Galot2.mat');
numSkip = 3;
for fieldname = fieldnames(track)'
    fieldname = fieldname{1};
    track.(fieldname) = track.(fieldname)(1:numSkip:end);
end
% track.phi = -track.phi;
track.z = smooth(track.z,5,'sgolay')';
dyaw = track.dyaw;
% r   = track.r;
% phi = track.phi;
s = track.s;
x   = track.x;
y   = -track.y;
z   = track.z;
% v = v(1:numSkip:end);
totalS = track.totalS;

[ds, slope] = DispIncr(s, z);
track.ds = ds;
track.slope = slope;

%% Prepare function for optimization

% Mass of car
mCar = 21;
mDriver = 50;
m = mCar + mDriver;

% Energy function
Energy = @(V) Eval(track, m, V);

% Initial guess
N = length(x);

% Constraints
Lower = zeros(1, N);
Upper = 10*ones(1, N);
tMax = 295; % Maximum time
pMax = 144; % Maximum power to/from motor

% Options
Options = optimoptions('fmincon', ...
                       'MaxFunctionEvaluations', 1e6, ...
                       'MaxIterations', 1e3, ...
                       'ConstraintTolerance', 0, ...
                       'Display', 'off');

%% Optimize
iters = 10;
allV = zeros([iters,length(x)]);
allE = zeros([iters,1]);
for restarts = 1:iters
    v = rand(size(x)).*2-1 + 6.7;
%     v = ones(size(x)) * 6.7;
    v0 = v;
    tic
    for i = 1:1
        A = [];
        b = [];
%         A = zeros(size(v));
%         A(1) = 1;
%         A(end) = -1;
%         b = 0;
        v = fmincon(Energy, v, [], [], A, b, ...
                    Lower, Upper, @(V) Constr(V, tMax, pMax), ...
                    Options);
    end
    toc
    plotStuff;
    allV(restarts,:) = v;
    allE(restarts) = Eval(track, m, v);
    fprintf('Iter %d: Lap energy = %.2fJ -> %.2f mi/kWh\n', ...
        restarts, allE(restarts), totalS*100/2.54/12/5280 / (allE(restarts)/3600/1000));
    
    pause(1);
end

%%
allE;
bestScore = 1.21 / (min(allE) / 3600 / 1000);
[v,ind] = min(allE);
v = allV(ind,:);
fprintf('best iteration: Lap energy = %.2fJ -> %.2f mi/kWh\n', ...
        allE(ind), bestScore);
plotStuff;

%% Constraint function

function [C, CEQ] = Constr(V, TMAX, PMAX)
    global track m
    
    % Time < TMAX,
    C = [Eval(track, m, V, 'time') - TMAX];
    CEQ = [];
end