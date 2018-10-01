clc;
clear;
close all;

%% Parameters
S0 = 14;
K = 15; % Strike price
r = 0.1;% rate
sigma = 0.25; 
timesteps = 1000;
T = 0.5;
simulations = 100000;

%Fixed strike price
[V0fixed, errorFixedV, Z0fixed, errorFixedZ] = asianOptionMCfixed(S0,K,sigma,r,T,timesteps,simulations);
%disp(V0fixed)
%disp(errorFixedV)
%disp(Z0fixed)
%disp(errorFixedZ)

% %Floating strike price
%[V0floating, errorFloating] = asianOptionMCfloating(S0,sigma,r,T,timesteps,simulations);

% %Delta using fixed strike
deltaS = 0.001;
[V0fixedDelta, errorFixedVDelta, Z0fixedDelta, errorFixedZDelta] = ...
    asianOptionMCfixed(S0-deltaS,K,sigma,r,T,timesteps,simulations);

%DeltaV = (V0fixed - V0fixedDelta)/deltaS;
DeltaZ = (Z0fixed - Z0fixedDelta)/deltaS;
disp(DeltaZ)

