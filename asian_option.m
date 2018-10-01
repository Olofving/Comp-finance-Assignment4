% Computational Finance
% Monte Carlo simulation of European call option

% Group 4
clc 
clear

%% Parameter Values
s = 14;
K = 15; 
r = 0.1;
sig = 0.25;
T = 0.5;
gamma = 1;%0.8;

%% Exact - bsexact 
bsexa = bsexact(sig, r, K, T, s);

%% Antithetic variables
samples = 1000;              % Varying no. samples
iterations = 10000;                  % Fix itereations
dt = (T/iterations);
E = zeros(length(samples),2);

asian_exact = asianOptionexact(sig, r, iterations, K, T, s)

for k=1:samples               % Sample trajectories
    S0 = s;
    %Sp0 = s;
    %Sn0 = s;
    for t=1:iterations      % Ito time discretisation
        q = randn;
        %Sp1 = Sp0 + r*Sp0*dt + sig*(Sp0^gamma)*sqrt(dt)*q;
        %Sn1 = Sn0 + r*Sn0*dt + sig*(Sn0^gamma)*sqrt(dt)*(-q);
        S1 = S0 + r*S0*dt + sig*(S0^gamma)*sqrt(dt)*q;
    
        S(t) = S1;
        
        %Sp0 = Sp1;
        %Sn0 = Sn1;
        S0 = S1;
    end
    %SpT(k) = Sp1;
    %SnT(k) = Sn1;
    ST(k) = S1;
    S_avg(k) = (1/iterations)*sum(S);
    V_asian_fix(k) = max(S_avg(k) - K,0);
    Vgeo(k) = max(exp(mean(log(S)))-K,0);
    V_asian_float(k) = max(S1 - S_avg(k),0);
    %VpT = max(Sp1-K,0);
    %VnT = max(Sn1-K,0);
    %VaT(k) = (VpT+VnT)/2;
    VT(k) = max(S1-K,0);
end
%E_Va = mean(VaT);
E_Vgeo = exp(-r*T)*mean(Vgeo)
E_V = exp(-r*T)*mean(VT)
E_V_asian_fix = exp(-r*T)*mean(V_asian_fix)
E_V_asian_float = exp(-r*T)*mean(V_asian_float)
Err = abs(bsexa - E_V);
%Err(i,2) = abs(bsexa - E_Va);

% figure
% plot(samples,avgErr(:,1),samples,avgErr(:,2))
% legend('std. var','antithetic var')
% title('Error of "normal"- vs antithetic variables over 20 simulations')
% xlabel('samples')
% ylabel('Average Error')

function sol = asianOptionexact(sigma, r, timesteps, E, T, s)

sigsqT = sigma^2*T*(2*timesteps + 1)/(6*timesteps + 6);
muT = 0.5*sigsqT + 0.5*(r - 0.5*sigma^2)*T;

d1 = (log(s/E) + (muT + 0.5*sigsqT))/(sqrt(sigsqT));
d2 = d1 - sqrt(sigsqT);

sol = exp(-r*T)*(s*exp(muT)*normcdf(d1) - E*normcdf(d2)); 
end