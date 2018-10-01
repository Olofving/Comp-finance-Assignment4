%Function to calculate the option value for a asian option
%with a fixed strike

function [V0avg,stdErrorV,Z0avg,stdErrorZ] = asianOptionMCfixed(S0,K,sigma,r,T,timesteps,simulations)
    dt = T/timesteps;
    sqrtTime = sqrt(dt);

    V = zeros(simulations + 1,1); 
    Vgeo = zeros(simulations + 1,1);
    S = zeros(timesteps + 1,1);
    Sgeo = zeros(timesteps + 1,1);
    
    for i = 1:simulations
        %Euler method
        S(1,1) = S0;
        Sgeo(1,1) = S0;
        for t = 1:timesteps
           S(t+1,1) = S(t,1) + r*S(t,1)*dt + sigma*S(t,1)*sqrtTime*randn;
           Sgeo(t+1,1) = Sgeo(t,1) + r*Sgeo(t,1)*dt + sigma*Sgeo(t,1)*sqrtTime*randn;
        end
        V(i,1) = max(mean(S)-K,0);
        Vgeo(i,1) = max(exp(mean(log(S)))-K,0);
    end
    V = exp(-r*T)*V;
    Vgeo = exp(-r*T)*Vgeo;
    
    covGeo = cov(V,Vgeo);%Covariance matrix, diagonal is the normal variance
    c = -covGeo(1,2)/covGeo(2,2);
    exactSol = asianOptionexact(sigma, r, timesteps, K, T, S0);
    Z = V + c*(Vgeo - exactSol);

    V0avg = mean(V);
    stdErrorV = std(V);
    
    Z0avg = mean(Z);
    stdErrorZ = std(Z);
end