%Function to calculate the option value for a asian option
%with a floating strike

function [V0,stdError] = asianOptionMCfloating(S0,sigma,r,T,timesteps,simulations)
    dt = T/timesteps;
    sqrtTime = sqrt(dt);

    V = zeros(simulations,1); 
    S = zeros(timesteps,1);
    
    for i = 1:simulations
        %Euler method
        S(1,1) = S0;
        for t = 1:timesteps
           randomnumber = randn;
           S(t+1,1) = S(t,1) + r*S(t,1)*dt + sigma*S(t,1)*sqrtTime*randomnumber;
        end
        V(i,1) = max(S(end,1)-mean(S),0);
    end

    EQ = mean(V);
    stdError = std(V)/sqrt(simulations);
    V0 = exp(-r*T)*EQ;
end