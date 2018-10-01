function sol = asianOptionexact(sigma, r, timesteps, E, T, s)

sigsqT = sigma^2*T*(2*timesteps + 1)/(6*timesteps + 6);
muT = 0.5*sigsqT + 0.5*(r - 0.5*sigma^2)*T;

d1 = (log(s/E) + (muT + 0.5*sigsqT))/(sqrt(sigsqT));
d2 = d1 - sqrt(sigsqT);

sol = exp(-r*T)*(s*exp(muT)*normcdf(d1) - E*normcdf(d2)); 

