% European Option pricing by Black-Scholes equations (benchmark value)

function [type, BS_price] = BS_option_price(S0,k,sigma,r,T,type)

    d1 = (log(S0/k) + (r + sigma^2/2)*T)/(sigma*sqrt(T));
    d2 = d1 - sigma*sqrt(T);
    if strcmp(type, 'call')
        BS_price = S0*normcdf(d1) - k*exp(-r*T)*normcdf(d2);
    else
        BS_price = k*exp(-r*T)*normcdf(-d2) - S0*normcdf(-d1);
    end