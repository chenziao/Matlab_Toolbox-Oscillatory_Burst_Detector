function y = maxwboltzcdf( x, a )
% Maxwell-Boltzmann cumulative distribution function.
y = erf(x/a/sqrt(2))-sqrt(2/pi)/a*x.*exp(-x.^2/(2*a^2));
y(x<0) = 0;
y(x==inf) = 1;
end

