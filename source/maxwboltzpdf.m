function y = maxwboltzpdf( x, a )
% Maxwell-Boltzmann probability density function.
y = sqrt(2/pi)/a^3*x.^2.*exp(-x.^2/(2*a^2));
y(x<0|x==inf) = 0;
end

