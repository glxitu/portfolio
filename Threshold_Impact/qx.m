function [ y ] = qx( x )

y = (1/2)*( 1 - erf( x/sqrt(2) ) );

end

