function [prefactor,exponent] = getproductivity(magnitude,N,varargin)
% measure the asftershock productivity based on the number of earthquakes
% 'N'  of magnitude 'magnitude'. TO deal with background seismicity a fit a
% peicewies function

% calls on 'makepeicewisefit' and 'lmsengine'

badInd = isnan(magnitude) | isnan(N) | N==0;
magnitude = magnitude(~badInd);
N         = N(~badInd);

logN            = log10(N);

numberOfKnots   = 3; % 2 ends on cusp

if strcmp(varargin,'peicewise')
    [a,b]           = makepeicewisefit(magnitude,logN,numberOfKnots);
    a               = a(2);
    b               = b(2);
else
    d               = polyfit(magnitude,logN,1);
    a               = d(1);
    b               = d(2);
end

% productivity: N = 10^b*10^(a*mag) = profactor*10^(exponent*mag)
prefactor       = 10^b;
exponent        = a;


end