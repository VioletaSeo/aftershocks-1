function [xq,yq,varargout] = mov_mean(x,y,N,r)

sX = size(x); sY = size(y);
if length(x)~= sX(1); x = x';  end
if any(sY ~= size(x)); y = y'; end
if nargin < 4; r = range(x')/N;end

Inan = isnan(x) | isnan(y);
x = x(~Inan);
y = y(~Inan);

minmaxX = minmax(x');
xq      = linspace(minmaxX(1),minmaxX(2),N);
yq      = zeros(N,1);
yerr    = zeros(N,2);
for n = 1:N
    I = x>(xq(n)-r) & x<(xq(n)+r);
    if sum(I) > 3
        yq(n)   = median(y(I));
        yerr(n,:) = prctile(y(I),[25,75]);
    else
        yq(n) = nan;
        yerr(n,:) = [nan,nan];
    end
end
varargout{1} = yerr;
end
