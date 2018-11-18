function [segmentSlopes,segmentIntercepts] = makepiecewisefit(x,y,N)
% Small function that make a peicewise linear fit with N knows (N = ends (2) +
% breakpoints) to a set of points p = (x,y). The function also plots the
% outputed fit

% uses amazing toolbox lmsengine

% fit peice wise linear model (degree 1), make the plot ('plot', 'on')
% with N knots, ('knots',N), with the knots locations free to move
% ('InteriorKnots','free')

slm = slmengine(x,y,'degree',1,'plot','on','knots',N,'InteriorKnots','free');

% get knots
knotX = slm.knots;
knotY = slm.coef;

segmentSlopes       = (knotY(2:end)-knotY(1:end-1))-(knotX(2:end)-knotX(1:end-1));
segmentIntercepts   = knotY(1:end-1)-segmentSlopes.*knotX(1:end-1);

end

