function [N] = eq_in_radius(lat,lon,depth, radius)

tic 
numPts = numel(lat);
wgs84 = wgs84Ellipsoid('meters');
[X,Y,Z] = geodetic2ecef(wgs84, lat', lon', -depth'*1000);
OT = OcTree([X,Y,Z],'minSize',2*radius,'binCapacity',100);
toc

tic
eqInRadius = zeros(numPts,1);
for n = 1:numPts
    
    pt          = [X(n),Y(n),Z(n)];
    
%     % octree
%     tic
%     numNeighbours = 6;
%     
%     bounds      = [eye(3,3) * radius; ....
%                    eye(3,3) *-radius]; 
%     querries    = repmat(pt,numNeighbours,1) + bounds;
%     
%     binQ    = OT.query(querries);
%     
%     
%     binQ            = unique(binQ);
%     
%     inBin           = ismember(OT.PointBins, binQ);
%     PT              = repmat(pt,sum(inBin),1);
%     I               = radius > sqrt(sum((PT-OT.Points(inBin)).^2,2));
%     eqInRadius(n)   = sum(I);
%     disp(eqInRadius(n))
%     disp('octree');toc
%     
%     tic
    
    % brute foce
%     tic
    I               = radius > sqrt(sum((pt-[X,Y,Z]).^2,2));
    eqInRadius(n)   = sum(I)-1;
%     disp(eqInRadius(n))
%     disp('coord');toc
    
%     % original method
%     tic
%     I               = radius/1000 >  ...
%         sqrt(sum([depth(n)-depth' , deg2km(distance(lat(n),lon(n),lat,lon))'].^2,2));
%     eqInRadius(n) = sum(I);
%     disp(eqInRadius(n))
%     disp('original');toc
end
toc

N = eqInRadius;


end