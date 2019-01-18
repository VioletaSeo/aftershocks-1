load eq_catalog.mat
numEq = length(lat);

wgs84 = wgs84Ellipsoid('meters');
[X,Y,Z] = geodetic2ecef(wgs84, lat', lon', -depth'*1000);

figure
scatter3(X,Y,Z,'.')
axis equal
hold on

OT = OcTree([X,Y,Z],'minSize',50*1000,'binCapacity',100);

boxH = OT.plot;
cols = lines(OT.BinCount);
doplot3 = @(p,varargin)plot3(p(:,1),p(:,2),p(:,3),varargin{:});
for i = 1:OT.BinCount
    set(boxH(i),'Color',cols(i,:),'LineWidth', 1+OT.BinDepths(i))
    doplot3(pts(OT.PointBins==i,:),'.','Color',cols(i,:))
end
axis image, view(3)