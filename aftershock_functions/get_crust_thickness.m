function AVAILABLECRUST = get_crust_thickness(lat,lon,varargin)

if nargin >= 3
    fms = varargin{1};
end

crustThicknessFileName  = 'crust1thickness.xyz';    % file name of the file containing xyz data of crust thickness
availableCrust          = 'whole';                  % appraoch used to measure available crust

% load in data%         %testing if there where no correlation (shuffled vector)
%         for n = 1:10
%             availableCrust = availableCrust(randperm(length(availableCrust)));
%             correlate_productivity((availableCrust),50,'Available Crust')
%         end

crustalThicknessData = load(crustThicknessFileName);

% Get the thickness of the crust for each earthquake
% regrid the xyz data
xDim =  sum(crustalThicknessData(:,1)==crustalThicknessData(1,1));
yDim =  length(crustalThicknessData)/xDim;

crustLon        = crustalThicknessData(:,1);
crustLat        = crustalThicknessData(:,2);
crustThickness  = crustalThicknessData(:,3);

interpolant = scatteredInterpolant(crustLat,crustLon,crustThickness,...
    'linear','nearest');
eqCrustThickness = interpolant(lat,lon);


% potential approaches to determining available crust
% 1) whole:             thickness of the crust
% 2) alongPlane:        thickness of the crust along fault plane
% 3) afterShockZone:    available 'aftershock zone' (N(r) = kr^(-1.3).

if      strcmp(availableCrust, 'whole')
    availableCrust = eqCrustThickness;
elseif  strcmp(availableCrust, 'alongPlane')
    % assume andersonian fault and get plane dip
    faultDip = zeros(1,length(eqCrustThickness));
    faultDip(fms == 1) = 90;
    faultDip(fms == 2) = 60;
    faultDip(fms == 3) = 15;
    
    % correct crustal thickness
    availableCrust = eqCrustThickness./sind(faultDip);
elseif  strcmp(availableCrust, 'afterShockZone')
    
    % the idea hare will be to use the distance decay equation to
    % essentially look at the fraction (?) of earthquakes that are "eaten
    % up" by non-seismogenic
end

AVAILABLECRUST = availableCrust;

end % get the crustal thickness of the crust