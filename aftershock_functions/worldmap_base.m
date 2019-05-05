function worldmap_base

figure
hold on
worldmap({'world'})
load coastlines
%plotm(coastlat, coastlon)
geoshow('landareas.shp', 'FaceColor', [0.9 0.9 0.9])

fn      = 'PB2002_steps.dat';
fileID  = fopen(fn);
C       = textscan(fileID,'%f %s %f %f %f %f %f %f %f %f %f %f %f %f %s');
PBclass = C{end};
expectedPBClass             = {'OSR','OTF','OCB','CRB','CTF','CCB','SUB'};
NClass  = length(expectedPBClass);
h = cell(NClass,1);
colors     = get(gca, 'ColorOrder');

for iClass = 1:NClass
    PBI     = contains(PBclass,expectedPBClass{iClass});
    [latPB,lonPB] = goodind(PBI,C{[4,3]});
    [latE, lonE]  = goodind(PBI,C{[6,5]});
    
    lat = merge_vec_nan(latPB,latE);
    lon = merge_vec_nan(lonPB,lonE);
    
    h{iClass} = plotm(lat,lon,'color',colors(mod(iClass-1,length(colors))+1,:),'LineWidth',2);
end

%varargout{1} = legend([h{:}],expectedPBClass);
    
end

function d = merge_vec_nan(a,b)

d = nan(3*length(a),1);
d(1:3:end-2) = a;
d(2:3:end-1) = b;

end

function [varargout]    = goodind(goodInd,varargin)

for n=1:length(varargin)
    varargout{n} = varargin{n}(goodInd);  %#ok<AGROW>
end

end