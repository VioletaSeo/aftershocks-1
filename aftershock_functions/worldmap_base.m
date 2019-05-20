function worldmap_base(region,newplot)

if nargin == 2 && ~strcmp(newplot,'no')
    figure
end

hold on
if nargin == 0
    h = worldmap([-90 90], [-360 0]);
else
    h = worldmap(region);
end

set(findall(h,'Tag','PLabel'),'visible','off')
set(findall(h,'Tag','MLabel'),'visible','off')
    
load coastlines
%plotm(coastlat, coastlon)
geoshow('landareas.shp', 'FaceColor', [0.7 0.7 0.7],'EdgeColor','none')

fn      = 'PB2002_steps.dat';
fileID  = fopen(fn);
C       = textscan(fileID,'%f %s %f %f %f %f %f %f %f %f %f %f %f %f %s');
PBclass = C{end};
expectedPBClass             = {'OSR','OTF','OCB','CRB','CTF','CCB','SUB'};
NClass  = length(expectedPBClass);
h = cell(NClass,1);
colors     = [0.4660    0.6740    0.1880; ...
                   0    0.4470    0.7410; ...
              0.6350    0.0780    0.1840; ...
              0.4660    0.6740    0.1880; ...
                   0    0.4470    0.7410; ...
              0.6350    0.0780    0.1840; ...
              0.6350    0.0780    0.1840];
                   
              
%colors     = get(gca, 'ColorOrder');

for iClass = 1:NClass
    PBI     = contains(PBclass,expectedPBClass{iClass});
    [latPB,lonPB] = goodind(PBI,C{[4,3]});
    [latE, lonE]  = goodind(PBI,C{[6,5]});
    
    lat = merge_vec_nan(latPB,latE);
    lon = merge_vec_nan(lonPB,lonE);
    h{iClass} = plotm(lat,lon,'color',[0.5350    0.0780    0.1840],'LineWidth',0.8);
    %h{iClass} = plotm(lat,lon,'color',colors(mod(iClass-1,length(colors))+1,:),'LineWidth',2);
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