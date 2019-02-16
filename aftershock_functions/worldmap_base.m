function worldmap_base

figure
hold on
worldmap({'World'})
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
    h{iClass} = scatterm(latPB,lonPB,[],colors(mod(iClass-1,length(colors))+1,:),'.');
end

legend([h{:}],expectedPBClass)
    
end

function [varargout]    = goodind(goodInd,varargin)

for n=1:length(varargin)
    varargout{n} = varargin{n}(goodInd);  %#ok<AGROW>
end

end