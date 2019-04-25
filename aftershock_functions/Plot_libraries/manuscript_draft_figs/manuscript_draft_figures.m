% This file (hopefully) contains manuscript-ready figures for the
% aftershock productivity project. The structure should be very similar to
% that of the proposal script
% 

%% preamble:
% clear
% close all
SAVEFIG = 'no';
p = mfilename('fullpath');
[filepath,~,~] = fileparts(p);
ftsz    = @(fh,fontSize) set(findall(fh,'-property','FontSize'),'FontSize',fontSize);
setsize = @(fh,dim1,dim2) set(fh,...
    'Units',        'Inches', ...
    'Position',     [0,0,dim1,dim2],...
    'PaperUnits',   'Inches',...
    'PaperSize',    [dim1,dim2]);

%% Load data and format
load('IRIS_DMC_with_FMS_and_energy.mat')
CAT     = iris_dmc_cat_with_fms_and_energy;

% MAGIC NUMBERS                                                 Justification
% -------------                                                 -----------
minMag      = 6.5; % Minimum magnituce of mainshocks considered (Mc + 2)
maxDepth    = 55;  % Maximum mainshock depth considered         (No deep focus earthquakes)
completeness= 4.5; % Completeness of the catalog                (from KS test for Mid ocean rideges)
distance2pb = 400; % Maximum distance to plate                  ~ maxDepth/tand(10) 


% aftershock statistics
[ASinfo,k,alpha] = aftershock_productivity_kernel(...
    CAT.time, ...
    CAT.lat, ...
    CAT.lon, ...
    CAT.depth, ...
    CAT.M, ...
    CAT.fms, ...
    'MinMainshockMag',minMag, ...
    'DepthRange',[0,maxDepth], ...
    'ReturnCatalog', 'yes', ...
    'SaveCatalog', 'no', ...
    'PlotYN','no', ...
    'Completeness',completeness);
MSCat = CAT(ASinfo.ID,:);
MSCat.MSres = ASinfo.MSres;
MSCat.MSprod= ASinfo.MSprod;

% shorter variable names
M       = MSCat.M;
lat     = MSCat.lat;
lon     = MSCat.lon;
depth   = MSCat.depth;
t       = MSCat.time;

% aftershock stats
res     = MSCat.MSres;
prod    = MSCat.MSprod;

% source properties
fms     = MSCat.fms;
E       = MSCat.BBtotalEnergy;
Enorm   = E./mw2mo(M);
dur     = MSCat.RuptureDuration;
durNorm = dur./(mw2mo(M)).^(1/3);
RT      = MSCat.RuptureDuration;

% assign ages
age = get_crust_age(MSCat.lat, MSCat.lon, MSCat.depth);

% plate boundaries
%                        see above   (enforce fms)                             
assign_PB_class(lat,lon, distance2pb,'yes',fms);


% plotting colors (be focal mechanism)
colors = {[0.7      0.7       0.7], ...
          [0        0.4470    0.7410],...
          [0.4660   0.6740    0.1880], ...
          [0.6350   0.0780    0.1840]}; % colors
c = assign_color(fms+1,colors);


%% main document figures:
% 1) World map:
lgd = worldmap_base();
scatterm(lat,lon,[],res, 'filled')
h = colorbar;
lgd.NumColumns = 4;
ylabel(h,'Relative Productivity')
setsize(gcf,8,3.5)
savefigure(gcf,'worldmap_res',SAVEFIG)

% 2) Productivity by plate boundary

%% Supplement:



%% functions:

function C = assign_color(X,c)
nX = length(X);
C = zeros(nX,3);
for n = 1:nX
    C(n,:) = c{X(n)};
end
end

function fade_plot(X,yPos,c)
sz     = 200;
scatter(X,repmat(yPos,1,length(X)),sz,'filled','MarkerFaceAlpha',0.03,'MarkerFaceColor',c)
plot(prctile(X,[25,75]),[yPos,yPos],'LineWidth',3,'Color',c/1.4)  
scatter(median(X),yPos,sz,'filled','MarkerFaceColor',c/1.4)
end

function savefigure(fh,figureName,plotYN)
if strcmp(plotYN,'yes')
saveas(fh,figureName);
print(fh,[figureName,'.png'],'-dpng','-r300');
end
end