function [CAT, varargout] = make_catalog_FSP(quakes,NEWCATALOGNAME)
% badly written
% HAS PATH DEPENDENCIES!!! 
% needs quakes

figure
path2data =  '/auto/home/kdascher/Documents/UCSC/projects/aftershock_productivity/aftershock_data/GCMT/eq_catalog.mat';
% aftershock_productivity_kernel(path2data,'DepthRange', [0,55]);
aftershock_productivity_kernel(quakes.otime, quakes.lat, quakes.lon, quakes.depth, quakes.mag, ...
    'DepthRange', [0,55],...
    'MinMainshockMag',6);



%%%% build comprehensive catalog:
% productivity: (here) - gcmt
% source complexity (here)/Hayes

% Load in each catalog and specify what the format of the Magnitude, time
% and spatial info is in (for merging purposes). Everything then gets
% appended into one big gnarly catalog. 

%% MainShock (MS) catalog 
MScat           = load('mainshock_info.mat');                 
MScat           = struct2table(MScat);
MSformat        = {'MSt','MSlat','MSlon','MSmag'};
MScat.SeismoWidth=get_seismo_width();
MScat.Background= get_background_rate();

%% source attribute catalog
dir = '/auto/home/kdascher/Documents/UCSC/general_data/Earthquake_data/hayes2017kinematic/FSPs/Inversion_files';
addpath(dir);
addpath('/auto/home/kdascher/Documents/UCSC/projects/stress_drop');

try 
    FSPs = load('fsp.mat','fspInvObjArray');
    FSPs = FSPs.fspInvObjArray;
catch 
    FSPs  = loop_over_dir(dir);     % get fsp_inversion object array and compute stress drops
end

FSPcat = flat(FSPs(1));         % turn fsp obj to flat table
for n = 2:length(FSPs)
    FSPcat = [FSPcat;flat(FSPs(n))]; %#ok<AGROW>
end
FSPformat = FSPcat.Properties.VariableNames([1:3,5]); 

% build table wit t, lat, lon, depth, Mag, 
    

%% merge catalogs

[mergedCat,I] = merge_eq_catalog(FSPcat, FSPformat, ...
                                [1,1,1,1],...
                                 MScat,MSformat);
distilledCat = mergedCat(I,:);
%% catalog subset
subCat = mergedCat(I,:);
distilledCat = subCat;
% %% distill catalog to desired source parameters
% 
% distilledCat = subCat(:,{'MSdepth', ...
%     'MSmag', ...
%     'MSres', ...
%     'MSprod', ...
%     'StressDrop',...
%     'aspectRatio',...
%     'heterogeneity'});
% 
% 
% % alt: distilledCat= table(Background, SeismoWidth, aspectRatio, perimeter, heterogenity);

%% OUT:
fms = subCat.MSfms; %#ok<NASGU>
CAT = distilledCat;
res = subCat.MSres; 
save(NEWCATALOGNAME,'CAT','res','fms')

%% varargout:
varargout{1} = res;
varargout{2} = subCat.MSfms;

end
