function [CAT, varargout] = make_catalog_FSP_2()
%%
load('FSPs.mat')
FSPcat = flat(FSPs(1));         % turn fsp obj to flat table
for n = 2:length(FSPs)
    FSPcat = [FSPcat;flat(FSPs(n))]; %#ok<AGROW>
end
FSPformat = FSPcat.Properties.VariableNames([2:3,5]); 
 
%%
load('IRIS_DMC_with_FMS_and_energy.mat');
CAT = iris_dmc_cat_with_fms_and_energy;
minMag  = 6;
[ASinfo] = aftershock_productivity_kernel(...
    CAT.time, ...
    CAT.lat, ...
    CAT.lon, ...
    CAT.depth, ...
    CAT.M, ...
    CAT.fms, ...
    'MinMainshockMag',minMag, ...
    'DepthRange',[0,55], ...
    'ReturnCatalog', 'yes', ...
    'SaveCatalog', 'no', ...
    'PlotYN','no');
MSCat = CAT(ASinfo.ID,:);
MSCat.MSres = ASinfo.MSres;
MSCat.MSprod= ASinfo.MSprod;
MSCat.age = get_crust_age(MSCat.lat, MSCat.lon, MSCat.depth);
MSformat = MSCat.Properties.VariableNames([3:4,6]);

%% merge catalogs

[mergedCat,I] = merge_eq_catalog(FSPcat, FSPformat, [1,1,1.5], ...
                                 MSCat,MSformat);
distilledCat = mergedCat(I,:);
%% catalog subset

%% OUT:
fms = distilledCat.fms_appended_cat1;
CAT = distilledCat;
res = distilledCat.MSres_appended_cat1; 

%% varargout:
varargout{1} = res;
varargout{2} = fms;

end
