%% make the declustered mainshock catalog (CAT)
load('IRIS_DMC_with_FMS.mat')
tempCatFN= 'temp_earthquake_cat_merge.mat';
minMag = 6;
delete_temp         = @() system(sprintf('rm %s',tempCatFN));
load('IRIS_DMC_with_FMS.mat')
aftershock_productivity_kernel(...
    CAT.time, ...
    CAT.lat, ...
    CAT.lon, ...
    CAT.depth, ...
    CAT.M, ...
    CAT.fms, ...
    'SaveCatalog',tempCatFN, ...
    'MinMainshockMag',minMag, ...
    'PlotYN', 'no');
clear CAT
CAT = load(tempCatFN); delete_temp();

%% merge the enrgy estimates if available:
load eq_energy.mat;
eqenergy.Time = datenum(eqenergy.Time); 

%% need to get files from school computer

%% merge

% CAT:
MSformat = {'MSt','MSlat','MSlon','MSmag'};
REformat = {'Time','Latidude','Longitude','Magnitude'};
SDformat = {'Time','Lat','Lon','Mw'};

[C,~] = merge_eq_catalog(CAT, MSformat, ...
                         eqenergy, REformat);
mergedCat = merge_eq_catalog()

