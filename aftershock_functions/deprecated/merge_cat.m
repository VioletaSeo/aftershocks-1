function [cat1] = merge_cat(fn1,fn2)
% this function merges seismic catalogs based on occurence time.

% load catalogs 

% load('eq_catalog.mat'); % 
% load('eq_catalog.mat');


cat1 = load(fn1); % 
cat2 = load(fn2);

% parse catalogs/get t-m-x-y-z pts matrix

%% parse catalog 1
% this may need some editing depending on the catalog that you are using
p1 = [cat1.MSt',   cat1.MSlat',   cat1.MSlon',   cat1.MSdepth',     cat1.MSmag'];


%% parse catalog 2
% this may need some editing depending on the catalog that you are using
cat2 = cat2.catalog;
cat2 = cat2(2:end,:);

timeString  = cat2.Date;
dateTime    = datetime(timeString,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS'); 
t           = datenum(dateTime);

nanInd      = isnan(t);
t(nanInd)   = datenum(datetime(timeString(nanInd),'InputFormat','yyyy-MM-dd''T''HH:mm:ss')); 

mag         = log10(cat2.m0*10^7)/1.5-10.7;


p2          = [t,cat2.Lat, cat2.Lon, cat2.Dep, mag];


%% merge

[Idx,~]                 = knnsearch(p1,p2,'Distance','seuclidean');

cat1.stress_drop        = nan(1,length(cat1.MSt));
% cat1.stress_drop(Idx)   = get_stressdrop(cat2.m0, cat2.rlen*1000,cat2.rwid*1000);
cat1.stress_drop(Idx)   = get_stressdrop(cat2.mu, cat2.av_slip, cat2.rlen*1000,cat2.rwid*1000);


cat1.area               = nan(1,length(cat1.MSt));
cat1.area(Idx)          = cat2.rlen*1000.*cat2.rwid*1000;

cat1.rlen               = nan(1,length(cat1.MSt));
cat1.rwid               = nan(1,length(cat1.MSt));
cat1.rlen(Idx)               = cat2.rlen*1000;
cat1.rwid(Idx)               = cat2.rwid*1000;




end

% function DELTASIGMA = get_stressdrop(m0,a1,a2)
% 
% DELTASIGMA = 7/16.*m0./((a1+a2)/2).^3;
% 
% end

function DELTASIGMA = get_stressdrop(mu,D,a1,a2)

DELTASIGMA = 7*pi/16 * mu .* D ./ ((a1.*a2).^(1/2));

end