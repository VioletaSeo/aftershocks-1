
function [CAT, varargout] = make_catalog(NEWCATALOGNAME)
figure
[MSprefactor,MSalpha] = aftershock_productivity_kernel('/auto/home/kdascher/Documents/UCSC/general_data/Earthquake_data/NCEDC/NCEDC_CAT.mat','DepthRange', [0,55]);

% build comprehensive catalog:
% dimensions: Gavin Hayes, 2016
% Energy    : Andy Newman
% productivity: (here) - gcmt
% source complexity (here)/Hayes

% Load in each catalog and specify what the format of the Magnitude, time
% and spatial info is in (for merging purposes). Everything then gets
% appended into one big gnarly catalog. 

% MainShock (MS) catalog 
MScat = load('mainshock_info.mat');                 
    MScat           = struct2table(MScat);
    % removed this calculation and folded it into the kernel function
    %MScat.MSres     = getMSres(MScat.MSmag,MScat.MSprod,MSprefactor,MSalpha
    MSformat        = {'MSt','MSlat','MSlon','MSdepth','MSmag'};
    MScat.SeismoWidth=get_seismo_width();
    MScat.Background= get_background_rate();
    

% Radiated Energy (RE) catalog
REcat = load('EQ_Energy_cat.mat');                  
    REcat           = REcat.cat;
    REformat        = {'t','lat','lon','dep','mag'};

% Source Dimensions (SD) catalog
SDcat = load('eq_catalog_moderate_to_large_1FMSP.mat'); % ARBITRARILY CHOOSING FOCAL MECHANISM SOLUTION !!!   
    SDcat           = SDcat.cat;
    SDcat.mag       = 2/3*log10(SDcat.m0)-6.07;
    SDformat        = {'t','Lat','Lon','Dep','mag'};
    timeString      = SDcat.Date;
    dateTime        = datetime(timeString,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS'); 
    t               = datenum(dateTime);
    nanInd          = isnan(t);
    t(nanInd)       = datenum(datetime(timeString(nanInd),'InputFormat','yyyy-MM-dd''T''HH:mm:ss'));
    SDcat.t         = t;
    
% Source Complexity (SC) catalog get tacked onto the Source dimensions
SCcat       = store_slipmaps;
    iedList1        = SDcat.EID; % tack on IED specific values: (this is a bit ugly)
    iedList2        = {SCcat.EID}; 
    [idx1,idx2]     = merge_EID(iedList1,iedList2); % lookup corresponding locations
   
    SDcat.EllipseErrorNorm       = nan(height(SDcat),1);
    SDcat.EllipseErrorNorm(idx1) = [SCcat(idx2).EllipseErrorNorm]';

    SDcat.Dip                    = nan(height(SDcat),1);
    SDcat.Dip(idx1)              = [SCcat(idx2).Dip]';
    
% 	SDcat.StressDrop             = nan(height(SDcat),1);
%     SDcat.StressDrop(idx1)       = [SCcat(idx2).StressDrop]';
    
    SDcat.Perimeter              = nan(height(SDcat),1);
    for n = 1:length(idx1)   
        SDcat.Perimeter(idx1(n))        = faultperimeter(SCcat(idx2(n)),0.15);
    end
    
    % ... add more things here as needed
    
    
    % remove nans
    SDcat       = SDcat(~isnan(SDcat.EllipseErrorNorm),:); 
    
    


%% merge catalogs
%[mergedCat,I] = merge_eq_catalog(MScat,MSformat, ...
%                                  REcat,REformat, ...
%                                  SDcat, SDformat);

[mergedCat,I] = merge_eq_catalog(MScat,MSformat, ...
                                 SDcat, SDformat);

%% catalog subset
subCat = mergedCat(I,:);

option = 1;
switch option
    case 1
        %% distill catalog to desired source parameters
%         distilledCat = subCat(:,{'MSdepth', ...
%             'MSmag', ...
%             'MSres', ...
%             'MSprod', ...
%             'MeBB_appended_cat1', ...
%             'MeHF_appended_cat1', ...
%             'av_slip_appended_cat2', ...
%             'rlen_appended_cat2', ...
%             'rwid_appended_cat2', ...
%             'mu_appended_cat2'  , ...
%             'EllipseErrorNorm_appended_cat2', ...
%             'StressDrop_appended_cat2',...
%             'Perimeter_appended_cat2'});

        distilledCat = subCat(:,{'MSdepth', ...
            'MSmag', ...
            'MSres', ...
            'MSprod', ...
            'av_slip_appended_cat1', ...
            'rlen_appended_cat1', ...
            'rwid_appended_cat1', ...
            'mu_appended_cat1'  , ...
            'EllipseErrorNorm_appended_cat1', ...
            'Perimeter_appended_cat1',...
            'SeismoWidth', ...
            'Background'});
        
    case 2
        %% alternative set of values:
        % - depth
        % depth       = subCat.MSdepth*10^3; 
        
        % - stress drop (approx)
%         stressDrop  = subCat.Mo_appended_cat1./(subCat.rlen_appended_cat2*10^3.*subCat.rwid_appended_cat2*10^3).^(3/2);
%        stressDrop  = subCat.StressDrop_appended_cat2;

        Background = subCat.Background;
        SeismoWidth= subCat.SeismoWidth;
        
        
        % - aspect ratio
        aspectRatio = subCat.rlen_appended_cat2./subCat.rwid_appended_cat2;
        aspectRatio(aspectRatio<1) = 1./aspectRatio(aspectRatio<1);
        
        % perimeter
        perimeter  = subCat.Perimeter_appended_cat2;
        
        % - heterogeneity
        heterogenity= subCat.EllipseErrorNorm_appended_cat2;
        
        magnitude 
        
%         % - material
%         mu          = subCat.mu_appended_cat2;
%         
%         % - "missing volume"
%         missingLen = subCat.rlen_appended_cat2*10^3./(depth*10^3./sind([subCat.Dip_appended_cat2])); 

%         distilledCat= table(depth, stressDrop, aspectRatio, perimeter, heterogenity);
        distilledCat= table(Background, SeismoWidth, aspectRatio, perimeter, heterogenity);

end

%% OUT:
fms = subCat.MSfms;
CAT = distilledCat;
res = subCat.MSres; 
save(NEWCATALOGNAME,'CAT','res','fms')

%% varargout:
varargout{1} = res;
varargout{2} = subCat.MSfms;
end

function [IDX1,IDX2] = merge_EID(IED1,IED2)

nIn = length(IED2);
IDX1 = zeros(nIn,1);
IDX2 = 1:nIn;

for n = 1:nIn
    logLoc = strcmp(IED2{n},IED1);
    if any(logLoc)
        IDX1(n) = find(logLoc);
    else
        IDX1(n) = nan;
        IDX2(n) = nan;
        
    end    
end

IDX1 = IDX1(~isnan(IDX1));
IDX2 = IDX2(~isnan(IDX2));

end    


% removed this and folded it into the kernel to avoid redundant
% calculations
%     function MSRES = getMSres(MSMAG,MSPROD,PREFACTOR,ALPHA)
%         
%         allAS       = @(x,y) log10(y) - log10(PREFACTOR*10.^(ALPHA*x));
%         MSRES       = allAS(MSMAG,MSPROD);
%         
%     end

