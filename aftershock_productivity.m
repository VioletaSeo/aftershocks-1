function aftershock_productivity(gCMT)

% calculate earthquake density

% code obtained from Emily Brodsky
% /auto/home/brodsky/Foreshock

% Usage
% drag and drop catalog into workspace

% function calls:

% get_prod_parameters
%   getproductivity
%       makepeicewise (optional)
%   parsetime

% set the following to mat        goodInd  = lat >= latRange(1) & lat <latRange(2);

% % set time wind of the search
% startTime   = datenum('2014-01-01 00:00:00');
% endTime     = datenum('2015-01-01 00:00:00');
%
% %set location window of search
% minLat      =-90;
% maxLat      = 90;
%
% minLong     =-180;
% maxLong     = 180;

% % get data
% [year, month, day, hour, sec, lat, long, depth, mag, magType] ...
%     = LoadComCat(startTime,endTime,minMagnitude); %,minLat,maxLat,minLong,maxLong);
% clear;

analysisType= 'fms';% type of analysis to perform ('spatial','fms')
catType = 'GCMT';   % type of catalalog - this is far from perfect and not raw data


if strcmp(catType, 'comcat')
    load MyCatalog;     % this should be the name of the .mat catalog getting imported
    C       = mag25plus20102017;
    t       = parsetime(C.time)';% time of the earthquake
    lat     = C.latitude';                          % latitude
    lon     = C.longitude';                         % longitude
    depth   = C.depth';                             % depth (km)
    M       = C.mag';                               % ? does it matter if the catalog is mixed
    
elseif strcmp(catType, 'GCMT')
    %load gCMT                                       % gCMT must be on the working path
    t       = gCMT.datenum';
    lat     = gCMT.latitude';
    lon     = gCMT.longitude';
    depth   = gCMT.depth';
    M       = gCMT.mb';
        
    fms     = gCMT.FMS_Type';                        % the focal mechanism (1: strike 2: normal 3: thrust)
    
end

mth = get_completeness(M);

if strcmp(analysisType, 'spatial');     make_crosssection(mth,t,lat,lon,depth,M);
elseif strcmp(analysisType,'fms');      compare_fms;
    
    
        
end % strcmp(analysisType,

% main blocks:

    function make_crosssection
% select data:

numberOfBoxes = 2;
minmaxLat   = minmax(lat);
latBracket  = linspace(minmaxLat(1), minmaxLat(2), numberOfBoxes+1);

% initialize paramaters:
alphaArray      = zeros(1,numberOfBoxes);
prefactorArray  = zeros(1,numberOfBoxes); 



    for latBlock = 1:numberOfBoxes
        
        latRange = latBracket(latBlock:(latBlock+1));
        goodInd  = lat >= latRange(1) & lat <latRange(2);
        
        selectedT       = t(goodInd);
        selectedLat     = lat(goodInd);
        selectedlon     = lon(goodInd);
        selectedDepth   = depth(goodInd);
        selectedM       = M(goodInd);
        
        [prefactorArray(latBlock), alphaArray(latBlock)] = ...
            get_prod_parameters(selectedT, selectedLat, selectedlon, selectedDepth, selectedM, mth);
    end
    
    figure
    
    % here you can choose to plot 'prefactor' or 'alpha'
    plotOption = 'prefactor';
    
    yyaxis left
    scatter(lat,lon,M,depth)
    axis equal
    
    yyaxis right
    if strcmp(plotOption, 'prefactor')
        plot((latBracket(1:end-1)+latBracket(2:end))/2, prefactorArray)
        ylabel('prefactor')
    elseif strcmp(plotOption,'alpha')
        plot((latBracket(1:end-1)+latBracket(2:end))/2, alphaArray)
        ylabel('\alpha')
    end


end
    function compare_fms
        
        % initialize paramaters:
        alphaArray      = zeros(1,3);
        prefactorArray  = zeros(1,3);
        
        for iFms = 1:3
            
        goodInd  = fms == iFms;
        selectedT       = t(goodInd);
        selectedLat     = lat(goodInd);
        selectedlon     = lon(goodInd);
        selectedDepth   = depth(goodInd);
        selectedM       = M(goodInd);    
        [prefactorArray(iFms), alphaArray(iFms)] = ...
            get_prod_parameters(selectedT, selectedLat, selectedlon, selectedDepth, selectedM, mth);
                
        end
        
        disp('              Prefactor   Alpha')
        disp('Strike slip:  ',num2str(prefactorArray(1)),'     ',num2str(alphaArray(1)));
        disp('Noraml slip:  ',num2str(prefactorArray(2)),'     ',num2str(alphaArray(2)));
        disp('Thrust slip:  ',num2str(prefactorArray(3)),'     ',num2str(alphaArray(3)));
        
    end

end

function mth = get_completeness(M)

[N,magBin] = hist(M);
mth=magBin(N==max(N));

end  

    





