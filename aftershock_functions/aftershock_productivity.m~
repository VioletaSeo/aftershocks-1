function aftershock_productivity

% calculate earthquake density and perform various correlations between
% productivity parameters and physical parameters

% code obtained from Emily Brodsky
% /auto/home/brodsky/Foreshock

% inputs:
analysisType    = 'depth';  % type of analysis to perform ('spatial','fms','crustalThickness','depth', 'crust')
subsampling     = {'plate_boundaries'};
catType         = 'GCMT';   % type of catalalog - this is far from perfect and not raw data ('CGMT','comcat','ANSS')
catFileName     = 'GCMT';   % 'california_ANSS.csv';


% load and get data from catalog
try     [t,lat,lon,depth,M, numEq,fms]  = readcatalog(catType,catFileName);
catch;  [t,lat,lon,depth,M, numEq]      = readcatalog(catType,catFileName);
end

% make an overview of the catalog 
catalog_overview(t,lat,lon,depth,M,numEq)

% subsample data
subsample_catalag(subsampling)


    function subsample_catalog(varargin)
        
        if      any(strcmp(varargin, 'plate_boundaries')); get_PB_eq
            
        end
        
        
        function get_PB_eq(lat, lon)
            %  load file with the plate boundary information - right now
            %  this is just the surface expression of plate boundaries,
            %  still (could still be intraplate)
            fn  = 'PB2002_steps.dat';
            fileID = fopen(fn);
            C = textscan(fileID,'%f %s %f %f %f %f %f %f %f %f %f %f %f %f %s');
            lonPB = C{3};
            latPB  = C{4};
            
            distance2PB = get_dist(lat, lon, latPB, lonPB);
            
            
            
            
            function DIST = get_dist(lat1, lon1, lat2. lon2)
                
                num1 = length(lat1);
                num2 = length(lat2); 
                
                dist1to2  = zeros(num1,num2);
                
                for n = 1:num1
                    for m = 1:num2
                        dist1to2(n,m) = distance('gc',lat1(n), lon1(n), lat2(m), lon2(m));
                    end
                end
                
                dists = min(dist1to2, [], 2);
                                
                % convert to metric distances 
                DISTS = deg2km(dists);
                
            end

        end
        
    end
% pipeline to function
if     strcmp(analysisType,	'spatial'           ); 	correlate_productivity(lat,3,'latitude');
elseif strcmp(analysisType, 'fms'               );	compare_fms;
elseif strcmp(analysisType, 'crustalThickness'  );  compare_thickness; % refer to function for specific input paramaters
elseif strcmp(analysisType, 'depth'             ); 	compare_depth('density','in crust');
elseif strcmp(analysisType, 'crust'             ); 	compare_dist_crust_edge;
end % strcmp(analysisType,

% main blocks:

    function correlate_productivity(correlatingParameter, numberOfBlocks, parameterName)
          
        
        % select data:
        
        numberOfBoxes = numberOfBlocks;
        minmaxPar   = minmax(correlatingParameter);
        parBracket  = linspace(minmaxPar(1), minmaxPar(2), numberOfBoxes+1);
        
        % initialize paramaters:
        alphaArray      = zeros(1,numberOfBoxes);
        prefactorArray  = zeros(1,numberOfBoxes);
        figure
        
        subplot(1,2,1)
        
        for parBlock = 1:numberOfBoxes
            
            parRange = parBracket(parBlock:(parBlock+1));
            goodInd  = correlatingParameter >= parRange(1) & correlatingParameter <parRange(2);
            
            selectedT       = t(goodInd);
            selectedLat     = lat(goodInd);
            selectedlon     = lon(goodInd);
            selectedDepth   = depth(goodInd);
            selectedM       = M(goodInd);
            
            if sum(goodInd) > 20
            
            [prefactorArray(parBlock), alphaArray(parBlock), magnitudeArray, magnitudeCountArray] = ...
                get_prod_parameters(selectedT, selectedLat, selectedlon, selectedDepth, selectedM);
            
            % plot an aggregat plot of the counts
            hold on
            scatter(magnitudeArray,magnitudeCountArray,[], ...
                [parBlock/numberOfBoxes, 0, 1-parBlock/numberOfBoxes], ...
                'filled')
            end
            
            
        end
        
        set(gca,'XScale', 'log', ...
            'YScale', 'log')
        xlabel('Magnitude')
        ylabel('Number Of Aftershocks') % anotating the aftershock prod plot
        
        % plot the evolution of the fit parameters with respect to the the
        % parameter of interest
        subplot(1,2,2)
        
        yyaxis left
        plot((parBracket(1:end-1)+parBracket(2:end))/2, prefactorArray)
        ylabel('prefactor')
        set(gca,'YScale','log')
        
        yyaxis right
        plot((parBracket(1:end-1)+parBracket(2:end))/2, alphaArray)
        ylabel('\alpha')
        
        xlabel(parameterName)

    
    end % correlate productivity

    function compare_fms

        % initialize paramaters:
        alphaArray      = zeros(1,3);
        prefactorArray  = zeros(1,3);
        colorArray      = {'r', 'g','b'};

        figure
        hold on


        for iFms = 1:3
            goodInd  = fms == iFms;
            selectedT       = t(goodInd);
            selectedLat     = lat(goodInd);
            selectedlon     = lon(goodInd);
            selectedDepth   = depth(goodInd);
            selectedM       = M(goodInd);
            
            [prefactorArray(iFms), alphaArray(iFms),magnitudeInc,magCount] = ...
                get_prod_parameters(selectedT, selectedLat, selectedlon, selectedDepth, selectedM);
            
            goodInd         = magCount ~= 0 & ~isnan(magCount);
            magCount        = magCount(goodInd);
            magnitudeInc    = magnitudeInc(goodInd);
            
            scatter(magnitudeInc,magCount, [colorArray{iFms},'o'])
            
            minmaxMag = minmax(magnitudeInc);
            plot(minmaxMag,prefactorArray(iFms)*10.^(alphaArray(iFms)*minmaxMag),[colorArray{iFms},'--']);
            

        end

        set(gca,'YScale','log')
        title('Aftershock Productivity')
        xlabel('Magnitude')
        ylabel('Average Number of Aftershocks')

        disp('                  Prefactor       Alpha')
        disp(['Strike faults:  ',num2str(prefactorArray(1)),'     ',num2str(alphaArray(1))]);
        disp(['Normal faults:  ',num2str(prefactorArray(2)),'     ',num2str(alphaArray(2))]);
        disp(['Thrust faults:  ',num2str(prefactorArray(3)),'     ',num2str(alphaArray(3))]);

    end

    function compare_thickness
        % this function look at potential correlations between crustal thickness and aftershock productivity
        availableCrust = get_crust_thickness(lat,lon,fms);

        %%%%%%
%         %testing if there where no correlation (shuffled vector)
%         for n = 1:10
%             availableCrust = availableCrust(randperm(length(availableCrust)));
%             correlate_productivity((availableCrust),50,'Available Crust')
%         end
        %%%%%%
        
        correlate_productivity((availableCrust),5,'Available Crust')
    end

    function compare_depth(varargin)
        
        metric = varargin{1};
        method = varargin{2};
        
        % thought process:
        % are aftershock just sampling the "local" productiviy of
        % earquakes
        
        if strcmp(method,'in crust')
            dist2crust = depth-get_crust_thickness(lat, lon);
            goodInd = dist2crust<0;       
            [t,lat,lon,depth,M] = goodind(goodInd,t,lat,lon,depth,M);
        end
        
        if strcmp(metric,'absolute')
            % quick and dirty
            correlate_productivity(depth, 10, 'depth')
            
        elseif strcmp(metric,'density')
            % normalized to the pre-existing distribution of earthquakes:
            [N,edges] = histcounts(depth,100);
            Ntot      = length(depth);
            localDepthDensity = interp1((edges(1:(end-1))+edges(2:end))/2,N,depth)/Ntot;
            correlate_productivity(localDepthDensity, 5, 'depth density')
        end
        
    end

    function compare_dist_crust_edge
        
        top     = 0;
        bottom  = get_crust_thickness(lat, lon);
        
        dist2crust = min([depth-top; depth-bottom]);
        
        goodInd = dist2crust>0;       
        [dist2crust,t,lat,lon,depth,M] = goodind(goodInd,dist2crust,t,lat,lon,depth,M);
        
        correlate_productivity(dist2crust,5,'Distance to edge of crust')
              
        
    end
end

function [t,lat,lon,depth,M,numEq, varargout] = readcatalog(catType,fileName)
% read the earthquake catalogs

%% input:
% catType:  string specifying the type of catalog used
%     'comcat'
%     'GCMT'
%     'ANSS'
%   *note that these are a bit hacked together so may nit work for any catalog
%   of this type
% fileName:     'path/to/the/fileName' of the catalog (string)

%% output: (truncated to magnitude completeness)
% t:      time (datetime)
% lat:    latitude (deg)
% lon:    longitude (deg)
% depth:  depth (km)
% M:      Moment magnitude (Ms)
% numEq:  number of eqrthquakes in the catalog

%%  varargout:
% varartout{1}: focal mechanism

%%

if strcmp(catType, 'comcat')
    load(fileName);     % this should be the name of the .mat catalog getting imported
    C       = mag25plus20102017;
    
    mth     = get_completeness(C.mag);
    C       = C(C.mag>mth,:);

    t       = parsetime(C.time)';% time of the earthquake
    lat     = C.latitude';                          % latitude
    lon     = C.longitude';                         % longitude
    depth   = C.depth';                             % depth (km)
    M       = C.mag';                               % ? does it matter if the catalog is mixed
    
    [t,lat,lon,depth,M] = goodind(M>mth,t,lat,lon,depth,M);
    
elseif strcmp(catType, 'GCMT')
    gCMT    = load(fileName);                                   % gCMT must be on the working path 'GCMT.mat'
    gCMT    = gCMT.gCMT; % yuck
    
    mth     = get_completeness(gCMT.Mw);
    goodInd= gCMT.Mw > mth & (gCMT.FMS_Type == 1 | gCMT.FMS_Type == 2 | gCMT.FMS_Type == 3); % remove magnitude less than threshold, and without focal mechanism solution
    gCMT    = gCMT(goodInd,:);
    
    t       = gCMT.datenum';
    lat     = gCMT.latitude';
    lon     = gCMT.longitude';
    depth   = gCMT.depth';
    M       = gCMT.Mw';
    fms     = gCMT.FMS_Type';                        % the focal mechanism (1: strike 2: normal 3: thrust)
    varargout{1} = fms;
    
elseif strcmp(catType, 'ANSS')
    ANSS    = importdata(fileName);
    mth     = 2.5; %get_completeness(ANSS.data(:,4));
    ANSS.data=ANSS.data(ANSS.data(:,4)>mth,:); 
    
    t       = ANSS.data(:,5)';
    lat     = ANSS.data(:,1)';
    lon     = ANSS.data(:,2)';
    depth   = ANSS.data(:,3)';
    M       = ANSS.data(:,4)';
end

numEq = length(t);

end

function catalog_overview(t,lat,lon,depth,M,numEq,varargin)
% make an overview of the catalog 
figure
figureDim   = [2,3];
figCount    = 1;

%% 1. histogram
subplot(figureDim(1),figureDim(2),figCount); figCount = figCount+1;

histogram(M)
title(['Histogram, n= ',num2str(numEq)])
ylabel('Number of Earthquakes') 
xlabel('Magnitude')

set(gca,'YScale','log')

%% 2. GRM
subplot(figureDim(1),figureDim(2),figCount); figCount = figCount+1;

numDivisions    = 50;
magInc          = linspace(min(M),max(M),numDivisions);
MCount          = zeros(1,numDivisions);

for n = 1:numDivisions
    MCount(n) = sum(M>magInc(n));
end

scatter(magInc,MCount)
title('Gutenberg-Righter')
ylabel('N>M')
xlabel('Magnitude')

set(gca,'YScale','log')

%% 3. Occurence
subplot(figureDim(1),figureDim(2),figCount); figCount = figCount+1;

minPlotM = 4; % plot only magnitudes larger than this
goodInd = M>minPlotM;
stem(t(goodInd),M(goodInd))

title('Occurence')
ylabel('Magnitude')
xlabel('Time (datetime')

%% 4. depth hist
subplot(figureDim(1),figureDim(2),figCount); figCount = figCount+1;
[N,edges] = histcounts(depth);
barh(-(edges(1:(end-1))+edges(2:end))/2,N,'hist')

title('Depth Distribution')
xlabel('Number of Earthquakes')
ylabel('Depth (km)')

%% 5. map
subplot(figureDim(1),figureDim(2),figCount); figCount = figCount+1;
scatter(lon,lat,M,-depth,'filled');

title('Map')
xlabel('longitude')
ylabel('latitude')

axis equal
colorbar

%% 6. Moment Release
subplot(figureDim(1),figureDim(2),figCount); figCount = figCount+1;
M0 = 10.^(3/2*(M + 10.7));
plot(t,cumsum(M0));
title('Cummulative moment release')
ylabel('Moment (erg)')
xlabel('Time (datetime)')


end

function mth            = get_completeness(M)

edges       = floor(min(M)):0.1:ceil(max(M));
[N,magBin]  = histcounts(M,edges);
mth=magBin(N==max(N));

end % obtain the completeness of the earthquake catalog

function AVAILABLECRUST = get_crust_thickness(lat,lon,varargin)

if nargin >= 3
    fms = varargin{1};
end

crustThicknessFileName  = 'crust1thickness.xyz';    % file name of the file containing xyz data of crust thickness
availableCrust          = 'whole';                  % appraoch used to measure available crust

% load in data%         %testing if there where no correlation (shuffled vector)
%         for n = 1:10
%             availableCrust = availableCrust(randperm(length(availableCrust)));
%             correlate_productivity((availableCrust),50,'Available Crust')
%         end

crustalThicknessData = load(crustThicknessFileName);

% Get the thickness of the crust for each earthquake
% regrid the xyz data
xDim =  sum(crustalThicknessData(:,1)==crustalThicknessData(1,1));
yDim =  length(crustalThicknessData)/xDim;

crustLon        = crustalThicknessData(:,1);
crustLat        = crustalThicknessData(:,2);
crustThickness  = crustalThicknessData(:,3);

interpolant = scatteredInterpolant(crustLat,crustLon,crustThickness,...
    'linear','nearest');
eqCrustThickness = interpolant(lat,lon);


% potential approaches to determining available crust
% 1) whole:             thickness of the crust
% 2) alongPlane:        thickness of the crust along fault plane
% 3) afterShockZone:    available 'aftershock zone' (N(r) = kr^(-1.3).

if      strcmp(availableCrust, 'whole')
    availableCrust = eqCrustThickness;
elseif  strcmp(availableCrust, 'alongPlane')
    % assume andersonian fault and get plane dip
    faultDip = zeros(1,length(eqCrustThickness));
    faultDip(fms == 1) = 90;
    faultDip(fms == 2) = 60;
    faultDip(fms == 3) = 15;
    
    % correct crustal thickness
    availableCrust = eqCrustThickness./sind(faultDip);
elseif  strcmp(availableCrust, 'afterShockZone')
    
    % the idea hare will be to use the distance decay equation to
    % essentially look at the fraction (?) of earthquakes that are "eaten
    % up" by non-seismogenic
end

AVAILABLECRUST = availableCrust;

end % get the crustal thickness of the crust

function [varargout]    = goodind(goodInd,varargin)

for n=1:length(varargin)
    varargout{n} = varargin{n}(goodInd);
end

end