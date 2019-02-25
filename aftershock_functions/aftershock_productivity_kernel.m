function varargout = aftershock_productivity_kernel(varargin)
% get productivity parameters for GCMT earthquake catalog

% input:

% aftershock_productivity_kernel(fn)
% specify a .mat filename, fn, or a gcmt catalog file name containing earthquake
% data.

% aftershock_productivity_kernel(t,lat,lon,depth,M)
% directly specify time, latitude, longitude, depth, magnitude, as N by 1
% arrays. Note that dummy value (0) is assigned to focal mechanism

% aftershock_productivity_kernel(t,lat,lon,depth,M,fms)
% directly specify time, latitude, longitude, depth, magnitude, as N by 1
% arrays.

% aftershock_productivity_kernel(...,'Spec',SpecVal)     
% see parseinput for possible specifciations

% output:

% [p,alpha, FSp, FS ] = aftershock_productivity_kernel(...)
% prefactor (p) and alpha prameter of the productivity law:
% N = k10^alpha*M
% where M is the mangitude


%% get data
[t,lat,lon,depth,M,fms,userSpec]= get_data(varargin);
ID = (1:length(t))';

%% parse input specifications
[p,notUsingDefault] = parse_input(userSpec);

%% initial subsample data
I = subsample_catalog(t,lat,lon,depth,M,fms,parse_input({})); % completeness
[ID,t,lat,lon,depth,M,fms] = goodind(I,ID,t,lat,lon,depth,M,fms);

%% shuffle/synthetic time (probably not doing anthing unless specified)
t = redefine_time(t,p);

%% analyse data
[mainshockIndex,numberOfAftershocks,numberOfForeshocks] = ...
    get_aftershock_count(t,lat,lon,depth,M,fms,p); % input must be passed up

%% subsample data
[ID,MSt,MSlat,MSlon,MSdepth,MSmag,MSfms] = goodind(mainshockIndex,ID,t',lat',lon',depth',M',fms');
I       = subsample_catalog(MSt',MSlat',MSlon',MSdepth',MSmag',MSfms',p);
[ID,MSt,MSlat,MSlon,MSdepth,MSmag,MSfms,MSprod,FSprod] = goodind(I,ID,MSt,MSlat,MSlon,MSdepth,MSmag,MSfms, numberOfAftershocks',numberOfForeshocks'); %#ok<ASGLU>

%% get productivity law
[MAG,MAGCOUNT,prefactor,alpha]          = productivity_law(MSmag,MSprod); % aftershock

%% get residual productivity
[MSres]     = getMSres(MSmag,MSprod,prefactor,alpha);
%% output
plot_output(MAG,MAGCOUNT,prefactor,alpha,p,notUsingDefault);   
varargout   = create_output(prefactor,alpha,p, notUsingDefault, MSmag, FSprod,M,MAG,MAGCOUNT);
if ~strcmp(p.Results.SaveCatalog,   'no'); save(p.Results.SaveCatalog,'MSt','MSlat','MSlon','MSdepth','MSmag','MSfms','MSprod','MSres','FSprod'); end
if ~strcmp(p.Results.ReturnCatalog, 'no')
    CAT = table(ID,MSprod,MSres,FSprod);
    CAT.Properties.UserData = p; % append the run information for reference
    varargout = {CAT,varargout{:}};    %#ok<CCAT>
end
%% nested functions
function    [p,notUsingDefault] = parse_input(input)

%% parse input

% 1) location
%     - lat lon
%     - geological environment
%         * type
%         * distance
%     - depth
% 2) time
% 3) magnitude
% 4) focal mechanism
% More: 
% - The code needs to be able to take on a combination of the above.
% - The needs to be good defaults

p =  inputParser;
p.FunctionName = 'aftershock_productivity_kernel';

%% default settings

% 1) subsampling

defaultLatRange         	= minmax(lat);  % everything
defaultLonRange             = minmax(lon);  % everything
defaultTimeRange            = minmax(t);    % everything
defaultDepthRange           = minmax(depth);% everything
defaultCrustThicknessRange  = [0, 500];  % conservative bounds
defaultMagRange             = 'completeness';% completeness of the catalog
defaultFocalMechanism       = 'all';        % everything  (options below)
expectedFocalMechanism      = {'all','normal','reverse','strike slip'};

defaultPlateBoundary        = 'all';        % everything  (options below)
expectedPlateBoundary       = {'all','plate boundaries','divergent','convergent','transform','none'};
defaultPlateBoundaryClass   = 'all';
defaultPlateBoundaryDist    = 100; % km     
expectedPBClass             = {'all','OSR','OTF','OCB','CRB','CTF','CCB','SUB'}; 
% Boundary types: CCB continental convergent boundary, CTF continental
% transform fault, CRB continental rift boundary, OSR oceanic spreading
% ridge, OTF oceanic transform fault, OCB oceanic convergent boundary, SUB
% subduction zone

% 2) analysis
defaultTimeWindow           = 100;      % days 
defaultTimeSelectionWindow  = 60;      % days
defaultForeshockTimeWindow  = 1;
defaultSearchRadius         = 5;        % times source dimensions
defaultSelectionRadius      = 4;        % times source dimensions
defaultMinMainshockMag      = 'min';   % minimum magnitude of mainshovk to analyse


% 3) options

defaultSaveData         = 'mainshock_info.mat'; % alternatively an actual file name
defaultReturnCatalog  = 'no';

% 4) plotting output
defaultPlotYN           = 'yes';        % plot a figure of the output
defaultNewPlotYN        = 'no';         % plot output on a new figure (without hold on)
defaultShowHistogramYN  = 'no';         % show a residuals histogram (to do)
defaultSaveFigureYN     = 'no';         % save the figure
defaultPrintFigureYN    = 'no';         % print a paper version of the figure
defaultShowOverviewYN   = 'no';         % show a figure that shows an overview of the catalog

% 5) error analysis

defaultShuffleYN        = 'no';         % shuffle the time vector
defaultSyntheticTimeYN  = 'no';

expectedYN              = {'yes','no'}; % expected input for the above inputs

%% error messages

% check if inputs are numbers and positive
numericErrorMessage     = @(x) sprintf('%s input value must be positive, scalar, and numerical.',class(x));
numericValidationFcn    = @(x) assert(isnumeric(x), numericErrorMessage(x));

% 'N'  of magnitude 'magnitude'. TO deal with background seismicity a fit a
% peicewies function
% check if inputs are numbers and positive and length 2
numeric2ErrorMessage     = @(x) sprintf('%s input value must be scalar, numerical and of length 2.',class(x));
numeric2ValidationFcn    = @(x) assert(isnumeric(x) && length(x) == 2, numeric2ErrorMessage(x));

% check if inpput is one of a string
stringErrorMessage      = @(x) ['Input must be one of:', sprintf(' ''%s''',x{:})];
stringValidationFcn     = @(x,y) assert(any(strcmp(x,y)),stringErrorMessage(y));

%% add parameters (in respective order)
addParameter(p,'LatRange',          defaultLatRange,    numeric2ValidationFcn);
addParameter(p,'LonRange',          defaultLonRange,    numeric2ValidationFcn);
addParameter(p,'timeRange',         defaultTimeRange,   numeric2ValidationFcn);
addParameter(p,'MagRange',          defaultMagRange,    numeric2ValidationFcn);
addParameter(p,'DepthRange',        defaultDepthRange,  numeric2ValidationFcn);
addParameter(p,'CrustThicknessRange',defaultCrustThicknessRange,numeric2ValidationFcn);
addParameter(p,'FocalMechanism',    defaultFocalMechanism,      @(x) stringValidationFcn(x,expectedFocalMechanism));
addParameter(p,'PlateBoundary',     defaultPlateBoundary,       @(x) stringValidationFcn(x,expectedPlateBoundary)); 
addParameter(p,'PlateBoundaryDist', defaultPlateBoundaryDist,   numericValidationFcn);
addParameter(p,'PlateBoundaryClass',defaultPlateBoundaryClass,  @(x) stringValidationFcn(x,expectedPBClass));

addParameter(p,'SaveCatalog',       defaultSaveData)

addParameter(p,'TimeWindow',        defaultTimeWindow,          numericValidationFcn);
addParameter(p,'TimeSelectionWindow',defaultTimeSelectionWindow,numericValidationFcn);
addParameter(p,'ForeshockTimeWindow',defaultForeshockTimeWindow,numericValidationFcn);
addParameter(p,'SearchRadius',   	defaultSearchRadius,        numericValidationFcn);
addParameter(p,'SelectionRadius',  	defaultSelectionRadius,     numericValidationFcn);
addParameter(p,'MinMainshockMag',  	defaultMinMainshockMag,     numericValidationFcn);

addParameter(p,'PlotYN',            defaultPlotYN,          @(x) stringValidationFcn(x,expectedYN));
addParameter(p,'NewPlotYN',       	defaultNewPlotYN,       @(x) stringValidationFcn(x,expectedYN));
addParameter(p,'ShowHistogram',     defaultShowHistogramYN,	@(x) stringValidationFcn(x,expectedYN));
addParameter(p,'SaveFigureYN',    	defaultSaveFigureYN,    @(x) stringValidationFcn(x,expectedYN));
addParameter(p,'PrintFigureYN',     defaultPrintFigureYN,   @(x) stringValidationFcn(x,expectedYN));
addParameter(p,'ShowOverviewYN',    defaultShowOverviewYN,  @(x) stringValidationFcn(x,expectedYN));
addParameter(p,'ReturnCatalog',     defaultReturnCatalog,   @(x) stringValidationFcn(x,expectedYN));


addParameter(p,'ShuffleYN',         defaultShuffleYN,       @(x) stringValidationFcn(x,expectedYN));
addParameter(p,'SyntheticTimeYN',   defaultSyntheticTimeYN, @(x) stringValidationFcn(x,expectedYN));

parse(p,input{:});

if ~(length(p.Parameters) == length(p.UsingDefaults))
    
    notUsingDefault = cell(1,length(p.Parameters)-length(p.UsingDefaults));
    for n = 1:length(notUsingDefault)
        ind = ~any(strcmp(p.Parameters{n},p.UsingDefaults));
        notUsingDefault{n} = string(p.Parameters{ind});
    end
    
else
    notUsingDefault = {};
end 
end
end

%% main blocks
function    [t,lat,lon,depth,M,fms,userSpec]    = get_data(userInput)
if isa(userInput{1},'char')
    fileName        = userInput{1};
    % fileName        = 'eq_catalog.mat';
    catType         = 'GCMT';
    catFileName     = 'GCMT';
    [t,lat,lon,depth,M,fms] = load_data(catType, catFileName, fileName);
    userSpec   = userInput(2:end);

else % a bit convoluted here
    try numericInput    = isa([userInput{1:5}],'double'); catch; numericInput 	 = 0; end
    try numericInputwFMS= isa([userInput{1:6}],'double'); catch; numericInputwFMS= 0; end 
    
    if numericInput
        t   = userInput{1}';
        lat = userInput{2}';
        lon = userInput{3}';
        depth=userInput{4}';
        M   = userInput{5}';
        fms = zeros(length(t),1)';
        
        userSpec = userInput(6:end);
    end
    
    if numericInputwFMS
        fms = userInput{6}'; 
        userSpec = userInput(7:end);
    end
    
    if ~numericInput 
    error(['input structure not recognized.', ...
           'Input must either specify a filename, or include', ...
           'time, lat, lon, depth, Magnitude, and (optionally) fms']);
    end
end

end 
function    [t,lat,lon,depth,M,fms]             = load_data(catType,catFileName,fileName)

if isfile(fileName)
    savedData = load(fileName);
    t       = savedData.t;
    lat     = savedData.lat;
    lon     = savedData.lon;
    depth   = savedData.depth;
    M       = savedData.M;
    fms     = savedData.fms;
    
else
    % load and get data from catalog
    try     [t,lat,lon,depth,M, ~,fms]  = readcatalog(catType,catFileName);
    catch;  [t,lat,lon,depth,M, ~]      = readcatalog(catType,catFileName);
    end
    
    save(fileName,'t','lat','lon','depth','M','fms');
    
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
        
        if strcmp(catType, 'GCMT')
            gCMT    = load(fileName);                                   % gCMT must be on the working path 'GCMT.mat'
            gCMT    = gCMT.gCMT; % yuck
            
            
            goodInd =  (gCMT.FMS_Type == 1 | gCMT.FMS_Type == 2 | gCMT.FMS_Type == 3); % remove magnitude less than threshold, and without focal mechanism solution
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
        varargout{1} = fms; 
        
    end
end
function  	I                                   = subsample_catalog(t,lat,lon,depth,M,fms, userInput)

    subsample = userInput.Results;
    
    % space-time-magnitude window
    ISTM      = [	lat > subsample.LatRange(1)     ;... 
                    lat < subsample.LatRange(2)     ;... 
                    lon > subsample.LonRange(1)     ;...
                    lon < subsample.LonRange(2)     ;...
                    t   > subsample.timeRange(1)    ;...
                    t   < subsample.timeRange(2)    ;...
                    depth>subsample.DepthRange(1)   ;...
                    depth<subsample.DepthRange(2)   ];
    % completeness
    if strcmp(subsample.MagRange, 'completeness')
        Icomplete   = M >= get_completeness(M);
    else 
        Icomplete   = [M > subsample.MagRange(1)    ;...
                       M < subsample.MagRange(2)];
    end
    
    % focal mechanism
    if ~strcmp(subsample.FocalMechanism,'all')
        Ifms        = get_fms(fms, subsample.FocalMechanism);
    else
        Ifms        = ones(1,length(fms));
    end
    
    % plate boundary (would be worth merging the below into one - some redundancy)
    if ~strcmp(subsample.PlateBoundary,'all')
        IplateBoundary = get_PB_eq(lat, lon, subsample.PlateBoundaryDist , subsample.PlateBoundary);
    else
        IplateBoundary = ones(1,length(lat));
    end
    
    % plate boundary class
    if ~strcmp(subsample.PlateBoundaryClass,'all')
        IPBclass        = get_PB_class(lat,lon, subsample.PlateBoundaryDist , subsample.PlateBoundaryClass);
    else     
        IPBclass = ones(1,length(lat));
    end
    % crustal thickness
    crustThickness = get_crust_thickness(lat,lon,fms);
    IcrustThickness= [subsample.CrustThicknessRange(1)< crustThickness ; ...
                      subsample.CrustThicknessRange(2)> crustThickness ];
                  
    % collapse all indices:
    Iaggregate = [ ISTM             ;...
                   Icomplete        ;...
                   Ifms             ;...
                   IplateBoundary   ;...
                   IPBclass         ;...
                   IcrustThickness  ];
    
    I = ~any(~Iaggregate);
    
    % display an overview of the catalog
    if strcmp(subsample.ShowOverviewYN, 'yes')
        catalog_overview(t,lat,lon,depth,M);
    end             
     
    function ind = get_PB_eq(lat, lon, criticalDistance, plateBoudaryType)
        %  load file with the plate boundary information - right now
        %  this is just the surface expression of plate boundaries,
        %  still (could still be intraplate)
        
        if      strcmp(plateBoudaryType,'convergent');  PBType = 1;
        elseif  strcmp(plateBoudaryType,'divergent');   PBType = 2;
        elseif  strcmp(plateBoudaryType,'transform');   PBType = 3;
        else;                                           PBType = 0;
        end
        
        %  An updated digital model of plate boundaries
        %  Authors
        %  Peter Bird
        
        fn  = 'PB2002_steps.dat';
        fileID = fopen(fn);
        C = textscan(fileID,'%f %s %f %f %f %f %f %f %f %f %f %f %f %f %s');
        lonPB = C{3};
        latPB = C{4};
        divPB       = C{11};
        lateralPB   = C{12};
        
        if PBType ~= 0
            
            numData = length(lonPB);
            
            PBMotion    = [divPB,lateralPB];
            absPBMotion = abs(PBMotion);
            
            PBClassification = zeros(numData,1);
            
            PBClassification(absPBMotion(:,1)>absPBMotion(:,2) & PBMotion(:,1)<0)= 1; % convergent
            PBClassification(absPBMotion(:,1)>absPBMotion(:,2) & PBMotion(:,1)>0)= 2; % divergent
            PBClassification(absPBMotion(:,1)<absPBMotion(:,2))                  = 3; % transform
            
            [latPB, lonPB] = goodind(PBClassification == PBType, latPB, lonPB);
            
        end
        
        distance2PB = get_dist(lat, lon, latPB, lonPB);
        ind = distance2PB<criticalDistance;
        
        if strcmp(plateBoudaryType,'none')
            ind = ~ind;
        end
    end
    function ind = get_PB_class(lat, lon, criticalDistance, plateBoudaryClass)
        fn      = 'PB2002_steps.dat';
        fileID  = fopen(fn);
        C       = textscan(fileID,'%f %s %f %f %f %f %f %f %f %f %f %f %f %f %s');
        PBclass = C{end};
        %expectedPBClass             = {'all','OSR','OTF','OCB','CRB','CTF','CCB','SUB'}; 
        % Boundary types: CCB continental convergent boundary, CTF continental
        % transform fault, CRB continental rift boundary, OSR oceanic spreading
        % ridge, OTF oceanic transform fault, OCB oceanic convergent boundary, SUB
        % subduction zone 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%% refine classification %%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        PBI     = contains(PBclass,plateBoudaryClass);
        [latPB,lonPB] = goodind(PBI,C{[4,3]});
        
        ind     = get_dist(lat,lon,latPB,lonPB) < criticalDistance;
    end
    function DISTS = get_dist(lat1, lon1, lat2, lon2)
        
        num1 = length(lat1);
        num2 = length(lat2);
        
        dists     = zeros(1,num1);
        
        for n = 1:num1
            LAT1 = repmat(lat1(n),num2,1);
            LON1 = repmat(lon1(n),num2,1);
            dist1to2 = distance('gc',LAT1, LON1, lat2, lon2);
            dists(n) = min(dist1to2);
        end
        
        % convert to metric distances
        DISTS = deg2km(dists);
    end
    function ind = get_fms(fms,fmsChoice)
        if      strcmp(fmsChoice,'strike slip');        fmsChoice = 1;
        elseif  strcmp(fmsChoice,'normal');             fmsChoice = 2;
        elseif  strcmp(fmsChoice,'reverse');            fmsChoice = 3;
        end
        ind = fms == fmsChoice;    
    end

end
function    TIMEOUT                             = redefine_time(TIMEIN,USERSPEC)

    % shuffle catalog
    if strcmp(USERSPEC.Results.ShuffleYN,'yes')
        TIMEOUT = TIMEIN(randperm(length(TIMEIN)));    
    % generate synthetic time vector
    elseif strcmp(USERSPEC.Results.SyntheticTimeYN,'yes')
        TIMEOUT = rand(1,length(TIMEIN))*range(TIMEIN);       
    else
        TIMEOUT = TIMEIN;        
    end

end
function    [mainShockIndex,numberOfAftershock, numberOfForeshocks] = ...
                                                  get_aftershock_count(t,lat,lon,depth,M,fms,userInput) %#ok<INUSL>
% get the aftershock productivity of earthquakes
% given a catalog with time (t), location (lat, lon), depth and magnitude
% (M)

% out:
% out{1}: indexes of the mainshocks in the catalog
% out{2}: number of aftershocks of the earthquakes
% out{3}: number of foreshocks

x=lat;
y=lon;
z=depth;

mag=M; tsort=t;

% sort earthquakes with respect to size
[mag,I]=sort(mag,'descend');
tsort=tsort(I);
x   = x(I);
y   = y(I);
z   = z(I);
%fms = fms(I);

% convert coordinates into earth centric coordinate
wgs84       = wgs84Ellipsoid('kilometer');
[x,y,z]     = geodetic2ecef(wgs84,x,y,-z);

%%
if strcmp(userInput.Results.MinMainshockMag,'min')
        minmag=min(mag); 
else
        minmag = userInput.Results.MinMainshockMag;
end

maxmag=max(mag);

flag                =zeros(size(tsort));
ijk                 =0;

mainShockIndex      = zeros(1,length(tsort));
numberOfAftershock  = zeros(1,length(tsort));
numberOfForeshocks  = zeros(1,length(tsort));
iMainShock = 0;

Dt              = userInput.Results.TimeSelectionWindow; % time window in days 
DtFS            = userInput.Results.ForeshockTimeWindow;
Dtbig           = userInput.Results.TimeWindow;          % this is the time search window (eq flagged);
DmaxScaling     = userInput.Results.SearchRadius;        % factor of source dimension
DScaling        = userInput.Results.SelectionRadius;     % factor of source dimension

% iterate from largest magnitude down in a hiearchical way
for mm=maxmag:-0.1:minmag
    
    ijk        	= ijk+1;
    NumMain    	= 0;
    I1       	= find(abs(mag-mm)<0.05);
    Nparfor    	= length(I1);
    
    mo          = (10^((mm+10.73)*1.5))*10^-7;
    stressDrop  = 10^5*10;
    sourceDim   = 2 * 0.001 * (mo*(7/16) / stressDrop) ^ (1/3); % in km tbl 9.1 modern global seismo lay
    Dmaxbig     = DmaxScaling*sourceDim;
    DistCutoff  = DScaling*sourceDim;  

    
    if (Nparfor>0)
        for i =1:Nparfor
            if (flag(I1(i))==0)

                t0    = tsort(I1(i));
                x0  = x(I1(i));
                y0  = y(I1(i));
                z0  = z(I1(i));
                
                trel  = tsort-t0;

                Dsquared = (x-x0).^2+(y-y0).^2+(z-z0).^2;
                
                Ibig  = ((trel>-DtFS)   &   (trel<Dtbig)) ...
                    & ...
                    Dsquared < Dmaxbig^2;

                IAS   = (trel>-DtFS)  &   (trel<Dt)   ...
                    & ...
                    Dsquared < DistCutoff^2;

                IAS(I1(i))      = 0; % guard against including mainshock
                NumMain         = NumMain+1;
                
                % record the index of the mainshock in the original input array
                eqSeqt                          = trel(IAS);
                iMainShock                      = iMainShock + 1;
                mainShockIndex(iMainShock)      = I(I1(i)); % (I is the hash map for sorting and I1(i) is the main shock index in the sorted list of eq)
                numberOfAftershock(iMainShock)  = sum(eqSeqt > 0);
                numberOfForeshocks(iMainShock)  = length(eqSeqt) - numberOfAftershock(iMainShock);
                flag = flag|Ibig;
                
            end
        end
    end    
end

mainShockIndex      = mainShockIndex(1:iMainShock);
numberOfAftershock  = numberOfAftershock(1:iMainShock);
numberOfForeshocks  = numberOfForeshocks(1:iMainShock);

end
function    [magArray,prodArray,K,A]            = productivity_law(MAGNITUDE,AFTERSHOCKS)

increment = 0.1; 
minmaxMag = minmax(MAGNITUDE');

magArray    = minmaxMag(1):increment:minmaxMag(2);
numMag      = length(magArray);
prodArray   = zeros(1,numMag);

for iM = 1:(numMag-1)
    % fit a normal distribution parameters mu and sigma to the
    % distribution of earthquake productivities:
    asSub    = AFTERSHOCKS(MAGNITUDE>=magArray(iM) & MAGNITUDE<magArray(iM+1));
    prodArray(iM) = median(asSub);
end

[K,A] = getproductivity(magArray,prodArray);

end
function                                          plot_output(MAG,MAGCOUNT,prefactor,alpha, userInput,~)

userSpec = userInput.Results;

if strcmp(userSpec.PlotYN,'no');        return; end
if strcmp(userSpec.NewPlotYN,'yes');    figure; end
fig = gcf;
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for the time being this will be done manually
% legendEntry = make_legend_entry(userInput,notUsingDefault);

% legendEntry  = sprintf('Distance: %f km',userInput.Results.PlateBoundaryDist);
% legendEntry  = sprintf('Depth: %0.2g to %f km',userInput.Results.DepthRange(1),userInput.Results.DepthRange(2));
 legendEntry  = sprintf('Time window is %0.1g',userInput.Results.TimeWindow);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scatterHandle   = scatter(MAG, MAGCOUNT,'filled', ...
    'DisplayName',legendEntry1);

lineLegend  = sprintf('N_{eq} = %0.1g * 10^{%0.2g M}',prefactor, alpha);
lineColor   = get(scatterHandle,'CData');

lineHandle = plot(minmax(MAG),prefactor*10.^(alpha*minmax(MAG)), ...
    'DisplayName',lineLegend);
set(lineHandle, 'Color',    lineColor   , ...
                'LineWidth',2           );

title('Mean number of aftershock as function of magnitude');

xlabel('Magnitude')
ylabel('Average Number of Eq')
set(gca, 'YScale', 'log')
orient(fig,'landscape')
grid on

% save the figure
if strcmp(userSpec.SaveFigureYN,'yes')
    % addParameter(p,'SaveFigureYN',    	defaultSaveFigureYN,    @(x) stringValidationFcn(x,expectedYN));
    print(sprintf('mean_number_of_eq_%s',date), '-dpdf')
end

% print the figure
if strcmp(userSpec.PrintFigureYN,'yes')
    % addParameter(p,'PrintFigureYN',     defaultPrintFigureYN,   @(x) stringValidationFcn(x,expectedYN));
    print('-dprnc')
end

legendHandle = legend('show');
set(legendHandle,'Location','eastoutside')

end 
function    MSRES                               = getMSres(MSMAG,MSPROD,PREFACTOR,ALPHA)
% get the residual productivity of each mainshock         
        allAS       = @(x,y) log10(y) - log10(PREFACTOR*10.^(ALPHA*x));
        MSRES       = allAS(MSMAG,MSPROD);
        
end
function OUT = create_output(prefactor,alpha,p, notUsingDefault,MSmag, FSprod,M,MAG,MAGCOUNT)

OUT{1} = prefactor; % aftershocks
OUT{2} = alpha;

%% *** optional (foreshocks)
if p.Results.ForeshockTimeWindow ~= 0 && ~isempty(FSprod)
    [FSMAG,FSMAGCOUNT,FSprefactor,FSalpha]  = productivity_law(MSmag,FSprod); % foreshocks
    plot_output(FSMAG,FSMAGCOUNT,FSprefactor,FSalpha,p,notUsingDefault);
    OUT{3} = FSprefactor; % foreshocks
    OUT{4} = FSalpha;
else
    OUT{3} = nan; % foreshocks
    OUT{4} = nan;
end

OUT{5} = min(M); % completeness
OUT{6} = MAG;
OUT{7} = MAGCOUNT;

end
%% accessory functions

function [varargout]    = goodind(goodInd,varargin)

for n=1:length(varargin)
    varargout{n} = varargin{n}(goodInd);  %#ok<AGROW>
end

end
function mth            = get_completeness(M)

edges       = floor(min(M)):0.1:ceil(max(M));
[N,magBin]  = histcounts(M,edges);
mth=magBin(find(N==max(N),1))+0.3;

% if magBin(find(N==max(N),1)) > 4
% mth = 5.0;
% else 
% mth = 2.5;
% end

end % obtain the completeness of the earthquake catalog
% function [Il]           = ind2logind(I,L)
% Il = zeros(1,L);
% Il(I) = 1;
% end