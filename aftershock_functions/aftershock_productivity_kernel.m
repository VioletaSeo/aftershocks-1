function aftershock_productivity_kernel(varargin)

catType         = 'GCMT';
catFileName     = 'GCMT';
fileName        = 'eq_catalog.mat';

%% get data
[t,lat,lon,depth,M,fms] = load_data(catType, catFileName, fileName);

%% parse input
[p,notUsingDefault] = parse_input(varargin);

%% subsample data
[t,lat,lon,depth,M,fms] = subsample_catalog(t,lat,lon,depth,M,fms,p);

%% analyse data
[prefactor,alpha,MAG, MAGCOUNT,mainshockIndex,numberOfAftershocks]=get_prod_parameters(t,lat,lon,depth,M,p); % inpute must be passed up

%% save mainshock catalog
save_mainshock_info(mainshockIndex,numberOfAftershocks,t,lat,lon,depth,M,fms)

%% plot output
plot_output(MAG,MAGCOUNT,prefactor,alpha,p,notUsingDefault)         

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

defaultLatRange         = minmax(lat);  % everything
defaultLonRange         = minmax(lon);  % everything
defaultTimeRange        = minmax(t);    % everything
defaultDepthRange       = minmax(depth);% everything
defaultCrustThicknessRange = [0, 500];  % conservative bounds
defaultMagRange         = 'completeness';% completeness of the catalog
defaultFocalMechanism   = 'all';        % everything  (options below)
expectedFocalMechanism  = {'all','normal','reverse','strike slip'};

defaultPlateBoundary    = 'all';        % everything  (options below)
expectedPlateBoundary   = {'all','plate boundaries','divergent','convergent','transform','none'};
defaultPlateBoundaryDist= 100;          % km

% 2) analysis
defaultTimeWindow           = 800;      % days (for mag 9 event)
defaultTimeSelectionWindow  = 600;      % days (for mag 9 event)
defaultSearchRadius         = 5;        % times source dimensions
defaultSelectionRadius      = 3;        % times source dimensions

% 3) options

defaultSaveData         = 'no';         % alternatively an actual file name

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
numericValidationFcn    = @(x) assert(isnumeric(x) && ~any(x<0),numericErrorMessage(x));

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
addParameter(p,'FocalMechanism',    defaultFocalMechanism,  @(x) stringValidationFcn(x,expectedFocalMechanism));
addParameter(p,'PlateBoundary',     defaultPlateBoundary,   @(x) stringValidationFcn(x,expectedPlateBoundary)); 
addParameter(p,'PlateBoundaryDist', defaultPlateBoundaryDist,   numericValidationFcn);

addParameter(p,'SaveCatalog',       defaultSaveData)

addParameter(p,'TimeWindow',        defaultTimeWindow,          numericValidationFcn);
addParameter(p,'TimeSelectionWindow',defaultTimeSelectionWindow,numericValidationFcn);
addParameter(p,'SearchRadius',   	defaultSearchRadius,        numericValidationFcn);
addParameter(p,'SelectionRadius',  	defaultSelectionRadius,     numericValidationFcn);

addParameter(p,'PlotYN',            defaultPlotYN,          @(x) stringValidationFcn(x,expectedYN));
addParameter(p,'NewPlotYN',       	defaultNewPlotYN,       @(x) stringValidationFcn(x,expectedYN));
addParameter(p,'ShowHistogram',     defaultShowHistogramYN,	@(x) stringValidationFcn(x,expectedYN));
addParameter(p,'SaveFigureYN',    	defaultSaveFigureYN,    @(x) stringValidationFcn(x,expectedYN));
addParameter(p,'PrintFigureYN',     defaultPrintFigureYN,   @(x) stringValidationFcn(x,expectedYN));
addParameter(p,'ShowOverviewYN',    defaultShowOverviewYN,  @(x) stringValidationFcn(x,expectedYN));

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
function    [t,lat,lon,depth,M,fms]                	= load_data(catType,catFileName,fileName)

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

function  	[t,lat,lon,depth,M,fms]                     = subsample_catalog(t,lat,lon,depth,M,fms, userInput)
    
    subsample = userInput.Results;
    
    % catalog parameters
    [t,lat,lon,depth,M,fms] = goodind(lat > subsample.LatRange(1)  & lat < subsample.LatRange(2),  t,lat,lon,depth,M,fms);
    [t,lat,lon,depth,M,fms] = goodind(lon > subsample.LonRange(1)  & lon < subsample.LonRange(2),  t,lat,lon,depth,M,fms);
    [t,lat,lon,depth,M,fms] = goodind(t   > subsample.timeRange(1) & t   < subsample.timeRange(2), t,lat,lon,depth,M,fms);
    [t,lat,lon,depth,M,fms] = goodind(depth>subsample.DepthRange(1)& depth<subsample.DepthRange(2),t,lat,lon,depth,M,fms);
    
    if strcmp(subsample.MagRange, 'completeness')
        [t,lat,lon,depth,M,fms] = goodind(M>get_completeness(M),t,lat,lon,depth,M,fms);
    else 
        [t,lat,lon,depth,M,fms] = goodind(M > subsample.MagRange(1) & M < subsample.MagRange(2),t,lat,lon,depth,M,fms); 
    end
    
    % focal mechanism
    if ~strcmp(subsample.FocalMechanism,'all')
        [t,lat,lon,depth,M,fms] = goodind(get_fms(fms, subsample.FocalMechanism),t,lat,lon,depth,M,fms);
    end
    
    % plate boundary
    if ~strcmp(subsample.PlateBoundary,'all')
        [t,lat,lon,depth,M,fms] = goodind(get_PB_eq(lat, lon, subsample.PlateBoundaryDist , subsample.PlateBoundary),t,lat,lon,depth,M,fms);
    end
    
    % crustal thickness
    crustThickness = get_crust_thickness(lat,lon,fms);
    [t,lat,lon,depth,M,fms] = goodind(subsample.CrustThicknessRange(1)< crustThickness  & ...
                                      subsample.CrustThicknessRange(2)> crustThickness  , ...
                                      t,lat,lon,depth,M,fms);
                                  
    % shuffle catalog
    if strcmp(subsample.ShuffleYN,'yes')
        t = t(randperm(length(t)));
    end
    
    % generate synthetic time vector
    if strcmp(subsample.SyntheticTimeYN,'yes')
        t = rand(1,length(t))*range(t);
    end
    
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
            
        function DISTS = get_dist(lat1, lon1, lat2, lon2)
            
            num1 = length(lat1);
            num2 = length(lat2);
            
            dists     = zeros(1,num1);
            
            parfor n = 1:num1
                LAT1 = repmat(lat1(n),num2,1);
                LON1 = repmat(lon1(n),num2,1);
                dist1to2 = distance('gc',LAT1, LON1, lat2, lon2);
                dists(n) = min(dist1to2);
            end
            
            % convert to metric distances
            DISTS = deg2km(dists);
        end
    end
    function ind = get_fms(fms,fmsChoice)
        if      strcmp(fmsChoice,'strike slip');        fmsChoice = 1;
        elseif  strcmp(fmsChoice,'normal');             fmsChoice = 2;
        elseif  strcmp(fmsChoice,'reverse');            fmsChoice = 3;
        end
        ind = fms == fmsChoice;    
    end

end
    
function    [prefactor,alpha,MAG,MAGCOUNT,varargout]= get_prod_parameters(t,lat,lon,depth,M, userInput)
% get parameter k (the prefactor), aplha (the slope - often 1), and the
% average number of aftershock (MAGCOUNT) as a function of magnitude (MAG)
% associtated to the aftershock productivity of earthquakes
% given a catalog with time (t), location (lat, lon), depth and magnitude
% (M)

% varargout:
%   varargout{1}: indexes of the mainshocks in the catalog
%   varargout{2}: number of aftershocks of the earthquakes 



x=lat;
y=lon;
z=depth;

mag=M; tsort=t;

% sort earthquakes with respect to size
%[tsort,I]=sort(tsort);
%mag=mag(I);
[mag,I]=sort(mag,'descend');
tsort=tsort(I);
x=x(I);
y=y(I);

Dtfor=0;            % foreshock time window

%% 
minmag=min(mag); %just search for all magnitudes

flag=tsort*0;
ijk=0;

% define structure
Rec(1).DD=[];
Rec(1).Mag=[];
Rec(1).Dfor=[];% aftershock_productivity_kernel('ShowOverviewYN','yes',...
%                                'NewPlotYN','yes');

Rec(1).C=[];
Rec(1).NumMain=[];
Nmag=ceil((max(mag)-minmag)*10);
Rec(Nmag)=Rec(1);
maxmag=max(mag);

mainShockIndex      = zeros(1,length(tsort));
numberOfAftershock  = zeros(1,length(tsort));
iMainShock = 0;


dt9max          = userInput.Results.TimeWindow;          % this is the time search window (eq flagged);
dt9             = userInput.Results.TimeSelectionWindow; % time window in days
DmaxScaling     = userInput.Results.SearchRadius;        % factor of source dimension
DScaling        = userInput.Results.SelectionRadius;     % factor of source dimension

% iterate from largest magnitude down in a hiearchical way
for mm=maxmag:-0.1:minmag
    
    ijk             = ijk+1;
    Rec(ijk).Mag    = mm;
    NumMain         = 0;
    I1              = find(abs(mag-mm)<0.05);
    Nparfor         = length(I1);
    
    if mm == 7.5
        
    end
    
    if (Nparfor>0)
        clear Recm;
        Recm(Nparfor).DD    =[];
        Recm(Nparfor).Dfor  =[];
        Recm(Nparfor).C     =[];
        
        for i =1:Nparfor
            
            t0    = tsort(I1(i));
            lat0  = x(I1(i));
            lon0  = y(I1(i));
            
            trel  = tsort-t0;
            
            D     =distance(x,y,lat0,lon0);
            D     =deg2km(D);            
            
            % scales with source dimension:
            mo          = (10^((mm+10.73)*1.5))*10^-7;
            stressDrop  = 10^5*30;
            sourceDim   = 0.001*(mo*(7/16)/(stressDrop))^(1/3); % in km tbl 9.1 modern global seismo lay
            Dmaxbig     = DmaxScaling*sourceDim; 
            DistCutoff  = DScaling*sourceDim;
            
            Dt          = dt9*10^(mm-9);
            Dtbig       = dt9max*10^(mm-9);
            
            Itbig =((trel>0)&(trel<Dtbig));
            Idbig =(abs(D)<Dmaxbig);
            Ibig  =and(Itbig,Idbig);
            
    
    
            if (flag(I1(i))==0)
                IAS             =(trel>0) &   (trel<Dt)   &   (D<DistCutoff);
                Ifor          =(trel<0)&(trel>Dtfor)&(D<10);
                IAS(I1(i))      =0; % guard against including mainshock
                Ifor(I1(i))   =0;
                %    if ((sum([flag(I) flag(Ifor)]))<1)   % none of the aftershocks or foreshocks are flagged
                NumMain=NumMain+1;
                Recm(i).DD=D(IAS);
                Recm(i).Dfor=D(Ifor);
                
                % record the index of the mainshock in the original input array
                iMainShock = iMainShock + 1;
                mainShockIndex(iMainShock)      = I(I1(i)); % (I is the hash map for sorting and I1(i) is the main shock index in the sorted list of eq) 
                numberOfAftershock(iMainShock) 	= sum(IAS);
            end
            
            
            flag(Ibig)=1;
        end
        
        Rec(ijk).DD=[Recm.DD];
        Rec(ijk).Dfor=[Recm.Dfor];
        
    end
    
    Rec(ijk).C          = length(Rec(ijk).DD)/Nparfor; % averaging
    Rec(ijk).NumMain    = NumMain;               % recorn the number of main shocks
end

% trim outputs
ijklast=ijk;
Rec=Rec(1:ijklast);
mainShockIndex      = mainShockIndex(1:iMainShock);
numberOfAftershock  = numberOfAftershock(1:iMainShock);

varargout{1} = mainShockIndex;
varargout{2} = numberOfAftershock;

[prefactor,alpha] = getproductivity([Rec.Mag], [Rec.C]); 
% 'threshold1' indicates to tuncate the model fit to earthqakes that
% generate, on average, more than 1 earthquake

MAG         = [Rec.Mag];
MAGCOUNT    = [Rec.C];

end

function    save_mainshock_info(mainshockIndex,numberOfAftershocks,t,lat,lon,depth,M,fms)
% save the mainshock data for later call (this simplfies the bookkeeping,
% though it is not ideal)

% input: 
% mainshockIndex, binary array of indices indicating the location if the
% mainshocks withing the space-time-magnitude information of the earthquake
% catalog

% numberOfAftershocks, corresponding productivity of each mainshock

% space-time-mag-fms infomation of all earthquakes:
%t,
%lat,
%lon,
%depth,
%M,
%fms


[MSt,MSlat,MSlon,MSdepth,MSmag,MSfms] = goodind(mainshockIndex,t,lat,lon,depth,M,fms);
MSprod = numberOfAftershocks;
save('mainshock_info.mat','MSt','MSlat','MSlon','MSdepth','MSmag','MSfms','MSprod');

end

function    plot_output(MAG,MAGCOUNT,prefactor,alpha, userInput,notUsingDefault)

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
    'DisplayName',legendEntry);

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

%     function legendEntry = make_legend_entry(userInput,notUsingDefault)
%         
%         if isempty(notUsingDefault) == 0
%             legendEntry = 'Entire Catalog';
%         elseif length(notUsingDefault) == 1
%             legendEntry = '';
%             parameterVal = userInput.Results.(notUsingDefault{:});
%             if any(strcmp(notUsingDefault{:},{'DepthRange','timeRange','LatRange','LonRange','MagRange'}))
%                 legendEntry = [legendEntry,sprintf('%s between %f and %f;', ...
%                     parameter,parameterVal(1), parameterVal(2))];
%             elseif any(strcmp(notUsingDefault{:},{'FocalMechanism','PlateBoundary'}))
%                 legendEntry = [legendEntry,sprintf('%s is %s', ...
%                     parameter,parameterVal)];
%             elseif any(strcmp(notUsingDefault{:},{'PlateBoundaryDist', 'SearchRadius', 'SelectionRadius', 'TimeSelectionWindow','TimeWindow'}))
%                 legendEntry = [legendEntry,sprintf('%s is %f', ...
%                     parameter,parameterVal)];
%             end
%         else
%             legendEntry = '';
%             for n = 1:length(notUsingDefault)
%                 parameterInd = zeros(1,length(userInput.Parameters));
%                 for m = 1:length(userInput.Parameters)
%                     par = userInput.Parameters{m};
%                     parameterInd(m) = strcmp(notUsingDefault{n},par);
%                 end
%                 
%                 parameterInd = find(parameterInd);
%                 parameter = userInput.Parameters{parameterInd};
%                 parameterVal = userInput.(parameter);
%                 
%                 if any(strcmp(parameter,{'DepthRange','timeRange','LatRange','LonRange','MagRange'}))
%                     legendEntry = [legendEntry,sprintf('%s between %f and %f;', ...
%                         parameter,parameterVal(1), parameterVal(2))];
%                 elseif any(strcmp(parameter,{'FocalMechanism','PlateBoundary'}))
%                     legendEntry = [legendEntry,sprintf('%s is %s', ...
%                         parameter,parameterVal)];
%                 elseif any(strcmp(parameter,{'PlateBoundaryDist', 'SearchRadius', 'SelectionRadius', 'TimeSelectionWindow','TimeWindow'}))
%                     legendEntry = [legendEntry,sprintf('%s is %f', ...
%                         parameter,parameterVal)];
%                 end
%             end
%         end
%         
%     end

end 

%% accessory functions

function [varargout]    = goodind(goodInd,varargin)

for n=1:length(varargin)
    varargout{n} = varargin{n}(goodInd);
end

end

function mth            = get_completeness(M)

edges       = floor(min(M)):0.1:ceil(max(M));
[N,magBin]  = histcounts(M,edges);
mth=magBin(N==max(N));

end % obtain the completeness of the earthquake catalog
