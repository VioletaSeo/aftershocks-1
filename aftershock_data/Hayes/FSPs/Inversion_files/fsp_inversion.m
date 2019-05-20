classdef fsp_inversion < handle
    %fsp_inversion parse and assembles finite fault inversion parameters
    %   Detailed explanation goes here
    
    properties
        FileName
        Event
        Date % the constructor do not get the time it just gets the date
        Tag
        Location
        Size % NOT RUPTURE DIMENSIONS!
        Mechanism
        Rupture
        Parameters
        VelocityStructure
        Multisegment
        FiniteFault
    end    
    properties (Dependent)
        ElasticProperties
        FlatFiniteFault
        AspectRatio
        StressDrop
        DateNum
        RuptureDimensions
    end    
    properties (SetAccess = private, Hidden = true)
        pStressDrop = [];       
    end   
    properties (Hidden = true)
        fitType = 'spheroid';
    end
    
    methods
        % Constructor
        function FSP = fsp_inversion(FN)
            %fsp_inversion Construct an instance of this class
            %   Parse .fsp file and store finite fault inversion parameters
            
            % Read commented section starting with  % sign
            %   Add in parameters as I go:
            % 1) Event info
            % 2) Inversion related parameters
            % 3) Velocity structure
            
            if nargin == 0
                return
            end
            
            fid = fopen(FN);
            
            %initialize variables
            aftercolon  = @(in) in((find(in == ':')+2):end);
            
            % 'read_nums' looks for the numbers 2 entries down from the specifies input
            fgetl(fid);                                                             %1
            fgetl(fid);                                                             %2
            Event       = aftercolon(fgetl(fid));                                   %3
            key = find(Event == '/'); Date = Event((1:10)+(key(1)-5)); % a bit (lot) of a hack
            Tag         = aftercolon(fgetl(fid));                                   %4
            fgetl(fid);                                                             %5
            
            Loc         = read_nums(fgetl(fid), 'LAT' , 'LON', 'DEP');              %6
            Size        = read_nums(fgetl(fid), 'LEN' , 'WID', 'Mw', 'Mo');         %7
            Mech        = read_nums(fgetl(fid), 'STRK', 'DIP', 'RAKE', 'Htop');     %8
            Rupt        = read_nums(fgetl(fid), 'HypX', 'Hypz', 'avTr', 'avVr');    %9
            
            fgetl(fid);                                                             %10
            fgetl(fid);                                                             %11
            fgetl(fid);                                                             %12
            
            Invs1       = read_nums(fgetl(fid), 'Nx', 'Nz', 'Fmin', 'Fmax');        %13
            Invs2       = read_nums(fgetl(fid), 'Dx', 'Dz');                        %14
            Invs3       = read_nums(fgetl(fid), 'Ntw', 'Nsg');                      %15
            Invs4       = read_nums(fgetl(fid), 'LEN', 'SHF');                      %16
            Inv         = aggregate_struct(Invs1,Invs2,Invs3,Invs4);
            
            skipl(10,fid)
            
            numberOf    = read_nums(fgetl(fid), 'layers');
            
            skipl(3,fid)
            
            varNames = {'Depth','PVEL','SVEL','DENS','QP','QS'};
            velRaw = zeros(numberOf.layers,length(varNames));
            for n = 1:numberOf.layers
                line = fgetl(fid);
                line = line(2:end);
                splitLine = split(line);
                splitLine = splitLine((end-(length(varNames)-1)):end); % HACK
                velRaw(n,:) = str2double(splitLine)';
            end
            VelocityStruct = array2table(velRaw,'VariableNames',varNames);
            
            % skip lines till start of inversion
            
            
            if Inv.Nsg == 1; Multisegment = false;
            while 1
                line = fgetl(fid);
                spltLine = split(line);
                if ismember('Nsbfs', spltLine); Nsbfs = num2down(spltLine,'Nsbfs'); end
                
                % start inversion
                if line(1) ~= '%'
                    varNames = {'LAT','LON','X','Y','Z','SLIP','RAKE','TRUP','RISE','SF_MOMENT'};
                    finiteFaultArray = zeros(Nsbfs,numel(varNames));
                    for n = 1:Nsbfs
                        if n~=1; line = fgetl(fid); end
                        line = strtrim(line);
                        finiteFaultArray(n,:) = str2double(split(line))';
                    end
                    finiteFault = array2table(finiteFaultArray,...
                        'VariableNames',varNames);
                    break
                end % of inversion
            end % while
            else; Multisegment = true; finiteFault = []; % giving up
            end % if Inv.Nsg == 1
            
            % assign a focal mechanism
            rakeWin = 45;
            if      Mech.RAKE > (90-rakeWin) && Mech.RAKE  <= (90+rakeWin)
                        Mech.FocalMechanism = 'Reverse';
            elseif  Mech.RAKE > (270-rakeWin) && Mech.RAKE <= (270+rakeWin)
                        Mech.FocalMechanism = 'Normal';
            elseif  Mech.RAKE >=(0)             && Mech.RAKE  <=(rakeWin)       || ...
                    Mech.RAKE > (180-rakeWin)   && Mech.RAKE  <=(180+rakeWin)   || ...
                    Mech.RAKE > (360-rakeWin)   && Mech.RAKE  <=(360)
                        Mech.FocalMechanism = 'StikeSlip';
            else
                        Mech.FocalMechanism = 'Oblique';
            end
            
            % output assingment:
            FSP.FileName            = FN;
            FSP.Event               = Event;
            FSP.Date                = Date; 
            FSP.Tag                 = Tag;
            FSP.Location            = Loc;
            FSP.Size                = Size;
            FSP.Mechanism           = Mech;
            FSP.Rupture             = Rupt;
            FSP.Parameters          = Inv;
            FSP.VelocityStructure   = VelocityStruct;
            FSP.Multisegment        = Multisegment;
            FSP.FiniteFault         = finiteFault;            
        end
        function out = get.ElasticProperties(FSP)
            numGC = length(FSP.FiniteFault.Z);
            layerRef = discretize(FSP.FiniteFault.Z, ...
                [FSP.VelocityStructure.Depth',inf]);
            
            
            meanProp = [];
            propNames= FSP.VelocityStructure.Properties.VariableNames;
            numProp  = length(propNames);
            
            for n = 1:numProp
                
                propN       = FSP.VelocityStructure.(propNames{n});
                propGrid = zeros(1,numGC);
                for m = 1:numGC
                    propGrid(m)    = propN(layerRef(m));
                end
                weight      = FSP.FiniteFault.SLIP/sum(FSP.FiniteFault.SLIP);
                
                meanProp.(propNames{n}) = propGrid*weight;
            end
            
            vp        = 1000*meanProp.PVEL;
            vs        = 1000*meanProp.SVEL;
            p         = 1000*meanProp.DENS;
            
            out = [];
            out.Vp      = meanProp.PVEL*1000;
            out.Vs      = meanProp.SVEL*1000;
            out.Dens    = meanProp.DENS*1000;
            out.QP      = meanProp.QP;
            out.QS      = meanProp.QS;
            out.Lame    = p*(vp^2+2*vs^2);
            out.Young   = (p*vs^2*(3*vp^2-4*vs^2))/(vp^2-vs^2);
            out.Shear   = p*vs^2;
            out.Bulk    = p*(vp^2-4/3*vs^2);
            out.Poisson = (vp^2-2*vs^2)/(2*(vp^2-vs^2));
            out.M       = p*vp^2;
            
        end
        function out = get.FlatFiniteFault(FSP)
            
            [xGrid,yGrid, slipMap, rakeMap, moMap] = grid_slip_inv(FSP);
            X       = xGrid(:);
            Y       = yGrid(:);
            SLIP    = slipMap(:);
            RAKE    = rakeMap(:);
            MOMENT  = moMap(:);
            
            out     = table(X,Y,SLIP,RAKE,MOMENT);
            
        end
        function out = get.StressDrop(FSP) 
            if isempty(FSP.pStressDrop)
                FSP = FSP.computeStressDrop;
            end
            out = FSP.pStressDrop;
            
        end
        function out = get.DateNum(FSP)
            dateTime= datetime(FSP.Date,'InputFormat','yyyy/MM/dd');
            out     = datenum(dateTime); 
        end
        function out = get.AspectRatio(FSP)
            out     = FSP.RuptureDimensions.Length/FSP.RuptureDimensions.Width;
        end 
        function out = get.RuptureDimensions(FSP)
            % obtain rupture dimensions following autocorelation method
            [xGrid,yGrid, slipMap, ~, ~] = grid_slip_inv(FSP);
            vertStack = sum(slipMap,2);
            horzStack = sum(slipMap,1);
            
            out.Width = autocorrelation_width(yGrid(:,1),vertStack);
            out.Length= autocorrelation_width(xGrid(1,:),horzStack);
        end
        %other
        function out = SlipHeterogeneity(FSP,fitType)

%%          make padded grid to ensure that the slip is captured with model
            % get_slip_heterogeneity(paddedXGrid, paddedYGrid, paddedSlipMap, xF, yF, slip, xPtSp, yPtSp, mo, fitType)
            % IDEA: thorne style
            % produce the slip distribution for an idealized slip ditribution
            if nargin == 1
                fitType = 'ellipsoid';
            end
            
            [~,~, slipMap, ~, ~] = grid_slip_inv(FSP);
                        
            paddedSlipMap   = pad_grid(slipMap); %ADD PADDED SLIPMAP
            paddedSlipMap   = inpaint_nans(paddedSlipMap);
            paddedDim       = size(paddedSlipMap);
            
            X = FSP.FlatFiniteFault.X;
            Y = FSP.FlatFiniteFault.Y;
            
            xSprd = diff(minmax(X'));
            ySprd = diff(minmax(Y'));
            
            [paddedXGrid , paddedYGrid]     = meshgrid(linspace(min(X)-xSprd, max(X)+xSprd,paddedDim(1)),...
                linspace(min(Y)-ySprd, max(Y)+ySprd, paddedDim(2)));

            % one of the below fit routines
            
            switch fitType               
                case 'spheroid'
                    %% oblate spheroid
                    fitFun = fittype('oblate_spheroid(a,c,x0,y0,x,y)', ...
                        'dependent',{'z'}, ...
                        'independent', {'x','y'}, ...
                        'coefficient', {'a','c','x0','y0'});
                    
                    startPoint = [max(X)-min(X), max(slipMap(:)), mean(X), mean(Y)];
                    idealSlip  = fit([paddedXGrid(:),paddedYGrid(:)], paddedSlipMap(:), fitFun, ...
                        'StartPoint', startPoint);
                    
                case 'ellipsoid'
                    %% ellipsoid                    
                    fitFun = fittype('positive_ellipsoid(a,b,c,x0,y0,x,y)', ...
                        'dependent',{'z'}, ...
                        'independent', {'x','y'}, ...
                        'coefficient', {'a','b','c','x0','y0'});
                    startPoint = [max(X)-min(X), max(Y)-min(Y) max(slipMap(:)), mean(X), mean(Y)];
                    idealSlip  = fit([paddedXGrid(:),paddedYGrid(:)], paddedSlipMap(:), fitFun, ...
                        'StartPoint', startPoint);
                    
                case 'ellipsoid_fixed_moment' %  DOES NOT WORK AT THE MOMENT
                    fitFun = fittype('positive_ellipsoid_fixed_moment(a,b,mo,x0,y0,x,y)', ...
                        'dependent',{'z'}, ...
                        'independent', {'x','y','mo'}, ...
                        'coefficient', {'a','b','x0','y0'});
                    startPoint = [max(X)-min(X), max(Y)-min(Y) max(slipMap), mean(X), mean(Y)];
                    idealSlip  = fit([paddedXGrid(:),paddedYGrid(:), mo], paddedSlipMap(:), fitFun, ...
                        'StartPoint', startPoint);
                    
                otherwise
                    error('''fitType'' must be one of ''spheroid'', ''ellipsoid'', ''ellipsoid_fixed_moment''')
            end
            
            cellArea        = xSprd*ySprd/numel(slipMap);
            
            
            elpseError      = sum(abs((paddedSlipMap(:) - feval(idealSlip,[paddedXGrid(:),paddedYGrid(:)])))*cellArea);
            elpseErrorNorm  = elpseError./(xSprd*ySprd)*mean(slipMap(:)); % normalized by the area and mean slip of the fault
            %idealSlipGrid = reshape(feval(idealSlip,[paddedXGrid(:),paddedYGrid(:)]),paddedDim(1), paddedDim(2));
            
            out = elpseErrorNorm;


        end
        function out = flat(FSP)
           % create a vector with time, lat, lon, depth, Magnitude
           % infomations
           
           Time             = FSP.DateNum;
           Lat              = FSP.Location.LAT;
           Lon              = FSP.Location.LON;
           Depth            = FSP.Location.DEP;
           Mw               = FSP.Size.Mw;
           Mo               = FSP.Size.Mo;
           Width            = FSP.RuptureDimensions.Width;
           Length           = FSP.RuptureDimensions.Length;
           AspectRatio      = FSP.AspectRatio;
           Strike           = FSP.Mechanism.STRK;
           Dip              = FSP.Mechanism.DIP;
           Rake             = FSP.Mechanism.RAKE;
           RuptureVelocity  = FSP.Rupture.avVr;
           RuptureTime      = FSP.Rupture.avTr;
           StressDrop       = FSP.StressDrop;
           Heterogeneity    = FSP.SlipHeterogeneity; 
           Young = FSP.ElasticProperties.Young;
           Poisson= FSP.ElasticProperties.Poisson;
           
           
           out  = table(Time,Lat,Lon,Depth,Mw,Mo,Width, ...
                        Length,AspectRatio,Strike,Dip, ...
                        Rake,RuptureVelocity, RuptureTime, ...
                        StressDrop, Heterogeneity, Young, Poisson);           
        end
        
        % overloads
        function plot(FSP)
            [xGrid,yGrid, slipMap, rakeMap, moMap] = grid_slip_inv(FSP);
            figure; hold on;
            imagesc(xGrid(1,:),-yGrid(:,1),moMap)
            quiver(xGrid(:),-yGrid(:),slipMap(:).*cosd(rakeMap(:)),slipMap(:).*sind(rakeMap(:)),'w')
            xlabel('Strike Direction (km)')
            ylabel('Up-Dip Direction (km)')
            colorbar
            set(gca,'Ylim',[-max(yGrid(:)),-min(yGrid(:))],'DataAspectRatio',[1,1,1])
            
            title(sprintf('%s, M0=%g, Stress Drop = %g',FSP.Event, FSP.Size.Mw, FSP.StressDrop));
            
            set(findall(gcf, 'Type', 'Text'),'FontWeight', 'Normal','FontSize',12)
        end
        
    end    
    methods(Access = private)
          function FSP = computeStressDrop(FSP)
               % functionName = mfilename('fullpath');
            functionDir = '/home/kdascher/Documents/UCSC/projects/stress_drop';
            % [functionDir, ~] = fileparts(functionName);
            
            %% file location
            % keep the file location local to this function
            
            %% read in file
            numSub  = height(FSP.FlatFiniteFault);
            fn = FSP.FileName;
            %% build new fort 37 file
            newFN   = [functionDir,'/',fn,'_fort.37']; % !!! It might be worth thinking a bit more about the path here
            FID     = fopen(newFN,'w');
            dummyCount = numSub;
            fprintf(FID,'%s\n',                    fn);
            fprintf(FID,'%s\n\n',                  'Gavin_Hayes_fsp_input_file_adapted');
            fprintf(FID,'   nt    dt   h0    Ne   beta0 \n');
            fprintf(FID,'%i \t %.2f \t %g \t %i \t %.2f \n\n', ...
                [dummyCount, ...            % nt (dummy)
                dummyCount, ...             % ne (dummy)
                min(FSP.FiniteFault.Z), ...  % depth km OF TOP OF THE INVERSION *HACK*
                dummyCount, ...             % Number of subevents  (dummy)
                0.5]);                      % rupture velocity? km/s
            
            %%% a bit pointless but necessary for the parsing in stresses.f
            momentScale = 10^-18; % because of normalization in the stresses.f script
            
            fprintf(FID,'< Subevent sequence >\n');
            fprintf(FID,'time \t x \t\t\t y \t\t\t Slip-angle\t Mo \t\t slip(m)\n');
            fprintf(FID,'%i \t\t %2.2e \t %2.2e \t %.2e \t %.2e \t %.3e\n', ...
                [zeros(numSub,1), ...                      	% slip in (m)
                FSP.FlatFiniteFault.X, ...                 	% x coord (km) (in fault plane)
                FSP.FlatFiniteFault.Y, ...               	% y coord (km) (in fault plane)
                FSP.FlatFiniteFault.RAKE, ...              	% rake
                momentScale * FSP.FlatFiniteFault.MOMENT, ...% Mo (norm) (Nm) = 10^-18 * Mo
                FSP.FlatFiniteFault.SLIP].');                 	% slip in (m)
            
            % ************************************* %
            
            %%% begin actual fault info here
            fprintf(FID,'X \t\t\t Y \t\t\t Mo \t\t Rake \t\t slip_m\n');
            fprintf(FID,'%2.2e \t %2.2e \t %.3e \t %2.2e \t %.3e\n', ...
                [FSP.FlatFiniteFault.X, ...                 	% x coord (km) (in fault plane)
                FSP.FlatFiniteFault.Y, ...                 	% y coord (km) (in fault plane)
                momentScale * FSP.FlatFiniteFault.MOMENT,...% Mo (Nm) = mu*A*D
                FSP.FlatFiniteFault.RAKE, ...              	% rake
                FSP.FlatFiniteFault.SLIP].');                 	% slip in (m)
            
            fprintf(FID,'Mo [Nm], Mw, rake =   %.4e %.2g %3.2g', ...
                FSP.Size.Mo, ...
                FSP.Size.Mw, ...
                FSP.Mechanism.RAKE);
            
            
            % there is additional stuff the the fort37 files but I assume I wont need
            % them
            fclose(FID);
            
            
            %% build i_stresses file
            FID     = fopen('i_stresses','w'); % rewrite
            fprintf(FID,'%s \n',                        newFN);
            fprintf(FID,'%.2e \t %2.2e \n',             FSP.Size.Mo, FSP.Mechanism.DIP); % !!! the dip calculation is not ideal
            fprintf(FID,'%.3e \t %.3e \t %s\n\n',       FSP.ElasticProperties.Lame,  ...
                FSP.ElasticProperties.Shear, ...
                'PREM'); %  fprintf(FID,'4.573e10 6.86e10   PREM'); % alternate numbers used in Chiaps inversion
            fprintf(FID,'%s\n\n',                   '------------------------------');
            fprintf(FID,'%s\n',                     '1.  fort.37 filename');
            fprintf(FID,'%s\n',                     '2.  reference moment (Nm),   dip angle');
            fprintf(FID,'%s\n',                     '3.  Lambda, Mu       (Pa)');
            fclose(FID);
            
            %% run stresses.f MUST BE COMPILED
            [~,txtOut]  = system([functionDir,'/stresses']);
            
            %% read output
            
            C = textscan(txtOut, 'Sigma_E:  %3.2f MPa');
            STRESSDROP = C{1}; % MPa
            
            system('rm i_stresses');
            system(['rm ', newFN]);    
            
            FSP.pStressDrop = STRESSDROP;
            
          end
    end
end

function [xGrid,yGrid, slipMap, rakeMap, moMap] = grid_slip_inv(FSP)

            xyz = flatten_XYZ(...
                [FSP.FiniteFault.X,...
                FSP.FiniteFault.Y,...
                -FSP.FiniteFault.Z]);
            rotX  = xyz(:,1);
            rotY  = xyz(:,2);
            
            ang     = atan((rotY(2)-rotY(1))/(rotX(2)-rotX(1)));
            rotMat  = [cos(ang) , sin(ang); ...
                -sin(ang), cos(ang)];
            
            XPYP    = rotMat*[rotX';rotY'];
            xF      = XPYP(1,:);
            yF      = XPYP(2,:);
            
            xF      = xF - min(xF);
            yF      = yF - min(yF);
            
            xPtSp   = FSP.Parameters.Dx;
            yPtSp   = FSP.Parameters.Dz;
            
            xAxis           = min(xF): xPtSp : max(xF);
            yAxis           = min(yF): yPtSp : max(yF);
            [xGrid,yGrid]   = meshgrid(xAxis, yAxis);
            
            [slipMap, rakeMap, moMap] = interp_maps(xF,yF,xGrid,yGrid, ...
                FSP.FiniteFault.SLIP, ...
                FSP.FiniteFault.RAKE, ...
                FSP.FiniteFault.SF_MOMENT);
            
            slipMap(slipMap<0) = 0; % set extrapolations that have negative values to zero
            moMap(moMap<0) = 0;

end

function out            = num2down(cellStr,str)
out = str2double(cellStr(find(strcmp(cellStr,str))+2));
end
function structOut      = read_nums(line,varargin)
structOut = [];
for n = 1:length(varargin)
    structOut.(varargin{n}) = num2down(split(line),varargin{n});
end
end
function skipl(N,fid)
for n = 1:N
    fgetl(fid);
end
end
function structOut      = aggregate_struct(varargin)
structOut = [];
for n = 1:length(varargin)
    structN = varargin{n};
    FNS = fieldnames(structN);
    for m = 1:length(FNS)
        structOut.(FNS{m}) = structN.(FNS{m});
    end
end
end
function varargout      = interp_maps(x,y,xg,yg,varargin)

for n = 1:length(varargin)
    map = griddata(x,y,varargin{n},xg,yg);
    varargout{n} = inpaint_nans(map);
end

end
function B              = inpaint_nans(A,method)
% INPAINT_NANS: in-paints over nans in an array
% usage: B=INPAINT_NANS(A)          % default method
% usage: B=INPAINT_NANS(A,method)   % specify method used
%
% Solves approximation to one of several pdes to
% interpolate and extrapolate holes in an array
%
% arguments (input):
%   A - nxm array with some NaNs to be filled in
%
%   method - (OPTIONAL) scalar numeric flag - specifies
%       which approach (or physical metaphor to use
%       for the interpolation.) All methods are capable
%       of extrapolation, some are better than others.
%       There are also speed differences, as well as
%       accuracy differences for smooth surfaces.
%
%       methods {0,1,2} use a simple plate metaphor.
%       method  3 uses a better plate equation,
%                 but may be much slower and uses
%                 more memory.
%       method  4 uses a spring metaphor.
%       method  5 is an 8 neighbor average, with no
%                 rationale behind it compared to the
%                 other methods. I do not recommend
%                 its use.
%
%       method == 0 --> (DEFAULT) see method 1, but
%         this method does not build as large of a
%         linear system in the case of only a few
%         NaNs in a large array.
%         Extrapolation behavior is linear.
%         
%       method == 1 --> simple approach, applies del^2
%         over the entire array, then drops those parts
%         of the array which do not have any contact with
%         NaNs. Uses a least squares approach, but it
%         does not modify known values.
%         In the case of small arrays, this method is
%         quite fast as it does very little extra work.
%         Extrapolation behavior is linear.
%         
%       method == 2 --> uses del^2, but solving a direct
%         linear system of equations for nan elements.
%         This method will be the fastest possible for
%         large systems since it uses the sparsest
%         possible system of equations. Not a least
%         squares approach, so it may be least robust
%         to noise on the boundaries of any holes.
%         This method will also be least able to
%         interpolate accurately for smooth surfaces.
%         Extrapolation behavior is linear.
%
%         Note: method 2 has problems in 1-d, so this
%         method is disabled for vector inputs.
%         
%       method == 3 --+ See method 0, but uses del^4 for
%         the interpolating operator. This may result
%         in more accurate interpolations, at some cost
%         in speed.
%         
%       method == 4 --+ Uses a spring metaphor. Assumes
%         springs (with a nominal length of zero)
%         connect each node with every neighbor
%         (horizontally, vertically and diagonally)
%         Since each node tries to be like its neighbors,
%         extrapolation is as a constant function where
%         this is consistent with the neighboring nodes.
%
%       method == 5 --+ See method 2, but use an average
%         of the 8 nearest neighbors to any element.
%         This method is NOT recommended for use.
%
%
% arguments (output):
%   B - nxm array with NaNs replaced
%
%
% Example:
%  [x,y] = meshgrid(0:.01:1);
%  z0 = exp(x+y);
%  znan = z0;
%  znan(20:50,40:70) = NaN;
%  znan(30:90,5:10) = NaN;
%  znan(70:75,40:90) = NaN;
%
%  z = inpaint_nans(znan);
%
%
% See also: griddata, interp1
%
% Author: John D'Errico
% e-mail address: woodchips@rochester.rr.com
% Release: 2
% Release date: 4/15/06


% I always need to know which elements are NaN,
% and what size the array is for any method
[n,m]=size(A);
A=A(:);
nm=n*m;
k=isnan(A(:));

% list the nodes which are known, and which will
% be interpolated
nan_list=find(k);
known_list=find(~k);

% how many nans overall
nan_count=length(nan_list);

% convert NaN indices to (r,c) form
% nan_list==find(k) are the unrolled (linear) indices
% (row,column) form
[nr,nc]=ind2sub([n,m],nan_list);

% both forms of index in one array:
% column 1 == unrolled index
% column 2 == row index
% column 3 == column index
nan_list=[nan_list,nr,nc];

% supply default method
if (nargin<2) || isempty(method)
  method = 0;
elseif ~ismember(method,0:5)
  error 'If supplied, method must be one of: {0,1,2,3,4,5}.'
end

% for different methods
switch method
 case 0
  % The same as method == 1, except only work on those
  % elements which are NaN, or at least touch a NaN.
  
  % is it 1-d or 2-d?
  if (m == 1) || (n == 1)
    % really a 1-d case
    work_list = nan_list(:,1);
    work_list = unique([work_list;work_list - 1;work_list + 1]);
    work_list(work_list <= 1) = [];
    work_list(work_list >= nm) = [];
    nw = numel(work_list);
    
    u = (1:nw)';
    fda = sparse(repmat(u,1,3),bsxfun(@plus,work_list,-1:1), ...
      repmat([1 -2 1],nw,1),nw,nm);
  else
    % a 2-d case
    
    % horizontal and vertical neighbors only
    talks_to = [-1 0;0 -1;1 0;0 1];
    neighbors_list=identify_neighbors(n,m,nan_list,talks_to);
    
    % list of all nodes we have identified
    all_list=[nan_list;neighbors_list];
    
    % generate sparse array with second partials on row
    % variable for each element in either list, but only
    % for those nodes which have a row index > 1 or < n
    L = find((all_list(:,2) > 1) & (all_list(:,2) < n));
    nl=length(L);
    if nl>0
      fda=sparse(repmat(all_list(L,1),1,3), ...
        repmat(all_list(L,1),1,3)+repmat([-1 0 1],nl,1), ...
        repmat([1 -2 1],nl,1),nm,nm);
    else
      fda=spalloc(n*m,n*m,size(all_list,1)*5);
    end
    
    % 2nd partials on column index
    L = find((all_list(:,3) > 1) & (all_list(:,3) < m));
    nl=length(L);
    if nl>0
      fda=fda+sparse(repmat(all_list(L,1),1,3), ...
        repmat(all_list(L,1),1,3)+repmat([-n 0 n],nl,1), ...
        repmat([1 -2 1],nl,1),nm,nm);
    end
  end
  
  % eliminate knowns
  rhs=-fda(:,known_list)*A(known_list);
  k=find(any(fda(:,nan_list(:,1)),2));
  
  % and solve...
  B=A;
  B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);
  
 case 1
  % least squares approach with del^2. Build system
  % for every array element as an unknown, and then
  % eliminate those which are knowns.

  % Build sparse matrix approximating del^2 for
  % every element in A.
  
  % is it 1-d or 2-d?
  if (m == 1) || (n == 1)
    % a 1-d case
    u = (1:(nm-2))';
    fda = sparse(repmat(u,1,3),bsxfun(@plus,u,0:2), ...
      repmat([1 -2 1],nm-2,1),nm-2,nm);
  else
    % a 2-d case
    
    % Compute finite difference for second partials
    % on row variable first
    [i,j]=ndgrid(2:(n-1),1:m);
    ind=i(:)+(j(:)-1)*n;
    np=(n-2)*m;
    fda=sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
      repmat([1 -2 1],np,1),n*m,n*m);
    
    % now second partials on column variable
    [i,j]=ndgrid(1:n,2:(m-1));
    ind=i(:)+(j(:)-1)*n;
    np=n*(m-2);
    fda=fda+sparse(repmat(ind,1,3),[ind-n,ind,ind+n], ...
      repmat([1 -2 1],np,1),nm,nm);
  end
  
  % eliminate knowns
  rhs=-fda(:,known_list)*A(known_list);
  k=find(any(fda(:,nan_list),2));
  
  % and solve...
  B=A;
  B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);
  
 case 2
  % Direct solve for del^2 BVP across holes

  % generate sparse array with second partials on row
  % variable for each nan element, only for those nodes
  % which have a row index > 1 or < n
  
  % is it 1-d or 2-d?
  if (m == 1) || (n == 1)
    % really just a 1-d case
    error('Method 2 has problems for vector input. Please use another method.')
    
  else
    % a 2-d case
    L = find((nan_list(:,2) > 1) & (nan_list(:,2) < n));
    nl=length(L);
    if nl>0
      fda=sparse(repmat(nan_list(L,1),1,3), ...
        repmat(nan_list(L,1),1,3)+repmat([-1 0 1],nl,1), ...
        repmat([1 -2 1],nl,1),n*m,n*m);
    else
      fda=spalloc(n*m,n*m,size(nan_list,1)*5);
    end
    
    % 2nd partials on column index
    L = find((nan_list(:,3) > 1) & (nan_list(:,3) < m));
    nl=length(L);
    if nl>0
      fda=fda+sparse(repmat(nan_list(L,1),1,3), ...
        repmat(nan_list(L,1),1,3)+repmat([-n 0 n],nl,1), ...
        repmat([1 -2 1],nl,1),n*m,n*m);
    end
    
    % fix boundary conditions at extreme corners
    % of the array in case there were nans there
    if ismember(1,nan_list(:,1))
      fda(1,[1 2 n+1])=[-2 1 1];
    end
    if ismember(n,nan_list(:,1))
      fda(n,[n, n-1,n+n])=[-2 1 1];
    end
    if ismember(nm-n+1,nan_list(:,1))
      fda(nm-n+1,[nm-n+1,nm-n+2,nm-n])=[-2 1 1];
    end
    if ismember(nm,nan_list(:,1))
      fda(nm,[nm,nm-1,nm-n])=[-2 1 1];
    end
    
    % eliminate knowns
    rhs=-fda(:,known_list)*A(known_list);
    
    % and solve...
    B=A;
    k=nan_list(:,1);
    B(k)=fda(k,k)\rhs(k);
    
  end
  
 case 3
  % The same as method == 0, except uses del^4 as the
  % interpolating operator.
  
  % del^4 template of neighbors
  talks_to = [-2 0;-1 -1;-1 0;-1 1;0 -2;0 -1; ...
      0 1;0 2;1 -1;1 0;1 1;2 0];
  neighbors_list=identify_neighbors(n,m,nan_list,talks_to);
  
  % list of all nodes we have identified
  all_list=[nan_list;neighbors_list];
  
  % generate sparse array with del^4, but only
  % for those nodes which have a row & column index
  % >= 3 or <= n-2
  L = find( (all_list(:,2) >= 3) & ...
            (all_list(:,2) <= (n-2)) & ...
            (all_list(:,3) >= 3) & ...
            (all_list(:,3) <= (m-2)));
  nl=length(L);
  if nl>0
    % do the entire template at once
    fda=sparse(repmat(all_list(L,1),1,13), ...
        repmat(all_list(L,1),1,13) + ...
        repmat([-2*n,-n-1,-n,-n+1,-2,-1,0,1,2,n-1,n,n+1,2*n],nl,1), ...
        repmat([1 2 -8 2 1 -8 20 -8 1 2 -8 2 1],nl,1),nm,nm);
  else
    fda=spalloc(n*m,n*m,size(all_list,1)*5);
  end
  
  % on the boundaries, reduce the order around the edges
  L = find((((all_list(:,2) == 2) | ...
             (all_list(:,2) == (n-1))) & ...
            (all_list(:,3) >= 2) & ...
            (all_list(:,3) <= (m-1))) | ...
           (((all_list(:,3) == 2) | ...
             (all_list(:,3) == (m-1))) & ...
            (all_list(:,2) >= 2) & ...
            (all_list(:,2) <= (n-1))));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(all_list(L,1),1,5), ...
      repmat(all_list(L,1),1,5) + ...
        repmat([-n,-1,0,+1,n],nl,1), ...
      repmat([1 1 -4 1 1],nl,1),nm,nm);
  end
  
  L = find( ((all_list(:,2) == 1) | ...
             (all_list(:,2) == n)) & ...
            (all_list(:,3) >= 2) & ...
            (all_list(:,3) <= (m-1)));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(all_list(L,1),1,3), ...
      repmat(all_list(L,1),1,3) + ...
        repmat([-n,0,n],nl,1), ...
      repmat([1 -2 1],nl,1),nm,nm);
  end
  
  L = find( ((all_list(:,3) == 1) | ...
             (all_list(:,3) == m)) & ...
            (all_list(:,2) >= 2) & ...
            (all_list(:,2) <= (n-1)));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(all_list(L,1),1,3), ...
      repmat(all_list(L,1),1,3) + ...
        repmat([-1,0,1],nl,1), ...
      repmat([1 -2 1],nl,1),nm,nm);
  end
  
  % eliminate knowns
  rhs=-fda(:,known_list)*A(known_list);
  k=find(any(fda(:,nan_list(:,1)),2));
  
  % and solve...
  B=A;
  B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);
  
 case 4
  % Spring analogy
  % interpolating operator.
  
  % list of all springs between a node and a horizontal
  % or vertical neighbor
  hv_list=[-1 -1 0;1 1 0;-n 0 -1;n 0 1];
  hv_springs=[];
  for i=1:4
    hvs=nan_list+repmat(hv_list(i,:),nan_count,1);
    k=(hvs(:,2)>=1) & (hvs(:,2)<=n) & (hvs(:,3)>=1) & (hvs(:,3)<=m);
    hv_springs=[hv_springs;[nan_list(k,1),hvs(k,1)]];
  end

  % delete replicate springs
  hv_springs=unique(sort(hv_springs,2),'rows');
  
  % build sparse matrix of connections, springs
  % connecting diagonal neighbors are weaker than
  % the horizontal and vertical springs
  nhv=size(hv_springs,1);
  springs=sparse(repmat((1:nhv)',1,2),hv_springs, ...
     repmat([1 -1],nhv,1),nhv,nm);
  
  % eliminate knowns
  rhs=-springs(:,known_list)*A(known_list);
  
  % and solve...
  B=A;
  B(nan_list(:,1))=springs(:,nan_list(:,1))\rhs;
  
 case 5
  % Average of 8 nearest neighbors
  
  % generate sparse array to average 8 nearest neighbors
  % for each nan element, be careful around edges
  fda=spalloc(n*m,n*m,size(nan_list,1)*9);
  
  % -1,-1
  L = find((nan_list(:,2) > 1) & (nan_list(:,3) > 1)); 
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([-n-1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end
  
  % 0,-1
  L = find(nan_list(:,3) > 1);
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([-n, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % +1,-1
  L = find((nan_list(:,2) < n) & (nan_list(:,3) > 1));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([-n+1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % -1,0
  L = find(nan_list(:,2) > 1);
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([-1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % +1,0
  L = find(nan_list(:,2) < n);
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % -1,+1
  L = find((nan_list(:,2) > 1) & (nan_list(:,3) < m)); 
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([n-1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end
  
  % 0,+1
  L = find(nan_list(:,3) < m);
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([n, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % +1,+1
  L = find((nan_list(:,2) < n) & (nan_list(:,3) < m));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([n+1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end
  
  % eliminate knowns
  rhs=-fda(:,known_list)*A(known_list);
  
  % and solve...
  B=A;
  k=nan_list(:,1);
  B(k)=fda(k,k)\rhs(k);
  
end

% all done, make sure that B is the same shape as
% A was when we came in.
B=reshape(B,n,m);
end
function neighbors_list = identify_neighbors(n,m,nan_list,talks_to)
% identify_neighbors: identifies all the neighbors of
%   those nodes in nan_list, not including the nans
%   themselves
%
% arguments (input):
%  n,m - scalar - [n,m]=size(A), where A is the
%      array to be interpolated
%  nan_list - array - list of every nan element in A
%      nan_list(i,1) == linear index of i'th nan element
%      nan_list(i,2) == row index of i'th nan element
%      nan_list(i,3) == column index of i'th nan element
%  talks_to - px2 array - defines which nodes communicate
%      with each other, i.e., which nodes are neighbors.
%
%      talks_to(i,1) - defines the offset in the row
%                      dimension of a neighbor
%      talks_to(i,2) - defines the offset in the column
%                      dimension of a neighbor
%      
%      For example, talks_to = [-1 0;0 -1;1 0;0 1]
%      means that each node talks only to its immediate
%      neighbors horizontally and vertically.
% 
% arguments(output):
%  neighbors_list - array - list of all neighbors of
%      all the nodes in nan_list

if ~isempty(nan_list)
  % use the definition of a neighbor in talks_to
  nan_count=size(nan_list,1);
  talk_count=size(talks_to,1);
  
  nn=zeros(nan_count*talk_count,2);
  j=[1,nan_count];
  for i=1:talk_count
    nn(j(1):j(2),:)=nan_list(:,2:3) + ...
        repmat(talks_to(i,:),nan_count,1);
    j=j+nan_count;
  end
  
  % drop those nodes which fall outside the bounds of the
  % original array
  L = (nn(:,1)<1)|(nn(:,1)>n)|(nn(:,2)<1)|(nn(:,2)>m); 
  nn(L,:)=[];
  
  % form the same format 3 column array as nan_list
  neighbors_list=[sub2ind([n,m],nn(:,1),nn(:,2)),nn];
  
  % delete replicates in the neighbors list
  neighbors_list=unique(neighbors_list,'rows');
  
  % and delete those which are also in the list of NaNs.
  neighbors_list=setdiff(neighbors_list,nan_list,'rows');
  
else
  neighbors_list=[];
end
end
function [pointCloudOut]= flatten_XYZ( XYZ )
%flatten_point_cloud: transforms the input pointcloud by a rotation to a flat
%plane and a translation to the zero plane

%   input: point cloud object
%   output: xyz point cloud flattened according to best fit horizontal
%   plane

% fit planar model to point cloud
[n, ~, ~]           = affine_fit(XYZ);

% create rotation matrix
normalToFlat        = [0 0 1];
rotationAxis        = cross(normalToFlat, n);
rotationAngle       = acos(dot(normalToFlat, n) / ...
                      (norm(normalToFlat) * norm(n))); 
                  
u                   = rotationAxis/norm(rotationAxis);                  
I                   = eye(3);
uCross              = [0, -u(3), u(2); u(3), 0, -u(1),; -u(2), u(1), 0];
uTens               = kron(u,u');    
R                   = cos(rotationAngle)*I + sin(rotationAngle)*uCross ...
                      + (1-cos(rotationAngle))*uTens;
                                  

% rotate the point cloud around the vector defined by the cross product of
% the normal verctor and a vertical vertor by the angle defined by their
% dot product

flatXYZ             = XYZ*R;

% set the data in reference to the average heigth
meanHeight          = nanmean(XYZ(:,3));
zpts                = flatXYZ(:,3)-meanHeight;
xpts                = flatXYZ(:,1)-min(flatXYZ(:,1));
ypts                = flatXYZ(:,2)-min(flatXYZ(:,2));
pointCloudOut       = [xpts,ypts,zpts];

end
function [n,V,p]        = affine_fit(X)
    %Computes the plane that fits best (lest square of the normal distance
    %to the plane) a set of sample points.
    %INPUTS:
    %
    %X: a N by 3 matrix where each line is a sample point
    %
    %OUTPUTS:
    %
    %n : a unit (column) vector normal to the plane
    %V : a 3 by 2 matrix. The columns of V form an orthonormal basis of the
    %plane
    %p : a point belonging to the plane
    %
    %NB: this code actually works in any dimension (2,3,4,...)
    %Author: Adrien Leygue
    %Date: August 30 2013
    
    %the mean of the samples belongs to the plane
    p = mean(X,1);
    
    %The samples are reduced:
    R = bsxfun(@minus,X,p);
    %Computation of the principal directions if the samples cloud
    [V,D] = eig(R'*R);
    %Extract the output from the eigenvectors
    n = V(:,1);
    V = V(:,2:end);
end
function [PADDEDGRID]   = pad_grid(GRID)
% add a zeros of width/length of grid around grid

gridSize = size(GRID);
zeroGrid = zeros(gridSize(1),gridSize(2));

PADDEDGRID = [repmat(zeroGrid,1,3)      ;...
    zeroGrid, GRID, zeroGrid  ;...
    repmat(zeroGrid,1,3)      ];

end
function width          = autocorrelation_width(X,Y)
    dx = (X(2)-X(1)); 
    [C,lags]    =xcorr(Y);
    Area1       =trapz(lags*dx,C);
    width       =Area1/max(C);    
end


















