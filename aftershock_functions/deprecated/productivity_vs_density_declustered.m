function productivity_vs_density_declustered()
% correlate productivity to local eq density without double counting
% earthquakes

subsampledDataFN = 'data_subsampled.mat';

savedData = load(subsampledDataFN);
    t       = savedData.t;
    lat     = savedData.lat;
    lon     = savedData.lon;
    depth   = savedData.depth;
    M       = savedData.M;
    fms     = savedData.fms;
    numEq   = length(t);
    
    
% input: 
numBlock    = 50;    % segmentation of the data into individual block (-or windows) of density ranges 
numberOfBoxes = numBlock; % i dont want to talk about it
numMag      = 20;  % number of magnitudes for the the plot of number of earthquakes vs mag
numDepthBin = 20;  % number of bins used to split up earthqukes

% editing
windowSize = 0.3; % proportion of the range in values

% feed all eqs into get_prod_parameters - will need to add the indices of
% the main shocks to the function

[~,~,~,~,MSInd, ASCount] = get_prod_parameters(t,lat,lon,depth,M);
[MSInd,I]=sort(MSInd,'descend');
ASCount = ASCount(I);


% subsample data
[~,~,~,depthMS,MMS] = goodind(MSInd,t,lat,lon,depth,M);

% make histogram of mainshocks 
[N,edges] = histcounts(depthMS,numDepthBin);
Ntot      = length(depthMS);
localDepthDensity = interp1((edges(1:(end-1))+edges(2:end))/2,N,depthMS)/Ntot;

% initializing things for first for loop (iterating of depth density
minmaxDD   = minmax(localDepthDensity);
absoluteWindowSize = windowSize*diff(minmaxDD);
startDD    = minmaxDD(1) + absoluteWindowSize;
endDD      = minmaxDD(2) - absoluteWindowSize;
DDArray    = linspace(startDD, endDD, numBlock);

prefactor  = zeros(1,numBlock);

% initialize variable regarding the magnitude ranges
magnitudeBracketArray = linspace(min(MMS),max(MMS),numMag); % center the magnitudes
magnitudeArray = (magnitudeBracketArray(1:end-1) + magnitudeBracketArray(2:end))/2;

figure

for n = 1:numBlock
    
    % select eqs that are in depth density range
    ind = localDepthDensity >= DDArray(n)-absoluteWindowSize & ...
          localDepthDensity <  DDArray(n)+absoluteWindowSize; 
      
    selectedMMS     = MMS(ind);
    selectedASCount = ASCount(ind);
    
    % intialize for second for loop
    magnitudeCountArray = zeros(1,numMag-1);
    for iMagInc = 1:numMag-1  
        ind         = selectedMMS >= magnitudeBracketArray(iMagInc) & ...
                      selectedMMS <  magnitudeBracketArray(iMagInc+1);
        numMS                           = sum(ind);
        magnitudeCountArray(iMagInc)    = sum(selectedASCount(ind))/numMS;
    end
    
    % plot an aggregate plot of the counts
    hold on
    parBlock = n;
    
    scatter(magnitudeArray,magnitudeCountArray,[], ...
        [parBlock/numberOfBoxes, 0, 1-parBlock/numberOfBoxes], ...
        'filled')
    
    % get productivity parameters
    [prefactor(n),~] = getproductivity(magnitudeArray,magnitudeCountArray,'alpha1');
    
end

set(gca,'XScale', 'log', ...
    'YScale', 'log')
xlabel('Magnitude')
ylabel('Number Of Aftershocks') % anotating the aftershock prod plot

% plot the change in in k as a function of depth
figure
plot(DDArray,prefactor)
xlabel('Seismic index')
ylabel('k''')
title('Productivity as a function seismic density')


end


function [varargout]    = goodind(goodInd,varargin)

for n=1:length(varargin)
    varargout{n} = varargin{n}(goodInd);
end

end