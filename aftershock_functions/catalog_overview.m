function catalog_overview(t,lat,lon,depth,M)
% make an overview of the catalog 
figure
figureDim   = [2,3];
figCount    = 1;
numEq = length(t); 

magInc          = min(M):0.1:max(M);

%% 1. histogram
subplot(figureDim(1),figureDim(2),figCount); figCount = figCount+1;

histogram(M,magInc)
title(['Histogram, n= ',num2str(numEq)])
ylabel('Number of Earthquakes') 
xlabel('Magnitude')

set(gca,'YScale','log')

%% 2. GRM
subplot(figureDim(1),figureDim(2),figCount); figCount = figCount+1;

numDivisions    = numel(magInc);
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