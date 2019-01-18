% quick script to determine how anomalous an earthquake is based 

%% inputs
eq.M        = 6.4;  % magnitude of the mainshock of interest
eq.timeWin  = 30;   % days
eq.spaceWin = 50;   % km
researchWin = [-180, 180; ...
               -90,  90 ];  % size of the catalog consisdered for comparison
depthRange  = [0,55];  % depth cutoff
catFN       = '/home/kdascher/Downloads/query (9).csv';

%% analysis 

% scale time and space windows
% time
scaledWindow    = eq.timeWin/10^(eq.M-9);

% space
stressDrop      = 10^5*10;
mo              = (10^((eq.M+10.73)*1.5))*10^-7;
sourceDim       = 2 * 0.001 * (mo*(7/16) / stressDrop) ^ (1/3); % in km tbl 9.1 modern global seismo lay
scaledRadius    = eq.spaceWin./sourceDim; 

% get reference
figure
[k,a,~,~,mth] = aftershock_productivity_kernel(...
    'DepthRange',        [depthRange], ...
    'LatRange',          researchWin(2,:), ...
    'LonRange',          researchWin(1,:), ...
    'TimeWindow',        scaledWindow*1.5, ...
    'TimeSelectionWindow',scaledWindow, ...
    'SelectionRadius',   scaledRadius);

load mainshock_info.mat
hold on
h = scatter(MSmag,MSprod,[],[0.3 0.3 0.3],'filled','MarkerFaceAlpha',0.4);
uistack(h,'down')

% count number of aftershocks for the eq of interest
data        = readtable(catFN);
timeS       = datetime(data.time, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z');
time        = datenum(timeS);

tol         = 0.001; % avoid roundoff error issues
mainshockInd= find(abs(data.mag - eq.M)<tol,1);

IAS         = sqrt(distance(data.latitude(mainshockInd),data.longitude(mainshockInd),data.latitude,data.longitude).^2 ...
    + (data.depth(mainshockInd)-data.depth).^2) < eq.spaceWin...
    & ...
    (time - time(mainshockInd)) < eq.timeWin ...
    & ...
    (time - time(mainshockInd)) > 0 ...
    & ...
    data.mag >= mth;
numAS       = sum(IAS);

hold on
scatter(eq.M,numAS, 200,'r','filled', 'MarkerFaceAlpha', 0.7)
legend({'N_{AS}','N_{AS}(M_i)','Productivity Law','Query'})

figure
hold on
stem(time,data.mag,'filled')
hold on
stem(time(mainshockInd),data.mag(mainshockInd),'r','filled')
plot(minmax(time'),[mth,mth],'r','LineWidth',2)
stem(time(IAS),data.mag(IAS),'r')

xlabel('Time')
ylabel('Magnitude')
legend({'','MainShock','M_{th}'},'Location','southeast')
%set(gca,'YLim',[min(data.mag)-1,max(data.mag)+0.5])








