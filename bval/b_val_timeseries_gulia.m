function b_val_timeseries_gulia(T,M,I0,Npre,Npost,Nmin,plotTrueFalse)

% please for the love of god put the data in matlab datetime

% for ComCat:
% t = datetime(query21.time,'InputFormat','uuuu-MM-dd''T''HH:mm:ss.SSS''Z');
% if using decimal year, convert to datetime using the following:
% dt = datetime(floor(x), 1, 1) + years(x-floor(x))

[T,I] = sort(T);
M = M(I);

Tin = T; 
Min = M;

% Magic numbers:
inc                 = 0.01; % mininum data increment in the data
corr                = 0.2;  % buffer past maximum curvature
DeltaM              = 2.5;    % what consists of a major event
majorEventMinMag    = max(M)-DeltaM;  % minimum magnitde to be considered to significantly affect catalog completess 
time2crop           = @(m) sqrt(10.^(3/2*m)/10.^(3/2*7)); % (days) currently scaled to be one day for an Mw7


% The only inputs required are T and M, the start of the sequence is
% identified as the first event before the mainshock (the largest event)
% that is 1 magnitude unit smaller than the mainshock.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parsing input %%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==2 % bare minimum catalog
    % determine mainshock:
    Imain = M==max(M);

    % determined potential foreshocks as events within 6 months prior to the
    % sequence with magnitude within 1 magnitude units of the mainshock mag.
    I0 = find(M>(max(M)-DeltaM),1,'first');
end

if nargin <= 3 %  and start time - default values
    Npre = 50;
    Npost= 400;
    Nmin = 10;
    plotTrueFalse = true;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sequenceStart       = T(I0);

% remove earthquakes after significant events:
majorEventId = find(M>=majorEventMinMag);
screen = ones(size(M));
for n = 1:length(majorEventId) 
    iI = majorEventId(n);
    screen = screen & (T<T(iI) | T>T(iI)+time2crop(M(iI)));
end
[M, T] = goodind(screen, M,T); 

% ok now we split the sequences
Ipre = T<sequenceStart;
[Mpre, Tpre ] = goodind( Ipre,M, T); 
[Mpost,Tpost] = goodind(~Ipre,M, T);
    
% run a sliding window for both the subcatalogs. This is where the b-value calculated
[apre, bpre, bstdpre, mcpre, tpre ]    = sliding_window(Mpre, Tpre, Npre, Nmin,inc,corr);
[apost,bpost,bstdpost,mcpost,tpost]    = sliding_window(Mpost,Tpost,Npost,Nmin,inc,corr);

% Concateneate the time series before (pre) and after (post) the sequence
A = [apre;      apost];
B = [bpre;      bpost];
Err=[bstdpre;   bstdpost];
Tc =[tpre;      tpost];
MC =[mcpre;     mcpost];

% Define the background b-value as the median of the background time-series
Bback= nanmedian(bpre); 

% revert to alternate definition of background seismicity insufficient events
if isnan(Bback) 
    error('Not enough events in the catalog to determine a stable value for the background')
end 


%% Plotting the results
if plotTrueFalse
    figure;

    XLim = [sequenceStart-2,sequenceStart+20];
    YLim = [0,max(Min)+0.2];
    
    subplot(2,1,1)
    yyaxis left; hold on
    scatter(Tin,Min, (Min).^2,'r','filled','MarkerFaceAlpha',0.2)
    scatter(T  ,M  , (M).^2,'k','filled','MarkerFaceAlpha',0.8)
    ylabel('Magnitude')
    ylim(YLim); 
    xlim(XLim)
    %xticks([])
    
    yyaxis right
    plot(Tc,MC,'LineWidth',2)
    ylim(YLim);
    yticks([])
    ylabel('M_c')
    
    ax = gca;
    ax.YAxis(1).Color = 'k';

    subplot(2,1,2)
    
    yyaxis right
    lN = plot(Tc,A,'color',[0.7 0.7 0.7]);
    ylabel('Number of events')
    set(gca,'yscale','log')
    
    yyaxis left; hold on
    lb = plot(Tc,B/Bback*100,'k','LineWidth',3);
    plot(Tc,(B-Err)/Bback*100,'--k');
    plot(Tc,(B+Err)/Bback*100,'--k');
    lback = plot([T(1),T(end)],[1,1]*100,'-','LineWidth',2);
    ylabel('relative b-values (%)')
    YLim = [0,150];
    
    I = Min>majorEventMinMag;
    stem(Tin(I),YLim(2)*ones(size(Tin(I))), ...
        'LineStyle','--', ...
        'Marker','none', ...
        'color','r', ...
        'Linewidth',2) 
    
    ylim(YLim)
    xlabel('Time')
    xlim(XLim)
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = [0.6 0.6 0.6];
       
    legend([lb,lback,lN], ...
        {'(b\pm1\sigma)/b_{pre}',sprintf('b_{pre} = %0.1f',Bback), sprintf('N_{min} = %i',Nmin)}, ...
        'Location','southwest');
    
    set(findall(gcf,'-property','FontSize'),'FontSize',10);
    set(gcf,...
    'Units',        'Inches', ...
    'Position',     [0,0,7,7],...
    'PaperUnits',   'Inches',...
    'PaperSize',    [7,5]);
    figureName = 'Gulia_Puerto_Rico';
    saveas(gcf,figureName);
    print(gcf,[figureName,'.jpg'],'-djpeg','-r300');
end

end

function [a,b,bstd,mc,t] = sliding_window(M,T,Nsub,Nmin,inc,corr)

% initial screen: 
[M,t] = goodind(M>=max_curvature(M, 0.1, corr), M, T);

Ntot = length(M);

[a,b,bstd,mc] = deal(nan(Ntot,1));

% maximum likelyhood b and one standard deviation range (shi and bolt)
b_val  = @(m,mc)     1/log(10) * 1/mean(m-mc + inc/2);
sigma_b= @(a,b,m)    2.3*b^2 * sum((m-mean(m).^2))/(a*(a-1));


% sliding window:[B,I] = sort(___)
for n = (Nsub+1):Ntot
    % get group of earthquakes
    Mg      = M(n-Nsub:n);
    mc(n)   = max_curvature(Mg, inc, corr);
    Mg      = Mg(Mg>=mc(n));
    a(n)    = length(Mg);
    % im sorry im sorry im sorry im sorry im sorry im sorry im sorry im sorry
    use_Gulia_estimation = ~true;
    if use_Gulia_estimation
        dummyCat        = nan(length(Mg),6);
        dummyCat(:,6)   = Mg;
        [~,b(n),flag,bstd(n)] = NLIndex_version1_sigma(dummyCat,min(Mg),'PreDefinedMc');
        if ~(flag == 2 || flag == 3); [b(n),bstd(n)] = deal(nan); end
    else

        b(n)    = b_val(Mg,mc(n));
        bstd(n) = sigma_b(a(n),b(n),Mg);
    end
    
end

% screen for Nmin
b(a<Nmin) = nan;

end


function Mc = max_curvature(M, inc, corr) 

[hc,E] = histcounts(M,min(M):inc:max(M));
I = find(hc==max(hc));
Mc = E(I(1))+corr;

end

function [varargout]    = goodind(goodInd,varargin)

for n=1:length(varargin)
    varargout{n} = varargin{n}(goodInd);  %#ok<AGROW>
end

end

