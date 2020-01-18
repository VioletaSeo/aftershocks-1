function b_val_timeseries_gulia(T,M,I0,Npre,Npost,Nmin,plotTrueFalse)

% please for the love of god put the data in matlab datetime
% for ComCat:
% t = datetime(query21.time,'InputFormat','uuuu-MM-dd''T''HH:mm:ss.SSS''Z');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parsing input %%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==2 % bare minimum catalog
    % determine mainshock:
    Imain = M==max(M);

    % determined potential foreshocks as events within 6 months prior to the
    % sequence with magnitude within 2 magnitude units of the mainshock mag.
    Ifore = (M>(max(M)-2)) & T<T(Imain);
end

if nargin <= 3 %  and start time - default values
    Npre = 250;
    Npost= 400;
    Nmin = 50;
    plotTrueFalse = true;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inc                 = 0.01; % mininum data increment in the data
corr                = 0.2;  % buffer past maximum curvature
bufferAfterEvent    = 1; % days after mainshock

% remove events near major earthquakes:
for iT = T([Imain,Ifor])
    T
end

% split the catalog into pre- and post- sequence
Ipre = T<min(T(I0));
[Mpre, Tpre ] = goodind( Ipre,M,T); 
[Mpost,Tpost] = goodind(~Ipre,M,T);
    
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

if plotTrueFalse
    figure;

    XLim = [T0-5*1/365,max(T)];
    
    subplot(2,1,1)
    yyaxis left
    scatter(T,M,(M+1).^2,'filled','MarkerFaceAlpha',0.5)
    ylabel('Magnitude')
    YLim = ylim; 
    xlim(XLim)
    xticks([])
    
    yyaxis right
    plot(Tc,MC,'--')
    ylim(YLim);
    yticks([])
    ylabel('M_c')

    subplot(2,1,2)
    
    yyaxis right
    lN = plot(Tc,A,'color',[0.8 0.8 0.8]);
    ylabel('Number of events')
    set(gca,'yscale','log')
    
    yyaxis left; hold on
    lb = plot(Tc,B/Bback*100,'k','LineWidth',3);
    plot(Tc,(B-Err)/Bback*100,'--k');
    plot(Tc,(B+Err)/Bback*100,'--k');
    lback = plot([T(1),T(end)],[1,1]*100,'-r','LineWidth',2);
    ylabel('relative b-values (%)')
    YLim = [0,150];
    ylim(YLim)
    xlabel('Time')
    xlim(XLim)
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = [0.8 0.8 0.8];
       
    legend([lb,lback,lN], ...
        {'(b\pm1\sigma)/b_{pre}',sprintf('b_{pre} = %0.1f',Bback), sprintf('N_{min} = %i',Nmin)}, ...
        'Location','southwest');
    
    set(findall(gcf,'-property','FontSize'),'FontSize',10);
    set(gcf,...
    'Units',        'Inches', ...
    'Position',     [0,0,7,7],...
    'PaperUnits',   'Inches',...
    'PaperSize',    [7,5]);
    figureName = 'Gulia_Ridgecrest';
    saveas(gcf,figureName);
    print(gcf,[figureName,'.jpg'],'-djpeg','-r300');
end

end

function [a,b,bstd,mc,t] = sliding_window(M,T,Nsub,Nmin,inc,corr)

% initial screen: 
[M,t] = goodind(M>=max_curvature(M, inc, corr), M, T);

Ntot = length(M);

[a,b,bstd,mc] = deal(nan(Ntot,1));

% maximum likelyhood b and one standard deviation range (shi and bolt)
b_val  = @(m,mc)     1/log(10) * 1/mean(m-mc + inc/2);
sigma_b= @(a,b,m)    2.3*b^2 * sum((m-mean(m).^2))/(a*(a-1));

% sliding window:
for n = (Nsub+1):Ntot
    % get group of earthquakes
    Mg      = M(n-Nsub:n);
    mc(n)   = max_curvature(Mg, inc, corr);
    Mg      = Mg(Mg>=mc(n));
    a(n)    = length(Mg);
    b(n)    = b_val(Mg,mc);
    bstd(n) = sigma_b(a(n),b(n),Mg);
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
