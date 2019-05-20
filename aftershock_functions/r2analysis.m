function [Rvalues] = r2analysis(X,Y,SourceParameters,mit)    

% store variables int
badI = any([isnan(X),isinf(X),isnan(Y),isinf(Y)],2);
X = X(~badI,:);
Y = Y(~badI);
SX= size(X);
X = X-mean(X,1);
X = X./std(X);
Y = Y - mean(Y);
Y = Y./std(Y); % ?

%% whole modelRrand
% BETA = (X'*X)\X'*Y;
% R2         = 1  -  sum((X*BETA - Y).^2)/sum(Y.^2);
% R2adjusted = 1  -  (SX(1)-1)/(SX(1)-SX(2)) * sum((X*BETA - Y).^2)/sum(Y.^2);
% % remove outlier (rinse and repeat)
% var = (X*BETA - Y).^2;
% outlierInd = var > 4*mean(var);
% X = X(~outlierInd,:);
% Y = Y(~outlierInd);
% SX= size(X);
% X = X - repmat(mean(X,1),SX(1),1); % remove average value
% Y = Y - mean(Y);
% BETA = (X'*X)\X'*Y;
% R2         = 1  -  sum((X*BETA - Y).^2)/sum(Y.^2);
% R2adjusted = 1  -  (SX(1)-1)/(SX(1)-SX(2)) * sum((X*BETA - Y).^2)/sum(Y.^2);
% 
% %% Random whole model
%  for m = 1:mit
%         Yshfl  = Y(randperm(length(Y)));
%         BETArand= (X'*X)\X'*Yshfl;
%         Rrand(m)= 1  -  sum((X*BETArand - Yshfl).^2)/sum(Y.^2);
%  end
%  

 %% individual models
 [corrSign,R] = deal(zeros(SX(2),1));
for n = 1:SX(2)
    Xcol    = X(:,n);
    BETA    = (Xcol'*Xcol)\Xcol'*Y;
    R(n)    = 1  -  sum((Xcol*BETA - Y).^2)/sum(Y.^2);
    corrSign(n) = BETA > 0;
end
Rvalues = R; %table(SourceParameters',R2');

%% random shuffled data
maxRshuffle = zeros(mit,1);
for m = 1:mit
    Rshuffle = zeros(mit,1);
    for n = 1:SX(2)
        Xcol    = X(:,n);
        Yshuffle= Y(randperm(length(Y)));
        BETAshuffle= (Xcol'*Xcol)\Xcol'*Yshuffle;
        Rshuffle(n)= 1  -  sum((Xcol*BETAshuffle - Yshuffle).^2)/sum(Yshuffle.^2);
    end
    maxRshuffle(m) = max(Rshuffle);
end

%% plots

h = histogram(maxRshuffle,'FaceColor', [0.8 0.8 0.8],'EdgeColor',[1 1 1]);
hold on

xlabel('Variance Reduction')
set(gca,'Xlim',[0,0.2])
ah = gca;
yLim = ah.YLim;
recX = linspace(0,max(maxRshuffle),20);
recX = recX(2:end);

% plot(repmat(prctile(maxRshuffle,95),1,2),yLim,'LineWidth',3,'Color',[1 0 0 0.2])
% text(prctile(maxRshuffle,96),0.9*yLim(2), 'p = 0.05','Color','r')

% for iDim = recX
%     rectangle('Position',[0,0,iDim,yLim(2)],'FaceColor',[1 0 0 0.005],'EdgeColor','none')
% end


corrSign = logical(corrSign);
SPh = stem(R(corrSign),2.5*sqrt(mit)*ones(length(R(corrSign)), 1), ...
    'Color', [1, 0, 0, 0.1], ...
    'LineWidth',2);

SNh = stem(R(~corrSign),2.5*sqrt(mit)*ones(length(R(~corrSign)), 1), ...
    'Color', [0, 0, 1, 0.1], ...
    'LineWidth',2);

legend([h,SPh,SNh],{'Max by permution','Positive corr.','Negative corr.'},'Location','northeast')

texth = text(R, 3*sqrt(mit)*ones(length(R),1), SourceParameters);
set(texth, 'Rotation', 45)
set(gca,'YLim',yLim,'Ytick',[])
set(findall(gcf,'-property','FontSize'),'FontSize',12)
ax1 = gca;
ax2 = axes('Position',ax1.Position,...
    'Ytick',[],...
    'XAxisLocation','top',...
    'Color','none', ...
    'Xlim',[0,0.2], ...
    'Xtick',0:0.02:2);
xlabel('Family-wise p-value')
xt = get(ax2,'Xtick');
for n = 1:length(xt)
   xtl{n} = sprintf('%0.1g',(sum(maxRshuffle > xt(n))/length(maxRshuffle)));
end
set(ax2,'Xtick',xt,'XtickLabel',xtl)

end





































