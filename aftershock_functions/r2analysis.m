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
BETA = (X'*X)\X'*Y;
R2         = 1  -  sum((X*BETA - Y).^2)/sum(Y.^2);
R2adjusted = 1  -  (SX(1)-1)/(SX(1)-SX(2)) * sum((X*BETA - Y).^2)/sum(Y.^2);
% remove outlier (rinse and repeat)
var = (X*BETA - Y).^2;
outlierInd = var > 4*mean(var);
X = X(~outlierInd,:);
Y = Y(~outlierInd);
SX= size(X);
X = X - repmat(mean(X,1),SX(1),1); % remove average value
Y = Y - mean(Y);
BETA = (X'*X)\X'*Y;
R2         = 1  -  sum((X*BETA - Y).^2)/sum(Y.^2);
R2adjusted = 1  -  (SX(1)-1)/(SX(1)-SX(2)) * sum((X*BETA - Y).^2)/sum(Y.^2);

%% Random whole model
 for m = 1:mit
        Yshfl  = Y(randperm(length(Y)));
        BETArand= (X'*X)\X'*Yshfl;
        Rrand(m)= 1  -  sum((X*BETArand - Yshfl).^2)/sum(Y.^2);
 end
 

 %% individual models
for n = 1:SX(2)
    Xcol    = X(:,n);
    BETA    = (Xcol'*Xcol)\Xcol'*Y;
    R(n)    = 1  -  sum((Xcol*BETA - Y).^2)/sum(Y.^2);
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

figure 
h = histogram(maxRshuffle,'FaceColor', [0.8 0.8 0.8],'EdgeColor',[1 1 1]);
hold on


xlabel('Variance Reduction')
ylabel('N')
set(gca,'Xlim',[0,0.3])
ah = gca;
yLim = ah.YLim;
recX = linspace(0,max(maxRshuffle),20);
recX = recX(2:end);
for iDim = recX
    rectangle('Position',[0,0,iDim,yLim(2)],'FaceColor',[1 0 0 0.005],'EdgeColor','none')
end
stem(R,3*sqrt(mit)*ones(length(R), 1), ...
    'Color', [0.4660    0.6740    0.1880], ...
    'LineWidth',2)
texth = text(R, 3.3*sqrt(mit)*ones(length(R),1), SourceParameters);
set(texth, 'Rotation', 45)
set(gca,'YLim',yLim)
set(findall(gcf,'-property','FontSize'),'FontSize',12)

end





































