function [MERGEDCAT,IDX] = merge_eq_catalog(catIn, pointStr, varargin)
% merge earthquake clogs. 

%e.g.
% pointStr = {'t','lat','lon','dep','mag'};

check_input(varargin);
%% get t-m-x-y-z pts matrix for the input catalogeq

p0 = get_point(catIn, pointStr);


if isa(varargin{1},'double')
    th = varargin{1};
    varargin = varargin(2:end);
else
    th = max(p0);
end
%% get str format
numIn = length(varargin)/2;
pointFormat = cell(1,numIn);
for n = 1:numIn
    pointFormat{n} = varargin{2*n};   
end

%% order catalogs

catalogs        = cell(1,numIn);
orderedCat      = cell(1,numIn);
orderedFormat   = cell(1,numIn);

for n = 1:numIn
    catalogs{n} = varargin{2*n-1};
end

tableHeightArray = zeros(1, numIn);
for n = 1:numIn
   tableHeightArray(n) = height(catalogs{n}); 
end

[~,I] = sort(tableHeightArray,'descend');

for n = 1:numIn
    orderedCat{n}       = catalogs{I(n)};
    orderedFormat{n}    = pointFormat{I(n)};
end

%% iteratively append catalogs 

% loop over catalogs to append
newCat = catIn;
IDX    = ones(1,height(catIn));


%% 
% p0 = (p0-repmat(min(p0,[],1),height(catIn),1));
% p0 = p0./repmat(max(p0,[],1),height(catIn),1);
%%

for n = 1:numIn
    
    catN    = orderedCat{n};
    catFormat = orderedFormat{n};
    pN      = get_point(catN,catFormat);
    sizePn  = size(pN);
    nPoints = sizePn(1);
    
    % standerdize coordinates    
%     pN = (pN-repmat(min(pN,[],1),nPoints,1));
%     pN = pN./repmat(max(pN,[],1),nPoints,1);
    
    % standardized eucledian search of nearest neighbour
    [Idx,eDistance] = knnsearch(p0,pN,'Distance','seuclidean',...
        'Scale',[std(p0(:,1:2)),max(p0(:,3))]);
    
    % deal with exceptions: - this gets ugly
        % find indexes that are assigned twice
        exception = zeros(1,nPoints);
        
        for iInd = 1:nPoints
            duplicates = find(Idx == Idx(iInd));
            if any(duplicates)
                badInd  = eDistance(duplicates) ~= min(eDistance(duplicates)); % only choose the closest point
                exception(duplicates(badInd)) = true; % flag data that cannot be sorted appropriately
            end
        end
        
        % deal with clear outliers
        outliers = eDistance > 3*geomean(eDistance(~exception));
        outliers = outliers | any(abs(p0(Idx,:)-pN) > repmat(th,nPoints,1),2);
       
        % discard the earthquakes that don't fit
        goodEq = ~(exception' | outliers);
        Idx = Idx(goodEq);
    
    fieldName = catN.Properties.VariableNames;
    for k = 1:width(catN)
        newFieldName            = sprintf('%s_appended_cat%g',fieldName{k},I(n)); % avoid redundance
        dataColumn              = catN.(fieldName{k});
        if isa(dataColumn,'double')
            newCat.(newFieldName)       = nan(height(newCat),1);
            newCat.(newFieldName)(Idx)  = dataColumn(goodEq);
        end
    end
    
    tempInd = zeros(1,length(IDX));
    tempInd(Idx) = true;
    IDX     = IDX & tempInd;    

end

% output:
MERGEDCAT  = newCat;
    
end

function check_input(input)

if length(input) < 1
    error('There must be at least two inputs')
end 

end

function P = get_point(cat,pointStr)

P = zeros(height(cat),length(pointStr));
    for n = 1:length(pointStr)
        P(:,n) = cat.(pointStr{n});
    end
end
