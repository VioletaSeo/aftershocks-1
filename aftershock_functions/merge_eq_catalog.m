function [MERGEDCAT,IDX] = merge_eq_catalog(catIn, pointStr, varargin)
% merge earthquake catalogs. 

%e.g.
% pointStr = {'t','lat','lon','dep','mag'};


check_input(varargin);

numIn = length(varargin)/2;

%% get str format
pointFormat = cell(1,numIn);
for n = 1:numIn
    pointFormat{n} = varargin{2*n};   
end

%% get t-m-x-y-z pts matrix for the input catalog

p0 = get_point(catIn, pointStr);

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

for n = 1:numIn
    
    catN    = orderedCat{n};
    catFormat = orderedFormat{n};
    pN      = get_point(catN,catFormat);
    sizePn  = size(pN);
    nPoints = sizePn(1);
    
    % standardized eucledian search of nearest neighbour
    [Idx,eDistance] = knnsearch(p0,pN,'Distance','seuclidean'); 
    
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
        
        % show what earthquakes couldnt be matched
        if any(exception)
            disp(sprintf('The following %g earthquakes could not be merged:', ...
                sum(exception)));
            for iInd = 1:nPoints
                if exception(iInd)
                    disp(sprintf('Time: %g       Lat: %g     Lon: %g     Depth: %g      Mag: %g', ...
                        catN.(catFormat{1})(iInd), ...
                        catN.(catFormat{2})(iInd), ...
                        catN.(catFormat{3})(iInd), ...
                        catN.(catFormat{4})(iInd), ...
                        catN.(catFormat{5})(iInd)))
                      % part of the debuggin
                      %scatterm(catN.(catFormat{2})(iInd), ...
                      %  catN.(catFormat{3})(iInd),'r')
                end
            end           
        end
        % discard the earthquakes that don't fit
        
        Idx = Idx(~exception);
        
        % look at whether the bad earthquakes can be attributed to another
        % earthquake
    
    % append data (notation indicates the originalname_appended_cat#)
    
    fieldName = catN.Properties.VariableNames;
    for k = 1:width(catN) 
            newFieldName            = sprintf('%s_appended_cat%g',fieldName{k},I(n)); % avoid redundance
            dataColumn              = catN.(fieldName{k});
            newCat.(newFieldName)       = nan(height(newCat),1);
            newCat.(newFieldName)(Idx)  = dataColumn(~exception);
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
