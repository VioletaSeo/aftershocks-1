function compare_nearby(X1,X2,y1,y2,R)

Idx = rangesearch(X2,X1,R);
get_med = @(y,I) median(y(I));
yq = cellfun(@(IDX) get_med(y2,IDX), Idx);

scatter(y1,yq)
    
end