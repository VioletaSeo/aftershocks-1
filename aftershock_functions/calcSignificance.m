function [p] = calcSignificance(data,model)

n = length(data);
xi = data(:,1);
yi = data(:,2);

a = model(1);
b = model(2);

yiHat = a*xi+b;

xBar = mean(xi);

SE = sqrt(sum((yi-yiHat).^2)/(n-2))/sqrt(sum((xi-xBar).^2));

t = a/SE;

p = 1-tcdf(t,n-2);

end
