clear;


iM=0;
Mvect=5:0.1:9.2;

for M=Mvect,
iM=iM+1;
    lambda=(6.2e-6)*10.^(M);
lambdavect(iM)=lambda;

k_vect=[0:3*lambda];
for k=k_vect
    if isfinite(prod(1:k))
        p(k+1)=(lambda.^k.*exp(-lambda))./(prod(1:k));
    else
        p(k+1)=(lambda.^k.*exp(-lambda))./(realmax/1e10);
    end
end
%plot(k_vect,p,lambda,0.1,'r*')

CDF=cumsum(p(isfinite(p)));
I2=find((CDF<=0.99),1,'last');
I1=find((CDF>=0.01),1,'first');
DeltaN(iM,:)=[k_vect(I1)/lambda, k_vect(I2)/lambda];
end


plot(Mvect,log10(lambdavect./lambdavect));
hold on
plot(Mvect,log10(DeltaN));
hold off
%sigma=sqrt(lambda);
%I_CI=(abs(k_vect-lambda)<=sigma);
%sum(p(I_CI))


%log10(k_vect(I)/lambda)
