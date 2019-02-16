clear;
%AS prod & attributes
% test of implementation by predicting Focal Mechanism from Rake
% example: mdl=fitcnb(xt,Y,'CategoricalPredictors',[1],'Distribution','mvmn')
load everything_merged.mat;

%FitCat=mergedCat(:,[10 12 13 16 17  28 39:40]);
FitCat=mergedCat(:,[ 13  16]); 
% dip, longitude and erbb/m0 and stress drop are all helpful
%[h,e]=hist(mergedCat.MSres_appended_cat1,4);
e=[-1.2 -0.2 0.2 1.2];
FitCat.D=discretize(mergedCat.MSres_appended_cat1,e,'categorical',{'low','med','high'});
%FitCat.ErM0=mergedCat.EBB_appended_cat1./mergedCat.Mo;

%mdl=fitcnb(mergedCat,'MSres_appended_cat1','CategoricalPredictors',[7],'Distribution','mvmn')
% need to discretize residuals
% I=isnan(FitCat.EBB_appended_cat1);
% FitCat=FitCat(~I,:);
%

Ntrain=50;
Itrain=[1:Ntrain];
Itest=[Ntrain+1:height(FitCat)];
FitCatTrain=FitCat(Itrain,:);
FitCatTest=FitCat(Itest,:);
% FitCatTrain=FitCat;
% FitCatTest=FitCat;
%mdl=fitcnb(FitCatTrain,'D');
mdl=fitcknn(FitCatTrain,'D');


% FitCat=table(mergedCat.FocalMechanism, mergedCat.Rake,'VariableNames',{'FocalMechanism', 'Rake'});
% I=(eq(FitCat.FocalMechanism,{'Oblique'}));
% FitCat=FitCat(~I,:);
% 
% %FitCat=table(FitCat,);
% mdl=fitcnb(FitCat,'FocalMechanism');%,'Distribution','mvmn');

P=predict(mdl,FitCat);
Correct=eq(FitCat.D,P);
FracCorrect=sum(Correct)/length(Correct)
