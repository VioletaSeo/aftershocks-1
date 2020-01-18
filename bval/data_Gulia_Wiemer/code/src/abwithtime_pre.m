function [pr_median,magco_median, bv_median, av_ann_median, timev, bprevalue,Mcgeneric,timewindow,result_flag_median,bpreplus,bpreminus]=abwithtime_pre(catSpre,N,mbin,Mms,corr)

% This function computes the time series of a and b-values and related uncertainties in the
% pre-mainshco time period. 
% Authors: Laura Gulia and Thessa Tormann, 2016-2018
% Modified by Laura Gulia, 2018-2019


Mtarg=Mms;
newt2=catSpre;
Nmin=50; % Minum number of events above completeness, otherwise no value is computed
step=1; % Move forward in overlapping windows, one event at a time to have the maximum resolution
Mrange=min(newt2(:,6)):mbin:max(newt2(:,6));

timev=newt2(:,3); %  event-by-event
disp('pre: event-by-event')

no_estimates_pre=length(timev);


%% 
% prepare output variables (nan-vectors with length of the timev)
Mcgeneric=ones(no_estimates_pre,1)*nan;
magco_median=ones(no_estimates_pre,1)*nan;

bv_median=ones(no_estimates_pre,1)*nan;
bv_minussig=ones(no_estimates_pre,1)*nan;
bv_plussig=ones(no_estimates_pre,1)*nan;

av_ann_median=ones(no_estimates_pre,1)*nan;
pr_median=ones(no_estimates_pre,1)*nan;

result_flag_median=ones(no_estimates_pre,1)*nan;

timewindow=ones(no_estimates_pre,1)*nan;

sig1_median=ones(no_estimates_pre,1)*nan;
bpreplus=ones(no_estimates_pre,1)*nan;
bpreminus=ones(no_estimates_pre,1)*nan;
%% Now compute the values at each time step, 

% for each time step
for j = 1:step:no_estimates_pre
    if rem(j,100)==0
        disp(['processing pre-time step ',num2str(j),' of ',num2str(length(timev))])
    end
    ind1=find(newt2(:,3)<=timev(j),1,'last');
    if ind1<Nmin
        b_orig=[];
    elseif and(ind1>=Nmin,ind1<N)
        b_orig=newt2(1:ind1,:);
    else  
        % select last ni events
        b_orig=newt2(ind1-N+1:ind1,:);
    end
    
    if ~isempty(b_orig)
        if length(b_orig(:,1))>=Nmin
            
            bv=nan; av_ann=nan; pr=nan; result_flag=nan;
            b=b_orig;
            F=1;          %*% % change between annual (F=1) versus daily (F=365) rates
            T=F*(timev(j)-min(b(:,3)));
            
            % estimate for this sampke max curv Mc and cut catalog
            cntnumb=hist(b(:,6),Mrange);
            [~,ind]=max(cntnumb);
            magco=Mrange(ind)+corr;
            ll=b(:,6)>=magco;
            b=b(ll,:);
            N3=sum(ll);
            
            % if enough events, estimate parameters and uncertainties, also check for
            % linearity 
            if N3>=Nmin
                [bestmc,bv,result_flag,sig1]=NLIndex_version1_sigma(b,magco,'PreDefinedMc');
                 if or(result_flag==3,result_flag==2)  
                    
                    magco=bestmc;
                    ll=b(:,6)>=magco;
                    N_b=sum(ll);
                    
                    % annualise the a-value and compute probabilities
                    av_ann=log10(N_b/T)+ bv*magco;
                    tr=1/10.^(av_ann-bv*Mtarg);
                    pr=1-exp(-1/tr);
                 else % no good answer found, set nan
                    bv=nan;
                    magco=bestmc;
                    av_ann=nan;
                    pr=nan;
                    sig1=nan;
                end
                
                magco_median(j)=magco;
                bv_median(j)=bv;
                sig1_median(j)=sig1;
                pr_median(j)=pr;
                av_ann_median(j)=av_ann;
                result_flag_median(j)=result_flag;
            end
        end
    end
end


% flag vector indicating for each time whether b is value or nan
bprevalue=isfinite(bv_median);
bpreplus=[bv_median+sig1_median];
bpreminus=[bv_median-sig1_median];
