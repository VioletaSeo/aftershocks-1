% ------------------- NLIndex NonLinearityIndex ---------------------------
% -------------------------------------------------------------------------

% Measure to detect whether the extrapolation of a given FMD to high
% magnitudes is applicable or whether the shape of the FMD suggests that
% an extrapolation would over- or underestimate the probably true rates for
% large events
%
% author: Thessa Tormann
% date: October 2012
%
%**************************************************************************
%
% Input parameters:
%
% 1) cat = catalog to test
% 2) Mcmin = suggested magnitude of completeness, depending on the mode, this
%    Mc estimate can be corrected to higher values in order to optimize the
%    linear fit
% 3) Nmin = minimum number of events required to attempt a b-value estimate
%    commonly Nmin=50
% 4) Mtarg = target magnitude of large events, commonly Mtarg = 6
% 5) mode = sets the mode in which to run the filter:
%       'PreDefinedMc' = uses the given Mcmin as magnitude of completeness
%       'OptimizeMc' = checks whether a higher Mc could optimize the
%                      NLIndex
% 6) binnumb = minimum number of different Mc for which N>=Nmin, i.e. b can
%        be estimated
% 7) sigMcDif = difference between Mmin and Mc that is regarded
%        significant, e.g. sigMcDif=1
% 8) slope_pos = minimum slope of b-value trend that indicates an
%        overestimation of large magnitude rates, e.g. slope_pos = 0.05
% 9) slope_neg = maximum slope of b-value trend that indicates an
%        underestimation of large magnitude rates, e.g. slope_neg = -0.05
%
%**************************************************************************
%
% Algorithm:
%
% 1) calculate b-value for all Mcut from Mcmin to the highest Mcut for
%       which still Nmin events are sampled
% 2) mode: 'PreDefinedMc'
%       --> calculate mean and standard deviation of b-values estimated in
%        step 1), do this if at least 5 estimates were possible
%       --> divide standard deviation by the largest individual b-value
%        uncertainty, this value is the NLIndex
%       --> if NLIndex is <=1, FMD is linear
%       --> if NLIndex is >1, and slope of b(mcut) is clearly positive,
%           FMD overestimates large M rates
%       --> if NLIndex is >1, and slope of b(mcut) is clearly negative,
%           FMD underestimates large M rates
% 3) mode: 'OptimizeMc'
%       --> calculate mean and standard deviation of b-values estimated in
%        step 1), do this if at least 5 estimates were possible
%       --> divide standard deviation by the largest individual b-value
%        uncertainty, do this for each possible mcut, and for each mcut
%        this value is the NLIndex
%       --> divide the NLIndex for each mcut by the number of estimated
%        b-values to weigh the result by data density
%       --> find the minimum weighted NLIndex to find the bestmc that
%        produces the most linear FMD fit
%       --> find the bestb by taking the mean of all b-value of those mcut
%        that produced linear estimates
%       --> if NLIndex is <=1, FMD is linear
%       --> if NLIndex is >1, and slope of b(mcut) is clearly positive,
%           FMD overestimates large M rates
%       --> if NLIndex is >1, and slope of b(mcut) is clearly negative,
%           FMD underestimates large M rates
%
%**************************************************************************
%
% Output parameters:
%
% 1) bestmc = best value for Mc (=Mcmin, if mode='PreDefinedMc')
% 2) bestb = best b-value estimate, i.e.:
%       b for all M>=Mcmin if mode='PreDefinedMc'
%       median of all b values calculated for M>=Mc for acceptable Mc>=Mcmin
% 3) result_flag = flag describing the outcome of the analysis:
%           1: catalog with N<Nmin events, no b-value estimate
%           2: catalog with N>=Nmin events, but not enough to compute NLI
%           3: FMD is linear with Mcmin<=Mc<=Mcmin+sigMcDif
%           4: FMD is linear for significantly increased Mc>Mcmin+0.5
%           5: FMD is unstable
%           6: FMD underestimates Mtarg rates
%           7: FMD overestimates Mtarg rates
%
%**************************************************************************

function [bestmc,bestb,result_flag,sig1]=NLIndex_version1_sigma(cat,Mcmin,mode)

clear marker mark uomark uomark2 b db
bestmc=[]; bestb=[];

% -------------------------------------------------------------------------
% SET PARAMETERS:
% -------------------------------------------------------------------------

% 1) CHOOSE MODE (comment one)

%mode='PreDefinedMc';
%mode='OptimizeMc';

% 2) SET NUMBERS

Nmin=50;
Mtarg=6;
binnumb=5;
sigMcDif=0.5;
slope_pos=0.05;
slope_neg=-0.15;

%--------------------------------------------------------------------------
% Catalog preparation
%--------------------------------------------------------------------------



if length(cat(:,1))<Nmin
    
    result_flag=1;
    bestmc=nan;
    bestb=nan;
%     disp('not enough events to calculate b')
    
end

% determine magnitude range in catalog from Mcmin to Mmax
Mrange=Mcmin:0.1:max(cat(:,6));

% calculate FMD (cumulative number of events for M>=Mcmin)
Numb=hist(cat(:,6),Mrange);
Numbh=Numb(end:-1:1);
Ncumh=cumsum(Numbh);
Ncum=Ncumh(end:-1:1);

% disp(['Catalog has ', num2str(Ncum(1)),' events'])

% determine magnitude range for which Ncum(M) >= Nmin (i.e. max mag for
% which >=50 events exist)
xx=Ncum>=Nmin;
Mcrange=Mrange(xx);

if length(Mcrange)<binnumb
    
    result_flag=2;
bestmc=Mcmin;
bestb=(1/(mean(cat(:,6))-(bestmc-0.05)))*log10(exp(1));
        
        sig1 = (sum((cat(:,6)-mean(cat(:,6))).^2))/(length(cat(:,6))*(length(cat(:,6))-1));
        sig1 = sqrt(sig1);
        sig1 = 2.30*sig1*bestb^2;
    
else
    
    % ---------------------------------------------------------------------
    % Calculation of
    %     b - maximum likelihood b-value (Aki 1965), and
    %     std(b) - standard deviation (Shi & Bold 1982)
    % for each cutoff-magnitude of Mcrange
    % (save results in matrix 'b')
    % ---------------------------------------------------------------------
    
    b=ones(length(Mcrange),3)*nan;
    for i=1:length(Mcrange)
        
        ll=cat(:,6)>=Mcrange(i);
        
        b1=(1/(mean(cat(ll,6))-(Mcrange(i)-0.05)))*log10(exp(1));
        
        sig1 = (sum((cat(ll,6)-mean(cat(ll,6))).^2))/(sum(ll)*(sum(ll)-1));
        sig1 = sqrt(sig1);
        sig1 = 2.30*sig1*b1^2;
        
        b(i,:)=[Mcrange(i),b1,sig1];
    end
    
    [histb,xout]=hist(b(:,2));
    histb=histb./length(b(:,1));
    
    % ---------------------------------------------------------------------
    % Calculation of NLIndex:
    %      1) for each cutoffmagnitude, for which at least 'binnumb' b-values
    %         from higher cutoff magnitudes have been calculated,
    %      2) calculate the standard variation in the b-values,
    %         and divide by the largest individual std(b) (mark)
    %      3) divide by number of mcutbins underlying the estimate (markw)
    %      4) determine slope and intercept of b(Mc) fit (slope,intercept)
    % (save in matrix 'marker')
    % ---------------------------------------------------------------------
    
    k=0;
    
    for mcut=Mcmin:0.1:max(Mcrange)-(0.1*binnumb-0.2)
        
        k=k+1;
        
        NLIndex=std(b(k:end,2))/max(b(k:end,3));
        NLIndexw=1/(length(Mcrange)-k)*NLIndex;
        
        h1=robustfit(b(k:end,1),b(k:end,2));
        slope(k)=h1(2);
        intercept(k)=h1(1);
        
        marker(k,:)=[mcut,b(k,2),b(k,3),NLIndex,NLIndexw,slope(k)];
        
    end
    
    
    
    hhnonlin=marker(:,4)>1;
    numbnonlin=sum(hhnonlin);
    
    %----------------------------------------------------------------------
    % interpret marker values, depending on chosen 'mode'
    %----------------------------------------------------------------------
    
    switch mode
        
        %------------------------------------------------------------------
        % use the pre-defined completeness magnitude
        %   i.e. only interpret first line of marker matrix
        %------------------------------------------------------------------
        case 'PreDefinedMc'
            
%             disp('use pre-defined Mcmin')
            
            %--------------------------------------------------------------
            % use given Mc, and corresponding b-value, calculate a-value
            %   and rate of target magnitude events
            %--------------------------------------------------------------
            
            bestmc=marker(1,1);
            bestb=marker(1,2);
            besta=log10(Ncum(1))+bestb*bestmc;
            NMc=Ncum(1);
            NMtarg=10^(besta-bestb*Mtarg);
            
            %--------------------------------------------------------------
            % if marker for Mc is <=1, the fit is ok
            %--------------------------------------------------------------
            
            if marker(1,4)<=1;
                
                result_flag=3;
                
                
                %--------------------------------------------------------------
                % if marker for Mc is >1, the fit is not ok
                %--------------------------------------------------------------
                
            else
                
                %----------------------------------------------------------
                % if trend of single b-value estimates is increasing
                %   --> extrapolation overestimates
                %----------------------------------------------------------
                
                if marker(1,6)>slope_pos
                    
                    result_flag=7;
                    
                    
                    %----------------------------------------------------------
                    % if trend of single b-value estimates is decreasing
                    %   --> extrapolation underestimates
                    %----------------------------------------------------------
                    
                elseif marker(1,6)<slope_neg
                    
                    result_flag=6;
                    %----------------------------------------------------------
                    % if trend of single b-value estimates is flat
                    %   --> extrapolation unstable
                    %----------------------------------------------------------
                    
                else
                    
                    result_flag=5;
                    
                    
                end
                
                
                
                %----------------------------------------------------------
                % output figure: subplot lower right corner
                %   histogram of b-values for different cutoff magnitudes
                %       (blue)
                %   mean b-value (blue line) +/- stdvar (dotted blue)
                %   mean b-value +/- max individual error (dotted red)
                %----------------------------------------------------------
                
                
            end
            
            %------------------------------------------------------------------
            % optimize suggested completeness magnitude for linearity
            %------------------------------------------------------------------
            
        case 'OptimizeMc'
            
            
            hhlin=marker(:,4)<=1;
            
            %--------------------------------------------------------------
            % if at least half of the estimated b-values are based on a
            % linear fit
            %--------------------------------------------------------------
            
            if numbnonlin<=0.5*length(marker(:,1))
                
%                 disp('FMD extrapolation realistically estimates M6+ rates')
                result_flag=3;
                
                %----------------------------------------------------------
                % use only estimates with linear fit:
                %   estimate bestMc from the minimum marker value (scaled by
                %     number of estimates) --> most linear fit
                %   estimate bestb as the median b-value of all linear
                %     estimates
                %   calculate besta and the rate of target magnitudes
                %----------------------------------------------------------
                
                markeracc=marker(hhlin,:);
                [~,i2]=min(markeracc(:,5));
                markeracc=markeracc(i2:end,:);
                bestmc=markeracc(1,1);
                bestb=nanmedian(markeracc(:,2));
                
                i1=find(round(10*Mrange)==round(10*bestmc));
                besta=log10(Ncum(i1))+bestb*bestmc;
                NMc=Ncum(i1);
                NMtarg=10^(besta-bestb*Mtarg);
                
                %----------------------------------------------------------
                % report if bestMc is significantly larger than Mmin
                %----------------------------------------------------------
                
                if bestmc-Mrange(1)>=sigMcDif
                     disp('Optimized Mc significantly higher than suggested Mc')
                    result_flag=4;
                end
                
                
                %--------------------------------------------------------------
                % if more than half of the estimated b-values are based on a
                % non-linear fit
                %--------------------------------------------------------------
                
            else
                
                %----------------------------------------------------------
                % estimate bestMc as minimum Mc, bestb as median of all
                % estimates, and calculate besta and rate of target events
                %----------------------------------------------------------
                
                bestmc=Mrange(1);
                bestb=nanmedian(marker(:,2));
                besta=log10(Ncum(1))+bestb*bestmc;
                NMc=Ncum(1);
                NMtarg=10^(besta-bestb*Mtarg);
                
                %----------------------------------------------------------
                % if trend of single b-value estimates is increasing
                %   --> extrapolation overestimates
                %----------------------------------------------------------
                
                if marker(1,6)>slope_pos
                    
                    
                    result_flag=7;
                    
                    
                elseif marker(1,6)<slope_neg
                    
                    result_flag=6;
                    
                    
                else
                    
                    result_flag=5;
                    
                    
                end
                
               
            end
            
    end
    
end

