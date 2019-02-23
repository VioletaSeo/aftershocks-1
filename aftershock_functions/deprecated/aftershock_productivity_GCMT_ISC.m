function CAT = aftershock_productivity_GCMT_ISC(varargin)
% wrapper around aftershock productivity kernel to produce. Inputs are
% directly passed along to aftershock_productivity_kernel.m. See the
% associated documentation for information as to the format of the input.
% Cannot use the 'MinMainshockMag' and 'SaveCatalog' functionality

% RECOMMENDED USE:
% CAT = aftershock_productivity_GCMT_ISC('ForeshockTimeWindow', 10, 'DepthRange', [0,55]);

errorM   = @(str) sprintf('Do not use parameter ''%s'' in the merge process, functionality not allowed',str);
checkInp = @(par) assert(~any(strcmp(varargin,par)), errorM(par));

checkInp('MinMainShockMag');
checkInp('SaveCatalog');
checkInp('PlotYN');

% PATH DEPENDENCIES
addpath(genpath('/auto/home/kdascher/Documents/MATLAB/GISMO'));
load('ISC_1976-2018.mat'); % loads quakes 
MSformat            = {'MSt','MSlat','MSlon','MSmag'}; % standard format out of the analysis

% magic numbers
minMainshock        = 6.0; % minimum magnitude for the gCMT catalog
magnitudeBuffer     = 1; % magnitude error buffer allowed for catalog merge

% _________________________________________________________________________

tempCatFN           = 'temp_earthquake_cat_merge.mat';
delete_temp         = @() system(sprintf('rm %s',tempCatFN));

aftershock_productivity_kernel('eq_catalog.mat', ...
                               'SaveCatalog',tempCatFN,...
                               'MinMainshockMag',minMainshock, ...
                               'PlotYN', 'no', ...
                               varargin{:});
                                                     
GCMT_CAT = load(tempCatFN); delete_temp();

minMag = min(GCMT_CAT.MSmag);
aftershock_productivity_kernel( quakes.otime, ... 
                                quakes.lat, ...
                                quakes.lon, ...
                                quakes.depth, ...
                                quakes.mag, ...
                                'SaveCatalog',tempCatFN, ...
                                'MinMainshockMag',minMag-magnitudeBuffer, ... % ensure that mainshocks with different magnitudes get merged approriately
                                'PlotYN', 'no', ...
                                varargin{:});
ISC_CAT = load(tempCatFN); delete_temp();

[mergedCat,I] = merge_eq_catalog(struct2table(GCMT_CAT),MSformat, ...
                                 [1,1,1,1], ...
                                 struct2table(ISC_CAT),MSformat);
shCat = mergedCat(I,:); 
CAT = shCat;

end

