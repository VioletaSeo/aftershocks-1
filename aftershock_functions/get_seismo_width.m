function MSseismo = get_seismo_width()
%% get the local seismogenic thickness of an earthquake. 
% calculates the 10th to 90th percentile range of the declustered
% seismicity locally around every earthquake in the the mainshock_info.mat
% catalog. At the moment the search radius defining the local seismicity is
% set to be a factor (4) of the source radius assuming constant stress drop
% scaling and circular rupture. No explicit input but mainshock_info.mat
% (as derived from the function afterhsock_productivity_kernel.m) must be
% on the path

load mainshock_info.mat

% basic data
numMS           = length(MSt);
MSseismo        = zeros(numMS,1);

for n = 1:numMS
    mapD        = distance(MSlat,MSlon,MSlat(n),MSlon(n));
    D           = sqrt(deg2km(mapD).^2 + (MSdepth-MSdepth(n)).^2);
    spaceWindow = getspacewindow('dynamic',MSmag(n));
    Iin = 0;
    while sum(Iin) < 100 % ensure stable estimation of depth range
        Iin         = D < spaceWindow;
        spaceWindow = 2*spaceWindow;
    end
    MSdepthRange= prctile(MSdepth(Iin),[10,90]);
    MSseismo(n) = diff(MSdepthRange);        
end
end

%% accessory functions

function SPACEWINDOW = getspacewindow(varargin)
% define the radius of the space window.
% Usage: - one of -
% getspacewindow('dynamic',M,MUlT)
%   where 

% dynamic (magnitude dependent) space window
if strcmp(varargin{1},'dynamic') && (nargin == 2 || nargin == 3) 
    M       = varargin{2};
    if nargin == 3; sourceScale = varargin{3}; else; sourceScale = 4; end
    mo          = (10^((M+10.73)*1.5))*10^-7;
    stressDrop  = 10^5*10;
    sourceDim   = 2 * 0.001 * (mo*(7/16) / stressDrop) ^ (1/3); % in km tbl 9.1 modern global seismo lay
    SPACEWINDOW = sourceScale*sourceDim;                    % 4 as per default setting in

% fixed space window as specified by the single input
elseif isnumeric(varargin{1}) && nargin == 1
    SPACEWINDOW = varargin{1};

else
    error('input to getspacewindow is wrong')
end
end