function thick = get_litho(lat,lon,LAYER)
% get the depth of the a earth layer (in km) using the litho 1.0 model.
% Inputs are an array of latitudes and longitudes. (optional: specify layer
% name). Layer names include:

    %  LID-BOTTOM
    %  LID-TOP
    %  CRUST3-BOTTOM
    %  CRUST3-TOP
    %  CRUST2-BOTTOM
    %  CRUST2-TOP
    %  CRUST1-BOTTOM
    %  CRUST1-TOP
    %  SEDS2-BOTTOM
    %  SEDS2-TOP
    %  SEDS1-BOTTOM
    %  SEDS1-TOP
    %  WATER-BOTTOM
    %  WATER-TOP

    % default is LID-BOTTOM
    
% has path dependencies and requires a working and compiled version of
% litho1.0
path2litho = '/auto/home/kdascher/Documents/UCSC/general_data/Tectonic_data/LITHO1.0';
thick = zeros(size(lat));

if nargin == 2
    LAYER = 'LID-BOTTOM';
end

for n = 1:length(lat)
    [~,textOut,~] = jsystem(sprintf('%s/bin/access_litho -p %0.1f %0.1f', path2litho, lat(n), lon(n)));
    try
        C = textscan(textOut,'%f %f %f %f %f %f %f %f %f %s');
    
    % example output:
    %  142440.  3300.00  8015.90  4356.47    0.00   70.00  8015.90  4356.47 1.00000 ASTHENO-TOP
    %  142440.  3300.00  8179.49  4660.68    0.00  200.00  8179.49  4660.68 1.00000 LID-BOTTOM
    %   33717.  3300.00  8179.49  4660.68    0.00  200.00  8179.49  4660.68 1.00000 LID-TOP
    %   33717.  2935.92  6886.83  3936.30    0.00  600.00  6886.83  3936.30 1.00000 CRUST3-BOTTOM
    %   22317.  2935.92  6886.83  3936.30    0.00  600.00  6886.83  3936.30 1.00000 CRUST3-TOP
    %   22317.  2783.48  6370.11  3662.02    0.00  600.00  6370.11  3662.02 1.00000 CRUST2-BOTTOM
    %   12290.  2783.48  6370.11  3662.02    0.00  600.00  6370.11  3662.02 1.00000 CRUST2-TOP
    %   12290.  2654.57  5817.68  3306.18    0.00  600.00  5817.68  3306.18 1.00000 CRUST1-BOTTOM
    %    2456.  2654.57  5817.68  3306.18    0.00  600.00  5817.68  3306.18 1.00000 CRUST1-TOP
    %    2456.  2370.00  4000.00  2130.00    0.00  600.00  4000.00  2130.00 1.00000 SEDS2-BOTTOM
    %    1711.  2370.00  4000.00  2130.00    0.00  600.00  4000.00  2130.00 1.00000 SEDS2-TOP
    %    1711.  2051.29  2336.62   904.07    0.00  600.00  2336.62   904.07 1.00000 SEDS1-BOTTOM
    %     805.  2051.29  2336.62   904.07    0.00  600.00  2336.62   904.07 1.00000 SEDS1-TOP
    %     805.  1020.00  1500.00     0.00    0.00  600.00  1500.00     0.00 1.00000 WATER-BOTTOM
    %    -397.  1020.00  1500.00     0.00    0.00  600.00  1500.00     0.00 1.00000 WATER-TOP

    
        % The output is:
        %  depth(m) density(kg/m3) Vp(m/s) Vs(m/s) Qkappa Qmu Vp2(m/s) Vs2(m/s) eta layername
        depthAry = C{1};
        tagAry   = C{end};
        tagLoc   = strcmp(LAYER,tagAry);
        waterBotTag = strcmp('WATER-BOTTOM',tagAry);
        waterTopTag = strcmp('WATER-TOP',tagAry);
        thick(n) = (depthAry(tagLoc))/1000;

        if any(waterBotTag) 
            thick(n) = thick(n)-((depthAry(waterBotTag)-depthAry(waterTopTag))/1000);
        end
    catch
        thick(n) = nan;
    end
        
end

end

