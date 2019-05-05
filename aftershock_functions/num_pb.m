function NPB    = num_pb(lat, lon, criticalDistance)
        fn      = 'PB2002_steps.dat';
        fileID  = fopen(fn);
        C       = textscan(fileID,'%f %s %f %f %f %f %f %f %f %f %f %f %f %f %s');
        fclose(fileID);
        
        [latPB,lonPB]   = C{[4,3]};
        
        % turn coord into earth centric to accelerate computation
        wgs84               = wgs84Ellipsoid('kilometer');
        [x  ,y  ,z  ]       = geodetic2ecef(wgs84,lat  ,lon  ,0);
        [xPB,yPB,zPB]       = geodetic2ecef(wgs84,latPB,lonPB,0);
        
        N = length(x);
        NPB = zeros(N,1);
        
        for n = 1:N
            NPB(n) = sum(sqrt((x(n)-xPB).^2+(y(n)-yPB).^2+(z(n)-zPB).^2) < criticalDistance); 
        end           
 end
