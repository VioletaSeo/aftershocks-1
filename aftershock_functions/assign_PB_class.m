 function PBASSIGNMENT    = assign_PB_class(lat, lon, criticalDistance, enforceFMSyn, fms)
        fn      = 'PB2002_steps.dat';
        fileID  = fopen(fn);
        C       = textscan(fileID,'%f %s %f %f %f %f %f %f %f %f %f %f %f %f %s');
        fclose(fileID);
        PBclass = C{end};
        PBclass = erase(PBclass,':');
        PBclass = erase(PBclass,'*');
        %expectedPBClass             = {'all','OSR','OTF','OCB','CRB','CTF','CCB','SUB'}; 
        % Boundary types: CCB continental convergent boundary, CTF continental
        % transform fault, CRB continental rift boundary, OSR oceanic spreading
        % ridge, OTF oceanic transform fault, OCB oceanic convergent boundary, SUB
        % subduction zone 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%% refine classification %%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [latPB,lonPB]   = C{[4,3]};
        
        % turn coord into earth centric to accelerate computation
        wgs84               = wgs84Ellipsoid('kilometer');
        [x  ,y  ,z  ]       = geodetic2ecef(wgs84,lat  ,lon  ,0);
        [xPB,yPB,zPB]       = geodetic2ecef(wgs84,latPB,lonPB,0);
        [Idx, D]            = knnsearch([xPB,yPB,zPB],[x,y,z]);
        
        pbAssignment = cell(length(Idx),1);
        % assign nearest plate boundary
        for n = 1:length(Idx)
                pbAssignment{n} = PBclass{Idx(n)};
        end
        
        % remore assignment more than the critical distance away from a pb
        pbAssignment(D>criticalDistance) = {'none'};
                
        if nargin > 4
            if strcmp(enforceFMSyn,'yes')
                % turn incosistent assignment to INTRAPLATE? - might be
                % weird
                for n = 1:length(pbAssignment)
                   ipb = pbAssignment{n};
                   ifms= fms(n);
                   switch ipb
                       case {'OSR','CRB'}
                           if ifms ~= 2, pbAssignment{n} = 'intraplate'; end
                       case {'OTF','CTF'}
                           if ifms ~= 1, pbAssignment{n} = 'intraplate'; end
                       case {'CCB','SUB','OCB'}
                           if ifms ~= 3, pbAssignment{n} = 'intraplate'; end
                   end
                end
            end
        end            
           
        PBASSIGNMENT = pbAssignment;
 end
