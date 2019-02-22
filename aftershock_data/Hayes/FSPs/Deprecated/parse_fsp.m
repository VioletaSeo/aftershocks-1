function [Loc, Size, Mech, Rupt, Inv, VelocityStruct, Multisegment,finiteFault] =  ...
    parse_fsp(filename)

% Read commented section starting with  % sign
%   Add in parameters as I go:
% 1) Event info
% 2) Inversion related parameters
% 3) Velocity structure

fid = fopen(filename);

%initialize variables

aftercolon  = @(in) in((find(in == ':')+1):end);

% 'read_nums' looks for the numbers 2 entries down from the specifies input

fgetl(fid);                                                             %1
fgetl(fid);                                                             %2
Event       = aftercolon(fgetl(fid));                                   %3
Tag         = aftercolon(fgetl(fid));                                   %4
fgetl(fid);                                                             %5

Loc         = read_nums(fgetl(fid), 'LAT' , 'LON', 'DEP');              %6
Size        = read_nums(fgetl(fid), 'LEN' , 'WID', 'Mw', 'Mo');         %7
Mech        = read_nums(fgetl(fid), 'STRK', 'DIP', 'RAKE', 'Htop');     %8
Rupt        = read_nums(fgetl(fid), 'HypX', 'Hypz', 'avTr', 'avVr');    %9

fgetl(fid);                                                             %10
fgetl(fid);                                                             %11
fgetl(fid);                                                             %12

Invs1       = read_nums(fgetl(fid), 'Nx', 'Nz', 'Fmin', 'Fmax');        %13
Invs2       = read_nums(fgetl(fid), 'Dx', 'Dz');                        %14
Invs3       = read_nums(fgetl(fid), 'Ntw', 'Nsg');                      %15
Invs4       = read_nums(fgetl(fid), 'LEN', 'SHF');                      %16
Inv         = [Invs1,Invs2,Invs3,Invs4];

skipl(11)

numberOf    = read_nums(fgetl(fid), 'layers');

skipl(4)

varNames = {'Depth','PVEL','SVEL','DENS','QP','QS'};
velRaw = zeors(numberOf.layers,length(varNames));
for n = 1:numberOf.layers
    line = fgetl(fid);
    line = line(2:end);
    velRaw(:,n) = str2double(split(line))';
end
VelocityStruct = table(velRaw,'VariableNames',varNames);

% skip lines till start of inversion
Multisegment = false;
while 1
    line = fgetl(fid);
    spltLine = split(line);
    if ismenber('MULTISEGMEMT', spltLine) 
        Multisegment = true; 
        break % give up
    end
    
    if ismember('Nsbfs', spltLine); Nsbfs = num2down(spltLine,'Nsbfs'); end
    
    % start inversion
    if line(1) ~= '%'
        varNames = {'LAT','LON','X','Y','Z','SLIP','RAKE','TRUP','RISE','SF_MOMENT'};
        finiteFault = table('Size',          [Nsbfs,numel(varNames)] , ...
                               'VariableNames', varNames                );
        for n = 1:Nsbfs
            line = fgetl(fid);
            finiteFault(n,:) = str2double(split(line))';
        end
        break
    end
end
    
end


function out = num2down(cellStr,str)
out = str2double(cellStr(strcmp(cellStr,str)+2));
end

function structOut = read_nums(line,varargin)
structOut = [];
for n = 1:length(varargin)
    structOut.(varargin{m}) = num2down(split(line),varargin{n});
end
end

function skipl(N)
for n = 1:N
    fgetl(fid);
end
end
