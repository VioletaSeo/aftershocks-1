function [DATENUMBER] = parsetime(timeArray)
% This is a small function snippet to parse through time in a catalog, for
% example

% * this is a REALLY bad way of doing this... 

%  input should be an array of strings

%                 num      -     num     -  num      T    num      :       num     :    num        .   num  
expression = ['(?<year>\d+)-(?<month>\d+)-(?<day>\d+)T(?<hours>\d+):(?<minutes>\d+):(?<seconds>\d+).(?<miliseconds>\d+)'];
tokenNames = regexp(timeArray,expression,'names');

dateNumber = zeros(1,length(timeArray));

for i = 1:length(timeArray)
    timeCell  = tokenNames(i);
    timeStruct= timeCell{:};
    timeVec   = [str2num(char(timeStruct.year)),    ...
                           str2num(char(timeStruct.month)),    ...
                           str2num(char(timeStruct.day)),      ...
                           str2num(char(timeStruct.hours)),    ...
                           str2num(char(timeStruct.minutes)),  ...
                           str2num(char(timeStruct.seconds)),  ...
                           str2num(char(timeStruct.miliseconds))];
    
    dateNumber(i)= datenum(timeVec(1:end-1));
end

DATENUMBER = dateNumber;

end

