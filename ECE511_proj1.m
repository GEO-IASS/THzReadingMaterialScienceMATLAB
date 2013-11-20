function [t,E] = open_picotd(filename,rows)
% This function opens picotd files from the T-Ray 4000.
% e.g. [t,E]=open_picotd('6_8.picotd',4101)
% INPUTS
% 'filename' = filename to be opened, in single quotes
% rows = total number of rows in the data file
%
% OUTPUTS
% t = time delay, picoseconds 
% E = Electric Field Strength, arbritrary units
%
% Scott Schecklman
% Nov. 13, 2008

fid = fopen(filename);
% skip over header information
for n = 1:5
    stuff = fgetl(fid);
    if n == 2
        delta_t = strread(stuff, '%*s %*s %*s %*s %*s %8c %*s %*s %*s %*s %*s');
                % read time increment, (picoseconds), from the header file
        delta_t = str2num(delta_t); %convert variable type from char to number
    end
end

%create time axis
t_max = (rows-5)*delta_t;
t = 0:delta_t:t_max; 
t = t';

%read E field measurements
E = zeros(rows-5,1); %pre-allocate
for n = 1:rows-5
    stuff = fgetl(fid);
    data = str2num(stuff);
    E(n) = data(2);
end

fclose(fid);

end

[t,E]=open_picotd('Silicon.picotd',4096);