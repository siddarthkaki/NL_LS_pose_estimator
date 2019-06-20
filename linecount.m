function number_of_lines = linecount(filename)
%% open text file
fid = fopen(filename);

%% init variable(s)
res={};

%% iterate to count number of lines
while ~feof(fid)
    thisline = fgetl(fid);
    if ~ischar(thisline); break; end
    res{end+1,1} = thisline;
end

%% close text file
fclose(fid);

%% create output variable
number_of_lines = numel(res);

end