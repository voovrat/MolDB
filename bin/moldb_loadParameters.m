function parameters = moldb_loadParameters(fname)


f=fopen(fname);
data = textscan(f,'%s');
lines = data{1};
fclose(f);

N = length(lines);
parameters = struct;

for i=1:N

    NameValue = strread(lines{i},'%s','Delimiter','=');
    parameters.(NameValue{1}) = NameValue{2};

end


