function Refs=moldb_getRef(Chain,restrictions)


if nargin<2
    restrictions='';
else
    restrictions = [ '''' restrictions '''' ];
end

system([ 'moldb_getRef ' Chain ' ' restrictions ' > out.txt ']);

f=fopen('out.txt');
data = textscan(f,'%s');
folders = data{1};
fclose(f);

N = length(folders);

Refs = cell(N,1);

for i=1:N

    folder = folders{i};

    Refs{i} = struct;
    Refs{i}.folder = folder;

    if ~exist([folder '/PARAMETERS'])

       warning([ 'NOT exist:' folder '/PARAMETERS'])
       Refs{i}.parameters = struct;
    else

       parameters = moldb_loadParameters([ folder '/PARAMETERS' ]);
       Refs{i}.parameters = parameters;
   end
end

