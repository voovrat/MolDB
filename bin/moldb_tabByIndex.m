function [Tab,Keys]=moldb_tabByIndex(TabIndex,Values)
%
%  Function moldb_group returns the table: cell array of structs 
%  with values of some given parameter AND indexes of these groups
% 
%  This function allows you to get Values from these indeces
%  columns correspond to diffenent keys


Keys = fieldnames(TabIndex);

N = length(Keys);
M0 = 0;
rect=1;

CTab = cell(1,N);

for i=1:N

    key = Keys{i};

    C = getfield(TabIndex,key);

    M = length(C);

    if M~=M0 && M0~=0 
        %M 
        %M0
        rect = 0;
    end

    M0=M;

   % CTab{i} = mycellfun(fn,C);

    for j=1:M

        if iscell(Values)
            CTab{i}{j} = Values{ C(j)  };
        else
            CTab{i}{j} = Values( C(j) );
        end
    end

end

if rect

    Tab = cell(M,N);

    for i=1:M
        for j=1:N

            Tab{i,j} = CTab{j}{i};
        end
    end


else
    Tab = CTab;
end
