function I = moldb_filter(Refs,Cell)

Names = Cell(1:2:end);
Values = Cell(2:2:end);

Nargs = length(Names);
N = length(Refs);
I = zeros(N,1);

for i=1:N

    prm = Refs{i}.parameters;
    Ok = true;

    for j = 1:Nargs

        arg = Names{j};
        val = Values{j};

        S = sum(strcmp(prm.(arg),val) > 0);

        if ~S
            Ok = false; 
            break;
        end


    end
    
    if Ok 
        I(i) = 1;
    end

end

I = find(I);
