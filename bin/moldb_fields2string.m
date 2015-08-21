function str = moldb_fields2string(param_struct,field_list)
%
%  moldb_fields2string(param_struct,field_list)
%
%  Get all the fields in field_list and put them to the slash-separated list
%
%

if ischar(field_list)

    str = getfield(param_struct,field_list);

else

   N=length(field_list);
   str='';
   for i=1:N

      val = getfield(param_struct,field_list{i});
      str = [ str val ];
      if i<N
          str=[ str '/' ];
      end

   end


end
