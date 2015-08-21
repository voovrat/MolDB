function wait_command(fname)

while 1
   sleep(1)
   f=fopen(fname,'r');
   if f>0
      
      while ~feof(f)
         s=fgets(f);
         if sum(s<0)==1
            break
         end

         try
            eval(s)
         catch 
            fprintf('%s\n',[s ':error occured']);
            break
         end
      end
      fclose(f);
      f=fopen(fname,'w');
      fclose(f);
   end
end
