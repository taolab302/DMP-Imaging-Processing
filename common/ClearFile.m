function ClearFile(filename)
% Delete file contents

fid = fopen(filename,'w');
fprintf(fid,'');
fclose(fid);

end