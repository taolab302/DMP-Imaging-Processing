function save_neuron_pos(pos, filename)

fid = fopen(filename,'wt');
for j=1:size(pos,1)
    fprintf(fid,'%d    %d\n',pos(j,1),pos(j,2));
end
fclose(fid);
end

