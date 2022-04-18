function write_neuronpos(PosFolder,neuron_indice,xpos,ypos,write_op)
% write neuron positions to file
% write_op can be 'a' or 'w'(default)

if nargin == 4
    write_op = 'w';
end

if ~(strcmp(write_op,'a') == 1 || strcmp(write_op,'w') == 1)
    disp('Invalid write operation, only can be a/w');
    return;
end
if strcmp(neuron_indice, 'AVL') == 1 || strcmp(neuron_indice, 'DVB') == 1
    output_name = [PosFolder neuron_indice '.txt'];
        fid = fopen(output_name,write_op);
        for t = 1:length(xpos)
            fprintf(fid,'%d    %d\n',xpos(t),ypos(t));
        end
        fclose(fid);
else
    for n=1:length(neuron_indice)
        output_name = [PosFolder sprintf('neuron %02d',neuron_indice(n)),'.txt'];
        fid = fopen(output_name,write_op);
        for t = 1:size(xpos,1)
            fprintf(fid,'%d    %d\n',xpos(t,n),ypos(t,n));
        end
        fclose(fid);
    end
end