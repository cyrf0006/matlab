function sx_rawprof2matprof(filename, outfile)
      
[meta, data] = sx2mat(filename);

keyboard

varNames = [];
for i  = 1:numel(meta.variables)
    command = sprintf('%s = %d', meta.variables{i}, i);
    eval(command);
end

toSave = meta.variables;
toSave = [toSave; 'data'];
sensors_fich = meta.variables;
toSave = [toSave; 'sensors_fich'];



save(outfile, char(toSave'));
