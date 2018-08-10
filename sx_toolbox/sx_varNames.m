function sx_varnames(structFile)
    
% function sx_varnames(gliderFile)
%
% Function use to quickly check variables names in a structure file 
% 

load(structFile);

w = whos;
structName = w.name;
command = sprintf('s = %s ;',structName);
eval(command);

disp('Science data variables:')
s.datVarNames

disp('Log/pilot variables:')
s.logVarNames

