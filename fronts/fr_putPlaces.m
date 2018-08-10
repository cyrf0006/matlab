function fr_putPlaces(region)
    
    
filename = '/home/cyrf0006/research/fronts/matlab_workspace/placeID';
run(filename);


I = find(strcmp(placeStruct(:,1), region) == 1);
S = placeStruct(I,:);

for i = 1:size(S,1)
    m_text(cell2mat(S(i,3)), cell2mat(S(i,2)), char(S(i,5)), ...
           'rotation', cell2mat(S(i,4)), ... 
           'verticalAlignment', char(S(i,6)), ... 
           'horizontalAlignment', char(S(i,7)), ...
           'fontweight', 'bold', ...
           'fontsize', 10, ...
           'color', [.6 0 0]);
end

