function TFilled  =iow_fillNaN(T, n)

if size(n,1)~=1
    n = n';
end
if size(T,2)~=length(n)
    T = T';
end
    
TFilled = nan(size(T));
for i = 1:size(T,1)
    %    disp(sprintf('clean sensor %d/%d', i, size(T,1)))
    TVec = T(i,:);
    I = find(~isnan(TVec));
    if ~isempty(I)
        TFilled(i,:) = interp1(n(I), TVec(I), n); 
    end
end
   
