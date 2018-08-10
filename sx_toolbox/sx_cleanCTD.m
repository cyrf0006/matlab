function [Cc, Tc, Pc] = sx_cleanCTD(timeVec, zVec, C,T,P);



%% Despike (based on standard dev.)
% $$$ Cdspk = C;
% $$$ for i = 1:size(C,1)
% $$$     I = find( abs(C(i,:) - nanmean(C(i,:))) > 3*nanstd(C(i,:)));
% $$$     if ~isempty(I)
% $$$         Cdspk(i,I) = NaN;
% $$$     end
% $$$ end
% $$$ C = Cdspk;

%% Fill NaNs
for i = 1:size(C,2)
    I = find(~isnan(C(:,i)));
    if length(I)>2
        C(:,i) = interp1(zVec(I), C(I,i), zVec);
        T(:,i) = interp1(zVec(I), T(I,i), zVec);
        P(:,i) = interp1(zVec(I), P(I,i), zVec);
    end
end


Cc = C;
Tc = T;
Pc = P;
    
    
    