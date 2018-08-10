function Cc = sx_clean_missing(timeVec, zVec, C);

%% fill missing values because of new ALSEAMAR method to deal with
%% repetitive data.

Cc = nan(size(C));
%% Fill NaNs
for i = 1:size(C,2)
    I = find(~isnan(C(:,i)));
    if length(I)>1
        itp = interp1(zVec(I), C(I,i), zVec);
        Cc(:,i) = itp;
    end
end
