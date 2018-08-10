function S = hst_temp2sal(T, TSrelFile)

% function S = hst_temp2sal(T, TSrelFile)
%
% usage ex: S = hst_temp2sal(T, 'TSrel_cubic.mat')
%    OR
%           S = hst_temp2sal(T, TSrelFile) in hst_temp2turbulence.m
% 
% NOTE THAT IT CAN BE USED FOR ANY RELATIONSHIIP (E.G. Temp2rho...)
    
load(TSrelFile); % variable 'p' must exist in this file    


if length(p) == 2 % linear
    S = T*p(1) + p(2);    
elseif length(p) == 3
    disp(['Problem with TS relation file: quadratic rel. not implemented yet!'])
elseif length(p) == 4
    S = T.^3*p(1) + T.^2*p(2) + T*p(3) + p(4);
else
    disp('Problem with TS relation file')
end
    