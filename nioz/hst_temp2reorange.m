function hst_temp2reorange(tempFile, TSrelFile, prefix)

% function hst_temp2reorange(tempFile, prefix)
%
% Prepare files for reorange analysis.
% example to run in /home/cyrf0006/research/NIOZ/RockallBank
% usage ex:  > hst_temp2reorange('temp_wholeSerie_skip10.mat', 'TSrel_cubic', './reorange/fine_skip10')
%            > hst_temp2reorange('temp_wholeSerie_skip1000.mat', 'TSrel_cubic', './reorange/fine_skip1000')
%
% This function assume that the matrix 'Titp' exists in the .mat
% file.
%
%
% See also mat2reorange.m

% F. Cyr, April 2014
% --------------------------------------------------------- %


load(tempFile);
no_profiles = size(Titp, 2);

T = Titp;
S = hst_temp2sal(T, TSrelFile); 

SA = gsw_SA_from_SP(S, zVec, -15.1, 55.5);
CT = gsw_CT_from_t(SA,T,zVec);
SIG_T = gsw_sigma1(SA,CT);


for i = 1:no_profiles
    disp(sprintf('saving %d / %d', i, no_profiles));
    fname_out = [prefix '_' datestr(timeVec(i), 30) '.dat'];
    dlmwrite(fname_out, flipud([zVec, CT(:,i), SA(:,i), SIG_T(:,i)]), 'delimiter',' ','precision',10);
end


