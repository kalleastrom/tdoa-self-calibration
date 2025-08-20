function sol = misstdoa_remove_outliers(sol,threshold)
% MISSTDOA_BRIEF_VISUALIZATION Plot status report for current solution
%keyboard;
[zcalc,zok] = misstdoa_calc_z(sol);
zerr = zcalc-sol.z;
zinl = sol.inlmatrix;
sol.inlmatrix = (~isnan(sol.z)) & abs(zerr)<threshold;
