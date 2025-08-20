function [solout, inlier_count] = extend_vbo_ransac(sol, varargin)
% EXTEND_VBO_RANSAC Extend a relaxed TDOA solution over senders
%   [solout, inlier_count] = EXTEND_VBO_RANSAC(sol) extends a relaxed
%       solution with one column in v, b and o using robust methods.
%   [solout, inlier_count] = EXTEND_VBO_RANSAC(sol, ...) additional inputs:
%           iters - number of iterations in RANSAC loop.
%           min_inliers - minimum number of inliers in new column.
%           threshold - threshold for the absolute error in TDOA measurment
%               when classifying inliers/outliers.

% The problem involves estimating v with 2 or 3 elements depending on the
% rank of the problem as well as a with one element
% For TDOA - also estimate o
% For COTOA use the offset from other columns
% For TOA - o is zero.

% Parse inputs.
p = inputParser;
addParameter(p, 'iters', 100);
addParameter(p, 'min_inliers', 3);
addParameter(p, 'threshold', 0.2);
parse(p, varargin{:});
opts = p.Results;

% Default output.
solout = sol;
inlier_count = 0;

z = sol.z;
max_inliers = 0;
bestsol = sol;

possible_cols = 1:size(sol.z, 2);
possible_cols(sol.cols) = [];

switch sol.offset_type
    case 'tdoa'
        min_data = sol.rank+2;
    case {'toa','cotoa'}
        min_data = sol.rank+1;
end
% We need min_data rows for finding vnew, bnew and onew, and at least 1 for testing.
possible_cols(sum(isfinite(z(sol.rows, possible_cols)), 1) < min_data+1) = [];

if isempty(possible_cols)
    return;
end

for i = 1:opts.iters
    newcol = possible_cols(randi(length(possible_cols)));
    possible_rows = sol.rows(isfinite(z(sol.rows, newcol)));
    use_rows = possible_rows(randperm(length(possible_rows), min_data));
    test_rows = setdiff(possible_rows, use_rows);

    % Get indices in sol.rows for use_rows and test_rows.
    rowsorder = zeros(1, length(use_rows));
    for k = 1:length(use_rows)
        rowsorder(k) = find(sol.rows == use_rows(k));
    end
    rowsorder_test = zeros(1, length(test_rows));
    for k = 1:length(test_rows)
        rowsorder_test(k) = find(sol.rows == test_rows(k));
    end

    % Estimate new cols in v, b and o using z(use_rows, newcol).
    usub = sol.u(rowsorder, :);
    asub = sol.a(rowsorder);
    zsub = z(use_rows, newcol); % lurigt. U har r?tt ordning

    % Slightly different cases for different offset_types
    % AA = [-2 * zsub, ones(5, 1), -usub, -ones(5, 1)];
    % bb = asub - zsub.^2;
    % xx_part = AA \ bb;
    % xx_hom = [0, 1, 0, 0, 0, 1]';
    % onew = xx_part(1);
    % lamb = onew^2 - xx_part(2);
    % xx = xx_part + lamb * xx_hom;
    % vnew = xx(3:5) / (-2); % TODO: Add -2 to AA instead.
    % bnew = xx(6);
    switch sol.offset_type
        case 'tdoa'
            AAA = [-2 * zsub, ones(min_data, 1), (-usub), -ones(min_data, 1)];
            bbb = asub - zsub.^2;
            xx_part = AAA \ bbb;
            xx_hom = [0, 1, zeros(1,sol.rank), 1]';
            onew = xx_part(1);
            lamb = onew^2 - xx_part(2);
            xx = xx_part + lamb * xx_hom;
            vnew = xx(3:(3+sol.rank-1)) / (-2); % TODO: Add -2 to AA instead.
            bnew = xx(end);
        case 'cotoa'
            % (z(cc,j)-ony).^2  is equal to -2*(u*vny)+a+bny
            onew = sol.o(1);
            d2sub = (zsub-onew).^2;
            AAA = [(-2*usub),ones(min_data, 1)];
            bbb = d2sub - asub;
            xxx = AAA \ bbb;
            vnew = xxx(1:sol.rank); %
            bnew = xxx(end);
        case 'toa'
            % (z(cc,j)-ony).^2  is equal to -2*(u*vny)+a+bny
            onew = 0;
            d2sub = (zsub-onew).^2;
            AAA = [(-2*usub),ones(min_data, 1)];
            bbb = d2sub - asub;
            xxx = AAA \ bbb;
            vnew = xxx(1:sol.rank); %
            bnew = xxx(end);
    end;

    % Evaluate using rows in test_rows.
    utest = sol.u(rowsorder_test, :);
    atest = sol.a(rowsorder_test);
    err = z(test_rows, newcol) - (sqrt(-2 * utest * vnew + atest + bnew) + onew);

    inind = abs(err) < opts.threshold;
    nrin = sum(inind);

    if nrin > max_inliers
        max_inliers = nrin;
        bestsol = sol;
        bestsol.cols = [bestsol.cols, newcol];
        bestsol.v = [bestsol.v, vnew];
        bestsol.b = [bestsol.b, bnew];
        bestsol.o = [bestsol.o, onew];
        inlrows = [use_rows, test_rows(inind)];
        bestsol.inlmatrix(inlrows, newcol) = 1;
    end
end

if max_inliers >= opts.min_inliers
    solout = bestsol;
    inlier_count = max_inliers;
end
