function [con, lb, ub, props] = constraintmerge(allcon)
% [con, lb, ub] = constraintmerge(allcon)
%
% Rearranges constraints and bounds into separate structures (essentially
% flipping the nested structure of (name).{con, lb, ub} to {con, lb, ub}.(name).
%
% The single struct argument must have entries with fields 'con', 'lb', and
% 'ub'.
con = struct();
lb = struct();
ub = struct();
props = struct();

connames = fieldnames(allcon);
for i = 1:length(connames)
    n = connames{i};
    thiscon = allcon.(n);
    if ~isempty(thiscon.con)
        if ~checksize(thiscon.con, thiscon.lb, thiscon.ub);
            error('Inconsistent sizes for %s constraints!', n);
        end
        con.(n) = thiscon.con;
        lb.(n) = thiscon.lb;
        ub.(n) = thiscon.ub;
        if isfield(thiscon, 'props')
            props.(n) = thiscon.props;
        end
    end
end

end%function

function tf = checksize(con, lb, ub)
    % Checks that the given sizes are consistent.
    % TODO: check that this covers all edge cases.
    if isvector(con)
        con = numel(con);
    else
        con = size(con);
    end
    tf = isequal(con, boundsize(lb), boundsize(ub));
end%function

function s = boundsize(bound)
    % Returns a size vector that can be compared to the constraints.
    s = size(bound);
    s = s(2:end);
end%function

