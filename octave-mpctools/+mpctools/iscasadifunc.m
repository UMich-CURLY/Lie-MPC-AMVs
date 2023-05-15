function tf = iscasadifunc(x)
    % tf = iscasadifunc(x)
    %
    % Returns true or false whether x is a casadi.Function.
    tf = isa(x, 'casadi.Function') || isa(x, 'Function');
end%function
