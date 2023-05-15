function model = emptycon()
% model = emptycon()
%
% Returns a struct with field 'con' = {}, 'lb' = [], and 'ub' = []. Useful
% for omitting a certain type of constraint from the model.
model = struct();
model.con = {};
model.lb = zeros(0, 1);
model.ub = zeros(0, 1);
end%function

