function s = structcat(varargin)
% Function that concatenates structures and add empty fields for fields
% that are not defined in each structure.
% s = structcat(s_a,s_b,s_c,...)
% s_a,s_b,s_c = potentially dissimilar structures
%%
N_structs = numel(varargin);

% Get the fieldnames of each structure.
Fieldnames = cell(N_structs,1);
Fieldnames_all = {};
for i=1:N_structs
    Fieldnames{i} = fieldnames(varargin{i});
    Fieldnames_all = [Fieldnames_all; Fieldnames{i}];
end
Fieldnames_all = unique(Fieldnames_all,'stable');
N_fields = numel(Fieldnames_all);

% Add empty fields to structures that do not possess the fields of other
% structures
for i=1:N_structs
    for j=1:N_fields
        if ~isfield(varargin{i},Fieldnames_all{j})
            varargin{i}(1).(Fieldnames_all{j})=deal([]);
        end
    end
end

% Concatenate the similar structures.
s = [varargin{:}];
end