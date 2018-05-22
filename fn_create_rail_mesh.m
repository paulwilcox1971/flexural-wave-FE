function [nodes, elements, sleeper_nodes, forcing_node] = fn_create_rail_mesh(nominal_length, sleeper_spacing, forcing_position, max_node_spacing, length_option)
%create mesh to describe length of rail
no_sleepers = ceil(nominal_length / sleeper_spacing);
% nodes_per_sleeper = max(2, ceil(sleeper_spacing / wavelength));
switch length_option
    case 'exact'
        actual_length = nominal_length;
    case 'ends at sleepers'
        no_sleepers = no_sleepers + 1;
        actual_length = (no_sleepers - 1) * sleeper_spacing;
    case 'ends at midspans'
        actual_length = no_sleepers * sleeper_spacing;
end
offset = (actual_length - sleeper_spacing * (no_sleepers - 1)) / 2;

%define the key nodes
nodes = [0: no_sleepers - 1]' * sleeper_spacing + offset;
sleeper_nodes = ones(size(nodes));
if offset
    nodes = [0; nodes; actual_length];
    sleeper_nodes = [0; sleeper_nodes; 0];
end
if ~isempty(forcing_position)
    if min(abs(nodes - forcing_position)) > (sleeper_spacing / 100)
        nodes = [nodes; forcing_position];
        sleeper_nodes = [sleeper_nodes; 0];
    end
end

%add extra nodes to satisfy max node spacing
for ii = 1:length(nodes)-1
    gap = nodes(ii+1) - nodes(ii);
    if gap > max_node_spacing
        tmp = linspace(nodes(ii), nodes(ii + 1), ceil((nodes(ii+1) - nodes(ii)) / max_node_spacing) + 1)';
        tmp = tmp(2:end-1);
        nodes = [nodes; tmp];
        sleeper_nodes = [sleeper_nodes; zeros(size(tmp))];
    end
end

%sort nodes into order
[nodes, ii] = sort(nodes);
sleeper_nodes = sleeper_nodes(ii);


elements = [1:length(nodes) - 1; 2:length(nodes)]';
sleeper_nodes = find(sleeper_nodes);

if ~isempty(forcing_position)
    [~, forcing_node] = min(abs(forcing_position - nodes));
else
    forcing_node = [];
end

end



