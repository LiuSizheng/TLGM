function neighbor_of_i = neighbor_of_node_set( mixedsig, node_set )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
neighbor_of_i=[];
for i=1:length(node_set)
    neighbor_of_i=[neighbor_of_i,find(mixedsig(node_set(i),:))];
end
neighbor_of_i=unique(neighbor_of_i);
neighbor_of_i=setdiff(neighbor_of_i,node_set);
end

