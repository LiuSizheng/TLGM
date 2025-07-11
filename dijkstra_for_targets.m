function shortest_distances = dijkstra_for_targets(adj_matrix, target_nodes)
    N = size(adj_matrix, 1);
    shortest_distances = inf(length(target_nodes), N);
    
    for t = 1:length(target_nodes)
        source = target_nodes(t);
        distance = inf(1, N);
        distance(source) = 0;
        unvisited = 1:N;
        
        while ~isempty(unvisited)
            [~, idx] = min(distance(unvisited));
            current = unvisited(idx);
            unvisited(idx) = [];
            
            for neighbor = find(adj_matrix(current, :) > 0)
                alt = distance(current) + adj_matrix(current, neighbor);
                if alt < distance(neighbor)
                    distance(neighbor) = alt;
                end
            end
        end
        shortest_distances(t, :) = distance;
    end
end