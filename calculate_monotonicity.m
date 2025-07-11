function M = calculate_monotonicity(ranking)
    % CALCULATE_MONOTONICITY Computes the monotonicity (M) of a ranking list
    % Input:
    %   ranking - Vector of ranks assigned to nodes (e.g., [1, 2, 2, 3, 4, 5, 5])
    % Output:
    %   M - Monotonicity value (between 0 and 1)
    
    % Total number of nodes
    n = length(ranking);
    
    % Get unique ranks and count nodes per rank
    unique_ranks = unique(ranking);
    n_r = arrayfun(@(r) sum(ranking == r), unique_ranks); % Number of nodes per rank
    
    % Calculate the numerator: sum(n_r * (n_r - 1))
    numerator = sum(n_r .* (n_r - 1));
    
    % Calculate the denominator: n * (n - 1)
    denominator = n * (n - 1);
    
    % Compute monotonicity M
    if denominator == 0  % Edge case: if n <= 1
        M = 1;  % Undefined case, assume perfect monotonicity
    else
        M = (1 - numerator / denominator)^2;
    end
end