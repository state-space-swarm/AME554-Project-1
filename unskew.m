function vec = unskew(matrix)
    % Check if the input is a square matrix
    [rows, cols] = size(matrix);
    if rows ~= cols
        error('Input must be a square matrix');
    end
    
    % Initialize the vector
    vec = [-matrix(2,3); matrix(1,3); -matrix(1,2)];
end