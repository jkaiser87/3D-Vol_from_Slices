function smoothedVertices = laplacianSmooth(vertices, faces, lambda, iterations)
    % vertices: Nx3 matrix of vertex coordinates
    % faces: Mx3 matrix of indices into vertices
    % lambda: Smoothing factor, typical values are in the range 0.5 - 1
    % iterations: Number of times the smoothing operation is applied

    smoothedVertices = vertices;
    for iter = 1:iterations
        for i = 1:size(vertices, 1)
            % Find all faces that include this vertex
            [row, ~] = find(faces == i);
            % Get unique vertices connected to the current vertex
            neighborIdx = unique(faces(row, :));
            neighborIdx(neighborIdx == i) = [];  % Remove the vertex itself

            % Calculate the mean position of neighboring vertices
            meanPos = mean(smoothedVertices(neighborIdx, :), 1);

            % Update the vertex position
            smoothedVertices(i, :) = smoothedVertices(i, :) + lambda * (meanPos - smoothedVertices(i, :));
        end
    end
end