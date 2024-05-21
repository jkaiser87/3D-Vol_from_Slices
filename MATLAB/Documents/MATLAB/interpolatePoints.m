function new_points = interpolatePoints(points, factor)
    % Points is an Nx3 matrix where columns are x, y, z coordinates
    % Factor is the number of intermediate points to add between each original point

    if size(points, 1) < 2
        error('Not enough points to interpolate.');
    end

    % Pre-allocate space for new points
    new_points = zeros((size(points, 1) - 1) * (factor + 1) + 1, 3);
    index = 1;

    % Interpolate between each pair of points
    for i = 1:size(points, 1)-1
        for j = 0:factor
            t = j / (factor + 1);
            new_points(index, :) = (1 - t) * points(i, :) + t * points(i + 1, :);
            index = index + 1;
        end
    end
    new_points(end, :) = points(end, :); % Add the last point

    return;
end