function values = kinematic_path_values(alpha)
% Monthly state path -> x, v, a, rolling ten-year average acceleration.
% Columns use degrees C, degrees C/decade, and degrees C/decade squared.
values = [alpha(1,:)', 10*alpha(2,:)', 100*alpha(3,:)', nan(size(alpha,2),1)];
if size(values,1) > 120
    % The same joint path supplies both endpoints; elapsed time is one decade.
    values(121:end,4) = values(121:end,2) - values(1:end-120,2);
end
end
