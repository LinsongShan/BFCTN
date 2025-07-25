function relative_error_map = calculateRelativeError3D(y_true, y_pred)
      % 输入检查
    assert(all(size(y_true) == size(y_pred)), 'Data size mismatch!');
    assert(~any(isnan(y_true(:))) && ~any(isnan(y_pred(:))), 'NaN in data!');
    assert(~any(isinf(y_true(:))) && ~any(isinf(y_pred(:))), 'Inf in data!');
    
    % 动态epsilon（基于数据中位数）
    nonzero_values = y_true(y_true ~= 0);
    if isempty(nonzero_values)
        epsilon = 1e-10;
    else
        epsilon = 1e-6 * median(abs(nonzero_values));
    end
    
    % 对称相对误差
    denominator = abs(y_true) + abs(y_pred) + epsilon;
    relative_error = abs(y_pred - y_true) ./ denominator;
    relative_error_map = mean(relative_error, 'all');
end

