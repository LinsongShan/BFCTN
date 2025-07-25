function log_relative_error = computeRelativeError(prev_result, curr_result)
% 计算两个多维数组之间的相对误差
% 输入：
%   prev_result - 上一次的结果（多维数组）
%   curr_result - 这一次的结果（多维数组）
% 输出：
%   relative_error - 相对误差值

% 检查输入维度是否相同
if ~isequal(size(prev_result), size(curr_result))
    error('输入数组的维度必须相同');
end

% 计算Frobenius范数
norm_prev = norm(prev_result(:), 'fro');
norm_diff = norm(curr_result(:) - prev_result(:), 'fro');

% 计算相对误差
relative_error = norm_diff / norm_prev;


log_relative_error = log10(relative_error); % 常用对数(10为底)
    
% 显示结果（可选）
fprintf('对数相对误差: %.4f\n', log_relative_error);

end