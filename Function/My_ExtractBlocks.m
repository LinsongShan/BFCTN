function [blocks] = My_ExtractBlocks(img, params, overlap_sz)

sz = size(img);
block_sz = params.block_sz;
block_num = floor((sz(1:2) - overlap_sz)./(block_sz - overlap_sz));
blocks = zeros([block_sz, sz(3), prod(block_num)]);

for i = 1:block_num(1)
    for j = 1:block_num(2)
        ii = 1 + (i - 1)*(block_sz(1) - overlap_sz(1));
        jj = 1 + (j - 1)*(block_sz(2) - overlap_sz(2));
        idx = (j-1)*block_num(1) + i;
        blocks(:, :, :, idx) = ...
            img(ii:ii+block_sz(1)-1, jj:jj+block_sz(2)-1, :);
    end
end


