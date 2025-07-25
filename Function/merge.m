function img = merge(blocks, params,overlap_sz)

sz = params.sz;
block_num = params.block_num;
block_sz = params.block_sz;
img = zeros(sz(1),  sz(2), sz(3));
mult = zeros(sz(1), sz(2), sz(3));

for i = 1:block_num(1)
    for j = 1:block_num(2)
        ii = 1 + (i - 1)*(block_sz(1)-overlap_sz(1));
        jj = 1 + (j - 1)*(block_sz(2)-overlap_sz(2));
        idx = (j-1)*block_num(1) + i;
        mult(ii:ii+block_sz(1)-1, jj:jj+block_sz(2)-1, :) = mult(ii:ii+block_sz(1)-1, jj:jj+block_sz(2)-1, :) + 1;
        img(ii:ii+block_sz(1)-1, jj:jj+block_sz(2)-1, :) = img(ii:ii+block_sz(1)-1, jj:jj+block_sz(2)-1, :) +blocks(:, :, :, idx);
    end
end

img = img./mult;