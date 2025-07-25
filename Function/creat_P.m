function P  = creat_P(patchsize, space_ratio)
        P = zeros(patchsize/space_ratio, patchsize);
        for i = 1 : patchsize/space_ratio
                P(i, space_ratio*(i - 1)  + 1 : space_ratio*(i - 1)  + space_ratio) = 1/space_ratio;
        end
end