function res = diff_large(Org, fusion,x,y)

        diff_image = abs(double(Org) - double(fusion));
        color_map = jet(256); 
        gray_diff_image = mean(diff_image, 3); 
        scaled_diff_image = mat2gray(gray_diff_image, [min(gray_diff_image(:)) max(gray_diff_image(:))]); % 将差值图像缩放到[0,1]范围
        colored_diff_image = ind2rgb(uint8(scaled_diff_image * 255), color_map);
        
        res = enlarge(colored_diff_image,x,y);
end



