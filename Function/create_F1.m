function [U] = create_F1()
    F=load('R.mat');
    F=F.R;
    U=F(:,1:end-10);
    for band = 1:size(U,1)
        div = sum(U(band,:));
        for i = 1:size(U,2)
            U(band,i) = U(band,i)/div;
        end
    end
end