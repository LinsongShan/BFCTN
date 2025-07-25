function MSI = generate_HRMSI(orgTensor,F)
    MSI = ttm(tensor(orgTensor), F, 3);
    MSI = double(MSI);
end