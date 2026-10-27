function K = mode_freqs(gridsz, vox)
%MODE_FREQS  Angular frequency magnitude |k| at every Fourier mode of a
%     grid with the given size and voxel spacing.
kk = cell(1,3);
for i = 1:3
    n = gridsz(i);
    f = [0:floor(n/2), -(ceil(n/2)-1):-1] / (n * vox);
    kk{i} = 2 * pi * f;
end
[K1,K2,K3] = ndgrid(kk{1}, kk{2}, kk{3});
K = sqrt(K1.^2 + K2.^2 + K3.^2);
end
