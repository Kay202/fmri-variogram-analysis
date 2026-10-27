function o = om3(x)
%OM3  Omega_3(x) = sin(x)/x, the isotropic basis function in three
%     dimensions, with the removable singularity at zero filled in.
o = ones(size(x)); nz = x ~= 0; o(nz) = sin(x(nz))./x(nz);
end

