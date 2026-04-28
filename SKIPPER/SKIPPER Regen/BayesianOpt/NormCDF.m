function Value = NormCDF(Z)
    Value = 1/2 * (1 + erf(Z / sqrt(2)));
end