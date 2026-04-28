function Value = NormPDF(Z)
    Value = 1/sqrt(2 * pi) * exp( - Z.^2 / 2);
end