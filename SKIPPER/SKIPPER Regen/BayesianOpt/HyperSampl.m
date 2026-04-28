function LHS = HyperSampl(NumSamples, NumDimensions)
% Latin Hypercube Sampling for optimal coverage for the GP training 
LHS = zeros(NumSamples, NumDimensions);
for j = 1:NumDimensions
    P_ij = randperm(NumSamples)' - 1;

    U_j = rand(NumSamples, 1);

    LHS(:, j) = (P_ij + U_j) / NumSamples;
end
end