function dkern = gradkernCompute(GPprior, T, X)
 
dkern = zeros(GPprior.nParams,1);
switch GPprior.type 
 case 'se'
    T = GPprior.K.*T; 
    dkern(1) = sum(sum(T.*GPprior.XX))/2;
    dkern(2) = sum(sum(T));
end  