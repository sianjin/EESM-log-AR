function logSGNparam = sgnFitting(stdr)
n = length(stdr);
M = mean(stdr);
S = std(stdr);
logSGNparam = zeros(1,4);
fBest = -10e10;
randRange1 = 1;
randRange2 = 0.5;
logSGNparam0 = [M,S, -randRange1 + (2*randRange1)*rand, (2*randRange2)*rand];
A = [0, 0, 0, -1; 0, -1, 0, 0];
b = [0, 0];
for iter = 1: 100
    LikelihoodFun = @(x) -(-n/2* log(pi * x(2)^2/2) - 1/2/x(2)^2*sum((stdr-x(1)).^2) + sum(log(normcdf(x(3).*(stdr-x(1))./sqrt(x(2)^2+ x(4).*(stdr-x(1)).^2)))));
    [x, fval] = fmincon(LikelihoodFun,logSGNparam0,A,b);
    if fval>fBest
        fBest = fval;
        logSGNparam = x;
        logSGNparam0 = logSGNparam;
    end
end
end