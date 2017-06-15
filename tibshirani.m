function x= tibshirani(H,f,t,maxIter)
% X=TIBSHIRANI(H,F,T,MAXITER)
% perform L^1 constrained reciprocity using Tibshirani's (1996) method for
% L^1 constrained least squares
%
% H: quadratic term of cost function (see Dmochowski et al. (2017),
% NeuroImage)
% f: linear term of cost function
% t: L^1 constraint
% maxIter: maximum number of iterations in the search for solution
%
% x: L^1 constrained reciprocal tDCS currents
%
% Reference: Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. 
% Journal of the Royal Statistical Society. Series B (Methodological), 267-288.

f=f(:);
x=-H\f;
A=[]; b=[]; % initial constraints are empty
Aeq=ones(1,numel(x)); beq=0;
% run the loop
for i=1:maxIter
    if norm([x;-sum(x)],1)<t  % this is the only place where we add reference electrode
        break;
    else
        thisA=( 2*(x>0)-1 ).';
        thisb=t;
        A=cat(1,A,thisA);
        b=cat(1,b,thisb);
        x=quadprog(H,f,A,b,Aeq,beq);
    end
end
i % display iteration count at convergence
end
