function I = reciprocate(V,head,Imax)
%I=RECIPROCATE(R,V,IMAX) determine the TES montage that targets the sources expressed in
%scalp potentials V, generated as a linear combination of the columns of R,
%and constrained to a maximum delivered current of Imax
%   R: channels x (locs*3) lead-field matrix
%   V: channels x 1 vector of scalp potentials
%   Imax: total current delivered in A (defaults to 0.002) 

if nargin<3, Imax=0.002; end

RR=head.R*(head.R)';

rankRR=rank(RR);
if rankRR<size(head.R,1),
    warning('R is rank-deficient.  Results will likely be wrong');
    % try a regularized inverse
    warning('Trying to compute a regularized inverse...');
    I = regInv( RR,100 ) * V;
else
    I=RR\V;  % R is well-conditioned
end


Itotal=sum(abs([I; sum(I)]));
I=I/(Itotal/(2*Imax));


end

