function Dbval = fun_Dexplinear(t, bvec, Bbasisobj)
nbasis    = length(bvec);
if isa_fdPar(Bbasisobj)
    Bbasisobj = getbasis(Bbasisobj);
end
basismat  = eval_basis(t, Bbasisobj);
bval      = exp(basismat*bvec);
Dbval     = basismat.*repmat(bval,1,nbasis);
end