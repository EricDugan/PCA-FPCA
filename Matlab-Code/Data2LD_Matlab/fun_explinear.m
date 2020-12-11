function bval = fun_explinear(t, bvec, Bbasisobj)
if isa_fdPar(Bbasisobj)
    Bbasisobj = getbasis(Bbasisobj);
end
basismat = eval_basis(t,Bbasisobj);
bval     = exp(basismat*bvec);

