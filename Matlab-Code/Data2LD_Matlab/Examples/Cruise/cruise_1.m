function DSvec = cruise_1(t, Svec, SetPtfd)
DSvec    = zeros(2,1);
Uvec     = eval_fd(t, SetPtfd);
DSvec(1) =  -Svec(1) + Svec(2)/4;
DSvec(2) =  (Uvec - Svec(1));



