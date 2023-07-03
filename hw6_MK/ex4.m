M = [49 99 199]
for jtest=1:3
    m = M(jtest);
    %[h,k,error]=advection_LW_pbc(m);
    [h,k,error]=advection_up_pbc(m);
    hvals(jtest)=h;
    E(jtest)=error;
end 

error_table(hvals, E);   % print tables of errors and ratios
error_loglog(hvals, E);  % produce log-log plot of errors and least squares fit
