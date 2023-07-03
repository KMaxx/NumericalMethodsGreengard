CK = [0.4 0.2 0.1]
for jtest=1:3
    ck = CK(jtest);
    [h,k,error]=advection_up_pbc2(ck);
    hvals(jtest)=h;
    E(jtest)=error;
end 

error_table(hvals, E);   % print tables of errors and ratios
error_loglog(hvals, E);  % produce log-log plot of errors and least squares fit
