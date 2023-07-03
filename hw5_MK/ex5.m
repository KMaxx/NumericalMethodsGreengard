for jtest=1:10
    m = jtest*10-1;
    [h,k,error]=heat_FE(m);
    hvals(jtest)=h;
    E(jtest)=error;
end 

error_table(hvals, E);   % print tables of errors and ratios
error_loglog(hvals, E);  % produce log-log plot of errors and least squares fit
