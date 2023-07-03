om = 0.1:0.1:2;
for i=1:length(om)
    omega = om(i);
    run iter_bvp_Asplit;
    gw(i) = rhoG;
end
figure(3)
plot(om,gw)
title('g(omega) for omega in [0,2]','FontSize',15)
xlabel('omega')
ylabel('g(omega)')
