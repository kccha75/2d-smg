M1=linspace(1,4,4)+9;
M2=linspace(1,4,4)+9;

figure('Position',[300 300 600 300]); fsz=15; lw=2;

subplot(1,2,1)
semilogy(M1,t_cg511,'-x',M1,t_mg511(:,2),'-o',M1,t_mg511(:,1),'-o')
xlabel('$N$','interpreter','latex','fontsize',fsz)
ylabel('$t$','interpreter','latex','fontsize',fsz)
legend('CG','SMG coarse grid N=8','SMG coarse grid N=9','Location','SouthEast')
xticks([10 11 12 13])

subplot(1,2,2)
semilogy(M2,t_cg512,'-x',M2,t_mg512(:,3),'-o',M2,t_mg512(:,2),'-o',M2,t_mg512(:,1),'-o')
xlabel('$N$','interpreter','latex','fontsize',fsz)
ylabel('$t$','interpreter','latex','fontsize',fsz)
legend('CG','SMG coarse grid N=7','SMG coarse grid N=8','SMG coarse grid N=9','Location','SouthEast')
xticks([10 11 12 13])