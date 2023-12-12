% -------------------------------------------------------------------------
% Script to generate nice plots comparing DJL and 0th/1st order
% perturbation
% -------------------------------------------------------------------------

UU=linspace(DJL.C,DJL.C*1.2,1000);

delta=(U-DJL.C)/epsilon;
delta0=(UU-DJL.C)/epsilon;
delta1=(UU-DJL.C)/epsilon;
deltastar=delta1/(DJL.s*DJL.mu^2);

B=deltastar/2.*sech(sqrt(deltastar)/2.*domain.x{1}/pi*DJL.KAI).^2;
A=6*DJL.s*DJL.mu^2/DJL.r*B;

% -------------------------------------------------------------------------
% 0th order DJL approx
% -------------------------------------------------------------------------
A=reshape(A,[size(A,1),1,size(A,2)]);

% Find zeta from fkdv
zeta_0=pagemtimes(A,DJL.phi');

% -------------------------------------------------------------------------
% 1st order DJL approx
% -------------------------------------------------------------------------
% A_xx in x domain (KAI/mu)
A_xx=ifft(-(pi/(1/mu*KAI)*domain.k{1}).^2.*fft(A));

% an
an=pagemtimes(A_xx,DJL.a1)+pagemtimes(A.^2,DJL.a2);

% n=N case
an(:,mode,:)=0;

zeta_1=pagemtimes(an,DJL.phis');

% Perturbation solutions
zeta0=epsilon*zeta_0;
zeta1=epsilon*(zeta_0+epsilon*zeta_1);

% -------------------------------------------------------------------------
% Check overturning
% -------------------------------------------------------------------------
diffv_zeta0=zeros(length(UU),1);
diffv_zeta0(1:end)=max(max(2*real(ifct(chebdiff(fct(permute(zeta0(:,:,1:end),[2 1 3])),1)))));

diffv_zeta1=zeros(length(UU),1);
diffv_zeta1(1:end)=max(max(2*real(ifct(chebdiff(fct(permute(zeta1(:,:,1:end),[2 1 3])),1)))));

% make sure its monotonically increasing
if min(diff(diffv_zeta0))<0
    fprintf('Not monotonically increasing!\n');
end

if min(diff(diffv_zeta1))<0
    fprintf('Not monotonically increasing!\n');
end
% -------------------------------------------------------------------------
% Filter out overturning solutions
% -------------------------------------------------------------------------
% find length
i0=length(diffv_zeta0(diffv_zeta0<1));
i1=length(diffv_zeta1(diffv_zeta1<1));

% Filter to valid solution + first invalid solution
delta0=delta0(1:i0+1);
zeta0=zeta0(:,:,1:i0+1);

delta1=delta1(1:i1+1);
zeta1=zeta1(:,:,1:i1+1);

% Calculate momentum 0th order
P0=trapI((epsilon*zeta0).^2,domain.dx{1}); % Integrate x
P0=permute(P0,[2,1,3]);
P0=clenshaw_curtis(2*P0/pi*KAI/mu^2); % Integrate y
P0=permute(P0,[3,1,2]);

% Calculate momentum 1st order
P1=trapI(zeta1.^2,domain.dx{1}); % Integrate x
P1=permute(P1,[2,1,3]);
P1=clenshaw_curtis(2*P1/pi*KAI/mu^2); % Integrate y
P1=permute(P1,[3,1,2]);

figure;
plot(delta0,P0,delta1,P1,delta(1:end-1),P(1:end-1),delta(end),P(end),'o')
xlabel('\Delta')
ylabel('Momentum')
legend('0th order','1st order','DJL','Location','NorthWest')

figure;
hold on
plot(delta0(1:end-1),P0(1:end-1),'b-',delta0(end),P0(end),'bo',delta0(end-1:end),P0(end-1:end),'b:','LineWidth',1)
plot(delta1(1:end-1),P1(1:end-1),'r-',delta1(end),P1(end),'ro',delta1(end-1:end),P1(end-1:end),'r:','LineWidth',1)
plot(delta(1:end-1),P(1:end-1),'k-',delta(end),P(end),'ko',delta(end-1:end),P(end-1:end),'k:','LineWidth',1)
xlabel('\Delta');
ylabel('Momentum')
legend('0th order WNL solutions','0th order WNL overturning point','','1st order WNL solutions','1st order WNL overturning point','','Valid DJL solutions','First DJL invalid overturning solution','','Location','NorthWest')





