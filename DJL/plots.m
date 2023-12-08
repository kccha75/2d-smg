% Script to generate nice plots comparing DJL and 0th/1st order
% perturbation

delta=(U-DJL.C)/epsilon;
deltastar=delta/(DJL.s*DJL.mu^2);

B=deltastar/2.*sech(sqrt(deltastar)/2.*domain.x{1}/pi*DJL.KAI).^2;
A=6*DJL.s*DJL.mu^2/DJL.r*B;


% -------------------------------------------------------------------------
% 0th order DJL approx
% -------------------------------------------------------------------------
A=reshape(A,[size(A,1),1,size(A,2)]);

% Find zeta from fkdv
zeta0=pagemtimes(A,DJL.phi');

% -------------------------------------------------------------------------
% 1st order DJL approx
% -------------------------------------------------------------------------
% A_xx in x domain (KAI/mu)
A_xx=ifft(-(pi/(1/mu*KAI)*domain.k{1}).^2.*fft(A));

% an
an=pagemtimes(A_xx,DJL.a1)+pagemtimes(A.^2,DJL.a2);

% n=N case
an(:,mode,:)=0;

zeta1=pagemtimes(an,DJL.phis');

zeta=epsilon*(zeta0+epsilon*zeta1);

% Calculate momentum 0th order
P0=trapI((epsilon*zeta0).^2,domain.dx{1}); % Integrate x
P0=permute(P0,[2,1,3]);
P0=clenshaw_curtis(2*P0/pi*KAI/mu^2); % Integrate y
P0=permute(P0,[3,1,2]);

% Calculate momentum 1st order
P1=trapI(zeta.^2,domain.dx{1}); % Integrate x
P1=permute(P1,[2,1,3]);
P1=clenshaw_curtis(2*P1/pi*KAI/mu^2); % Integrate y
P1=permute(P1,[3,1,2]);

figure;
plot(delta,P0,delta,P1,delta,P)
xlabel('\Delta')
ylabel('Momentum')
legend('0th order','1st order','DJL')

% check overturning
diffv=2*ifct(chebdiff(fct(zeta(:,:,end)'),1));
if max(diffv(:))>1
    fprintf('Overturning detected!\n')
    return
end