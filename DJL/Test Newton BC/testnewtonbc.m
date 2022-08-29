% simple newton bc test 

clear;clc
N=64;
k = (0:N-1)';
x = cos(pi*k/(N(1)-1));
dx = x(2:end)-x(1:end-1);

a=1;
b=0;
c=-1;
d=3;

f=exp(-2*x).*(exp(x).*(2-4*x)+3.*(-1+x.^2).^2);

% Newton ...
e0=zeros(size(N));
v=(x.^2-1).*exp(-x)+0.2;

ue=(x.^2-1).*exp(-x);


% set bc here ...
v(1)=ue(1);v(end)=ue(end);

for i=1:10

    % Jacobian
    Ja=a;
    Jb=b;
    Jc=c+2*v*d;
    Jf=f-(a.*ifct(chebdiff(fct(v),2))+b.*ifct(chebdiff(fct(v),1))+c.*v+d.*v.^2);
  
    r=rms(Jf);
    fprintf('Residual Newton = %d\n',r)

    if r<=1e-12
        fprintf('Converged after %d Newton Iterations \n',i-1)
        break
    end

    % solve
    A=Ja.*ifct(chebdiff(fct(eye(N,N)),2))+ ...
        Jb.*ifct(chebdiff(fct(eye(N,N)),1))+ ...
        Jc.*speye(N);

    A(1,:)=0;A(1,1)=1;
    A(end,:)=0;A(end,end)=1;

    Jf(1)=0;
    Jf(end)=0;

    e=A\Jf;

    v=v+e;

end