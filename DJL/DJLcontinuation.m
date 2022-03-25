du=0.005;

global u
U(1)=u;
V(1)=max(abs(v(:)));

for ii=1:100
    fprintf('Step %d\n',ii)
    u=u+du;
    
    v0=v;
    DJLsolve

    if converge==false
        fprintf('OH NOOOO')
        break
    end

    U(ii+1)=u;
    V(ii+1)=max(abs(v(:)));

end

plot(U,V)
xlabel('u');ylabel('Amplitude');title('mode 1 DJL')