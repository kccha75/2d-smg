clear;close all;%clc

global u
u=0.235;

DJLinitialise

DJLv0

DJL_pde_initialise

DJLsolve

X2=X/pi*KAI*L/mu;
Y2=(Y+1)/2;

DJLcontinuation

% 
% surf(X2,Y2,v0);title('initial guess')
% figure;surf(X2,Y2,v);title('DJL solution');xlabel('x');ylabel('z')
% figure;surf(X2,Y2,v0-v);title('initial guess error');xlabel('x');ylabel('z')
% 
fprintf('Domain is %d\n',KAI*L/mu)