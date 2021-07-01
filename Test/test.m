% Testing dealias effects ... and order of operation

v1=sin(exp(-X.^2-Y.^2));

v2=v1.^2;
v2_dealias=dealias_2d(v1,v1);
surf(v2-v2_dealias)

% Coarse grid 

% restricted dealias
v2_dealias_c=fourier_restrict_2d_filtered(dealias_2d(v1,v1));

% dealiased restricted?
v2_c_dealias=dealias_2d(fourier_restrict_2d_filtered(v1),fourier_restrict_2d_filtered(v1));

z=fourier_restrict_2d_filtered(v1).^2;

figure;
surf(v2_dealias_c-v2_c_dealias)
% Quite different ...