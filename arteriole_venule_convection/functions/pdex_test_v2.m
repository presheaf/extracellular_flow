function u = pdex_test(r0,dmax,D,x,t)

m = 1; % cylinder symmetry in the diffusion equation
sol = pdepe(m,@(x,t,u,DuDx)pdex1pde(x,t,u,DuDx,D),@pdex1ic,@pdex1bc,x,t); 
% pdex1pde: assigns parameters to differential equation, se matlabs pdepe
% pdex1ic: assign initial condition (concentration = 1 all over)
% pdex1bc: assign boundary conditions (u=0 at x=r0, u=1 at x=xmax)

% Extract the first solution component as u.
u = sol(:,:,1);

% --------------------------------------------------------------
function [c,f,s] = pdex1pde(x,t,u,DuDx,D)
c = 1/D;
f = DuDx;
s = 0;
% --------------------------------------------------------------
function u0 = pdex1ic(x)
r0 = 30e-6;
%u0 = heaviside(x-r0); % initial condition: c = zero within paravascular space, 1 within interstitial space
u0 = ones(1,length(x));
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
pl = ul;
ql = 0;
pr = ur-1;
qr = 0;