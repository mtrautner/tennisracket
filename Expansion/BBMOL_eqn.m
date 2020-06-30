function [ut] = BBMOL_eqn(t,u,dx,n, EI, m, kb,mb,ni, br, L)
% this gives the rhs of the MOL 
% Since the equations are second order, all variables of x have a
% corresponding pair
ut = zeros(2*(n+4),1);

u_pos = u(n+5:end-1);
b_pos = u(end);

% u_pos = [u_pos; 0; 0];
u_vel = u(1:n+3);
b_vel = u(n+4);
ni_vel = u(ni);
ni_pos = u(ni+n+4);

rho = EI*dx/(m*L^3);
uxx = zeros(n+3,1);

uxx(3:n+1) = -rho*(u_pos(1:n-1) -4*u_pos(2:n) + 6*u_pos(3:n+1)-4*u_pos(4:n+2) + u_pos(5:n+3))/(dx^4) ;

ut(3:n+1) = uxx(3:n+1); %Populate second derivatives
ut(n+5:(end-1)) = u_vel; %Populate first derivatives
ut(end) = b_vel; %Ball velocity

ut(n+5) = 0; %velocity of clamped end is 0 

ut(2) = -rho/(dx^4)*(-4*u_pos(1) + 6*u_pos(2) -4*u_pos(3) + u_pos(4)); %second derivative of close to end needs separate claculation

% Boundary condition for clamped second end
% ut(n+4 + n+2) = 0; 
% ut(n+4 + n+3) = 0; 

% %Boundary conditions for free end: set velocities of n+2 and n+3
ut(n+4 + n+2 ) = 2*u_vel(n+1) - u_vel(n);
ut(n+4 + n+3 ) = 4*u_vel(n+1) -4*u_vel(n) +u_vel(n-1);

% Ball is affected by racket force only if hasn't left racket contact yet
% Assume contact lost after both pass through center again
if (u_pos(ni)== 0) || (b_pos < br)
    ut(n+4) = kb/mb*(br - (b_pos-u_pos(ni))); %ball acceleration
    ut(ni) = ut(ni) - kb/(m*L)*(br - (b_pos-u_pos(ni)));
    t
else
    ut(n+4) = 0; %ball acceleration is 0 
end

end