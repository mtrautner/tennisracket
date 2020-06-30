function v_out = beam_ball_MOL_func(ni)


%% beam + ball equation by the method of lines
% ut = uxx
% the boundary condition are: u(x=a) = ua , and u(x=b) = ub
% the initial condition is u(x, t=0) = u0(x)

% clear all; close all; clc;

% nomenclature:
% dx; % spatial discretisation 
% n;  % number of partition of the spatial interval

phi = zeros(1,10^5);
ct = 1;


n = 2^6;

L = 1.1; 
% Geometry
a = 0; b = 1;
h = 0.010;
w = 0.032;
dx = (b-a)/n;
x = a:dx:b; x = x';
EI = (w*(h^3)/12)*7*(10^10);
m = (h*w)*2700; %perhaps also times dx? 
kb = 2*10e4;
mb = 1;
%ni = 2^6-2^4;
br = 0.0325/L; %ball radius
% Initial Condition
v1 = -10; %Incident velocity of the ball
p1 = br; %Initial position of the ball
uinit = u0(x, v1, n, p1);

tspan = 0:0.0001:0.02;

[t, u] = ode15s(@(t,u) BBMOL_eqn(t,u,dx,n, EI, m, kb, mb, ni, br, L), tspan, uinit);

Z = u(:,n+5:end-3);
% surf(x,t,Z, 'EdgeColor', 'None')
% xlabel('x')
% ylabel('t')
% %% Plot results
% 
% figure(1)
% % 
% % obj = VideoWriter('Tennis_Ball_impact_2clamp.avi');
% % obj.Quality= 100;
% % obj.FrameRate = 20;
% % 
% % open(obj);
% 
% for ti = 1:length(tspan)
%     plot(x,Z(ti,:),'LineWidth',3,'Color','b');
%     hold on
%     scatter(dx*ni,u(ti,end),'r','LineWidth',3);
%     scatter(dx*ni,u(ti,end)- min([br, u(ti,end)-Z(ti,ni)]),'g','LineWidth',3);
%     xlim([0,b]);
%     ylim([-0.05,0.05]);
%     xlabel('Tennis racket length')
%     ylabel('Transverse Displacement')
%     title(sprintf('Tennis Racket on impact, t = %f', ti))
% %     sln = zeros(1,length(x));
% %     for n = 1:5
% %         sln = sln + y(t,x,n);
% %     end
% %     plot(x,sln,'LineWidth',3)
% %     xlim([0,1.2])
% %     ylim([-4,4])
% %     xlabel('x')
% %     ylabel('Transverse Displacement')
% %     title('Beam Vibrations: IC: f(x) = 0.25x^3, g(x) = -x')
%     hold off
%     pause(0.001);
% % %     pause(0.1)
% %     f = getframe(gcf);
% %     writeVideo(obj,f);
% end
% % 
% % obj.close()

v_out = u(end,n+4);


end
