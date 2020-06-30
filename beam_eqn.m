%Define constants
L = 1;
c1 = 1;
gamma = 3;

%Initial condition
f_t0 = @(x) 0.25*x.^3;
fp_t0 = @(x) -x;

% %% Solve for betas (resonance frequencies)
% betas = [];
% for i = 1:15
%     syms x
%     S = vpasolve(cos(x)*cosh(x) == -1, x,[i i+1]);
%     if ~isempty(S)
%         betas = [betas, double(S)];
%     end
% end
% 
% rhos = 1/L.*betas;

%% Evaluate v_n(x) over x = [0 1] for each n
x = linspace(0,1,100);
Vs = zeros(length(rhos), length(x));
for i = 1:length(rhos)
    rho = rhos(i);
    v_n = c1.*(cosh(rho*x)-cos(rho*x) - (cos(rho*L)+cosh(rho*L))/(sin(rho*L)+sinh(rho*L))*(sinh(rho*x)-sin(rho*x)));
    Vs(i,:) = v_n;
    plot(x,v_n,'LineWidth',2);
    hold on
end

xlabel('x: Beam Length = 1')
ylabel('Transverse Displacement')
title('Normal modes of vibration for a uniform beam')
legend('Mode 1: 1.87','Mode 2: 4.69','Mode 3: 7.85', 'Mode 4: 11.00', 'Mode 5: 14.14')
%% Solve for Fourier coefficients a_n
a_ns = zeros(1,length(rhos));
b_ns = zeros(1,length(rhos));

g = @(x, rho) f_t0(x).*c1.*(cosh(rho*x)-cos(rho*x) - (cos(rho*L)+cosh(rho*L))/(sin(rho*L)+sinh(rho*L))*(sinh(rho*x)-sin(rho*x)));
g2 = @(x, rho) fp_t0(x).*c1.*(cosh(rho*x)-cos(rho*x) - (cos(rho*L)+cosh(rho*L))/(sin(rho*L)+sinh(rho*L))*(sinh(rho*x)-sin(rho*x)));

for i = 1:length(rhos)
    rho = rhos(i);
    ub = pi/sqrt(gamma*(rho.^2));
    q = integral(@(x) g(x,rho),0,ub);
    p = integral(@(x) g2(x,rho),0,ub);
    a_ns(i) = 2*q/ub;
    b_ns(i) = 2*p/(ub*sqrt(gamma)*rho.^2);
end

%%
lambdas_sq = rhos.^(2);
t = linspace(0,1);
y = @(t,x,n) (a_ns(n)*cos(sqrt(gamma)*lambdas_sq(n).*t) + b_ns(n)*sin(sqrt(gamma)*lambdas_sq(n).*t))*c1*(cosh(rhos(n)*x)-cos(rhos(n)*x) - (cos(rhos(n)*L)+cosh(rhos(n)*L))/(sin(rhos(n)*L)+sinh(rhos(n)*L))*(sinh(rhos(n)*x)-sin(rhos(n)*x)));

figure(1)
% 
obj = VideoWriter('Clamped_beam.avi');
obj.Quality= 100;
obj.FrameRate = 20;

open(obj);

for t = linspace(0,20,301)
    sln = zeros(1,length(x));
    for n = 1:5
        sln = sln + y(t,x,n);
    end
    plot(x,sln,'LineWidth',3)
    xlim([0,1.2])
    ylim([-4,4])
    xlabel('x')
    ylabel('Transverse Displacement')
    title('Beam Vibrations: IC: f(x) = 0.25x^3, g(x) = -x')
    hold off
    pause( 0.01 );
   
%     pause(0.1)
    f = getframe(gcf);
    writeVideo(obj,f);
end

obj.close()



% xx = linspace(0,1);
% plot(xx,cos(xx).*cosh(xx))