

%% Starting with a rectangular grid
c =5;
xd= 41;
yd= 41;
x = linspace(0,1,xd);
y = linspace(0,1,yd);
[X,Y] = meshgrid(x,y);

n = xd*yd;
U = zeros(yd,xd);

%Initial Condition is Gaussian
U0 = 0.01*reshape(mvnpdf([X(:),Y(:)],[0.6,0.5],0.01*eye(2)),[yd,xd]);
U = U0;
% surf(X,Y,U0)


%Perform triangulation
nodes = [X(:),Y(:)]; % n rows with two dimensions (cols)
N = 2*(xd-1)*(yd-1);
elems = zeros(N,3); %N rows of elements with 3 indices in orientation order

%Make all the triangles
count = 1;
for xi = 1:(xd-1)
    for yi = 1:(yd-1)
        elems(count,:) = [yi + (xi-1)*yd, yi + xi*yd, yi + 1 + (xi-1)*yd];
        elems(count+1,:) = [yi+1+xi*yd, yi + 1 + (xi-1)*yd, yi + xi*yd];
        count = count + 2;
    end
end

%Set up calculations for S and T
%Integrations over the triangle for pairs 1:3 1:3
Tc = [[1/12, -1/24, -1/24];[-1/24, 1/4, 1/8]; [-1/24, 1/8, 1/12]];
Sc = [[1, -1/2, -1/2];[-1/2,1/2,0];[-1/2,0,1/2]];

S = zeros(n,n); 
T = zeros(n,n);
xx = nodes(:,1);
yy = nodes(:,2);

for k = 1:N %elements
    elem = elems(k,:);
    i1 = elem(1);
    i2 = elem(2);
    i3 = elem(3);
    J = (xx(i2)-xx(i1))*(yy(i3)-yy(i1))-(xx(i3)-xx(i1))*(yy(i2)-yy(i1));

    for i = 1:3
        for j = 1:3
            S(elem(i),elem(j)) = S(elem(i),elem(j)) + J*Sc(i,j);
            T(elem(i),elem(j)) = T(elem(i),elem(j)) + J*Tc(i,j);
        end
    end
end

T_inv = inv(T);
%Iterate in time
U1 = U0; %Initial velocity of 0
U_t0 = U0;
U_t1 = U1;
U_t2 = U;
dt = 0.1;
steps = 60;
M = cell(steps+1);
M{1}=U_t0;
M{2}=U_t1;

% obj = VideoWriter('rectangle_wave.avi');
% obj.Quality= 100;
% obj.FrameRate = 20;
% 
% open(obj);
for t = 2:steps
    U_t2 = reshape(-T_inv*S*c^2*U_t1(:)*dt^2 + 2*U_t1(:)-U_t0(:),[yd,xd]);
    %Impose boundary constraints
    U_t2([1,yd],:) = 0;
    U_t2(:,[1,xd]) = 0;
    
    %Plotting
    M{t+1} = U_t2;
    surf(X,Y,U_t2)
    caxis([-0.1 0.1])
    xlim([0 1])
    ylim([0 1])
    zlim([-0.2,0.2])
    xlabel('x')
    ylabel('y')
    zlabel('Transverse displacement')
    title(strcat('Wave equation on a rectangle FEM, ', num2str(t*dt,'time = %4.1f (sec)')));
    drawnow
    hold off
%     pause(0.1)
%     f = getframe(gcf);
%     writeVideo(obj,f);

    U_t0 = U_t1; 
    U_t1 = U_t2; 
end

% obj.close();

%% Solving the 2D Wave Equation on an "Egg Shape" 

