

%% Starting with a rectangular grid
c =5;
xd= 31;
yd= 31;
x = linspace(0,1,xd);
y = linspace(0,1,yd);
[X,Y] = meshgrid(x,y);

n = xd*yd;
U = zeros(yd,xd);

%Initial Condition is Gaussian
U0 = 0.005*reshape(mvnpdf([X(:),Y(:)],[0.6,0.4],0.01*eye(2)),[yd,xd]);
U = U0;
% surf(X,Y,U0)

%Identify boundary points
b=0.4;
a=0.5;
k=0.2;
yi1s = zeros(1,xd);
yi2s = zeros(1,xd);
boundary = zeros(yd,xd);
boundary_nodes = zeros(yd,xd);
nodes = [X(:),Y(:)]; % n rows with two dimensions (cols)

for i = 1:length(x)
    xi = x(i);
    yi1 = b*sqrt(a^2-(xi-a).^2)/(a*sqrt(1+k*(xi-a)))+b;
    yi2 = -b*sqrt(a^2-(xi-a).^2)/(a*sqrt(1+k*(xi-a)))+b;
    yi1s(i) = yi1;
    yi2s(i) = yi2;
    k1 = dsearchn(nodes,[xi,yi1]);
    node = nodes(k1,:);
    boundary(int8(node(2)*(yd-1)+1):end,int8(node(1)*(xd-1)+1))=1;
    boundary_nodes(int8(node(2)*(yd-1)+1),int8(node(1)*(xd-1)+1))=1;
    k2 = dsearchn(nodes,[xi,yi2]);
    node = nodes(k2,:);
    boundary(1:int8(node(2)*(yd-1)+1),int8(node(1)*(xd-1)+1))=1;
    boundary_nodes(int8(node(2)*(yd-1)+1),int8(node(1)*(xd-1)+1))=1;
end


%% Perform triangulation
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
steps = 100;

% obj = VideoWriter('Egg_wave1.avi');
% obj.Quality= 100;
% obj.FrameRate = 20;
% 
% open(obj);
figure
P = surf(X,Y,U0);
    hold on
plot3(x,yi1s,zeros(xd,1),'m-','LineWidth',2);
plot3(x,yi2s,zeros(xd,1),'m-','LineWidth',2);
for t = 2:steps
    U_t2 = reshape(-T_inv*S*c^2*U_t1(:)*dt^2 + 2*U_t1(:)-U_t0(:),[yd,xd]);

    %Impose boundary conditions- egg
    U_t2 = U_t2.*(~boundary);
    U_nan = U_t2;
    U_nan(((boundary-boundary_nodes)==1))=nan;
    %Plotting
%     M{t+1} = U_t2;
%     P = surf(X,Y,U_t2);
    set(P,'ZData',U_nan);
    caxis([-0.1 0.1])
    xlim([0 1])
    ylim([0 1])
    zlim([-0.2,0.2])
    xlabel('x')
    ylabel('y')
    zlabel('Transverse displacement')
    title(strcat('Wave equation on an egg shape, FEM, ', num2str(t*dt,'time = %4.1f (sec)')));
     drawnow
%     pause(0.1)
%     f = getframe(gcf);
%     writeVideo(obj,f);

    U_t0 = U_t1; 
    U_t1 = U_t2; 
end
% 
% obj.close();

%% Solving the 2D Wave Equation on an "Egg Shape" 

