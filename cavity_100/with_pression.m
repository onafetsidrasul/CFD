dx=0.01;
dy=0.01;

x=dx:dx:1;
y=dy:dy:1;
np=length(x);
nx=np;
nmt=6000;
dt=0.001;

d=0.01; %Re=100

u=zeros(np,np,nmt);
u_star=zeros(np,np);
v=zeros(np,np,nmt);
v_star=zeros(np,np);
S=zeros(np,np);
p=zeros(np,np,nmt);

% NOTE : the part below contains initialization of variables needed and
% coefficient matrix construction. To be add befor the time-loop

% Initial condition
niter     = 50;
dnx       = np*np;
pv        = zeros(dnx,1);
sv        = zeros(dnx,1);
res_story = zeros(niter,1);
ite_story = zeros(niter,1);


% build matrix A (laplacian operator) and N
A = zeros(dnx,dnx);
N = zeros(dnx,dnx);

for j=2:nx-1
    for i=2:nx-1
        index_ij  = nx*(j-1) + i;
        index_top = nx*(j-1) + i + nx;
        index_bot = nx*(j-1) + i - nx;
        index_lef = nx*(j-1) + i + 1;
        index_rig = nx*(j-1) + i - 1;
        %
        A(index_ij, index_ij)     = -4;    
        A(index_ij, index_rig) =  1;         % right coeff
        A(index_ij, index_lef) =  1;         % left coeff
        A(index_ij, index_bot) =  1;         % top coeff
        A(index_ij, index_top) =  1;         % bottom coeff
    end
end
% correct boundary conditions on laplacian operator
% FIRST/LAST nx rows - top/bottom boundary cond - 
% vertical derivative = 0 p_ij+1 = 0 and p_ij-2 x 2 and viceversa
% Then, left/right boundaries have horiz derivative = 0
for i=2:nx-1
        index_ij  =  i;
        index_bot =  nx + i;
        A(index_ij,  index_ij)     = -4;    
        A(index_ij,  index_ij + 1) =  1;         % right coeff
        A(index_ij,  index_ij - 1) =  1;         % left coeff
        A(index_ij,  index_bot)    =  2;         % bottom coeff - vertical derivative = 0
end
j=nx;
for i=2:nx-1
        index_ij  = nx*(j-1) + i;
        index_top = nx*(j-2) + i;
        A(index_ij, index_ij)     = -4;    
        A(index_ij, index_ij + 1) =  1;         % right coeff
        A(index_ij, index_ij - 1) =  1;         % left coeff
        A(index_ij, index_top)    =  2;         % top coeff - vertical derivative = 0
end
i = 1;
for j=2:nx-1
        index_ij  = nx*(j-1) + i;
        index_top = nx*(j-2) + i;
        index_bot = nx*(j  ) + i;
        A(index_ij, index_ij)     = -4;    
        A(index_ij, index_ij + 1) =  2;      % right coeff
        A(index_ij, index_top) =  1;         % top coeff
        A(index_ij, index_bot) =  1;         % bottom coeff
end
i = nx;
for j=2:nx-1
        index_ij  = nx*(j-1) + i;
        index_top = nx*(j-2) + i;
        index_bot = nx*(j  ) + i;
        A(index_ij, index_ij)     = -4;    
        A(index_ij, index_ij - 1) =  2;      % left coeff
        A(index_ij, index_top) =  1;         % top coeff
        A(index_ij, index_bot) =  1;         % bottom coeff
end
% corners
A(1,1)      = 1;
A(dnx,dnx)  = 1;
A(nx,nx)    = 1;
A(nx*(nx-1)+1,nx*(nx-1)+1) = 1;
%
for i=1:dnx
    N(i,i)   = A(i,i);
end
Nm1 = inv(N);

%ic
u(:,1,1)=1;
%shear on top

ss_x=d*dt/(dx^2);
ss_y=d*dt/(dy^2);
s_x=dt/dx;
s_y=dt/dy;

for n=2:nmt
    for i=2:np-1
        for j=2:np-1
            %u_star(i,j)=u(i,j,n-1)-s_x*(u(i,j,n-1)^2-u(i-1,j,n-1)^2)-s_y*(u(i,j,n-1)*v(i,j,n-1)-u(i,j-1,n-1)*v(i,j-1,n-1))...
                %+ss_x*(u(i+1,j,n-1)-2*u(i,j,n-1)+u(i-1,j,n-1))+ss_y*(u(i,j+1,n-1)-2*u(i,j,n-1)+u(i,j-1,n-1));
            %v_star(i,j)=v(i,j,n-1)-s_y*(v(i,j,n-1)^2-v(i,j-1,n-1)^2)-s_x*(v(i,j,n-1)*u(i,j,n-1)-v(i-1,j,n-1)*u(i-1,j,n-1))...
                %+ss_x*(v(i+1,j,n-1)-2*v(i,j,n-1)+v(i-1,j,n-1))+ss_y*(v(i,j+1,n-1)-2*v(i,j,n-1)+v(i,j-1,n-1));
            %u_star(i,j)=u(i,j,n-1)-s_x*u(i,j,n-1)*(u(i,j,n-1)-u(i-1,j,n-1))-s_y*v(i,j,n-1)*(u(i,j,n-1)-u(i,j-1,n-1))...
                %+ss_x*(u(i+1,j,n-1)-2*u(i,j,n-1)+u(i-1,j,n-1))+ss_y*(u(i,j+1,n-1)-2*u(i,j,n-1)+u(i,j-1,n-1));
            %v_star(i,j)=v(i,j,n-1)-s_y*v(i,j,n-1)*(v(i,j,n-1)-v(i-1,j,n-1))-s_x*u(i,j,n-1)*(v(i,j,n-1)-v(i-1,j,n-1))...
                %+ss_x*(v(i+1,j,n-1)-2*v(i,j,n-1)+v(i-1,j,n-1))+ss_y*(v(i,j+1,n-1)-2*v(i,j,n-1)+v(i,j-1,n-1));
            u_star(i,j) = u(i,j,n-1) ...
                - s_x * ( (u(i+1,j,n-1)^2 - u(i-1,j,n-1)^2) / 2 ) ...
                - s_y * ( (u(i,j+1,n-1)*v(i,j+1,n-1) - u(i,j-1,n-1)*v(i,j-1,n-1)) / 2 ) ...
                + ss_x * (u(i+1,j,n-1) - 2*u(i,j,n-1) + u(i-1,j,n-1)) ...
                + ss_y * (u(i,j+1,n-1) - 2*u(i,j,n-1) + u(i,j-1,n-1));

            v_star(i,j) = v(i,j,n-1) ...
                - s_y * ( (v(i,j+1,n-1)^2 - v(i,j-1,n-1)^2) / 2 ) ...
                - s_x * ( (v(i+1,j,n-1)*u(i+1,j,n-1) - v(i-1,j,n-1)*u(i-1,j,n-1)) / 2 ) ...
                + ss_x * (v(i+1,j,n-1) - 2*v(i,j,n-1) + v(i-1,j,n-1)) ...
                + ss_y * (v(i,j+1,n-1) - 2*v(i,j,n-1) + v(i,j-1,n-1));
        end
    end
    u_star(:,1)=1;
    u_star(1,:)=0;
    u_star(np,:)=0;
    u_star(:,np)=0;
    v_star(1,:)=0;
    v_star(:,1)=0;
    v_star(np,:)=0;
    v_star(:,np)=0;
    for i=2:np-1
        for j=2:np-1
            S(i,j)=((u_star(i+1,j)-u_star(i-1,j))/(2*dx*dt)+(v_star(i,j+1)-v_star(i,j-1))/(2*dy*dt))*dx^2;
        end
    end
    index=1;
    for i=1:np
        for j=1:np
            sv(index)=S(i,j);
            pv(index)=p(i,j,n-1);
            index=index+1;
        end
    end
    %solver
    eps        = 0.000001; % minimum residual
    nitermax   = 300;      % max iteration for the jacobi method
    res_mod    = 1; 
    exit_tol   = 1;
    count_iter = 1;
    xk         = pv;
    b          = sv;
    while count_iter < nitermax & exit_tol > eps
        res        = A*xk - b;
        %xkp1       = xk - Nm1*res;
        xkp1 = xk - res ./ diag(N); %CHANGED
        res_mod    = norm(A*xkp1 - b); % residual magnitude
        exit_tol   = norm(xkp1 - xk);
        xk         = xkp1;
        count_iter = count_iter + 1;
    end
    res_story(n,1) = res_mod;
    
    pv = xkp1;
    
    %back to matrix for p
    index=1;
    for i=1:np
        for j=1:np
            p(i,j,n)=pv(index);
            index=index+1;
        end
    end

    p(:,1,n) = p(:,2,n);
    p(:,nx,n) = p(:,nx-1,n);
    p(1,:,n) = p (2,:,n);
    p(nx,:,n) = p(nx-1,:,n);
    
    %corrector
    for i=2:np-1
        for j=2:np-1
            u(i,j,n)=u_star(i,j)-dt*(p(i+1,j,n)-p(i-1,j,n))/(2*dx);
            v(i,j,n)=v_star(i,j)-dt*(p(i,j+1,n)-p(i,j-1,n))/(2*dy);
        end
    end
    u(:,1,n)=1;
    u(1,:,n)=0;
    u(np,:,n)=0;
    u(:,np,n)=0;
    v(1,:,n)=0;
    v(:,1,n)=0;
    v(np,:,n)=0;
    v(:,np,n)=0;
end

figure(2);
for t = 1:nmt
    contourf(x, y, u(:,:,t));
    ylim([-0.1,1.1]);
    shading interp;
    colormap jet;
    colorbar;
    %clim([0 1]); % Normalizza il range
    title(['Time step: ', num2str(t)]);
    xlabel('X'); ylabel('Y'); zlabel('c(x,y,t)');
    pause(0.05);
end 
%momentum diffusion!