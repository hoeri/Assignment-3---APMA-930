function mit18086_navierstokes
%MIT18086_NAVIERSTOKES
%    Solves the incompressible Navier-Stokes equations in a
%    rectangular domain with prescribed velocities along the
%    boundary. The solution method is finite differencing on
%    a staggered grid with implicit diffusion and a Chorin
%    projection method for the pressure.
%    Visualization is done by a colormap-isoline plot for
%    pressure and normalized quiver and streamline plot for
%    the velocity field.
%    The standard setup solves a lid driven cavity problem.

% 07/2007 by Benjamin Seibold
%            http://www-math.mit.edu/~seibold/
% Feel free to modify for teaching and learning.

% 03/2017 by MC Kropinski - corrected BCs

close all; clear; clc;

Errors = [] ;

%for dt = [5e-1,1e-1,5e-2,1e-2,5e-3,1e-3]
Col = 0;

%[5e-1,1e-1,5e-2,1e-2,5e-3,1e-3]
for m = 50:50:500


%-----------------------------------------------------------------------
Re = 1e2;     % Reynolds number
dt = 1e-3;    % time step
tf = dt ;    % final time
lx = 1;       % width of box
ly = 1;       % height of box
dx = 1/m ;
%dx = 1e-2 ;   % Set delta x 
nx = lx/dx      % number of x-gridpoints
ny = nx;      % number of y-gridpoints
nsteps = 7;  % number of steps with graphic output


%-----------------------------------------------------------------------
%nt = ceil(tf/dt); dt = tf/nt;
x = linspace(0,lx,nx+1); hx = lx/nx;
y = linspace(0,ly,ny+1); hy = ly/ny;
[X,Y] = meshgrid(y,x);
%-----------------------------------------------------------------------
% initial conditions
V = zeros(nx,ny-1) ;

% Set the parabolic profile initial conditions. 
 U = zeros(nx-1,ny);
 
 for n = 1:nx-1 
     %U(n,:) = y(1:end-1).*(ly - y(2:end)) ;
     U(n,:) = avg(y).*(ly - avg(y)) ;
 end
 
% boundary conditions
 vN = avg(x)*0;
 vS = avg(x)*0;
 vW = y*0;
 vE = y*0;


% New boundary conditions for infinite channel like flow. The BCs for V
% don't cange and we set the BCs for U equal to the exact solution for
% channel flow on the boundaries. The top and bottom get nN and nS get no
% slip BCs. 

uW = avg(y).*(ly - avg(y)) ;
uE = avg(y).*(ly - avg(y)) ;
% uW = y(1:end-1).*(ly - y(2:end)) ;
% uE = y(1:end-1).*(ly - y(2:end)) ;

uN = x*0 ;
uS = x*0 ;


% uE(1) = 0;
% uE(end) = 0;
% uW(1) = 0;
% uE(end) = 0;


% Store the exact negative pressure gradient
DPexact = -2/Re*ones(1,numel(avg(y))) ;


%-----------------------------------------------------------------------
%Original BCs, with bug
%Ubc = dt/Re*([2*uS(2:end-1)' zeros(nx-1,ny-2) 2*uN(2:end-1)']/hy^2+...
%      [uW;zeros(nx-3,ny);uE]/hx^2);
%Vbc = dt/Re*([vS' zeros(nx,ny-3) vN']/hy^2+...
%      [2*vW(2:end-1);zeros(nx-2,ny-1);2*vE(2:end-1)]/hx^2);
%
% Corrected BCs:
Ubc = dt/Re*([2*uS(2:end-1)' zeros(nx-1,ny-2) 2*uN(2:end-1)']/hy^2+...
      [uW;zeros(nx-3,ny);uE]/hx^2);
 
Vbc = dt/Re*([vS' zeros(nx,ny-3) vN']/hy^2+...
      [2*vW(2:end-1);zeros(nx-2,ny-1);2*vE(2:end-1)]/hx^2);
  

  
fprintf('initialization')
Lp = kron(speye(ny),K1(nx,hx,1))+kron(K1(ny,hy,1),speye(nx));
Lp(1,1) = 3/2*Lp(1,1);
 
perp = symamd(Lp); Rp = chol(Lp(perp,perp)); Rpt = Rp';
Lu = speye((nx-1)*ny)+dt/Re*(kron(speye(ny),K1(nx-1,hx,2))+...
     kron(K1(ny,hy,3),speye(nx-1)));
peru = symamd(Lu); Ru = chol(Lu(peru,peru)); Rut = Ru';
Lv = speye(nx*(ny-1))+dt/Re*(kron(speye(ny-1),K1(nx,hx,3))+...
     kron(K1(ny-1,hy,2),speye(nx)));
perv = symamd(Lv); Rv = chol(Lv(perv,perv)); Rvt = Rv';
Lq = kron(speye(ny-1),K1(nx-1,hx,2))+kron(K1(ny-1,hy,2),speye(nx-1));
perq = symamd(Lq); Rq = chol(Lq(perq,perq)); Rqt = Rq';

fprintf(', time loop\n--20%%--40%%--60%%--80%%-100%%\n')
for k = 1:1
   % treat nonlinear terms
   gamma = min(1.2*dt*max(max(max(abs(U)))/hx,max(max(abs(V)))/hy),1);
   Ue = [uW;U;uE]; Ue = [2*uS'-Ue(:,1) Ue 2*uN'-Ue(:,end)];
   Ve = [vS' V vN']; Ve = [2*vW-Ve(1,:);Ve;2*vE-Ve(end,:)];
   Ua = avg(Ue')'; Ud = diff(Ue')'/2;
   Va = avg(Ve);   Vd = diff(Ve)/2;
   UVx = diff(Ua.*Va-gamma*abs(Ua).*Vd)/hx;
   UVy = diff((Ua.*Va-gamma*Ud.*abs(Va))')'/hy;
   Ua = avg(Ue(:,2:end-1));   Ud = diff(Ue(:,2:end-1))/2;
   Va = avg(Ve(2:end-1,:)')'; Vd = diff(Ve(2:end-1,:)')'/2;
   U2x = diff(Ua.^2-gamma*abs(Ua).*Ud)/hx;
   V2y = diff((Va.^2-gamma*abs(Va).*Vd)')'/hy;
   U = U-dt*(UVy(2:end-1,:)+U2x);
   V = V-dt*(UVx(:,2:end-1)+V2y);
   
   
   
   % implicit viscosity
   rhs = reshape(U+Ubc,[],1);
   u(peru) = Ru\(Rut\rhs(peru));
   U = reshape(u,nx-1,ny);
   rhs = reshape(V+Vbc,[],1);
   v(perv) = Rv\(Rvt\rhs(perv));
   V = reshape(v,nx,ny-1);
   
  
    
   
   % pressure correction
   rhs = reshape(diff([uW;U;uE])/hx+diff([vS' V vN']')'/hy,[],1);
   p(perp) = -Rp\(Rpt\rhs(perp));
   P = reshape(p,nx,ny);
   U = U-diff(P)/hx;
   V = V-diff(P')'/hy;
   
   
 
   
   % Extract and compute the pressure gradient
   DP = (1/dt)*diff(P)/hx ;
   %imagesc(DP')
   DP = DP(ceil(nx/2),:) ;
%    figure(1)
%    plot(avg(y),DP,'-','Color',[0.8, Col/9 ,0],'LineWidth',2, 'DisplayName', ['dx = ', num2str(dx)])
%      legend('-DynamicLegend', 'Location', 'north')
%    title(['Pressure Gradient at dt = ' num2str(dt)])
%    xlabel('y range')
%    ylabel('Pressure Gradient')
%    hold on
   
   Errors = [Errors norm(DP-DPexact,2)] ;
   
   Col = Col + 1; 
   
%    % visualization
%    if floor(25*k/nt)>floor(25*(k-1)/nt), fprintf('.'), end
%    if k==1|floor(nsteps*k/nt)>floor(nsteps*(k-1)/nt)
%       % stream function
%       rhs = reshape(diff(U')'/hy-diff(V)/hx,[],1);
%       q(perq) = Rq\(Rqt\rhs(perq));
%       Q = zeros(nx+1,ny+1);
%       Q(2:end-1,2:end-1) = reshape(q,nx-1,ny-1);
%       clf, 
%       %contourf(avg(x),avg(y),P',20,'w-'),
%       hold on
%       %contour(x,y,Q',20,'k-');
%       Ue = [uS' avg([uW;U;uE]')' uN'];
%       Ve = [vW;avg([vS' V vN']);vE];
%       Len = sqrt(Ue.^2+Ve.^2+eps);
%       %quiver(x,y,(Ue./Len)',(Ve./Len)',.4,'k-')
%       colorbar
%       contourf(x,y,Ve',20,'w-')
%       
%       hold off, axis equal, axis([0 lx 0 ly])
%       p = sort(p); caxis(p([8 end-7]))
%       title(sprintf('Re = %0.1g   t = %0.2g',Re,k*dt))
%       drawnow
%    end
end
fprintf('\n')

end

%plot(avg(y),DPexact,'b','LineWidth',2,'DisplayName','Exact Gradient')

hold off 

% Fit a line to the log of the data
Coeff = polyfit(log(1./(50:50:500)), log(Errors),1) ;
disp(Coeff(1))
clf
hold on
title(['LogLog plot of Error of Pressure Gradient vs dx, with slop = ' num2str(Coeff(1))])
plot(log(1./(50:50:500)),log(Errors),'b*')
plot(-7:0.1:-3.5, polyval(Coeff,-7:0.1:-3.5),'r')
xlabel('dx value')
ylabel('pressure gradient')
legend('Errors','Best fit line','Location','NorthWest')
hold off

% hold on
% loglog((1./(50:50:500)),[(Errors)],'b*')
% plot((1./(50:50:500)),(polyval(Coeff,1./(50:50:500))),'r')
% 
% hold off



% plot(avg(y),DPexact,'r','LineWidth',2,'Displayname','Exact Gradient')
% hold off
%=======================================================================

function B = avg(A,k)
if nargin<2, k = 1; end
if size(A,1)==1, A = A'; end
if k<2, B = (A(2:end,:)+A(1:end-1,:))/2; else, B = avg(A,k-1); end
if size(A,2)==1, B = B'; end

function A = K1(n,h,a11)
% a11: Neumann=1, Dirichlet=2, Dirichlet mid=3;
A = spdiags([-1 a11 0;ones(n-2,1)*[-1 2 -1];0 a11 -1],-1:1,n,n)'/h^2;


