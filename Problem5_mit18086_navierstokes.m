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

SteadyTime = [];


%-----------------------------------------------------------------------
Re = 3000;     % Reynolds number
dt = 1e-2;    % time step
tf = 13e0;    % final time
lx = 2;       % width of box
ly = 1;       % height of box
nx = 180;      % number of x-gridpoints
ny = 90;      % number of y-gridpoints
nsteps = 1000;  % number of steps with graphic output
%-----------------------------------------------------------------------
nt = ceil(tf/dt); dt = tf/nt;
x = linspace(0,lx,nx+1); hx = lx/nx;
y = linspace(0,ly,ny+1); hy = ly/ny;
[X,Y] = meshgrid(y,x);
%-----------------------------------------------------------------------


% initial conditions
 V = zeros(nx,ny-1);
 
 

 
 
 
 
% Set the parabolic profile initial conditions. 
 U = zeros(nx-1,ny);
 
 for n = 1:nx-1 
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
uN = x*0 ;
uS = x*0 ;
uW = avg(y).*(ly - avg(y)) ;
uE = avg(y).*(ly - avg(y)) ;




% Add a perturbation to V, close to the left boundary. 
Pert = 1*ones(10,10) ;
V(10+1:20,40+1:50) = Pert ;





% Initialize the Vectors to store the magnitude of the flow, to check if
% the disturbance is growning. 

% One array for 1 norm
GroOne = zeros(1,nt) ;

% One array for inf norm
GroInf = zeros(1,nt) ;





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
for k = 1:nt
   % treat nonlinear terms
   gamma = min(1.2*dt*max(max(max(abs(U)))/hx,max(max(abs(V)))/hy),1);
   Ue = [uW;U;uE]; 
   Ue = [2*uS'-Ue(:,1) Ue 2*uN'-Ue(:,end)];
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
   
   
   % Check for the relative error
%    tol = 1e-6 ;
%    if k>1
%        Error = norm(U - Uinter)/norm(U) ;
%        if Error < tol
%            disp(['Returned to steady state after t = ' num2str(k*dt)])
%            SteadyTime = [SteadyTime k*dt ];
%            break
%        end
%    end
%    Uinter = U ;
   
    rhs = reshape(diff(U')'/hy-diff(V)/hx,[],1);
      q(perq) = Rq\(Rqt\rhs(perq));
      Q = zeros(nx+1,ny+1);
      Q(2:end-1,2:end-1) = reshape(q,nx-1,ny-1);
      clf, contourf(avg(x),avg(y),P'/dt,20,'w-'), 
      GroOne(k) = norm(Q,1) ;
      GroInf(k) = norm(Q,inf) ;
   
   

   % visualization
   if floor(25*k/nt)>floor(25*(k-1)/nt), fprintf('.'), end
   if k==1|floor(nsteps*k/nt)>floor(nsteps*(k-1)/nt)
      % stream function
      rhs = reshape(diff(U')'/hy-diff(V)/hx,[],1);
      q(perq) = Rq\(Rqt\rhs(perq));
      Q = zeros(nx+1,ny+1);
      Q(2:end-1,2:end-1) = reshape(q,nx-1,ny-1);
      clf, contourf(avg(x),avg(y),P'/dt,20,'w-'), 
     
      hold on 
      contour(x,y,Q',20,'k-');
      Ue = [uS' avg([uW;U;uE]')' uN'];
      Ve = [vW;avg([vS' V vN']);vE];
      Len = sqrt(Ue.^2+Ve.^2+eps);
     % quiver(x,y,(Ue./Len)',(Ve./Len)',.4,'k-')
      hold off, axis equal, axis([0 lx 0 ly])
      p = sort(p); caxis(p([8 end-7]))
      title(sprintf('Re = %0.1g   t = %0.2g',Re,k*dt))
      drawnow
   end
end
fprintf('\n')


close
clf
hold on 
title('Change in the magnitude of the stream function')
plot(1:nt,GroInf,'r*') 
plot(1:nt,GroOne,'bo')
ylabel('Magnitude of stream function')
xlabel('time step')
legend('Infinity Norm','One Norm')
hold off


% clf
% hold on
% title('Time to reach steady state vs Re')
% ylabel('Time')
% xlabel('Reynolds number')
% plot([10 55 100 550 800 1000 3000 5500],SteadyTime,'b*','LineWidth',2)
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
