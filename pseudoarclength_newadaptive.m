% Function uses Adaptive step size Pseudo arclength continuation to trace out solution family
%
% Solves F(u,lambda,gamma)=0 and G(u,lambda,gamma)=0 then traces solution
% family with varying lambda gamma
%
% Preconditioned relaxations only!
%
% Uses adaptive step length ds, and if not converged, halves step size and
% tries again
%
% Continues until steps or if step size exceeds a minimum chosen
%
% Algorithm follows that of solution technique in nonlinear oceanography
%
% Inputs:
% v - initial solution          -u_xx+au_x+b*lambda*u+cu^2=gamma*f
% lambda - initial lambda       -u_xx+au_x+b*lambda*u+cu^2=gamma*f
% gamma - initiial gamma        -u_xx+au_x+b*lambda*u+cu^2=gamma*f
% dv - initial tangent of v (or approximation of)
% dlambda - initial tangent of lambda (or approximation of)
% dgamma - initial tangent of gamma (or approximation of)
% ds - step size (constant)
% steps - number of steps to take
% x - vector of x values
% pre_relaxation - iterative scheme used (see main)
% vcyclegrid - number of grids used in v-cycle
% cellN - cell of grid points
% cellk - cell of wave number
% Nd - number of iterations per grid down cycle
% Nu - number of iterations per grid up cycle
% cella - cell of function a(x) in -u_xx+au_x+b*lambda*u+cu^2=gamma*f
% cellb - cell of function b(x) in -u_xx+au_x+b*lambda*u+cu^2=gamma*f
% cellc - cell of function c(x) in -u_xx+au_x+b*lambda*u+cu^2=gamma*f
% cellRHS - cell of f(x) function
% preconditioner - preconditioner chosen
%
% Outputs:
% v - vector of solution curves
% lambda - vector of lambda values
% gamma - vector or gamma values

function [v,lambda,gamma]=pseudoarclength_newadaptive(v,lambda,gamma,dv,dlambda,dgamma,ds,steps,... 
    x,pre_relaxation,vcyclegrid,cellN,cellk,Nd,Nu,cella,cellb,cellc,cellRHS,preconditioner)

global L

psi1=1/cellN{1};
psi2=1;
psi3=1;

% nominal contraction rate
kappanom=0.0;
% nominal distance to curve
deltanom=0.5; % or max step size
% nominal angle
alphanom=0.001;

% Initial values (since need 2 steps before application)
kappa=kappanom;
delta=deltanom;
alpha=alphanom;
f=1;
dv_old=dv;dgamma_old=dgamma;dlambda_old=dlambda;

% Optimal number of iterations for Newton (roughly 5)
N_opt=3;

% Maximum number of newton loops before breaking and halving step size
max_newton_loops=2*N_opt;

% Minimum step length before breaking loop due to non convergence
ds_min=1e-6;

% Step counter
j=1;

% Save base variable value
b=cellb{1};
RHS=cellRHS{1};

% Set cells for z1,z2,z3,RHS1,RHS2,RHS3 (for multigrid solver)
cellz1=cell(vcyclegrid,1);
cellz2=cell(vcyclegrid,1);
cellz3=cell(vcyclegrid,1);
cellRHS1=cell(vcyclegrid,1);
cellRHS2=cell(vcyclegrid,1);
cellRHS3=cell(vcyclegrid,1);
cellbnew=cell(vcyclegrid,1);

cellz1{1}=zeros(cellN{1},1);
cellz2{1}=cellz1{1};
cellz3{1}=cellz1{1};

% Loops steps
while j<steps
    
    % Tangent predictor
    v(:,j+1)=v(:,j)+ds*dv;
    lambda(j+1)=lambda(j)+ds*dlambda;
    gamma(j+1)=gamma(j)+ds*dgamma;
    
    % find nominal contraction/distance to curve/angle
    kappa=(psi1*dot(dv,dv)+psi2*dgamma.^2+psi3*dlambda.^2)/(psi1*dot(dv_old,dv_old)+psi2*dgamma_old.^2+psi3*dlambda_old.^2);
    delta=(psi1*dot(dv,dv)+psi2*dgamma.^2+psi3*dlambda.^2);
    alpha=acos(psi1*dot(dv,dv_old)+psi2*dlambda*dlambda_old+psi3*dgamma*dgamma_old);
    
    % Update variable wrt parameter
    cellb{1}=lambda(j+1)*b;
    cellRHS{1}=gamma(j+1)*RHS;
    
    % New b(x) function in Newton
    bnew=cellb{1}+2*cellc{1}.*v(:,j+1);
    H=preconditioner(L,cellN{1},cella{1},bnew);
    
    % Correction loop using Newton iterations
    for i=1:max_newton_loops
        
        converged=false;
        
        if f==2
        
            fprintf('f too large\n')

            % Halve step size
            ds=ds/2;
            fprintf('Halving step size to %f\n',ds)
        
            % Break loop if minimum step size exceeded
            if ds<=ds_min
            	fprintf('Minimum step length of %f exceeded, breaking loop\n',ds_min)
            	return % end function
            end
        end
        
        % F and F derivatives
        F=NA(v(:,j+1),cellk{1},cella{1},cellb{1},cellc{1})-cellRHS{1};
        F_lambda=sign(b).*v(:,j+1);
        F_gamma=-RHS;
                
        % G_gamma calculation using trapezoidal integration
        RHSx=real(ifft(1i*cellk{1}'.*fft(RHS)));
        G_gamma=trapint(RHSx(1:cellN{1}/2+1).*v(1:cellN{1}/2+1,j+1),x(1:cellN{1}/2+1));
        
        % G calculation
        G=lambda(j+1)^3/54+trapint(gamma(j+1)*RHSx(1:cellN{1}/2+1).*v(1:cellN{1}/2+1,j+1),x(1:cellN{1}/2+1));
        
        % G_u calculation
        
        cn=zeros(cellN{1},1);
        cn(2:cellN{1}/2)=L/cellN{1};
        cn(1)=1/2*L/cellN{1};
        cn(cellN{1}/2+1)=1/2*L/cellN{1};
        G_u=gamma(j+1)*RHSx.*cn;
        
        % G_lambda
        G_lambda=lambda(j+1)^2/18;
        
        % Check convergence condition
        if rms(F)<=1e-10 && rms(G)<1e-10%v(cellN{1}/2+1,j+1)+1/3*lambda(j+1)<=1e-10%

            fprintf('Converged after %d Newton Iterations step = %d\n',i-1,j)
            j=j+1;
            converged=true;
            
            % Optimal step length control for next step
            f=max([sqrt(kappanom/kappa),sqrt(deltanom/delta),alphanom/alpha]);
            f=max(min(f,2),1/2);
            ds=ds/f;
            fprintf('New step size to %f\n',ds)
            
            % Max step size check
            if ds<ds_min
                fprintf('Minimum step length of %f exceeded, breaking loop\n',ds_min)
                return
            end
            
            break % End if converged to sufficient accuracy
        end

        % Arclength equation
        P=psi1*dot(dv,(v(:,j+1)-v(:,j)))+psi2*dlambda*(lambda(j+1)-lambda(j))+psi3*dgamma*(gamma(j+1)-gamma(j))-ds;
        
        % Initial RHS of linear equation
        RHS1=-F;
%         RHS2=[F_lambda F_gamma];
        RHS2=F_lambda;
        RHS3=F_gamma;
        
        % Set cell for bnew and step downs
        cellbnew=setcellsNewton(vcyclegrid,cellbnew,bnew);
        cellH=setcellsH(vcyclegrid,H,L,cellN,cella,cellbnew);
        
        % Update cell RHS with new values
        [cellRHS1,cellRHS2,cellRHS3]=setcellspseudo2(vcyclegrid,cellN,cellRHS1, ...
            cellRHS2,cellRHS3,RHS1,RHS2,RHS3); 
        
        % Solve F_u*z1=-F and F_u*z2=F_lambda (F_u is jacobian)
%         cellz1{1}=LA(eye(cellN{1},cellN{1}),cellk{1},cella{1},cellbnew{1})\RHS1;
%         cellz2{1}=LA(eye(cellN{1},cellN{1}),cellk{1},cella{1},cellbnew{1})\RHS2;
%         cellz3{1}=LA(eye(cellN{1},cellN{1}),cellk{1},cella{1},cellbnew{1})\RHS3;
%         [cellz1{1},r1]=vcycle_pre(pre_relaxation,1,vcyclegrid,cellN,cellk,1,Nd,Nu,cella,cellbnew,cellRHS1,cellH,cellz1);
%         [cellz2{1},r2]=vcycle_pre(pre_relaxation,1,vcyclegrid,cellN,cellk,1,Nd,Nu,cella,cellbnew,cellRHS2,cellH,cellz2);
%         [cellz3{1},r3]=vcycle_pre(pre_relaxation,1,vcyclegrid,cellN,cellk,1,Nd,Nu,cella,cellbnew,cellRHS3,cellH,cellz3);
%         [cellz1,r1]=FMG_pre(pre_relaxation,vcyclegrid,cellN,cellk,1,Nd,Nu,cella,cellbnew,cellRHS1,cellH,cellz1);
%         [cellz2,r2]=FMG_pre(pre_relaxation,vcyclegrid,cellN,cellk,1,Nd,Nu,cella,cellbnew,cellRHS2,cellH,cellz2);
%         [cellz3,r3]=FMG_pre(pre_relaxation,vcyclegrid,cellN,cellk,1,Nd,Nu,cella,cellbnew,cellRHS3,cellH,cellz3);
        [cellz1{1},r1]=Pcg(cellz1{1},cellk{1},cella{1},cellbnew{1},cellRHS1{1},H);
        [cellz2{1},r2]=Pcg(cellz2{1},cellk{1},cella{1},cellbnew{1},cellRHS2{1},H);
        [cellz3{1},r3]=Pcg(cellz3{1},cellk{1},cella{1},cellbnew{1},cellRHS3{1},H);

        % Solving for delta_v delta_lambda delta_gamma
        Y=([G_lambda G_gamma;psi2*dlambda psi3*dgamma]-[G_u';psi1*dv']*[cellz2{1} cellz3{1}])\([-G;-P]-[G_u';psi1*dv']*cellz1{1});
        delta_lambda=Y(1);
        delta_gamma=Y(2);
        delta_v=cellz1{1}-[cellz2{1} cellz3{1}]*Y;
        
        % Update variable for next Newton iteration
        v(:,j+1)=v(:,j+1)+delta_v;
        lambda(j+1)=lambda(j+1)+delta_lambda;
        gamma(j+1)=gamma(j+1)+delta_gamma;
        
        % Update variable wrt parameter
        cellb{1}=lambda(j+1)*b;
        cellRHS{1}=gamma(j+1)*RHS;
        
        % Update vector b for next iteration
        bnew=cellb{1}+2*cellc{1}.*v(:,j+1);
        H=preconditioner(L,cellN{1},cella{1},cellbnew{1});
        
    end

    if (i==max_newton_loops && converged~=true) || (f==2 && j>2)
        
        fprintf('Did not converge to required tolerance after %d Newton Iterations at step %d\n',i,j)
        % Do not save latest vectors if not converged
        v(:,j+1)=[];
        lambda(:,j+1)=[];
        gamma(:,j+1)=[];
        
        % Halve step size
        ds=ds/2;
        fprintf('Halving step size to %f\n',ds)
        
        % Break loop if minimum step size exceeded
        if ds<=ds_min
            fprintf('Minimum step length of %f exceeded, breaking loop\n',ds_min)
            return % end function
        end
        
    elseif converged==true

        % Solving for dlambda dgamma dv (new direction)
        Y=([G_lambda G_gamma;psi2*dlambda psi3*dgamma]-[G_u psi1*dv]'*[cellz2{1} cellz3{1}])\[0;1];
        
        dlambda_old=dlambda;
        dgamma_old=dgamma;
        dv_old=dv;
        
        dlambda=Y(1);
        dgamma=Y(2);
    
    %     dv=-cellz2{1}*Y;
        dv=-[cellz2{1} cellz3{1}]*Y;
    
        % Normalise
        mag=sqrt(psi1*dot(dv,dv)+psi2*dlambda^2+psi3*dgamma^2);
        dlambda=dlambda/mag;
        dgamma=dgamma/mag;
        dv=dv/mag;
        
        % check orientation and correct if necessary
        if psi1*dot(dv,dv_old)+psi2*dlambda*dlambda_old+psi3*dgamma*dgamma_old<0
            fprintf('Sign change detected!\n')
            disp(psi1*dot(dv,dv_old)+psi2*dlambda*dlambda_old+psi3*dgamma*dgamma_old)
            Y=-Y;
            dlambda=Y(1);
            dgamma=Y(2);
            dv=-[cellz2{1} cellz3{1}]*Y;
            
            % Normalise
            mag=sqrt(psi1*dot(dv,dv)+psi2*dlambda^2+psi3*dgamma^2);
            dlambda=dlambda/mag;
            dgamma=dgamma/mag;
            dv=dv/mag;
            disp(psi1*dot(dv,dv_old)+psi2*dlambda*dlambda_old+psi3*dgamma*dgamma_old)
        end
       
             
%         disp(psi1*dot(dv,dv1)+psi2*dlambda*dlambda1+psi3*dgamma*dgamma1)
        disp(psi1*dot(dv,dv_old)+psi2*dlambda*dlambda_old+psi3*dgamma*dgamma_old);
       
    end
    
end
end