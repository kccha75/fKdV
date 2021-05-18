% Function uses Pseudo arclength continuation to trace out solution family

% Preconditioned relaxations only!

% Uses fixed step length ds
% Continues until steps or if Newton does not converge (after 20
% iterations)

% Algorithm follows that of solution technique in nonlinear oceanography

% Inputs:
% v - initial solution 
% lambda - initial lambda
% dv - initial tangent of v (or approximation of)
% dlambda - initial tangent of lambda (or approximation of)
% ds - step size (constant)
% steps - number of steps to take
% x - vector of x values
% relaxation - iterative scheme used (see main)
% vcyclegrid - number of grids used in v-cycle
% cellN - cell of grid points
% cellk - cell of wave number
% Nd - number of iterations per grid down cycle
% Nu - number of iterations per grid up cycle
% cella - cell of function a(x) in -u_xx+au_x+bu+cu^2=f
% cellb - cell of function b(x) in -u_xx+au_x+bu+cu^2=f
% cellc - cell of function c(x) in -u_xx+au_x+bu+cu^2=f

% Outputs:
% v - vector of solution curves
% lambda - vector of lambda values

function [v,gamma]=pseudoarclength_pre(v,gamma,dv,dgamma,ds,steps,...
    pre_relaxation,vcyclegrid,cellN,cellk,Nd,Nu,cella,cellb,cellc,cellRHS,preconditioner)

global L

psi1=1;
psi2=0.1;

% Optimal number of iterations for Newton (roughly 5)
N_opt=3;

% Maximum number of newton loops before breaking and halving step size
max_newton_loops=2*N_opt;

% Minimum step length before breaking loop due to non convergence
ds_min=1e-6;

% Maximum step length
ds_max=0.5;

% Step counter
j=1;

% Save base variable value
b=cellb{1};
RHS=cellRHS{1};

% Set cells for z1,z2,RHS1,RHS2 (for FMG solver)
cellz1=cell(vcyclegrid,1);
cellz2=cell(vcyclegrid,1);
cellRHS1=cell(vcyclegrid,1);
cellRHS2=cell(vcyclegrid,1);
cellbnew=cell(vcyclegrid,1);

cellz1{1}=zeros(cellN{1},1);
cellz2{1}=cellz1{1};

% Loops steps
while (j<steps) && (abs(gamma(j))<500)

    % Tangent predictor
    v(:,j+1)=v(:,j)+ds*dv;
    gamma(j+1)=gamma(j)+ds*dgamma;
    
    % Update variable wrt parameter
    cellRHS{1}=gamma(j+1)*RHS;
    
    % New b(x) function in Newton
    bnew=cellb{1}+2*cellc{1}.*v(:,j+1);
    H=preconditioner(L,cellN{1},cella{1},bnew);
    
    % Correction loop using Newton iterations
    for i=1:max_newton_loops
        
        converged=false;
        
        % Initial RHS of linear equation
        F=(NA(v(:,j+1),cellk{1},cella{1},cellb{1},cellc{1})- ...
            cellRHS{1});
        F_gamma=-RHS;  
        
        % Check convergence condition
        if rms(F)<=1e-10%v(cellN{1}/2+1,j+1)-4/3*lambda(j+1)<=1e-10%
            fprintf('Converged after %d Newton Iterations step = %d\n',i-1,j)
            j=j+1;
            converged=true;
            
            % Optimal step length control for next step
            xi=N_opt/(i-1);
            
            if xi<0.5
                ds=ds*0.5;
            elseif xi>=0.5 && xi<=2
                ds=ds*xi;
            elseif xi>2
                ds=ds*2;
            end
            
            % Max step size check
            if ds>ds_max
                ds=ds_max;
            elseif ds<ds_min
                fprintf('Minimum step length of %f exceeded, breaking loop\n',ds_min)
                return
            end
            fprintf('New step size to %f\n',ds)
            
            break % End if converged to sufficient accuracy
        end
        
        % Arclength equation
        P=psi1*dot(dv,(v(:,j+1)-v(:,j)))+psi2*dgamma*(gamma(j+1)-gamma(j))-ds;
        
        % Initial RHS of linear equation
        RHS1=-F;
        RHS2=F_gamma;
        
        % Set cell for bnew and step downs
        cellbnew=setcellsNewton(vcyclegrid,cellbnew,bnew);
        cellH=setcellsH(vcyclegrid,H,L,cellN,cella,cellbnew);
        
        % Update cell RHS with new values
        [cellRHS1,cellRHS2]=setcellspseudo(vcyclegrid,cellN,cellRHS1, ...
            cellRHS2,RHS1,RHS2); 
        
        % FMG to solve F_x*z1=-F
%         [cellz1,r1]=FMG_pre(pre_relaxation,vcyclegrid,cellN,cellk,1,Nd,Nu,cella,cellbnew,cellRHS1,cellH,cellz1);
        [cellz1{1},r1]=vcycle_pre(pre_relaxation,1,vcyclegrid,cellN,cellk,1,Nd,Nu,cella,cellbnew,cellRHS1,cellH,cellz1);
%         cellz1{1}=LA(eye(cellN{1},cellN{1}),cellk{1},cella{1},cellbnew{1})\cellRHS1{1}; % Exact solver
%         [cellz1{1},r1]=Pcg(zeros(cellN{1},1),cellk{1},cella{1},cellbnew{1},cellRHS1{1},cellH{1}); % Pcg
%         r1=cellRHS1{1}-LA(cellz1{1},cellk{1},cella{1},cellbnew{1});
%           disp(rms(r1))
        
        % FMG to solve F_x*z2=F_lambda
%         [cellz2,r2]=FMG_pre(pre_relaxation,vcyclegrid,cellN,cellk,1,Nd,Nu,cella,cellbnew,cellRHS2,cellH,cellz2);
        [cellz2{1},r2]=vcycle_pre(pre_relaxation,1,vcyclegrid,cellN,cellk,1,Nd,Nu,cella,cellbnew,cellRHS2,cellH,cellz2);
%         cellz2{1}=LA(eye(cellN{1},cellN{1}),cellk{1},cella{1},cellbnew{1})\cellRHS2{1}; % Exact solver
%         [cellz2{1},r2]=Pcg(zeros(cellN{1},1),cellk{1},cella{1},cellbnew{1},cellRHS2{1},cellH{1}); % Pcg
%         r2=cellRHS2{1}-LA(cellz2{1},cellk{1},cella{1},cellbnew{1});
%           disp(rms(r2))

        % Solving for delta_v and delta_lambda (see nonlinear oceanography)
        delta_lambda=(ds-dot(dv,(v(:,j+1)-v(:,j)))-dgamma*(gamma(j+1)- ...
            gamma(j))-dot(dv,cellz1{1}))/(dgamma-dot(dv,cellz2{1}));
        delta_v=cellz1{1}-delta_lambda*cellz2{1};    

        % Update variable for next Newton iteration
        v(:,j+1)=v(:,j+1)+delta_v;
        gamma(j+1)=gamma(j+1)+delta_lambda;
        
        % Update variable wrt parameter
        cellRHS{1}=gamma(j+1)*RHS;
        
        % Update vector b for next iteration
        bnew=cellb{1}+2*cellc{1}.*v(:,j+1);
        H=preconditioner(L,cellN{1},cella{1},cellbnew{1});
        
    end

    if i==20
        fprintf('Did not converge to required tolerance after %d Newton Iterations at step %d\n',i,j)
        % Do not save latest vectors if not converged
        v(:,j+1)=[];
        gamma(:,j+1)=[];
        return % end function
    end

    if i==max_newton_loops && converged~=true
        
        fprintf('Did not converge to required tolerance after %d Newton Iterations at step %d\n',i,j)
        % Do not save latest vectors if not converged
        v(:,j+1)=[];
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

        % Solving for dlambda and dv
        dgamma1=1/(dgamma-dot(dv,cellz2{1}));
        dv1=-dgamma*cellz2{1};
    
        % Normalise
        mag=sqrt(psi1*dot(dv1,dv1)+psi2*dgamma1^2);
        dgamma1=dgamma1/mag;
        dv1=dv1/mag;
        
        if psi1*dot(dv,dv1)+psi2*dgamma*dgamma1<0
            fprintf('Sign change detected!\n')
            disp(psi1*dot(dv,dv1)+psi2*dgamma*dgamma1)
            dgamma1=-dgamma1;
            dv1=-dv1;
            
            % Normalise
            mag=sqrt(psi1*dot(dv1,dv1)+psi2*dgamma1^2);
            dgamma1=dgamma1/mag;
            dv1=dv1/mag;
            disp(psi1*dot(dv,dv1)+psi2*dgamma*dgamma1)
        end
             
%         disp(psi1*dot(dv,dv1)+psi2*dlambda*dlambda1+psi3*dgamma*dgamma1)
        disp(psi1*dot(dv,dv1)+psi2*dgamma*dgamma1);
        dgamma=dgamma1;
        dv=dv1;
       
    end
    
end