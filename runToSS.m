function [T,X,SUCCESS] = runToSS(model,dt,X0,maxIt,tol,parms)
%Runs ode15s in steps of 100*dt untill the model reaches steady state
%(|dxdt|<tol) 
%SUCCESS = bool: true if model converged

%model = @systemFunction
%dt = timestep; 1
%X0 = initial values; [10 10]
%maxIt = maximum iterations = timesteps*100; 1e2
%tol = tolerance for steady state; 1e-7
%parms = parameter set that modelFunction takes in {}; {ax,ay,b,c,d,eps,mu,tau}

warning('off','MATLAB:ode15s:IntegrationTolNotMet');

T = 0;
X = X0;
err = Inf;
for i = 1:maxIt
    if err < tol
        SUCCESS = true;
        warning('on','MATLAB:ode15s:IntegrationTolNotMet');
        return
    end
    
    [t,x] = ode15s(model,dt:dt:100*dt,X(end,:),[],parms{:});
    err = norm(model(t(end),x(end,:),parms{:}),1);

    T = [T;t+T(end)];
    X = [X;x];
end

SUCCESS = false;
if nargout < 3       %only complains if I didn't ask for flag 
    disp('function ran to maxIt without converging to <tol !')
end
warning('on','MATLAB:ode15s:IntegrationTolNotMet');