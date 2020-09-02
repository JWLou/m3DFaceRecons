function x = solverPCG(A, b, x0, maxit)
% Apply Preconditioned Conjugate Gradient method to obtain a numerical solution
% of the parameter update through an iterative optimization process

%% PCG process
diagA = diag(A); invM = diag(1./diagA); % Jacobi preconditioning
r0 = b-A*x0; z0 = invM*r0; p0 = z0;
rk = r0; pk = p0; zk = z0;
xk = x0;
for k = 0:maxit-1
    alphak = (rk'*zk)/(pk'*A*pk);
    xk1 = xk+alphak*pk;
    rk1 = rk-alphak*A*pk;
    zk1 = invM*rk1;
    betak = (zk1'*rk1)/(zk'*rk);
    pk1 = zk1+betak*pk;  
    rk = rk1; zk = zk1; pk = pk1; xk = xk1;
end
x = xk;

end