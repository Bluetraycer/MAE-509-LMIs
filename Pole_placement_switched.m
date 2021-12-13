As(:,:,1) = [0.4 -0.3;-0.5 0.5];
As(:,:,2) = [-0.3 -0.1;-0.2 -0.6];
Bs(:,:,1) = [0.1; 0.2];
Bs(:,:,2) = [-0.2; -0.1];

num_modes = size(As,3);
nx = size(As,1);
nu = size(Bs,2);

gamma = sdpvar(1);
P = sdpvar(nx);
Z = sdpvar(nu,nx,'full');
constraints = [];


alpha = -2;

for i = 1:num_modes
        A = As(:,:,i);
        B = Bs(:,:,i);
        
        constraints = P >= 0.001*eye(nx);
        constraints = [constraints, (A*P + B*Z)' + (A*P + B*Z) + 2*alpha*P <= 0];

end

options = sdpsettings('solver','mosek','verbose',0);
solution = optimize(constraints,[],options);

ZZ = value(Z);
YY = value(P);

K = ZZ/YY;