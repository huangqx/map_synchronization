function [X] = csm_sdp_opt(SDP, Para)
% this function solves the following optimization problem
% minimize    <SDP.C, X>
% subject to  A(X) ==b
%             X(SDP.ids_geq_0) >= 0;
%             X(SDP.ids_leq_1) <= 1;
%             X \succeq 0

% the variables to be optimized
dim = size(SDP.C, 1);
X = zeros(dim, dim);
S = zeros(dim, dim); % dual variable of X
z_leq_1 = zeros(length(SDP.ids_leq_1),1);

SDP.C = SDP.C/max(max(SDP.C));
mu = Para.mu_init;
for iter = 1:Para.nIterations
    TP = S + mu*X - SDP.C;
    TP(SDP.ids_leq_1) = TP(SDP.ids_leq_1) - z_leq_1;
    
    % optimize variable y
    y = SDP.A*reshape(TP, [dim*dim,1]) - mu*SDP.b;
    y = (SDP.A*SDP.A')\y;
    
    % optimize variable z_0
    z_geq_0 = -TP(SDP.ids_geq_0);
    ids = find(z_geq_0 < 0);
    z_geq_0(ids) = 0;
    
    % optimize variable z_1
    TP = S + mu*X - SDP.C - reshape(SDP.A'*y, [dim,dim]);
    z_leq_1 = TP(SDP.ids_leq_1) - mu*ones(length(SDP.ids_leq_1),1);
    ids = find(z_leq_1 < 0);
    z_leq_1(ids) = 0;
    
    
    % optimize S
    TP = reshape(SDP.A'*y, [dim, dim]);
    TP(SDP.ids_geq_0) = TP(SDP.ids_geq_0) - z_geq_0;
    TP(SDP.ids_leq_1) = TP(SDP.ids_leq_1) + z_leq_1;
    
    TP =(TP + TP') - diag(diag(TP));
    TP = TP + SDP.C -mu*X;
    TP = (TP+TP')/2;
    [U,V] = eig(TP);
    v = diag(V);
    ids = find(v > 0);
    S = U(:,ids)*V(ids,ids)*U(:,ids)';
    X = (S - TP)/mu;
    X = (X+X')/2;
    if mod(iter, 10) == 0
        r = norm(SDP.A*reshape(X, [dim*dim,1]) - SDP.b);
        fprintf(' iter = %d, e = %f, mu = %.3f, norm(A(X)-b) = %.4f\n', iter, sum(sum(SDP.C.*X)), mu, r);
    end
    mu = mu*Para.rho;
end
X = X(2:dim, 2:dim);
