function [angleDistCurve] = gt_eval_abs(Rgt, R_opt)
% This function evaluates optimized rotations with the ground-truth
% rotations

n = size(Rgt, 3);

subset = [];
if 1
    A = zeros(3*n,3);
    B = zeros(3*n,3);
    j = 0;
    for i = 1:n
        if(sum(sum((Rgt(:,:,i)==0)))==9||any(any(isnan(Rgt(:,:,i)))));
            continue;
        end
        j = j + 1;
        ids = (3*j-2):(3*j);
        A(ids,:) = R_opt(:,:,i);
        B(ids,:) = Rgt(:,:,i);
        subset = [subset, i];
    end
    A = A(1:(3*j),:);
    B = B(1:(3*j),:);

    A = A';
    B = B';

    [R] = rot_fitting(A, B);
    A = R*A;
else
    A = zeros(3, 3*n);
    B = zeros(3, 3*n);
    for i = 1:n
        ids = (3*i-2):(3*i);
        A(:, ids) = R_init(:,:,i)';
        B(:, ids) = Rgt(:,:,i)';
    end
    rootId = 3592;
    ids = (3*rootId-2):(3*rootId);
    R = B(:,ids)*A(:,ids)';
    A = R*A;
    j = n;
end

errors = zeros(1, j);
for i = 1:j
    ids = (3*i-2):(3*i);
    R1 = A(:, ids);
    R2 = B(:, ids);
    if(sum(sum((Rgt(:,:,i)==0)))==9||any(any(isnan(Rgt(:,:,i)))));
        errors(i) = 0;
    else
        errors(i) =...
            acos(max(min((R1(1,:)*R2(1,:)'+R1(2,:)*R2(2,:)'+R1(3,:)*R2(3,:)'-1)/2,1),-1));
    end
end
angleDevs = errors*360/2/pi;

stats = zeros(1, 181);
errors = zeros(1, 180);
for i = 1:181
    stats(i) = length(find(angleDevs <= i))/j;
    if i > 1
        errors(i-1) = stats(i) - stats(i-1);
    end
end
stats = stats(1:180);
angleDistCurve = stats;


function [R] = rot_fitting(A, B)

S = A*B';

N = zeros(4,4);
N(1,1) = S(1,1) + S(2,2) + S(3,3);
N(1,2) = S(2,3) - S(3,2);
N(1,3) = S(3,1) - S(1,3);
N(1,4) = S(1,2) - S(2,1);
N(2,1) = N(1,2);
N(3,1) = N(1,3);
N(4,1) = N(1,4);

N(2,2) = S(1,1) - S(2,2) - S(3,3);
N(2,3) = S(1,2) + S(2,1);
N(2,4) = S(3,1) + S(1,3);
N(3,2) = N(2,3);
N(4,2) = N(2,4);

N(3,3) = S(2,2) - S(1,1) - S(3,3);
N(3,4) = S(2,3) + S(3,2);
N(4,3) = N(3,4);

N(4,4) = S(3,3) - S(1,1) - S(2,2);

[u,v] = eig(N);

q = u(:,4);
q0 = q(1);
n = q(2:4);
cor = [0, -n(3), n(2);
    n(3), 0, -n(1);
    -n(2), n(1), 0];

R = (q0*q0-n'*n)*eye(3) + 2*n*n' + 2*q0*cor;
