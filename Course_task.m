s = struct;
s.A = [0 0; 0 0];
s.G = [3 1; 1 2];
B = [2 4; 3 1];
s.B = B - s.G;
s.H = [-1 0;0 -1];
s.x_0 = [2;3];
s.c_0 = [5;6];
s.T = 5;
s.h = 0.5;
s.n = 2;
s.g = [0;0];

% A = [3 1; 1 2];
% B = [2 4; 3 1];
% B = B - A;
% H = [-1 0;0 -1];
% x_0 = [2;3];
% c_0 = [5;6];
% T = 5;
% h = 0.5;
% n = 2;
% g = [0;0];

step = zeros(1,T/h);
u_1 = zeros(size(step));
u_2 = zeros(size(step));
i = 1;
for t = 0:h:T-h
    step(i) = t;
    i = i + 1;
end
[X, U] = Solve(0,s.x_0,s);

U_1 = U(1,:);
U_2 = U(2,:);
X_1 = X(1,:);
X_2 = X(2,:);

plot(X_1,X_2);
subplot(2,2,1);
title('The graphic of U_1(t)');
plot(step, U_1)

subplot(2,2,2); 
title('The graphic of U_2(t)');
plot(step, U_2)

subplot(2,2,3);
title('The graphic of X_1(t)');
plot(step, X_1)

subplot(2,2,4); 
title('The graphic of X_2(t)');
plot(step, X_2)


function [X, U] = Solve(r,z,s)% пока только для позиции (0,x_0) т.к не исправлял CreateLp
    [A_lp, b_lp, c_lp] = CreateLp(s.A, s.B, s.H, s.G, s.T, s.h, s.g, z, s.c_0);
    lb = zeros(size(c_lp));
    c_lp = c_lp.';
    u = linprog(-c_lp, A_lp, b_lp,[],[],lb);
    j = 1;
    u_transp = u.';
    for i = 1:2:size(u)
    u_1(j) = u_transp(i);
    u_2(j) = u_transp(i+1);
    j = j + 1;
    end
    U = cat(1,u_1,u_2); 
    %находим X
    X = zeros(s.n,1);
    h = s.h;
    T = s.T;
    for t = r:h:T-h
        count = t/h;
        f = quadv(@(x)FundMatrix(s.A, t - x) * s.B * U(:,count+1), r, t);
        X_Solve = FundMatrix(s.A,t - r)*z + f;
        if(t == r)
            X = X_Solve;
        else
            X = cat(2, X, X_Solve);
        end
    end
end

function [A_lp, b_lp, c_lp] = CreateLp(A, B, H, G, T, h, g, x_0, c_0)
    % t = 0
    Z = zeros(size(G));
    A_lp = G;
    b_lp = g - H*x_0; 
    Z_n = Z;
    c_lp = quadv(@(x)FundMatrix(A, T - x)*B, 0, h);
    c_lp = c_lp * c_0;
    % t\in[h;T-h]
    for t = h:h:T-h
        c = quadv(@(x)FundMatrix(A, T - x)*B, t, t + h);
        c = c * c_0;
        D = quadv(@(x) H*FundMatrix(A, t - x)*B, 0, h);
        b = g - H*FundMatrix(A, t)*x_0;
        A_lp = cat(2, A_lp, Z);
        for s = h:h:t
            D_s = quadv(@(x) H*FundMatrix(A, t - x)*B, s, s + h);
            if(s == t)
                D = cat(2, D, G);
            else
                D = cat(2, D, D_s);
            end
        end
        A_lp = cat(1, A_lp, D);
        Z = cat(1, Z, Z_n);
        b_lp = cat(1, b_lp, b);
        c_lp = cat(1, c_lp, c);
    end
end
function F = FundMatrix(A,t)
    F = expm(A*t);
end

