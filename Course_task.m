s = struct;
s.n = 2;
s.r = 2;
s.N = 30;
s.A = zeros(s.n,s.n);
s.G = [3 1; 1 2];
B = [2 4; 3 1];
s.B = B - s.G;
s.H = -eye(s.n);
s.x_0 = [15; 1];
s.c_0 = [5 2];
s.T = 3;
s.h = s.T/s.N;
s.g = [0;0];
step = 0:s.h:s.T-s.h;

[X, U] = Solve(0,s.x_0,s);

figure(1);
plot(X(1,:),X(2,:));

figure(2);
for i = 1:s.r
    subplot(2,2,i);
    plot(step, U(i,:))
    title("The graphic of U_" + i + "(t)");
end
for i = 1:s.r
    subplot(2,2,i+s.r);
    plot(step, X(i,:))
       title("The graphic of X_" + i + "(t)");
end


function [X, U] = Solve(r,z,s)% пока только для позиции (0,x_0) т.к не исправлял CreateLp
    [A_lp, b_lp, c_lp] = CreateLp(s);
    lb = zeros(size(b_lp));
    u = linprog(-c_lp, A_lp, b_lp,[],[],lb);
    U = reshape(u,[s.n,s.N]);
    X = [];
    for t = 0:s.h:s.T-s.h
        count = 1;  
        f_res = zeros(s.n,1);
        for q = s.h:s.h:t
        f = quadv(@(x)FundMatrix(s.A, t - x) * s.B, q, q + s.h)* U(:,count);
        count = count + 1;
        f_res = f_res + f;
        end
        X_Solve = FundMatrix(s.A,t - r)*z + f_res;
        X = [X X_Solve];
    end
end

function [A_lp, b_lp, c_lp] = CreateLp(s)
    A_lp = [];
    b_lp = []; 
    c_lp = [];
    for t = 0:s.h:s.T-s.h
        c_lp = [c_lp s.c_0 * quadv(@(tt)FundMatrix(s.A, s.T - tt)*s.B, t, t + s.h)];
        b_lp = [b_lp; s.g - s.H*FundMatrix(s.A,t)*s.x_0];
        D = [];
        for r = s.h:s.h:t
            D_s = quadv(@(x) s.H*FundMatrix(s.A, t - x)*s.B, r, r + s.h);
            D = [D D_s];
        end
        D = [D s.G];
        A_lp = [A_lp zeros(size(A_lp,1),s.r)];
        A_lp = [A_lp;D];
    end
end
function F = FundMatrix(A,t)
    F = expm(A*t);
end