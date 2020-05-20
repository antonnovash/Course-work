s = struct;
s.n = 2;
s.r = 2;
s.N = 100;
s.A = zeros(s.n,s.n);
s.G = [3 1; 1 2];
B = [2 4; 3 1];
s.B = B - s.G;
s.H = -eye(s.n);
s.M = [0.5 0; 0 0.3];
s.x_0 = [15; 1];
s.c_0 = [5 2];
s.T = 3;
s.h = s.T/s.N;
s.g = [0;0];
step = 0:s.h:s.T-s.h;

X_sol = [];
U_sol = [];
z = s.x_0;
for r=0:s.h:s.T-s.h
    [X_p, U_p] = Solve(s,r,z);
    X_sol = [X_sol X_p(:,1)];
    U_sol = [U_sol U_p(:,1)];
    [a,b] = size(X_p);
    if(b > 1)
    z = X_p(:,2);
    end
end
[X, U] = Solve(s,0,s.x_0);
figure(1);
p = plot(X(1,:),X(2,:),'K',X_sol(1,:),X_sol(2,:),'R');
p(1).LineWidth = 1.5;
p(2).LineWidth = 1.5;
xlabel('x_1');
ylabel('x_2');
grid on
figure(2);
for i = 1:s.r
    subplot(2,2,i);
    p = plot(step, U(i,:),'K', step, U_sol(i,:),'R');
    p(1).LineWidth = 1.3;
    p(2).LineWidth = 1.3;
    xlabel('t');
    ylabel("u_"+i);
end
for i = 1:s.r
    subplot(2,2,i+s.r);
    p = plot(step, X(i,:),'K',step, X_sol(i,:),'R');
    p(1).LineWidth = 1.3;
    p(2).LineWidth = 1.3;
    xlabel('t');
    ylabel("x_"+i);
end

function [X, U] = Solve(s,r,z)
    [A_lp, b_lp, c_lp] = CreateLp(s,r,z);
    lb = zeros(size(b_lp));
    u = linprog(-c_lp, A_lp, b_lp,[],[],lb);
    U = reshape(u,[s.n,length(A_lp)/2]);
    X = [];
    for t = r:s.h:s.T-s.h
        count = 1;  
        f_res = zeros(s.n,1);
        for q = r+s.h:s.h:t
        f = quadv(@(x)FundMatrix(s.A, t - x) * s.B, q, q + s.h)* U(:,count);
        count = count + 1;
        f_res = f_res + f;
        end
        if(r == 0)
        X_Solve = FundMatrix(s.A,t - r)*z + f_res;
        else
        w = -0.1*U(:,count);
        X_Solve = FundMatrix(s.A,t - r)*z + f_res + quadv(@(x)FundMatrix(s.A, r - x) * s.M * w, r, r + s.h);
        end
        X = [X X_Solve];
    end
end

function [A_lp, b_lp, c_lp] = CreateLp(s,r,z)
    A_lp = [];
    b_lp = []; 
    c_lp = [];
    for t = r:s.h:s.T-s.h
        c_lp = [c_lp s.c_0 * quadv(@(tt)FundMatrix(s.A, s.T - tt)*s.B, t, t + s.h)];
        b_lp = [b_lp; s.g - s.H*FundMatrix(s.A,t)*z];
        D = [];
        for q = r+s.h:s.h:t
            D_s = quadv(@(x) s.H*FundMatrix(s.A, t - x)*s.B, q, q + s.h);
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