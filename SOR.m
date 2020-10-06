function [x,iter] = SOR(A,b,x0,omega,tol)
% coded by Soo Xi Wen for ENG3456 lab 1
% 10th march 2019
%sizing matrices
L = tril(A,-1);
D = diag(A);
Dia = 1./D;
U = triu(A,1);
w = omega;
%counter
iter=0;
m = 1;
%initializing storage and useful parameters
max_diff = tol+10;
n0 = size(x0);
n = n0(1);
% creating storage array
X(:,1)=x0;
X(:,2)=zeros(size(x0));
max_diff_array = zeros(n,1);
x = zeros(size(x0));

% looping thru iterations, until the tolerance is met
while(max_diff>tol)
    % iteration counter
     iter = iter+1;
        
        % looping thru each variable
        for(i=1:n)
            
            % formula for SOR
            % x(m+1) = x(m) + w[ (1/D)[b-Lx(m+1)-Ux(m)]-x(m)]
            LX = (L(i,:)*X(:,m+1));
            UX = (U(i,:)*X(:,m));
            b_LU = (b(i)-LX-UX);
            X(i,m+1)= X(i,m) + w*( Dia(i,1)*b_LU - X(i,m));
           
        end
            % using if func to detect any division by zero and substitute 0
            % by a small number 1*10^-9
            if (X(:,m+1) == 0)
                X(:,m+1) = 10^-9;
            end
            % calculating max diff to check for tolerance
            max_diff_array(:,1) = abs((X(:,m+1)-X(:,m))./X(:,m+1));
            max_diff = max(max_diff_array);
            
            % allocating values into memory space after calc max diff
            X(:,m) = X(:,m+1);
end

% output val
x(:,1) = X(:,end);

