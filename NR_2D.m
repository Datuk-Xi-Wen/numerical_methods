function [x_root,y_root,iter] = NR_2D(x0,y0,f,g,dfdx,dfdy,dgdx,dgdy,tol)
% coded by Soo Xi Wen for ENG3456 comp lab 1 
% 13th march 2019
% input functions of f,g, and all their respective derivatives are assumed 
% to be in the form of F(x,y)

% initial val for iterations
m=1;
iter = 0;
conv1 = tol+1;
conv2 = conv1;
% xy iteration array
xy = [x0,0;y0,0];

if(iter < 1000)
   % calc loop
   while ((conv1 > tol) || (conv2  > tol))
        % counting iterations
        %(abs(conv1) > tol )|| (abs(conv2)  > tol)   
        iter = iter + 1;
        
        % NR-2D eqn calc
        fg_diff_i = [dfdx( xy(1,m) ),dfdy( xy(2,m) ) ; dgdx( xy(1,m) ),dgdy( xy(2,m) )];
        
        fg_i = [f( xy(1,m),xy(2,m) );g( xy(1,m),xy(2,m) )];
        xy(:,m+1) = xy(:,m) - ( inv(fg_diff_i)*fg_i );
        
        % allocating value into prev space
        xy(:,m) = xy(:,m+1);
        
        % convergence test for new values
        conv1 = abs(f(xy(1,m),xy(2,m)));
        conv2 = abs(g(xy(1,m),xy(2,m)));
   end
else
    [x_root,y_root] = ["Error", "No solution found"];
    fprintf("No solution found after 1000 cycles, this method does not converge")
end
% output values
x_root = xy(1,1);
y_root = xy(2,1);
