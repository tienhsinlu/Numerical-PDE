function p = cg(A,b,tol,k_max)
% Conjugate Gradient Method
    p = b;
    r = b - A*p;
    if norm(r) < tol
        return
    end
    y = r;
    w = A*y;
    s = y'*w;
    alpha = (r'*y)/s;
    p = p + alpha*y;
    for k = 1:k_max
       r = r - alpha*w;
       if(norm(r) < tol )
            return;
       end
       beta = (r'*w)/s;
       y = r + beta*y;
       w = A*y;
       s = y'*w;
       alpha = (r'*y)/s;
       p = p + alpha*y;
    end
 end