function [mu,A,B,D,obj] = Alg_SRRR_C(theta,gamma,xi,MaxItr,N,P,R,Q,S,Y,X,Z) 
  
    eps = 1e-5;
    One     = ones(N,1);
    OneT    = transpose(One);

    alpha  = 1/max(eig(X*transpose(X)))-eps;
    MaxItr_B = 1e3;
    
    lambda = 0;
    
    %%
    rng(3);
    mu      = randn(P,1);
    D       = randn(P,S);
    A       = randn(P,R);
    [U,~,V] = svd(A,'econ');   
    A       = U*transpose(V);
    B       = randn(Q,R);
    Reg_B = 0;
    for i = 1:Q
        bi_norm = norm(B(i,:),2);
        Reg_B = Reg_B+xi*bi_norm/(bi_norm+theta);
    end
    while Reg_B > gamma
        B = 1/2*B;
        Reg_B = 0;
        for i = 1:Q
            bi_norm = norm(B(i,:),2);
            Reg_B = Reg_B+xi*bi_norm/(bi_norm+theta);
        end
    end   
    
    itr = 1;
    obj(itr)= EvalObj(mu,A,B,D,Y,X,Z,theta,lambda,xi,0); 
  
    while itr <= MaxItr
        mu_old = mu;
        A_old = A;
        B_old = B;
        D_old = D;
        
        %% mu Step 
        mu       = (1/N)*(Y-A*transpose(B)*X-D*Z)*One;
        itr      = itr+1;
        obj(itr) = EvalObj(mu,A,B,D,Y,X,Z,theta,lambda,xi,1); 
        
        %% D Step
        D        = (Y-mu*OneT-A*transpose(B)*X)*transpose(Z)/(Z*transpose(Z));
        itr      = itr+1;
        obj(itr) = EvalObj(mu,A,B,D,Y,X,Z,theta,lambda,xi,2);        

        %% A Step
        Qa = (Y-mu*OneT-D*Z)*transpose(X)*B; 
        [U,~,V] = svd(Qa,'econ');  
        A = U*transpose(V);
        itr = itr+1;
        obj(itr) = EvalObj(mu,A,B,D,Y,X,Z,theta,lambda,xi,3);  

        %% B Step
        d = zeros(Q,1); 
        c = 0;
        for i = 1:Q
            bi_norm = norm(B(i,:),2);
            d(i) = theta/(bi_norm+theta)^2;
            c    = c + bi_norm/(bi_norm+theta)-d(i)*bi_norm;
        end

        for itr_B = 1:MaxItr_B
                B_in_old = B;
                G =  X*transpose(X)*B - X*transpose(Y-mu*OneT-D*Z)*A; 
                Qb =  B - alpha*G;
                B = Alg_GLC_Proj(Qb, xi.*d, gamma - c); 
                if norm(B-B_in_old,'fro')<=eps
                    itr = itr+1;
                    obj(itr) = EvalObj(mu,A,B,D,Y,X,Z,theta,lambda,xi,4);
                    break
                elseif itr_B == MaxItr_B
                    itr = itr+1;
                    obj(itr) = EvalObj(mu,A,B,D,Y,X,Z,theta,lambda,xi,4);
                    disp('Inner B loop attains MaxItr.')
                end                      
        end
        %%
        if norm(mu-mu_old,'fro')<=eps && norm(A-A_old,'fro')<=eps &&...
                norm(D-D_old,'fro')<=eps && norm(B-B_old,'fro')<=eps &&...
                norm(obj(itr)-obj(itr-1),'fro')<=eps
            break
        end
    end
