function [mu,A,B,D,obj] = Alg_SRRR_P_CF(theta,lambda,xi,MaxItr,N,P,R,Q,S,Y,X,Z) 
    
    eps = 1e-5;
    One     = ones(N,1);
    OneT    = transpose(One);
    
    Lambda_Max = max(eig(X*transpose(X))); %eigs(X*transpose(X),1);
    
    %%
    rng(3);
    mu      = randn(P,1);
    D       = randn(P,S);
    A       = randn(P,R); 
    [U,~,V] = svd(A,'econ');   
    A       = U*transpose(V);
    B       = randn(Q,R);
    
    itr = 1;
    obj(itr) = EvalObj(mu,A,B,D,Y,X,Z,theta,lambda,xi,0);
        
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
        d = randn(Q,1); 
        for i = 1:Q
            bi_norm = norm(B(i,:),2);
            d(i) = theta/(bi_norm+theta)^2;
        end
        First_Term  = eye(Q)-X*transpose(X)/Lambda_Max; % can be further simplified
        Second_Term =  X*transpose(Y-mu*OneT-D*Z)/Lambda_Max;
        Qb    =  First_Term*B + Second_Term*A;
        B = Alg_GST_Prox(Qb, xi.*d, lambda/Lambda_Max);
        itr = itr+1;
        obj(itr)= EvalObj(mu,A,B,D,Y,X,Z,theta,lambda,xi,4); 
        
        %%
        if norm(mu-mu_old,'fro')<=eps && norm(A-A_old,'fro')<=eps &&...
                norm(D-D_old,'fro')<=eps && norm(B-B_old,'fro')<=eps &&...
                norm(obj(itr)-obj(itr-1),'fro')<=eps
            break
        end  
    end 
end
