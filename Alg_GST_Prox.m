function [X] = Alg_GST_Prox(V, xi, lambda)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    Q = size(V,1);
    R = size(V,2);

    h = zeros(Q,1);
    X = zeros(Q,R);
    for  i =1:Q
        h(i) = norm(V(i,:),2);
        if h(i) > lambda*xi(i)
            X(i,:) = V(i,:)*(1-lambda*xi(i)/h(i));
        else
            X(i,:) = 0;
        end
    end
end

