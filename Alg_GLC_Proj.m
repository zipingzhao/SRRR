function [X] = Alg_GLC_Proj(V, xi, gamma)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%     V = 10*randn(8,5);
%     xi = ones(8,1);
%     gamma = 3;

    Q = size(V,1);
    R = size(V,2);

    h = zeros(Q,1);
    for i = 1:Q
        h(i) = norm(V(i,:),2);
    end

    h_xi = h./xi;
    [h_xi_sort, ind_sort] = sort(h_xi,'descend');
    h_xi_sort(Q+1) = 0;
    xi_sort = xi(ind_sort);
    if sum(xi.*h) <= gamma
        eta = 0;
    else               
        for I = 1:Q
            sum_h = -gamma;
            for i = 1:I % if use I-1 then an if-else clasue should be added before
                sum_h = sum_h + (xi_sort(i))^2*(h_xi_sort(i)-h_xi_sort(I));
            end
            if sum_h < 0
                I_star = I;
                continue
            else                 
                break
            end
        end
        sum_h = -gamma;
        for i = 1:I_star+1 - 1
            sum_h = sum_h + (xi_sort(i))^2*(h_xi_sort(i)-h_xi_sort(I_star+1));
        end
        if norm(sum_h) <= 1e-8
            eta = h_xi_sort(I_star+1);
        end
        lwb = h_xi_sort(I_star+1);
        upb = h_xi_sort(I_star);
        while norm(sum_h) > 1e-8
            eta = 1/2*(lwb + upb);
            sum_h = -gamma;
            for i = 1:I_star
                sum_h = sum_h + (xi_sort(i))^2*(h_xi_sort(i)-eta);
            end
            if sum_h >= 0
                lwb = eta;
            else
                upb = eta;
            end           
        end
    end
    X = zeros(Q, R);
    for i = 1:Q
        if h(i) <= eta*xi(i)
            X(i,:) = 0;
        else
            X(i,:) = V(i,:)*(1-eta*xi(i)/h(i));
        end
    end
% X
    %%
%     cvx_begin quiet
%         variable X_cvx(Q,R)
%         Reg = 0;
%         for i = 1:Q
%             Reg = Reg + xi(i)*norm(X_cvx(i,:));     
%         end 
%         minimize(.5*sum(sum_square_abs(X_cvx-V)));
%         subject to
%             Reg <= gamma
%     cvx_end
%     X_cvx
%     X = X_cvx;  
end



