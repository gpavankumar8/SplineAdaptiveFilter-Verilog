function [error, s_vec, U_vec, j_vec, s] = HSAFLMS_filter(x, s_vec, U_vec, j_vec, d, w, q, par)
    
    P = par.P;      
    Q = par.Q;
    delX = par.delX ;
    C = par.C;
    u_vec = zeros(P+1,1);
    u_vec(P+1) = 1;

    u = x/delX - floor(x/delX);
    j = floor(x/delX) + (Q-1)/2;
%     if(j<0)
%         disp(x)
%         j=0;
%         u_vec = zeros(P+1,1);
%     elseif(j > Q)
%         disp(x);
%         j = 9;
%         u_vec = zeros(P+1,1);
%     else
        j_vec = [j, j_vec(1:end-1)];
    
        %u_vec = [u,1]';
        u_vec(P) = u;
        for p_ind = 2:P
%             u_vec = [u^p_ind; u_vec];
            u_vec(P-p_ind+1) = u^p_ind;
        end
%     end

    U_vec = [u_vec, U_vec(:,1:end-1)];

    s = u_vec'*C*q(j+1:j+P+1);
    s_vec = [s, s_vec(1:end-1)];
    y = [0, s_vec(1:length(w)-1)] * w';
%     y = [1,s_vec(1:length(w)-1)] * w';
    error = d - y;

    % Mask out U_vec values which are not related span j
    
%     U_vec_masked = zeros(size(U_vec));
% 
%     for i=1:length(U_vec)
% %         for k=1:P+1
%             if( j_vec(i) == j )
% %             if( abs(j_vec(i) - j) < 2 )
%                U_vec_masked(:,i) = U_vec(:,i);
%             end
% %         end
%     end
