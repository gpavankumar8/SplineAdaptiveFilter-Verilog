function op = overflow_underflow_det(A, fraclen, wordlen)

%     global maxval;
%     global maxvar;
% 
% if(maxval < max(abs(A)))
%     maxval = max(abs(A));
%     maxvar = inputname(1);
% end
op = 0;

A_fi = A(:).*2^(fraclen);
rnd_A_fi = nearest(A_fi);

% Underflow detection
% for i=1:length(A_fi)
%     if(A_fi(i) ~= 0 && rnd_A_fi(i) == 0)
%         disp("Underflow");
% %         disp(A(i));
% %         disp(i);
%         op = 1;
%     end
% end

A_fi = rnd_A_fi;
cmp_arr = 2^(wordlen-1)*ones(size(A(:)));

cmp_vals = A_fi > cmp_arr;
sum = cmp_vals'*cmp_vals;
if( sum >  0 )
    disp('Overflow!');
    op = 1;
end