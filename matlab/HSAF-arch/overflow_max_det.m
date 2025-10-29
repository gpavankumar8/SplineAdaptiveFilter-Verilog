function op = overflow_max_det(A, fraclen, wordlen)

    global maxval;
    global maxvar;

if(maxval < max(abs(A)))
    maxval = max(abs(A));
    maxvar = inputname(1);
end

A_fi = A(:).*2^(fraclen);
A_fi = round(A_fi);
cmp_arr = 2^(wordlen-1)*ones(size(A(:)));
op = 0;

cmp_vals = A_fi > cmp_arr;
sum = cmp_vals'*cmp_vals;
if( sum >  0 )
    disp('Overflow!');
    op = 1;
end