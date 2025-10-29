function op = overflow_det(A, fraclen, wrdlen)

    A_fi = A(:).*2^(fraclen);
    A_fi = round(A_fi);
    
    max_val_arr = 2^(wrdlen-1)*ones(size(A(:)));   % (wrdlen-1) assuming MSB is sign bit
    op = 0;
    
    cmp_vals = abs(A_fi) > max_val_arr;
    cmp_true = cmp_vals'*cmp_vals;
    
    if( cmp_true >  0 )
        disp('Overflow!');
        op = 1;
    end