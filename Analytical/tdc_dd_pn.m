function pn = tdc_dd_pn(r_t,tdc_dead,n_tdc,n_pixels_to_tdc)
    t_dd = tdc_dead;
    n = n_tdc;
    n_pix = n_pixels_to_tdc;
    %%r_t_mean = mean(reshape(r_t,rat,size(r_t,2)/rat),1); <----------
    pn_nom = (n_pix.*r_t.*t_dd).^n ./ factorial(n);%% switch to r_t
 
    pn_denom = zeros(size(pn_nom),'single');

    for m=1:1:n
        pn_denom = pn_denom + (n_pix.*r_t.*t_dd) .^m ./ factorial(m);
%         pn_denom_cd = pn_denom_cd + (n_pix.*r_t_cd.*t_dd) .^m ./ factorial(m);
        
    end

    pn = (pn_nom./(pn_denom+1));
end