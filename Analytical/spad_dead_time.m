function pn = tdc_dd_pn(r_t,tdc_dead,n_tdc,n_pixels_to_tdc)
    t_dd = tdc_dead;
    n = n_tdc;
    n_pix = n_pixels_to_tdc;
    %%r_t_mean = mean(reshape(r_t,rat,size(r_t,2)/rat),1); <----------
    pn_nom = (n_pix.*r_t.*t_dd).^n ./ factorial(n);%% switch to r_t
%     pn_nom_cd = (n_pix.*r_t_cd.*t_dd).^n ./ factorial(n);
%     pn_spad_nom = (lambdas .* sim_res .* t_dd).^spad_per_pix ./ factorial(spad_per_pix);  
    pn_denom = zeros(size(pn_nom),'single');
%     pn_denom_cd = zeros(size(pn_nom),'single');
%     pn_spad_denom = zeros(size(pn_spad_nom));
    for m=1:1:n
        pn_denom = pn_denom + (n_pix.*r_t.*t_dd) .^m ./ factorial(m);
%         pn_denom_cd = pn_denom_cd + (n_pix.*r_t_cd.*t_dd) .^m ./ factorial(m);
        
    end
%     for s=1:1:spad_per_pix
%        pn_spad_denom = pn_spad_denom + (lambdas.*sim_res.*t_dd) .^s ./ factorial(s); 
%     end
%     pn_spad = pn_spad_nom ./ (pn_spad_denom+1);
    pn = (pn_nom./(pn_denom+1));
end