function pdf = coincidence_detection(lambda_t,c_f,c_t, n_spad_per_pix,sim_res)

    function prob = binomial(n,i,p)
        prob = nchoosek(n,i) .* p.^i.* (1-p).^(n-i);
    end

cd_w = zeros(size(lambda_t));
cd_w(1:uint16(c_t/sim_res)) = 1;

gamma_full = conv(1-exp(-lambda_t), cd_w).*sim_res;
prob_det = gamma_full(1:length(lambda_t));
p_det = lambda_t.*sim_res;
p_cd = zeros(size(lambda_t));
for i=1:length(lambda_t)
    p_cd(i) =  binomial(n_spad_per_pix-1,c_f-1, prob_det(i)) .* (1-p_det(i)).^(n_spad_per_pix-(c_f-1));
    if p_cd(i) < 1e-140
        p_cd(i) = 0;
    end
end

pdf = lambda_t .* p_cd;

end
