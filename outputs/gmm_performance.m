%diff_array_init
%diff_array = zeros(2,5*3);
%conf_array = zeros(2,5*3);
init = 5;
init_t = 110e-9;
sim_res = 10e-12;
for i=1:1:5
    pix_name = sprintf("pixel%d",i-1);
    gx = eval(pix_name+"_mu_convergance(end,:)");
    gmm_mus = eval(pix_name+"_mu_convergance(end,:)");
    gmm_pix = eval(pix_name+"_pi_convergance(end,:)");
    gmm_sigs = eval(pix_name+"_sig_convergance(end,:)");
    golden_mu = (init_t + 20e-9*(i-1))/sim_res;
    golden_mu2 = golden_mu + (20e-9)/sim_res;
    diff = abs(gmm_mus-golden_mu);
    diff2 = abs(golden_mu2 - gmm_mus);
    metric = gmm_pix ./ gmm_sigs;
    [out, index] = sort(metric,'descend');
    %[minimum_distance, min_idx] = min(diff);
    %[minimum_distance2, min_idx2] = min(diff2);
    minimum_distance=abs(gmm_mus(index(1))-golden_mu);
    minimum_distance2=abs(gmm_mus(index(2))-golden_mu2);
    %min_pi = gmm_pix(min_idx);
    min_sig = gmm_sigs(index(1));
    min_sig2 = gmm_sigs(index(2));
    diff_array(1,i+init) = minimum_distance;
    conf_array(1,i+init) = min_sig;
    diff_array(2,i+init) = minimum_distance2;
    conf_array(2,i+init) = min_sig2;
end

