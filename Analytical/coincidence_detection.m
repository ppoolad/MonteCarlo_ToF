function pdf = coincidence_detection(lambda_t,c_f,c_t, n_spad_per_pix,sim_res)
    %helper binomial function
    function prob = binomial(n,i,p)
        prob = double(nchoosek(n,i)) .* p.^i .* (1-p).^(n-i);
    end

cd_w = zeros(size(lambda_t));
cd_w(1:uint16(c_t/sim_res)) = 1;

%using convolution we can calculate probability of detection in cd windows
gamma_full = conv(1-exp(-lambda_t.*sim_res), cd_w);
%cut extra bits
prob_det = gamma_full(1:length(lambda_t));
%probability of detection
%p_det = lambda_t.*sim_res;
%p_cd = zeros(size(lambda_t));

%Probability of CD is to have at least c_f events in the preceding window
binoSum = zeros(size(lambda_t));
for c = c_f:n_spad_per_pix
    binoSum = binoSum + binomial(n_spad_per_pix,c, prob_det);
end


%each point of lambda_t (input rate) is scaled by the probability
%Also we have n i.i.d spads
pdf = n_spad_per_pix .* lambda_t .* (binoSum);


end
