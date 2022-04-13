function pdf = analytical_lidar(lambda_t,c_f,c_t, n_spad_per_pix,sim_res,adjust_cd)
%ANALYTICAL_LIDAR Calculates an approximate histogram analytically and very
% quickly
%   This functions gets calculated photon rate function \lambda(t) and
%   models the final histogram of the system
%   Args:
%   lambda_t: Photon rate function calculated based on the system metrics
%   and distance
%   c_f: Coincidence factor
%   c_t: coincidence Window
%   n_spad_per_pix: number of spads in a macro pixel
cd_w = zeros(size(lambda_t));
cd_w(1:uint16(adjust_cd.*c_t/sim_res)) = 1;
%cd_w = ones(1,uint16(adjust_cd.*c_t/sim_res));

%cd_w(1:uint16(c_t/sim_res)) = 1;
gamma_full = conv(n_spad_per_pix.* lambda_t, cd_w).*sim_res;
gamma = gamma_full(1:length(lambda_t));
inside = zeros(size(gamma));
for n=0:c_f-1
    temp  = (gamma .^ n) ./ factorial(n);
    inside = inside + temp;
end

outside = exp(-gamma);

cdf = 1 - outside .* inside; %<- should be multiplied by \lambda
pdf = zeros(size(cdf));
pdf = cdf;
%pdf(2:end) = diff(cdf); 
pdf = normalize(pdf,'norm',1);
pdf = lambda_t .* pdf;
end