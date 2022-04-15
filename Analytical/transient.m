function p_det_cut = transient(lambda_t,t_dead,sim_res)
T_dead = t_dead / sim_res; 
lambda_t_twice = [lambda_t, lambda_t, lambda_t, lambda_t];
p_0 = -lambda_t_twice(1) .* sim_res;
p_det = zeros(size(lambda_t_twice));
p_det(1) = p_0;
for i = 1:length(lambda_t_twice)
    p_ev = lambda_t_twice(i) .* sim_res;
    %dead_window = ones(1,uint16(T_dead));
    p_dead_window = sum(p_det(max(1,i-uint32(T_dead)):i)); %conv(P_ev,dead_window) .* sim_res;
    p_no_ev = 1 - p_dead_window;
    
    p_det(i) = (p_ev .* p_no_ev);
    
end
p_det_cut = p_det(3*length(lambda_t)+1:end);
%Let's do another one
% for i = 1:length(lambda_t)
%     p_ev = 1 - exp(-lambda_t(i) .* sim_res);
%     %dead_window = ones(1,uint16(T_dead));
%     det_index_start = i-uint16(T_dead);
%     if(det_index_start<1)
%         det_index_start = length(lambda_t) - det_index_start;
%         p_det_prev_it = sum(p_det(det_index_start:end)) + sum(p_det(1:i));
%     else
%         p_det_prev_it = sum(p_det(max(1,i-uint16(T_dead)):i));
%     end
%     p_dead_window = p_det_prev_it; %conv(P_ev,dead_window) .* sim_res;
%     p_no_ev = 1 - p_dead_window;
%     
%     p_det(i) = (p_ev .* p_no_ev);
%     
% end

end