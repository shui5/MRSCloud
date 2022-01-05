% Apply coherence filter
function d_out = sim_apply_pfilter(d_in,H,p)
for n=1:length(H)
    d_out{n} = d_in{n} .* (H(n).P_matrix==p);
end

end