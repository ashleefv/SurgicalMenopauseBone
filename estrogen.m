
% time dependent estrogen function
function e=estrogen(t,if_surgical,t_m,tau_e,k_deg,k_syn)

% t_m is the onset point of menopause

% if_surgical is boolean which is 0 for natural menopause, 1 for surgical
% menopause

% for natural menopause:
% tau_e is the decay rate of estrogen
% surgical menopause:
% rapid decline between e=1 and e=eovx

if if_surgical
    %  e=1.*(t<t_m)+eovx*(t>=t_m); % step function definition using embedded if statements (t<t_m)
    % e=1.*(t<t_m)+(t>=t_m).*((1-k_syn/k_deg).*exp(-k_deg.*(t-t_m))+k_syn/k_deg) % step function definition using embedded if statements (t<t_m)
    if t>=t_m
        e=((1-k_syn/k_deg).*exp(-k_deg.*(t-t_m))+k_syn/k_deg);
    else
        e=1;
    end
else
    e=1.*(t<t_m)+(t>=t_m)./(1+ (t-t_m)/tau_e);
end

end