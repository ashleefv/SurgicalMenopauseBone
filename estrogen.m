
% time dependent estrogen function
function E=estrogen(t,if_surgical,t_m,tau_E,kappa_E,k_syn)

% t_m is the onset point of menopause

% if_surgical is boolean which is 0 for natural menopause, 1 for surgical
% menopause

% for natural menopause:
% tau_e is the decay rate of estrogen
% surgical menopause:
% rapid decline between E=1 and E=Eovx

if if_surgical
    %  E=1.*(t<t_m)+Eovx*(t>=t_m); % step function definition using embedded if statements (t<t_m)
    % E=1.*(t<t_m)+(t>=t_m).*((1-k_syn/kappa_E).*exp(-kappa_E.*(t-t_m))+k_syn/kappa_E) % step function definition using embedded if statements (t<t_m)
    if t>=t_m
        E=((1-k_syn/kappa_E).*exp(-kappa_E.*(t-t_m))+k_syn/kappa_E);
    else
        E=1;
    end
else
    E=1.*(t<t_m)+(t>=t_m)./(1+ (t-t_m)/tau_E);
end

end