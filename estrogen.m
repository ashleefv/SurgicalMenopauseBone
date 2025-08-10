% % time dependent estrogen function
% function e=estrogen(t,if_surgical,t_m,tau_e,k_dec,k_syn)
% 
% % t_m is the onset point of menopause
% 
% % if_surgical is boolean which is 0 for natural menopause, 1 for surgical
% % menopause
% (t<t_m)
% % for natural menopause:
% % tau_e is the decay rate of estrogen
% % surgical menopause:is a step function between e=1 and e=eovx
% 
% if if_surgical==1
%     %  e=1.*(t<t_m)+eovx*(t>=t_m); % step function definition using embedded if statements (t<te)
% 
%     e=1.*(t<t_m)+(t>=t_m).*((1-k_syn/k_dec)*exp(-k_dec*(t-t_m))+k_syn/k_dec); % step function definition using embedded if statements (t<te)
% 
% else
%     'hre'
%     e=1.*(t<t_m)+(t>=t_m)/(1+ (t-t_m)/tau_e);
% end
% 
% end



% time dependent estrogen function
function e=estrogen(t,if_surgical,t_m,tau_e,k_dec,k_syn)

% t_m is the onset point of menopause

% if_surgical is boolean which is 0 for natural menopause, 1 for surgical
% menopause

% for natural menopause:
% tau_e is the decay rate of estrogen
% surgical menopause:is a step function between e=1 and e=eovx

if if_surgical
    %  e=1.*(t<t_m)+eovx*(t>=t_m); % step function definition using embedded if statements (t<te)
    % e=1.*(t<t_m)+(t>=t_m).*((1-k_syn/k_dec).*exp(-k_dec.*(t-t_m))+k_syn/k_dec) % step function definition using embedded if statements (t<t_m)
    if t>=t_m
        e=((1-k_syn/k_dec).*exp(-k_dec.*(t-t_m))+k_syn/k_dec);
    else
        e=1;
    end
else
    e=1.*(t<t_m)+(t>=t_m)./(1+ (t-t_m)/tau_e);
end

end