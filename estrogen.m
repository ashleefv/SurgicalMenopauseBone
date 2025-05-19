% % time dependent estrogen function
% function e=estrogen(t,if_surgical,t_e,tau_e,k_dec,k_syn)
% 
% % t_e is the onset point of menopause
% 
% % if_surgical is boolean which is 0 for natural menopause, 1 for surgical
% % menopause
% (t<t_e)
% % for natural menopause:
% % tau_e is the decay rate of estrogen
% % surgical menopause:is a step function between e=1 and e=eovx
% 
% if if_surgical==1
%     %  e=1.*(t<t_e)+eovx*(t>=t_e); % step function definition using embedded if statements (t<te)
% 
%     e=1.*(t<t_e)+(t>=t_e).*((1-k_syn/k_dec)*exp(-k_dec*(t-t_e))+k_syn/k_dec); % step function definition using embedded if statements (t<te)
% 
% else
%     'hre'
%     e=1.*(t<t_e)+(t>=t_e)/(1+ (t-t_e)/tau_e);
% end
% 
% end



% time dependent estrogen function
function e=estrogen(t,if_surgical,t_e,tau_e,k_dec,k_syn)

% t_e is the onset point of menopause

% if_surgical is boolean which is 0 for natural menopause, 1 for surgical
% menopause

% for natural menopause:
% tau_e is the decay rate of estrogen
% surgical menopause:is a step function between e=1 and e=eovx

if if_surgical
    %  e=1.*(t<t_e)+eovx*(t>=t_e); % step function definition using embedded if statements (t<te)
    % e=1.*(t<t_e)+(t>=t_e).*((1-k_syn/k_dec).*exp(-k_dec.*(t-t_e))+k_syn/k_dec) % step function definition using embedded if statements (t<te)
    if t>=t_e
        e=((1-k_syn/k_dec).*exp(-k_dec.*(t-t_e))+k_syn/k_dec);
    else
        e=1;
    end
else
    e=1.*(t<t_e)+(t>=t_e)./(1+ (t-t_e)/tau_e);
end

end