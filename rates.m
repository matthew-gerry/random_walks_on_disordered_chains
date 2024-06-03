function [k] = rates(sites,ga_a,ga_b,tau)
%Calculate the rate based on whether the site 
        if sites == 2
            k = tau^2/(ga_a);
        else
            k = tau^2/(ga_b);
        end

end