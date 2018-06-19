function death = vessel_death(dose)
% compute the probability of death of a vessel for a given dose
% Plos One Potiron
% data: dose = 7, 11, 13 15 17 25 death = 0 0.18 0.23 0.26 0.43 0.62
if dose<0.22/0.0344
    death = 0;
else if dose >1.22/0.0344
        death = 1;
    else
        death = 0.0344*dose - 0.22;
    end
end

