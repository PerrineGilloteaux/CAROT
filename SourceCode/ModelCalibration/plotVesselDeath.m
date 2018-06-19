% Figure 4 A Probability of vessel survival
figure
plot(0:35,1-arrayfun(@vessel_death,0:35));
xlabel('Dose in Gy');
ylabel ('1-vesseldeathproba');
title('Probability of Vessel Survival with one exposure');
