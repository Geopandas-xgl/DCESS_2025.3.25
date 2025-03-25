function total = input_total(time,rate)
idx_rate = ~isnan(rate);
total = round(trapz(time(idx_rate)*1000,rate(idx_rate)),0); % Pg C
end