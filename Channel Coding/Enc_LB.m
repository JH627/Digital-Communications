function resmbit = Enc_LB(mbit)

resmbit = zeros(1, (5 / 2) * length(mbit));

for i = 1:2:length(mbit)
    resmbit((i - 1) * 5 / 2 + 1) = mbit(i); % c1 = x1
    resmbit((i - 1) * 5 / 2 + 2) = mbit(i + 1); % c2 = x2
    resmbit((i - 1) * 5 / 2 + 3) = xor(mbit(i), mbit(i + 1)); % c3 = xor(x1, x2)
    resmbit((i - 1) * 5 / 2 + 4) = mbit(i + 1); % c4 = x2
    resmbit((i - 1) * 5 / 2 + 5) = mbit(i + 1); % c5 = x2
end