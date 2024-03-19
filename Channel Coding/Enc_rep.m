function resmbit = Enc_rep(mbit)

resmbit = zeros(1, 3*length(mbit));
for i = 1:length(mbit)
    if mbit(i) ==1
        resmbit(3*i-2:3*i) = 1;
    else
        resmbit(3*i-2:3*i) = 0;
    end
end