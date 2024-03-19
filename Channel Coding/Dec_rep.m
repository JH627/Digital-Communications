function mbit = Dec_rep(resmbit)

mbit = zeros(1, length(resmbit)/3);
for i = 1:length(mbit)
    if sum(resmbit(3*i-2:3*i))>1.5
        mbit(i) = 1;
    else
        mbit(i) = 0;
    end
end