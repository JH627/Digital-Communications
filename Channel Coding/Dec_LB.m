function mbit = Dec_LB(resmbit)

mbit = zeros(1, length(resmbit) / (5 / 2));
codeword = [0 0 0 0 0; 0 1 1 1 1; 1 0 1 0 0; 1 1 0 1 1];
origin = [0 0; 0 1; 1 0; 1 1];

for i = 1:5:length(resmbit)
    minDist = 6; % distance 초기화
    for j = 1:4
        xorResult = xor(resmbit(i:i+4), codeword(j, :)); % xor연산
        distance = sum(xorResult); % xor결과의 1의 개수를 더하여 distance 구함
        if distance <= minDist     % 만약 현재 최소거리보다 가까운 codeword를 발견할 경우
            if distance == minDist      % 현재의 최소 거리와 같다면
                if ((randi(2) - 1) == 1)  % 랜덤하게 0과 1중 뽑고 1이 나올 경우
                    continue;             % 기존값 유지
                end
            end
            minDist = distance;           % 가까운 값을 저장하는 로직
            idx = floor(i / 5) * 2 + 1;
            mbit(idx) = origin(j, 1);
            mbit(idx + 1) = origin(j, 2);
        end
    end
end