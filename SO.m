function [BestScore, BestPop, conv_curve] = EAO(N, MaxIter, Lb, Ub,dim, fojb)
% rng(0)
% rng(1)
% rng(42)
% rng(88)
% --- Initialization ---
Pops = Lb + (Ub - Lb) .* rand(N, dim);
fitness = arrayfun(@(i) fojb(Pops(i, :)), 1:N)'; 
[BestScore, idx] = min(fitness);
BestPop = Pops(idx, :); 
conv_curve = zeros(1, MaxIter);

%main loop
for It = 1:MaxIter
    for i = 1:N
        if rand < rand * (It / MaxIter)
            Candidate = (Ub + Lb + 2 .* BestPop - 2 .* Pops(i,:)) / 2 + ...
                        rand * abs(randn) * (BestPop - Pops(i, :)); 
        else
            if rand > rand
                S = randperm(N, 2);
                while any(S == i), S = randperm(N, 2); end
                R = Pops(S(1), :) - Pops(S(2), :);
                [~, index] = sort(fitness);
                T=[Pops(index(1:3),:)]; C=randperm(3,1);
                Vector = (T(C, :) - Pops(i, :));
                BinaryArray = randi([0 1], 1, dim);
                Oppose = ~ BinaryArray .* rand(1, dim);
                OffsetX = R .* ( cos(-100 * It / MaxIter * 2 * pi))  .* Oppose;
                OffsetY = -R   .* ( sin(-100 * It / MaxIter * 2 * pi)) .* BinaryArray;
                FullArray = (OffsetX + OffsetY);  FullArray(randi(dim)) = 0;
                Candidate = Pops(i, :)+(1 - exp(-5 * (It+20) / MaxIter)) .* Vector + FullArray;

            else
                dist = sqrt(sum((Pops - BestPop) .^ 2, 2));
                [~,index] = sort(dist);
                idx1 = [index(1: N/5, :); index(N*4/5: N, :)];
                if ismember(i, idx1)
                    Candidate = BestPop +  sign(rand - 0.5) * (normrnd(1, 0.03)) .* (Pops(i, :) - BestPop);
                else
                    S = randperm(N, 2);
                    while any(S == i), S = randperm(N, 2); end
                    if fitness(S(1)) < fitness(S(2))
                        Candidate = Pops(i, :) + (normrnd(0, 0.35)) .* (Pops(S(1), :)-Pops(S(2), :));
                    else
                        Candidate = Pops(i, :) + (normrnd(0, 0.35)) .* (Pops(S(2), :)-Pops(S(1), :));
                    end
                end
            end
        end
        
        %界限处理和适应度计算
        Candidate = Boundary_processing(Candidate, Ub, Lb);
        CandidateFitness = fojb(Candidate);  

        %贪婪选择最优
        if CandidateFitness < fitness(i)
            Pops(i, :) = Candidate;
            fitness(i) = CandidateFitness;
            if CandidateFitness < BestScore
                BestScore = CandidateFitness;
                BestPop    = Candidate;
            end
        end
    conv_curve(It) = BestScore;
    end

end
end

function New_pop = Boundary_processing( Pop, Ub, Lb )
    DimLocal = length(Pop);
    if rand < rand
        Over = (Pop > Ub)+(Pop < Lb);
        New_pop = (rand(1, DimLocal).*(Ub-Lb) + Lb).*Over + Pop.*(~Over);
    else
        New_pop = min(max(Pop, Lb), Ub);
    end
end