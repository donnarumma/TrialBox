function [Accuracy,kappaValue] = prob_nbpw(outNBPW,labs)

P_test = struct();
for n=1:length(labs)
    for i=1:length(outNBPW)
        P_test(n).prob(:,i) = outNBPW(i).test.prob(n,:)';
    end
end

p_pred = nan(length(labs),1);
lab_pred = nan(length(labs),1);
for ntr=1:length(labs)
    p_prod = prod(P_test(ntr).prob,2);
    p_final_prod = p_prod/sum(p_prod);
    % 
    % p_sum =  sum(P_test(ntr).prob,2);
    % p_final = p_sum/sum(p_sum);
    % p_max = max(P_test(ntr).prob,[],2);
    % p_final = p_max/sum(p_max);
    [p_pred(ntr), lab_pred(ntr)] = max(p_final_prod);
end

Accuracy = sum(lab_pred == labs)/length(labs)*100;
Cmatrxix = confusionmat(labs, lab_pred);
kappaValue  = kappaModel(Cmatrxix);