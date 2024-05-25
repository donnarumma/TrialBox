function out = NSAtimeEval(confError,pcaComp,pClasses)

[~,indMax] = max([confError.sum_ij]);

pij=pClasses(confError(indMax).class_i).comparisons(pcaComp,:,confError(indMax).class_j);
pij = squeeze(pij);
pij = smooth(pij);

time            = pClasses(1).timecomparisons;
[t,indT]= min(pij(3:end-2));

disp(t);
disp(time(indT+2));

out.t1 = t;
out.indT1 = indT;