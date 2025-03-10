
function [recall,precision,F1score,accuracyBalanced] = evaluateClassificationMetrics(confMatrix)

% confMatrix = confusionmat(y_true, y_pred);

% True Positives (TP) e False Negatives (FN)
TP1 = confMatrix(1,1); % Classe 1
FN1 = confMatrix(1,2);

TP2 = confMatrix(2,2); % Classe 2
FN2 = confMatrix(2,1);

% False Positives (FP) e True Negatives (TN)
FP1 = confMatrix(2,1); % Classe 1
TN1 = confMatrix(2,2);

FP2 = confMatrix(1,2); % Classe 2
TN2 = confMatrix(1,1);

% Recall (Sensibilità)
recall1 = TP1 / (TP1 + FN1);
recall2 = TP2 / (TP2 + FN2);
recall = [recall1;recall2];

% Precision
precision1 = TP1 / (TP1 + FP1);
precision2 = TP2 / (TP2 + FP2);
precision = [precision1;precision2];

% F1-score
F1score(1) = 2 * (precision1 * recall1) / (precision1 + recall1);
F1score(2) = 2 * (precision2 * recall2) / (precision2 + recall2);

% Balanced Accuracy
accuracyBalanced = (recall1 + recall2) / 2;