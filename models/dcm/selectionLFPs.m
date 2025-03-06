% function [selS,selK,depthS,depthK] = selectionLFPs(session_name,par)

function [selS,selK,depthS,depthK] = selectionLFPs(session_name,par)


behav_JOINT  = 'H';
% behav_EYE    = 'E';
signal_type  = 'Raw';
f_Joint      = load([session_name behav_JOINT '_' signal_type]);
Raw_H        = f_Joint.Trials;

%% select the good LFP indexes
Depth        = Raw_H(1).Depth;
iLFP         = false(1,length(Depth));
dLFP         = ~isnan(Depth);

vLFP_S       = iLFP;
vLFP_S(1:5)  = true;
vLFP_S       = vLFP_S & dLFP;
% iLFP_S       = find(vLFP_S);

vLFP_K       = iLFP;
vLFP_K(6:10) = true;
vLFP_K       = vLFP_K & dLFP;
% iLFP_K       = find(vLFP_K);

if strcmp(par.method,'Maximum')
    try
        [depthS,selS] = max(Depth(vLFP_S));
        [depthK,selK] = max(Depth(vLFP_K));
    catch
        selS          = [];
        selK          = [];
        return
    end
% elseif strcmp(par.method, 'Similar')
%     try
%         % Estraggo i valori validi (non NaN) per i due gruppi con i loro indici originali
%         validIdx_S = find(vLFP_S & dLFP); % Indici posizionali relativi a 1:5
%         validIdx_K = find(vLFP_K & dLFP); % Indici posizionali relativi a 6:10
% 
%         validDepth_S = Depth(validIdx_S); % Valori non NaN gruppo 1-5
%         validDepth_K = Depth(validIdx_K); % Valori non NaN gruppo 6-10
% 
%         % Se entrambi i gruppi hanno almeno un valore valido
%         if ~isempty(validDepth_S) && ~isempty(validDepth_K)
%             % Inizializzo la distanza minima ad un valore grande
%             minDist = Inf;
%             selS = [];
%             selK = [];
%             depthS = [];
%             depthK = [];
% 
%             % Ciclo su ogni valore del primo gruppo
%             for i = 1:length(validDepth_S)
%                 % Ciclo su ogni valore del secondo gruppo
%                 for j = 1:length(validDepth_K)
%                     % Calcolo la distanza assoluta tra i due valori
%                     dist = abs(validDepth_S(i) - validDepth_K(j));
% 
%                     % Se la distanza è la più piccola trovata, aggiorno gli output
%                     if dist < minDist
%                         minDist = dist;
%                         selS = i; % Indice relativo a 1:5
%                         selK = j; % Indice relativo a 1:5 (del gruppo 6:10)
%                         depthS = validDepth_S(i); % Valore corrispondente
%                         depthK = validDepth_K(j); % Valore corrispondente
%                     end
%                 end
%             end
% 
%             % Stampa i risultati
%             disp(['I due valori più vicini tra i gruppi sono: ', num2str(depthS), ' e ', num2str(depthK)]);
%             disp(['La distanza tra i due valori è: ', num2str(minDist)]);
%             disp(['Gli indici relativi sono: ', num2str(selS), ' e ', num2str(selK)]);
%         else
%             % Se uno dei due gruppi è vuoto, restituisci valori vuoti
%             selS = [];
%             selK = [];
%             depthS = [];
%             depthK = [];
%             return
%         end
%     catch
%         selS = [];
%         selK = [];
%         depthS = [];
%         depthK = [];
%         return
%     end
end