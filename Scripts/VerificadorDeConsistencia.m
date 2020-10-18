  
% SolOk = changeCobraSolver('gurobi','all');
% model = importExcelModel('MRR_draftGap3.xlsx',true);
% modelSce = importModel('yeastGEM840.xml',true);

model_cobra = ravenCobraWrapper(model);
SCM = cell(length(model.mets),1);
S2 = full(model.S);
for i= 1:length(model.mets)
    
    SCM{i,1} = char.empty(0);
    
    A = find(S2(i,:) == 0);
    if length(A) == length(model.rxns)- 1
       SCM(i,1)= model_cobra.metNames(i,1); 
       
    else
        continue
    end    
end

SCM = unique(SCM);

DEM = cell(length(model.mets),1);

R_f = find(model.rev == 1);
t_columna = size(S2,2);
nt_columna = t_columna + length(R_f);
S_a = S2;
S_a(:,t_columna + 1:nt_columna) = -S_a(:,R_f);

for i= 1:length(model.mets)
    
    DEM{i,1} = char.empty(0);
    
    A = find(S_a(i,:) > 0);
    B = find(S_a(i,:) < 0);
    if (length(S_a(i,:)) - length(find(S_a(i,:)==0))== length(A) ) || ...
       (length(S_a(i,:)) - length(find(S_a(i,:) == 0)) == length(B) )
        if length(A) > 1 || length(B) > 1
            DEM(i,1)= model_cobra.metNames(i,1); 
        end
        
    else
        continue
    end    
end

 DEM = unique(DEM);

[Vmin, Vmax] = fluxVariability(model,100,'max');
 ZFR = cell(length(model.rxns),1);

for i = 1:length(model.rxns)
    ZFR{i,1} = char.empty(0);
    
    if (Vmin(i,1) == 0) && (Vmax(i,1) == 0)
        ZFR(i,1) = model_cobra.rxns(i,1);
        
    else
        continue
    end
end

ZFR = unique(ZFR);
% Lista = cell(length(ZFR)-1,1);
% for i = 1:length(ZFR) -1
%     Lista{i,1} = char.empty(0);
%     
%     Lista(i,1) = model.subSystems{find(model_cobra.rxns == strrep(ZFR(i+1,1),"",''))}(1,1);
%     
%     
% end
% 
% Genes_ZFR = findGenesFromRxns(model_cobra, ZFR(2:753,1));
% Genes_ZFRO = cell(length(Genes_ZFR),1);
% for i = 1:length(Genes_ZFR)
%     Genes_ZFRO{i,1} = char.empty(0);
%     
%    Celda = length(Genes_ZFR{i,1});
%    if Celda == 1 || Celda == 0
%        Genes_ZFRO(i,1) = Genes_ZFR(i,1);
%        
%    elseif Celda > 1
%        a = [];
%        for j = 1:Celda
%            a = strcat(a, Genes_ZFR{i}(j,1),';');
%        end
%        Genes_ZFRO(i,1) = a(1,1);
%    end
%     
% end
% UR = cell(length(model.rxns),1);
% for i = 1:length(model.rxns)
%     UR{i,1} = char.empty(0);
%     
%     if model.rev(i,1) == 0
%         continue
%            
%     elseif (Vmin(i,1) >= 0) || (Vmax(i,1) <= 0)
%         UR(i,1) = model_cobra.rxnNames(i,1);
%     else
%         continue
%     
%     end
% end
%  
% UR = unique(UR);   
% 
% N = null(S_a,'r');
% nN_columnas = size(N,2);
% RR = find(model.rev ==0);
% CR = [];
% RCR = [];
% for i = 1:length(RR)
%     
%     for j = i+1:length(RR)
%         bFlag = true;
%         a = [];
%         for k = 1:nN_columnas
%             if (N(RR(i), k)~=0 && N(RR(j), k)~=0)
%                 a(k) = N(RR(i), k)/N(RR(j), k);
%             elseif (N(RR(i), k)==0 && N(RR(j), k)==0)
%                 a(k) = 0;
%             else
%                 bFlag = false;
%                 break
%             end
%         end
%         if (~isempty(a))
%             a(a==0) = [];
%         end
%         if (~isempty(a) && all(a==a(1)) && bFlag == true)
%                 CR = [CR; RR(i) RR(j)];
%             if a(1)<0
%                 RCR = [RCR; RR(i) RR(j)];
%             end
%         end
%     end
% end
%         
%         
