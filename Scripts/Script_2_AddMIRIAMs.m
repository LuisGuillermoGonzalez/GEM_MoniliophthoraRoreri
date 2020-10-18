% load('MRR_draft.mat')
% model2 = ravenCobraWrapper(MRR_draft);
% 
%  MIRIAMs_id = importExcelModel('MIRIAMS.xlsx',true);
% MIRIAMs_id2 = ravenCobraWrapper(MIRIAMs_id)
% 
% load('Rtho.mat')
% Rhto2 = ravenCobraWrapper(Rtho);



% model2.metBIGGID = cell(length(model2.metChEBIID),1);
% 
% 
% for i= 1:length(model2.metChEBIID)
%    model2.metBIGGID{i,1} = char.empty(0); 
%    A =  find(Rhto2.metNames == strrep(model2.metNames(i,1),'',""));
%   if A > 0
%       model2.metBIGGID(i,1) = Rhto2.metBIGGID(A,1);
%            
%   elseif A == 0        
%     continue
%       
%   end      
%      
% end
% 
%  [num,Metacyc_names] = xlsread('Metacyc_data.xlsx', 'A:B');
% 
% for i= transpose(find(model2.metBIGGID == ""))
%     
%    A =  find(Metacyc_names == strrep(model.metNames(i,1),'',""));
%   if A > 0 
%       model2.metChEBIID(i,1) = Metacyc_names(A,2);
%    
%   elseif A == 0
%     continue
%       
%   end      
%      
% end
% 
% [num2,Metacyc_synonyms] = xlsread('Metacyc_data','C:C');
% 
% 
% for i= transpose(find(model2.metBIGGID == ""))
%    
%     for j = 1:length(Metacyc_synonyms)
%       A = find(strsplit(Metacyc_synonyms{j,1},"//") == strrep(model.metNames(i,1),'',""));
%       if  A > 0
%           model2.metChEBIID(i,1) = Metacyc_names(j,2);
% 
%       elseif A == 0
%           continue
%           
%     
%       
%       end      
%     end  
% end

tic
ticBytes(gcp);
parfor i= transpose(find(model2.metBIGGID == ""))
    
  if model2.metChEBIID(i,1) == ""
     
     continue      
          
  elseif  model2.metChEBIID(i,1) ~=	""
      parfor j= 1:length(MIRIAMs_id2.metChEBIID)
         if  find(strsplit(MIRIAMs_id2.metChEBIID{j,1},";") == strrep(model2.metChEBIID{i,1},'',"")) > 0
             model2.metBIGGID{i,1} = MIRIAMs_id2.metBIGGID{j,1};
         elseif find(strsplit(MIRIAMs_id2.metChEBIID{j,1},";") == strrep(model2.metChEBIID{i,1},'',"")) == 0       
                continue
         end
         
      end
  end      
     
end
tocBytes(gcp)
toc
save('model2.mat')