% MRR = importExcelModel('MRR_draft.xlsx',true);
MRR = model;
MRR = replaceMets(MRR,"l-methionyl-trna","l-methionyl-trna(met)");
% MRR_cobra = ravenCobraWrapper(MRR);
metabolitos = cell(1774,1);
for i = 1:1774
    metabolitos{i,1} = char.empty(0);
    metabolitos(i,1) = MRR_cobra.metChEBIID(i,1); 
     
end



Mets =  cell(1138,1);
Replaces = cell(1138,1);


for i = 1775:2912
    
    Mets{i - 1774,1} = char.empty(0);
    Replaces{i - 1774,1} = char.empty(0);
    if MRR_cobra.metChEBIID(i,1) == ""
        continue
    
    elseif MRR_cobra.metChEBIID(i,1) ~= ""
        A = find(metabolitos == strrep(MRR_cobra.metChEBIID(i,1),'',""));
        A(2,1) = 1;
        if A(1,1) == 0
            continue
            
        elseif  A(1,1) > 0
            Replaces{i - 1774,1} = MRR.metNames{A(1,1),1};
            Mets{i - 1774,1} = MRR.metNames{i,1};
            
        end
    end
end
% 
Mets= unique(Mets,'stable');
Replaces = unique(Replaces,'stable');
Mets{148,1} = char.empty(0);
Replaces{148,1} = char.empty(0);
for i = 1:152
    if Mets(i,1) == ""
        continue
    elseif Mets(i,1) ~= ""    
        MRR = replaceMets(MRR,Mets(i,1),Replaces(i,1));
    end
end 

MRR = replaceMets(MRR,"nad(+)","NAD");
MRR = replaceMets(MRR,"nadp(+)","NADP(+)");
MRR = replaceMets(MRR,"nadh","NADH");
MRR = replaceMets(MRR,"nadph","NAD");
MRR = replaceMets(MRR,"oleoyl-coa","oleoyl-CoA");
MRR = replaceMets(MRR,"damp","dAMP");
MRR = replaceMets(MRR,"amp","AMP");
MRR = replaceMets(MRR,"ump","UMP");
MRR = replaceMets(MRR,"dump","dUMP");
MRR = replaceMets(MRR,"cmp","CMP");
MRR = replaceMets(MRR,"dcmp","dCMP");
MRR = replaceMets(MRR,"dgmp","dGMP");
MRR = replaceMets(MRR,"gmp","GMP");
MRR = replaceMets(MRR,"dgtp","dGTP");
MRR = replaceMets(MRR,"gtp","GTP");
MRR = replaceMets(MRR,"gdp","GDP");
MRR = replaceMets(MRR,"dgdp","dGDP");
MRR = replaceMets(MRR,"fad","FAD");
MRR = replaceMets(MRR,"fadh2","FADH2");
MRR = replaceMets(MRR,"fmn","FMN");
MRR = replaceMets(MRR,"cdp","CDP");
MRR = replaceMets(MRR,"dcdp","dCDP");
MRR = replaceMets(MRR,"datp","dATP");
MRR = replaceMets(MRR,"dctp","dCTP");
MRR = replaceMets(MRR,"ctp","CTP");
MRR = replaceMets(MRR,"dadp","dADP");
MRR = replaceMets(MRR,"dudp","dUDP");
MRR = replaceMets(MRR,"udp","UDP");
MRR = replaceMets(MRR,"dutp","dUTP");
MRR = replaceMets(MRR,"utp","UTP");
MRR = replaceMets(MRR,"itp","ITP");
MRR = replaceMets(MRR,"udp-n-acetyl-d-glucosamine","UDP-N-acetyl-alpha-D-glucosamine");
MRR = replaceMets(MRR,"udp-d-galactose","UDP-D-galactose");
MRR = replaceMets(MRR,"udp-glucose","UDP-D-glucose");
MRR = replaceMets(MRR,"udp-glucose","UDP-D-glucose");
MRR = replaceMets(MRR,"l-lysine","L-lysine");
MRR = replaceMets(MRR,"d-xylulose 5-phosphate","D-xylulose 5-phosphate");
MRR = replaceMets(MRR,"-5[(5-phospho-1-deoxy-d-ribulos-1-ylamino)methylideneamino]-1-(5-phospho-d-ribosyl)imidazole-4-carboxamide","5-[(5-phospho-1-deoxy-D-ribulos-1-ylamino)methylideneamino]-1-(5-phospho-D-ribosyl)imidazole-4-carboxamide [cytoplasm]");
MRR = replaceMets(MRR,"acetoacetyl-coa","acetoacetyl-CoA");
MRR = replaceMets(MRR,"acetyl-coa","acetyl-CoA");
MRR = replaceMets(MRR,"succinyl-coa","succinyl-CoA");
MRR = replaceMets(MRR,"decanoyl-coa","decanoyl-CoA");
MRR = replaceMets(MRR,"3-oxoicosanoyl-coa","3-oxoicosanoyl-CoA");
MRR = replaceMets(MRR,"lauroyl-coa","lauroyl-CoA");
MRR = replaceMets(MRR,"malonyl-coa","malonyl-CoA");
MRR = replaceMets(MRR,"myristoyl-coa","myristoyl-CoA");
MRR = replaceMets(MRR,"octanoyl-coa","octanoyl-CoA");
MRR = replaceMets(MRR,"palmitoyl-coa","palmitoyl-CoA");
MRR = replaceMets(MRR,"propionyl-coa","propionyl-CoA");
MRR = replaceMets(MRR,"stearoyl-coa","stearoyl-CoA");
MRR = replaceMets(MRR,"trans-octadec-2-enoyl-coa","trans-octadec-2-enoyl-CoA");
MRR = replaceMets(MRR,"adp","ADP");
MRR = replaceMets(MRR,"atp","ATP");
MRR = replaceMets(MRR,"dtdp","dTDP");
MRR = removeMets(MRR,'TRNA');
MRR = removeMets(MRR,'TRNAm');
MRR = replaceMets(MRR,"coenzyme a","coenzyme A");
MRR = replaceMets(MRR,"cytochrome c","Cytochrome c");
MRR = replaceMets(MRR,"apocytochrome c","Apocytochrome c");
MRR = replaceMets(MRR,"ubiquinone","ubiquinone-6");
MRR = replaceMets(MRR,"ubiquinol","ubiquinol-6");
MRR = replaceMets(MRR,"dihydrofolate","dihydrofolic acid");
MRR = removeReactions(MRR,'r0437');

