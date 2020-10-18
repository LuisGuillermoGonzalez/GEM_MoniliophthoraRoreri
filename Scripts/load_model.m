function model = load_model(metabolic_model)
% Load a genome scale model for the selected organism.
% Options:
% metabolic_model: 'MRR'
%
% Author: Victor LOpez, Guillermo González 2020 July.

%% MODEL PATH


    MRR2o2o = 'MRR_draftGap4.xlsx';
    
%% CHOICES    

    switch metabolic_model
    
                 
      case 'MRR'
            
                   MRR = importExcelModel(MRR2o2o,true);
                   model = MRR;
                   model.lb(model.lb < 0) = -1000;
                   model.ub(model.ub > 0) = 1000;
                   exchangeRxns  = model.rxns(findExcRxns(model)); 
                   


                   
                   model = changeRxnBounds(model,exchangeRxns,0,'l');
                  
                   IC = {'Ex_D-glucose','Ex_raffinose','Ex_L-asparagine','Ex_L-valine','Ex_d-mannitol','Ex_L-leucine'...
                         'Ex_L-phenylalanine','Ex_gamma-aminobutyrate','Ex_L-glutamine','Ex_myo-inositol','r_2056','r_1687'...
                         'r_1586','Ex_L-tyrosine','Ex_D-fructose','Ex_L-aspartate','Ex_glycogen','Ex_L-arginine'...
                           };


                   precursores_esenciales = {...
                                    'nh3IN','Ex_H2O','Ex_phosphate','Ex_potassium','Ex_Mg(2+)',...
                                    'Ex_sodium','Ex_oxygen','Ex_sulphate','Ex_iron(2+)','Ex_Mn(2+)','Ex_Zn(2+)'...
                                    'Ex_chloride','r_4600','r_4594'}; 
                                
                   model = changeRxnBounds(model,IC, -0.1,'l'); 
                   model = changeRxnBounds(model,precursores_esenciales, -1000,'l');             
      end         
    
end