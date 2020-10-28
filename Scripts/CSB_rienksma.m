%% Context-Specific-Biomass in sMtb2.0 Model - July 2020
% This script take into account the description of Rienksma et al., 2018 (https://doi.org/10.3389/fcimb.2018.00264) for generate context-specific-biomass reactions,
% and part of this implementation ("gene_to_reaction_levels.m") was modified from Machado et al., 2014  ( https://doi.org/10.1371/journal.pcbi.1003580) scripts.
% Author: Victor López, Guillermo González
% University of Antioquia - Colombia

% cd 'E:\FILES\Tesis-de-Grado-Guillermo\0Tesis Guillermo-DRIVE\Tesis Guillermo\Context-specific-BiomassReaction'
clear;
clc;

%% ADDING DIRECTORY PATH
% addpath ('model')
% addpath ('dataset')

%% PROBLEM SETUP
METABOLIC_MODEL = 'MRR';
TRANSCRIPTOMICS_CATIE_1000 = 'CATIE_1000_intracellular'; % Transcriptomics data from 10.1111/mpp.12134
TRANSCRIPTOMICS_CATIE_R4 = 'CATIE_R4_intracellular'; % 	   Transcriptomics data from 10.1111/mpp.12134
TRANSCRIPTOMICS_CATIE_R7 = 'CATIE_R7_intracellular'; %     Trancriptomics data from 10.1111/mpp.12134
TRANSCRIPTOMICS_POUND_7 = 'POUND_7_intracellular'; %       Transcriptomics data from 10.1111/mpp.12134

COBRA_SOLVER= 'gurobi';
solverOK = changeCobraSolver(COBRA_SOLVER,'all');

%% LOAD METABOLIC MODEL
disp('Integrating Normalized Base Mean into MRR model using the E-Flux approach');
time0 = cputime;

dispstr = sprintf('%5.1f second: reading Moniliophthora roreri network model constrained with in vivo medium...',cputime-time0);
disp(dispstr)

model = load_model(METABOLIC_MODEL);


%% READ GENE EXPRESSION DATA
dispstr = sprintf('%5.1f second: reading expression data...',cputime-time0);
disp(dispstr)

[expression_values_CATIE_1000,gene_id_CATIE_1000,mixed_table_CATIE_1000] = load_transcriptomics_data(TRANSCRIPTOMICS_CATIE_1000);
[expression_values_CATIE_R4,gene_id_CATIE_R4,mixed_table_CATIE_R4] = load_transcriptomics_data(TRANSCRIPTOMICS_CATIE_R4);
[expression_values_CATIE_R7,gene_id_CATIE_R7,mixed_table_CATIE_R7] = load_transcriptomics_data(TRANSCRIPTOMICS_CATIE_R7);
[expression_values_POUND_7,gene_id_POUND_7,mixed_table_POUND_7] = load_transcriptomics_data(TRANSCRIPTOMICS_POUND_7);


dispstr = sprintf('%5.1f second: EFLUX integration...',cputime-time0);
disp(dispstr)

%% MAPPING RXN LEVELS to GPR MRR Model
cutoff = 100; % A cutoff value of percentil 25 was used to identify lowly expressed genes
expression_values_CATIE_1000(expression_values_CATIE_1000 < cutoff) = 0; % , those who accomplishes < cutoff  were assigned a count value of zero
expression_values_CATIE_R4(expression_values_CATIE_R4 < cutoff) = 0; % , those who accomplishes < cutoff  were assigned a count value of zero
expression_values_CATIE_R7(expression_values_CATIE_R7 < cutoff) = 0; % , those who accomplishes < cutoff  were assigned a count value of zero
expression_values_POUND_7(expression_values_POUND_7 < cutoff) = 0; %, those who accomplishes < cutoff  were assigned a count value of zero

% CATIE_1000 infecting cacao fruit
model_CATIE_1000 = model;
levels_CATIE_1000 = gene_to_reaction_levels(model_CATIE_1000, gene_id_CATIE_1000,...
                            expression_values_CATIE_1000, @min, @(x,y)(x+y));
levels_CATIE_1000 (isnan(levels_CATIE_1000 )) = 1; % **** Reactions that received no counts using this method were not allowed to carry any flux
levels_CATIE_1000  = levels_CATIE_1000  / max(levels_CATIE_1000); % **** Rxn Levels Normalization

% CATIE_R4 infecting cacao fruit
model_CATIE_R4 = model;
levels_CATIE_R4 = gene_to_reaction_levels(model_CATIE_R4, gene_id_CATIE_R4,...
                            expression_values_CATIE_R4, @min, @(x,y)(x+y));
levels_CATIE_R4 (isnan(levels_CATIE_R4 )) = 1; %  Reactions that received no counts using this method were not allowed to carry any flux  
levels_CATIE_R4  = levels_CATIE_R4  / max(levels_CATIE_R4); % Rxn Levels Normalization

% CATIE_R7 infecting cacao fruit
model_CATIE_R7 = model;
levels_CATIE_R7 = gene_to_reaction_levels(model_CATIE_R7, gene_id_CATIE_R7,...
                            expression_values_CATIE_R7, @min, @(x,y)(x+y));
levels_CATIE_R7 (isnan(levels_CATIE_R7 )) = 1; %  Reactions that received no counts using this method were not allowed to carry any flux  
levels_CATIE_R7  = levels_CATIE_R7  / max(levels_CATIE_R7); % Rxn Levels Normalization

% POUND_7 infecting cacao fruit
model_POUND_7 = model;
levels_POUND_7 = gene_to_reaction_levels(model_POUND_7, gene_id_POUND_7,...
                            expression_values_POUND_7, @min, @(x,y)(x+y));
levels_POUND_7 (isnan(levels_POUND_7 )) = 1; %  Reactions that received no counts using this method were not allowed to carry any flux  
levels_POUND_7  = levels_POUND_7  / max(levels_POUND_7); % Rxn Levels Normalization


% Replacing LB and UBs into GSMNs
% identifying blocked reactions
blocked_lb = model.lb >= 0;
blocked_ub = model.ub <= 0;

% Replacing LB into CATIE_1000, CATIE_R4, CATIE_R7 and POUND_7
    model_CATIE_1000.lb = -levels_CATIE_1000; % Replacing expression levels.
    model_CATIE_1000.lb (blocked_lb) = 0;
    model_CATIE_R4.lb = -levels_CATIE_R4; % Replacing expression levels.
    model_CATIE_R4.lb (blocked_lb) = 0;
	model_CATIE_R7.lb = -levels_CATIE_R7; % Replacing expression levels.
    model_CATIE_R7.lb (blocked_lb) = 0;
    model_POUND_7.lb = -levels_POUND_7; % Replacing expression levels.
    model_POUND_7.lb (blocked_lb) = 0;
    
% Replacing UB into CATIE_1000, CATIE_R4, CATIE_R7 and POUND_7 
    model_CATIE_1000.ub = levels_CATIE_1000; % Replacing expression levels.
    model_CATIE_1000.ub (blocked_ub) = 0;
    model_CATIE_R4.ub = levels_CATIE_R4; % Replacing expression levels.
    model_CATIE_R4.ub (blocked_ub) = 0;
	model_CATIE_R7.ub = levels_CATIE_R7; % Replacing expression levels.
    model_CATIE_R7.ub (blocked_ub) = 0;
    model_POUND_7.ub = levels_POUND_7; % Replacing expression levels.
    model_POUND_7.ub (blocked_ub) = 0;
% ###########################################################################################################################
dispstr = sprintf('%5.1f second: Biomass Precursor Production from uptake nutrients representing in vivo conditions...',cputime-time0);
disp(dispstr)    
%% Producing Biomass Precursors of UT127
[BMprecursors_values,BMprecursors_names,BMmixed_table] = xlsread('ListBiomassPrecursors_mrr'); % read table of biomass precursors names and MW values
ListBMPs = BMprecursors_names(:,1);
Na_metabolites = BMprecursors_names(:,3);
ListBMPs(1) = [];
Na_metabolites(1) = [];

%% Add sink reaction of biomass precursors and maximizing flux
 
%% infecting cacao fruit
Flux_BMPs_CATIE_1000 = zeros(length(ListBMPs),1);
origStatus_CATIE_1000 = cell(length(ListBMPs),1);

for i = 1:length(ListBMPs)   
    model_CATIE_1000_sink = addSinkReactions(model_CATIE_1000, ListBMPs{i});
    model_CATIE_1000_sink = changeObjective(model_CATIE_1000_sink, strcat("sink_",cellstr(ListBMPs{i})));
    sol = optimizeCbModel(model_CATIE_1000_sink,'max');
    Flux_BMPs_CATIE_1000(i,1) = sol.f;
    origStatus_CATIE_1000{i,1} = sol.origStat;
end 

dispstr = sprintf('%5.1f second: Computing Constraint-Specific Biomass Reactions Coefficients for Mtb CATIE_1000...',cputime-time0);
disp(dispstr)  

MWs = BMprecursors_values(:,20); % Molecular Weights
BM_coeff_CATIE_1000 = (MWs.*Flux_BMPs_CATIE_1000)/sum(MWs.*Flux_BMPs_CATIE_1000); %******

Table_CATIE_1000 = table(ListBMPs,Na_metabolites,Flux_BMPs_CATIE_1000,origStatus_CATIE_1000,MWs,BM_coeff_CATIE_1000);
filename = 'results/CSB_CATIE_1000.xlsx';
writetable(Table_CATIE_1000,filename,'Sheet',1);

%% inserting Context-Specific Biomass Reaction 
CSI_coefficients = zeros(length(model.mets),1); % Context-Specific-Infection Biomass
% Add stoichiometric coefficient
for k = 1:length(ListBMPs)
    
    MetID = findMetIDs(model_CATIE_1000,ListBMPs(k)); % Metabolite IDs sMtb model
    CSI_coefficients(MetID) = -BM_coeff_CATIE_1000(k);
	
end

% ATP, ADP and PI coefficients in sMtb2.0
CSI_coefficients(findMetIDs(model_CATIE_1000,'s_0434')) = -186.2; % ATP coefficient
CSI_coefficients(findMetIDs(model_CATIE_1000,'s_0394')) = 186.2; % ADP coefficient
CSI_coefficients(findMetIDs(model_CATIE_1000,'s_1322')) = 186.2; % PI coefficient
CSI_coefficients(findMetIDs(model_CATIE_1000,'s_0803')) = -186.2; % H2O coefficient

% Creating Context-Specific Biomass Reaction
model_CSB_infection_CATIE_1000 = addReaction(model, 'CSB_infection_CATIE_1000', 'metaboliteList', model.mets ,...
'stoichCoeffList', CSI_coefficients); % insertion of Context-Specific biomass reaction


%% Solving Optimization CATIE_1000
model_CSB_infection_CATIE_1000  = changeObjective(model_CSB_infection_CATIE_1000 , 'CSB_infection_CATIE_1000');
sol_CATIE_1000 = optimizeCbModel(model_CSB_infection_CATIE_1000,'max');
 

%% Producing Biomass Precursors of CATIE_R4
% Add sink reaction of biomass precursors and maximize flux

% CATIE_R4 infecting cacao fruit
Flux_BMPs_CATIE_R4 = zeros(length(ListBMPs),1);
origStatus_CATIE_R4 = cell(length(ListBMPs),1);

for i = 1:length(ListBMPs)   
    model_CATIE_R4_sink = addSinkReactions(model_CATIE_R4, ListBMPs{i});
    model_CATIE_R4_sink = changeObjective(model_CATIE_R4_sink, strcat("sink_",cellstr(ListBMPs{i})));
    sol = optimizeCbModel(model_CATIE_R4_sink,'max');
    Flux_BMPs_CATIE_R4(i,1) = sol.f;
    origStatus_CATIE_R4{i,1} = sol.origStat;
end 

dispstr = sprintf('%5.1f second: Computing Constraint-Specific Biomass Reactions Coefficients for Mtb CATIE_R4...',cputime-time0);
disp(dispstr) 

MWs = BMprecursors_values(:,20); % Molecular Weights
BM_coeff_CATIE_R4 = (MWs.*Flux_BMPs_CATIE_R4)/sum(MWs.*Flux_BMPs_CATIE_R4);

Table_CATIE_R4 = table(ListBMPs,Na_metabolites,Flux_BMPs_CATIE_R4,origStatus_CATIE_R4,MWs,BM_coeff_CATIE_R4);
filename = 'results/CSB_CATIE_R4.xlsx';
writetable(Table_CATIE_R4,filename,'Sheet',1)

%% inserting Context-Specific Biomass Reaction of CATIE_R4
CSI_coefficients = zeros(length(model.mets),1);

% Add stoichiometric coefficient
for k = 1:length(ListBMPs)
    
    MetID = findMetIDs(model_CATIE_R4,ListBMPs(k)); % Metabolite IDs 
    CSI_coefficients(MetID) = -BM_coeff_CATIE_R4(k);
    
end 
% ATP, ADP, H2O and PI coefficients
CSI_coefficients(findMetIDs(model_CATIE_R4,'s_0434')) = -186.2; % ATP coefficient
CSI_coefficients(findMetIDs(model_CATIE_R4,'s_0394')) = 186.2; % ADP coefficient
CSI_coefficients(findMetIDs(model_CATIE_R4,'s_1322')) = 186.2; % PI coefficient
CSI_coefficients(findMetIDs(model_CATIE_R4,'s_0803')) = -186.2; % H2O coefficient
% Creating Context-Specific Biomass Reaction
model_CSB_infection_CATIE_R4 = addReaction(model, 'CSB_infection_CATIE_R4', 'metaboliteList', model.mets ,...
'stoichCoeffList', CSI_coefficients); % Biomass reaction inclusion

%% Solving Optimization CATIE_R4
model_CSB_infection_CATIE_R4  = changeObjective(model_CSB_infection_CATIE_R4 , 'CSB_infection_CATIE_R4');
sol_CATIE_R4 = optimizeCbModel(model_CSB_infection_CATIE_R4,'max');

%% Producing Biomass Precursors of CATIE_R7
% Add sink reaction of biomass precursors and maximize flux

% CATIE_R7 infecting cacao fruit
Flux_BMPs_CATIE_R7 = zeros(length(ListBMPs),1);
origStatus_CATIE_R7 = cell(length(ListBMPs),1);

for i = 1:length(ListBMPs)   
    model_CATIE_R7_sink = addSinkReactions(model_CATIE_R7, ListBMPs{i});
    model_CATIE_R7_sink = changeObjective(model_CATIE_R7_sink, strcat("sink_",cellstr(ListBMPs{i})));
    sol = optimizeCbModel(model_CATIE_R7_sink,'max');
    Flux_BMPs_CATIE_R7(i,1) = sol.f;
    origStatus_CATIE_R7{i,1} = sol.origStat;
end 

dispstr = sprintf('%5.1f second: Computing Constraint-Specific Biomass Reactions Coefficients for Mtb CATIE_R7...',cputime-time0);
disp(dispstr) 

MWs = BMprecursors_values(:,20); % Molecular Weights
BM_coeff_CATIE_R7 = (MWs.*Flux_BMPs_CATIE_R7)/sum(MWs.*Flux_BMPs_CATIE_R7);

Table_CATIE_R7 = table(ListBMPs,Na_metabolites,Flux_BMPs_CATIE_R7,origStatus_CATIE_R7,MWs,BM_coeff_CATIE_R7);
filename = 'results/CSB_CATIE_R7.xlsx';
writetable(Table_CATIE_R7,filename,'Sheet',1)

%% inserting Context-Specific Biomass Reaction of CATIE_R7
CSI_coefficients = zeros(length(model.mets),1);

% Add stoichiometric coefficient
for k = 1:length(ListBMPs)
    
    MetID = findMetIDs(model_CATIE_R7,ListBMPs(k)); % Metabolite IDs 
    CSI_coefficients(MetID) = -BM_coeff_CATIE_R7(k);
    
end 
% ATP, ADP, H2O and PI coefficients
CSI_coefficients(findMetIDs(model_CATIE_R7,'s_0434')) = -186.2; % ATP coefficient
CSI_coefficients(findMetIDs(model_CATIE_R7,'s_0394')) = 186.2; % ADP coefficient
CSI_coefficients(findMetIDs(model_CATIE_R7,'s_1322')) = 186.2; % PI coefficient
CSI_coefficients(findMetIDs(model_CATIE_R7,'s_0803')) = -186.2; % H2O coefficient
% Creating Context-Specific Biomass Reaction
model_CSB_infection_CATIE_R7 = addReaction(model, 'CSB_infection_CATIE_R7', 'metaboliteList', model.mets ,...
'stoichCoeffList', CSI_coefficients); % Biomass reaction inclusion

%% Solving Optimization CATIE_R7
model_CSB_infection_CATIE_R7  = changeObjective(model_CSB_infection_CATIE_R7 , 'CSB_infection_CATIE_R7');
sol_CATIE_R7 = optimizeCbModel(model_CSB_infection_CATIE_R7,'max');

%% Producing Biomass Precursors of POUND_7
% Add sink reaction of biomass precursors and maximize flux

% POUND_7 infecting cacao fruit
Flux_BMPs_POUND_7 = zeros(length(ListBMPs),1);
origStatus_POUND_7 = cell(length(ListBMPs),1);

for i = 1:length(ListBMPs)   
    model_POUND_7_sink = addSinkReactions(model_POUND_7, ListBMPs{i});
    model_POUND_7_sink = changeObjective(model_POUND_7_sink, strcat("sink_",cellstr(ListBMPs{i})));
    sol = optimizeCbModel(model_POUND_7_sink,'max');
    Flux_BMPs_POUND_7(i,1) = sol.f;
    origStatus_POUND_7{i,1} = sol.origStat;
end 

dispstr = sprintf('%5.1f second: Computing Constraint-Specific Biomass Reactions Coefficients for Mtb POUND_7...',cputime-time0);
disp(dispstr) 

MWs = BMprecursors_values(:,20); % Molecular Weights
BM_coeff_POUND_7 = (MWs.*Flux_BMPs_POUND_7)/sum(MWs.*Flux_BMPs_POUND_7);

Table_POUND_7 = table(ListBMPs,Na_metabolites,Flux_BMPs_POUND_7,origStatus_POUND_7,MWs,BM_coeff_POUND_7);
filename = 'results/CSB_POUND_7.xlsx';
writetable(Table_POUND_7,filename,'Sheet',1)

%% inserting Context-Specific Biomass Reaction of POUND_7
CSI_coefficients = zeros(length(model.mets),1);

% Add stoichiometric coefficient
for k = 1:length(ListBMPs)
    
    MetID = findMetIDs(model_POUND_7,ListBMPs(k)); % Metabolite IDs 
    CSI_coefficients(MetID) = -BM_coeff_POUND_7(k);
    
end 
% ATP, ADP, H2O and PI coefficients
CSI_coefficients(findMetIDs(model_POUND_7,'s_0434')) = -186.2; % ATP coefficient
CSI_coefficients(findMetIDs(model_POUND_7,'s_0394')) = 186.2; % ADP coefficient
CSI_coefficients(findMetIDs(model_POUND_7,'s_1322')) = 186.2; % PI coefficient
CSI_coefficients(findMetIDs(model_POUND_7,'s_0803')) = -186.2; % H2O coefficient
% Creating Context-Specific Biomass Reaction
model_CSB_infection_POUND_7 = addReaction(model, 'CSB_infection_POUND_7', 'metaboliteList', model.mets ,...
'stoichCoeffList', CSI_coefficients); % Biomass reaction inclusion

%% Solving Optimization POUND_7
model_CSB_infection_POUND_7  = changeObjective(model_CSB_infection_POUND_7 , 'CSB_infection_POUND_7');
sol_POUND_7 = optimizeCbModel(model_CSB_infection_POUND_7,'max');
%% Saving Models with CSB reaction
dispstr = sprintf('%5.1f second: Saving Models with Constraint-Specific Biomass Reaction...',cputime-time0);
disp(dispstr)  

save('results/model_CSB_infection_CATIE_1000.mat'); % saving model with context-specific biomass reaction
%writeCbModel(model_CSB_infection_CATIE_1000,'xlsx','sMtb2.0_CSB_UT127.xlsx');

save('results/model_CSB_infection_CATIE_R4.mat'); % saving model with context-specific biomass reaction
%writeCbModel(model_CSB_infection_CATIE_R4,'xlsx','sMtb2.0_CSB_CATIE_R4.xlsx');

save('results/model_CSB_infection_CATIE_R7.mat'); % saving model with context-specific biomass reaction
%writeCbModel(model_CSB_infection_CATIE_R7,'xlsx','sMtb2.0_CSB_UT127.xlsx');

save('results/model_CSB_infection_POUND_7.mat'); % saving model with context-specific biomass reaction
%writeCbModel(model_CSB_infection_POUND_7,'xlsx','sMtb2.0_CSB_CATIE_R4.xlsx');


%% Global results
disp('Context-Specific Biomass Reaction Flux')
fprintf('CATIE 1000 CSB flux: %.4f/h.\n', sol_CATIE_1000.f);
fprintf('CATIE R4 CSB flux: %.4f/h.\n', sol_CATIE_R4.f);
fprintf('CATIE R7 CSB flux: %.4f/h.\n', sol_CATIE_R7.f);
fprintf('POUND 7 CSB flux: %.4f/h.\n', sol_POUND_7.f);