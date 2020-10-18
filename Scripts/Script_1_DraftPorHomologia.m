% Este script permite generar un modelo metabólico a escala genómica(draft) 
% de Moniliophthora roreri usando como templetes dos modelos curados. El
% primero de S. Cerevisiae (metabolismo central) y el segundo de Cordyceps
% militaris (Hongo filogenéticamente relacionado)
% 

% Blastp por homología. Se descargaron los proteomas de EnsemblFungi
% Se modificó la nomenclatura de las proteínas para que coincidieran con
% los genes anotados CCM_0001... Archivo fasta: >CCM_0001...

% blastSce=getBlast('mrr', 'moro_mca2997.faa','sce','sce_s288c.faa');
% blastCmt=getBlast('mrr', 'moro_mca2997.faa','cmt','Cmt_01.fa');
% save( 'blastStruct.mat','blast*');

%load('blastStruct.mat');

% Carga el modelo de  S. cerevisiae. Descargado de: 
% https://github.com/SysBioChalmers/yeast-GEM/raw/master/ModelFiles/xml/yeastGEM.xml
% (version 8.4.0).

modelSce        = importModel('yeastGEM840.xml',true);
modelSce.id     = 'sce';


% Cordyceps 
modelCmt        = importExcelModel('model_.xlsx',true);
% % modelCmt        = importModel('CoMi_CM01.xml',true);
modelCmt.id     = 'cmt';

% Estandarizar GrRules ; por Or  : por and.

[modelCmt.grRules, modelCmt.rxnGeneMat]=standardizeGrRules(modelCmt);


save('modelTemplate.mat', 'model*');

%% Generate draft model por homología 

model   = getModelFromHomology(modelSce,blastSce,'mrr',{},1,false,10^-20,150,35);



%% Agregar las reacciones de Cordyceps militaris
model2    = getModelFromHomology(modelCmt,blastCmt,'mrr',{},1,false,10^-20,150,35);

% En caso de tener un id definida use este apartado para eliminar las
% reacciones duplicadas

% model2    = removeReactions(model2,contains(model2.rxns,model.rxns),true,true,true)

modelComb           = mergeModels({model,model2});
model               = contractModel(modelComb);

%% Agregar datos acerca de M. roreri

model.annotation.defaultLB    = -1000;
model.annotation.defaultUB    = +1000;
model.annotation.taxonomy     = 'taxonomy/1381753';
model.annotation.givenName    = 'Guillermo;Victor ';
model.annotation.familyName   = 'González;López';
model.annotation.email        = 'luisg.gonzalez@udea.edu.co,valonso.lopez@udea.edu.co';
model.annotation.organization = 'Universidad de Antioquia';
model.annotation.note         = 'Moniliophthora roreri - 2997';
model.id                      = 'mrr';
model.description             = 'Genome-scale metabolic model of Moniliophthora roreri';

save('MRR_draft.mat','model');

%Formato xml

% exportModel(model,'MMR_draft.xml')

%Formato xlsx
exportToExcelFormat(model,'MRR_draft.xlsx')

