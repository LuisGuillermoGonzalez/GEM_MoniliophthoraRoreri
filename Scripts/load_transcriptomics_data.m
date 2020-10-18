function [expression_values,gene_id,mixed_table] = load_transcriptomics_data(TRANSCRIPTOMICS)
% Load a transcriptomics data for the selected conditions.
% Options: 'CATIE_1000_intracellular', 'CATIE_R4_intracellular', 'CATIE_R7_intracellular', 'POUND_7_intracellular'
%
% Author: Victor A. Lopez-Agudelo, Guillermo González 2020 July.


    
    CATIE_1000_DATA = 'Datos_CATIE-1000_CMP.xlsx'; 
    CATIE_R4_DATA = 'Datos_CATIE-R4_CMP.xlsx';
	CATIE_R7_DATA = 'Datos_CATIE-R7_CMP.xlsx';
	POUND_7_DATA = 'Datos_Pound7_CMP.xlsx';

    switch TRANSCRIPTOMICS
        
                
        case 'CATIE_1000_intracellular'

            [expression_values,gene_id,mixed_table] = xlsread(CATIE_1000_DATA);
            
        case 'CATIE_R4_intracellular'

            [expression_values,gene_id,mixed_table] = xlsread(CATIE_R4_DATA);
		
		case 'CATIE_R7_intracellular'

            [expression_values,gene_id,mixed_table] = xlsread(CATIE_R7_DATA);
	
		case 'POUND_7_intracellular'

            [expression_values,gene_id,mixed_table] = xlsread(POUND_7_DATA);

    end
    
end