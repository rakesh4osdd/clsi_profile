# clsi_profile

 CLSI profile for ASIST using the CLSI standards 2020 by AST profiling for bacteria. It help to generate the input file for ASIST program based on MIC values provided by CLSI based AST profiling.

# ASIST: Antimicrobial Susceptibility standards

 In the first column strain names will be mentioned, the second column will be left blank for getting resistance phenotype, starting two rows will comprise antibiotic names (first row) and antibiotic classes (second row). The data for the resistance profile will then start from the C3 column of an excel file which can be converted into a .csv file. Since B3 is left blank, after implementing the algorithm, column B3 will be filled with the resistance phenotype (susceptible, MDR, XDR, PDR).

**Example Input CSV file:**

	Strain name,Resistance_phenotype,Antibiotic_A1,Antibiotic_A2,Antibiotic_A_N
	Strain_1,Phenotype_1,Resistant,Resistant,Resistant
	Strain_2,Phenotype_2,Resistant,Susceptible,Resistant

 Link to the code : https://github.com/rakesh4osdd/clsi_profile , https://github.com/rakesh4osdd/asist
 
 
# ASIST tool suite

 These tools can be install from Galaxy toolshed to a Galaxy installation https://galaxyproject.org/admin/get-galaxy/.
