# On-multifactorial-drivers-for-malaria-rebound-in-Brazil

You can replicate work results using the Rscripts from "Bayesian models and incidence prediction" folder. In this folder, you can reproduce figures and tables.

We also share a set of previous data files:

1. Deforestation data: we download deforestation indicators from PRODES website (http://www.dpi.inpe.br/prodesdigital/prodesmunicipal.php). This data was saved in "dados_desmatamento.Rdata". You can obtain the original data through PRODES website.

2. Incidence data: we requested incidence data to the Brazilian Health Ministry through SIVEP malaria system (http://200.214.130.44/sivep_malaria/). If you want to obtain the original data, you must create a formal request in SIVEP malaria website. Nevertheless, we debug the original data implementing "Depurate_data.R" script that generated the clean database ("BASE_ABA_042019.RData") from original data. 

3. INLA data: this folder contains two Rscripts for creating data bases in the correct spatial and temporal order to implement Bayesian model with INLA. "P_falciparum_data.R" script generates P. falciparum database and "P_vivax_data.R" script generates P. vivax database.

4. Maps: this folder contains a set of files to implement Brazilian map in R. These files were obtaining from the Brazilian Geography and statistic Institute IBGE through its website: https://www.ibge.gov.br/geociencias/downloads-geociencias.html. On the other hand, "Bayesian models and incidence prediction" folder contains the Brazilian maps per State, and we obtained these maps from IBGE website. If you wanto to use these files, you have to decompress the folder.

5. Population data: this folder contains population data for each municipality in Amazon basin from 2003 to 2018. We obtained these data base from Brazilian Geography and statistic Institute IBGE through the website: https://www.ibge.gov.br/estatisticas/sociais/populacao.html. We ordered data in "Population.xls" file.
