#importazione librerie
library(corrplot)
library(FEAR)
library(readxl)
data <- read_excel("DB.xlsx")

#conseguenza della prima corr_matrix(vedi dopo)
data$pop <- data$over65_pop + data$population

# Seleziona le colonne da convertire
cols_to_convert <- setdiff(names(data), "Country")

# Converte le colonne selezionate in numerico
data[cols_to_convert] <- lapply(data[cols_to_convert], function(x) as.numeric(as.character(x)))

# Seleziona solo le colonne numeriche (tranne la seconda che contiene gli anni)
data_numeric <- data[sapply(data, is.numeric)]
data_numeric <- data_numeric[, -1]
#Rimuovi le colonne contenenti dati non rilevanti nell'analisi(ex dati usati solo per scalare i valori)
data_numeric <- data_numeric[, -11]
data_numeric <- data_numeric[, -11]
data_numeric <- data_numeric[, -11]
data_numeric <- data_numeric[, -11]

#----------------------------------------------------------------------------
#normalizzazione per DEA

#definizione funzione normalizzazione
normalize_minmax <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))}

#normalizza i dati 
data_normalized <- as.data.frame(lapply(data_numeric, normalize_minmax))

standardize <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)}

#-----------------------------------------------------------------------
#controlla la correllazione tra tutte le variabili prese in considerazione
corr_matrix <- cor(data_normalized, use = "pairwise.complete.obs", method = "pearson")

corrplot(corr_matrix, method = "color", type = "upper",tl.col = "black", tl.srt = 45)

#outcome:uniamo population e pop_over65

#------------------------------------------------------------------------

#seleziona inputs e outputs INIZIALI
inputs <- as.matrix(data_normalized[, c("Patents_pc", "health_exp","avg_stay_length","n_phy","gdp_pc","tech_eq","tech_hexp_pc","pop")])
outputs <- as.matrix(data_normalized[, c("discharges","mort_inv","h_status")])

#TEST CORRELAZIONE
correlation_matrix <- cor(inputs, outputs, use = "pairwise.complete.obs")
print(correlation_matrix)

#DROP DISCHARGES, AVG_STAY,POP(negative correlation) N_PHY (weak correlation)

#--------------------------------------------------------------------------

#INPUT OUTPUT FINALI
inputs <- as.matrix(data_normalized[, c("Patents_pc", "health_exp","gdp_pc","tech_eq","tech_hexp_pc")])
outputs <- as.matrix(data_normalized[, c("mort_inv","h_status")])

#TEST CORRELAZIONE
correlation_matrix <- cor(inputs, outputs, use = "pairwise.complete.obs")
print(correlation_matrix)

#---------------------------------------------------------------------------
#test preliminari

#trasposte
x<-t(inputs)
y<-t(outputs)

#TEST RETURN TO SCALE
test.rts(x, y, ORIENTATION = 1, METRIC = 1, NSPLIT = 1,NREP = 2000)

#TEST CONVEXITY
test.convexity(x, y, ORIENTATION = 1, METRIC = 1, NSPLIT = 1, NREP = 2000)

#test separability
#library(nonparaeff)
z <- as.matrix(data_normalized$pop)
sep_res <- test.sep.cont(x, y, t(z), ESTIMATOR = 1, ORIENTATION = 1, METRIC = 1, NSPLIT = 1, NREP = 1000)

cat("Tau:", sep_res$tau, "\nP-value:",sep_res$pval, "\n")

#----------------------------------------------------------------------------
#applicazione DEA

detach("package:FEAR", unload = TRUE)
library(Benchmarking)

#DEA crs
dea_resc<- dea(inputs, outputs, RTS ="crs", ORIENTATION = "out", SLACK = TRUE,DUAL = TRUE)
summary(dea_resc)
dea.plot(inputs, outputs, RTS ="crs")

#DEA vrs
dea_resv<- dea(inputs, outputs, RTS ="vrs", ORIENTATION = "out", SLACK = TRUE,DUAL = TRUE)
summary(dea_resv)
dea.plot(inputs, outputs, RTS ="vrs")

#GRAFICI VARI

#----------------------------------------------------------------------------

#DEA bootstrap
detach("package:Benchmarking", unload = TRUE)
library(FEAR)

#evita valori nulli dovuti alla normalizzazione senza alterare la distribuzione
epsilon <- 1e-3
inputs_fix <- inputs + epsilon
outputs_fix <- outputs + epsilon

bootstrap_results <- boot.sw98(t(inputs_fix),t(outputs_fix),RTS = 2,ORIENTATION = 2)

print(bootstrap_results)

original_eff <- 1/bootstrap_results$dhat
bias_corrected_eff <- 1/bootstrap_results$dhat.bc
conf_int<-1/bootstrap_results$conf.int

#GRAFICI VARI

#--------------------------------------------------------------------------------
#SFA vs DEA

#DEA ----------------------------------------------------------------------
library(corrplot)
library(readxl)
data <- read_excel("DB.xlsx")

#conseguenza della prima corr_matrix(vedi dopo)
data$pop <- data$over65_pop + data$population

# Seleziona le colonne da convertire escludendo quelle non importanti
cols_to_convert <- setdiff(names(data), "Country")

# Converte le colonne selezionate in numerico
data[cols_to_convert] <- lapply(data[cols_to_convert], function(x) as.numeric(as.character(x)))

# Seleziona solo le colonne numeriche (tranne la seconda che contiene gli anni)
data_numeric <- data[sapply(data, is.numeric)]
data_numeric <- data_numeric[, -1]
data_numeric <- data_numeric[, -11]
data_numeric <- data_numeric[, -11]
data_numeric <- data_numeric[, -11]
data_numeric <- data_numeric[, -11]

#definizione funzione normalizzazione
normalize_minmax <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))}

#normalizza i dati 
data_normalized <- as.data.frame(lapply(data_numeric, normalize_minmax))

standardize <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)}

inputs <- as.matrix(data_normalized[, c("Patents_pc", "health_exp","gdp_pc","tech_eq","tech_hexp_pc")])
outputs <- as.matrix(data_normalized[, c("h_status")])

library(Benchmarking)

#DEA vrs
dea_resv<- dea(inputs, outputs, RTS ="vrs", ORIENTATION = "out", SLACK = TRUE,DUAL = TRUE)

#SFA (COBB)--------------------------------------------------------------------------
detach("package:Benchmarking", unload = TRUE)
library(frontier)

# Seleziona gli input e output
inputs_sfa <- data_normalized[, c("Patents_pc", "health_exp", "gdp_pc", "tech_eq", "tech_hexp_pc")]
output_sfa <- data_normalized$h_status  #SFA richiede un solo output (h_status o mort_inv)

inputs_sfa <-as.matrix(log(inputs_sfa))
output_sfa <- as.matrix(log(output_sfa))

data_sfa <- data.frame(output_sfa, inputs_sfa)

# Specifica il modello
formula_sfa <- output_sfa ~ .

# Stima del modello

sfa_model <- sfa(formula_sfa, data = data_sfa)


efficiency <- efficiencies(sfa_model)

# Riassunto del modello
summary(sfa_model)

#SFA TRANSLOG----------------------------------------------------------------------------
library(frontier)
inputs_sfa <- data_normalized[, c("Patents_pc", "health_exp", "gdp_pc", "tech_eq", "tech_hexp_pc")]
output_sfa <- data_normalized$h_status  # SFA richiede un solo output (h_status o mort_inv)

inputs_sfa <-as.matrix(log(inputs_sfa))
output_sfa <- as.matrix(log(output_sfa))

data_sfa <- data.frame(output_sfa, inputs_sfa)

formula_translog <- output_sfa ~ 
  Patents_pc + health_exp + gdp_pc + tech_eq + tech_hexp_pc + 
  I(Patents_pc^2) + I(health_exp^2) + I(gdp_pc^2) + I(tech_eq^2) + I(tech_hexp_pc^2) +
  Patents_pc:health_exp + Patents_pc:gdp_pc + Patents_pc:tech_eq + Patents_pc:tech_hexp_pc +
  health_exp:gdp_pc + health_exp:tech_eq + health_exp:tech_hexp_pc +
  gdp_pc:tech_eq + gdp_pc:tech_hexp_pc + tech_eq:tech_hexp_pc

# Fit the translog SFA model
sfa_model_translog <- sfa(formula_translog, data = data_sfa)

# Extract efficiencies
efficiencies_translog <- efficiencies(sfa_model_translog)

# Summary of the model
summary(sfa_model_translog)