#Project 3 - Investigating potential virulence markers in the Nonstructural protein 1 of Avian Influenza Virus

#Knowledge-Based Systems in Bioinformatics - 1MB465
#Uppsala University - Spring Semester 2024
#Miltiadis Kesidis
#Agapi Eleni Simaiaki

"""
Data origin: NCBI's Influenza resource - Until Janoury 2014

Type: A
Host: Avian
Country/Region: any
Protein: HA
Subtype: H: any N: any
Full length plus: checked
Collapse identical sequences: checked

The data used here were provided to us preprocessed and balanced for the decision class.
Amino acid sequences were aligned using MUSCLE (v3.8.31)
Annotations: 1-> High Pathogenicity (HP) / 0-> Low Pathogenicity (LP) / ? -> Gap in

The presented study wad based on the work: Khaliq, Z., Leijon, M., Belák, S., & Komorowski, J. (2015).
A complete map of potential pathogenicity markers of avian influenza virus subtype H5 predicted from 
11 expressed proteins. BMC microbiology, 15, 128. https://doi.org/10.1186/s12866-015-0465-x
"""


#Import data
init_data <- read.table("Project3.csv", header = TRUE, sep='\t', colClasses = 'character')
data <- init_data[,2:251] #keep only features/remove strain names

#Check if the dataset is balanced for decision class
table(init_data$Pathogenicity) 


################################# MCFS ###################################
library(rmcfs)

mcsf_fts <- mcfs(Pathogenicity~., data, projections=500, 
                projectionSize=50, splits=5, splitSetSize = 40,
                cutoffPermutations = 20, threadsNumber = 8)

#Inspect Results
mcsf_fts$RI[1:14,]
plot(mcsf_fts, type = "distances")
plot(mcsf_fts, type = "cmatrix")

#Plot MCFS network
gid <- build.idgraph(mcsf_fts)
plot.idgraph(gid, label_dist = 0.3)
plot(result4.1, type = "features")


################################ Rosetta #################################
library(devtools)
library(R.ROSETTA)


############## Johnson Reducer ###############

#Run Rosetta

dataJohnson <- rosetta(mcsf_fts$data,reducer = "Johnson", roc = TRUE, discrete = TRUE, clroc = "1")
rules_Johnson <- dataJohnson$main                              #extract rules
qual_Johnson <- dataJohnson$quality                            #extract quality metrics
table(rules_Johnson$decision)                                  #check the amount of rules per decision class
top_rules <- rules_Johnson[rules_Johnson$pValue < 0.05,]       #extract most significant rules

#ROC curve / Validation of method
plotMeanROC(dataJohnson)

#Check the amount of significant rules
tabS_Johnson <- table(rules_Johnson[rules_Johnson$pValue < 0.05,]$decision) #tested also for p-value < 0.1 and 0.2 
tabS_Johnson                                                                #small differences of additional 2 ruls

#View rules
viewRules(rules_Johnson[rules_Johnson$decision=="0",]) #for LP
viewRules(rules_Johnson[rules_Johnson$decision=="1",]) #for HP

#Top five rules
rules_Johnson[rules_Johnson$decision=="1",][1:5,]
rules_Johnson[rules_Johnson$decision=="0",][1:5,]


############# #Genetic reducer ##############

#Run Rosetta
dataGenetic <- rosetta(mcsf_fts$data, reducer = "Genetic", roc = TRUE, discrete = TRUE, clroc = "1")
rules_Genetic <- dataGenetic$main
qual_Genetic <- dataGenetic$quality

#ROC curve
plotMeanROC(dataGenetic)

#Check the amount of significant rules
tabS_Genetic <- table(rules_Genetic[rules_Genetic$pValue < 0.001,]$decision)
tabS_Genetic



############# Visualization with VisuNet ##############

library(VisuNet)
library(clusterProfiler)
library(org.Hs.eg.db)

vis_Johnson <- visunet(dataJohnson$main)
vis_Genetic <- visunet(dataGenetic$main)


#Additional exploration 0f the reducers
#Comparison of common nodes 

gf_J <- getFeatures(rules_Johnson,filter = T, filterType = "pvalue", thr = 0.05)
hp_J <- gf_J$features$"1"
lp_J <- gf_J$features$"0"


gf_G <- getFeatures(rules_Genetic, filter = T, filterType = "pvalue", thr = 0.05)
hp_G <- gf_G$features$"1"
lp_G <- gf_G$features$"0"


intersect(hp_J, hp_G)
intersect(lp_J, lp_G)
#Note: the nodes for both HP & LP of the Johnson are almost all includes in the corresponding
#      groups of the Genetic approach



'''
#Optimization of MCFS parameters
#projectionSize was tested for : 25, 50 and 75 based on the small # of features in the inital dataset
#The selection of the parameters was based on the distance plot and the confusion matrix results

#Setting the number of projection
#Test for 500
data <- init_data[,2:251]
result1 <- mcfs(Pathogenicity~., data, projections=500, 
               projectionSize=25, splits=5, splitSetSize=20,
               cutoffPermutations =20, threadsNumber=8)
head(result1$RI)
plot(result1, type = "distances")
plot(result1, type = "cmatrix")


result2 <- mcfs(Pathogenicity~., data, projections=500, 
                projectionSize=50, splits=5, splitSetSize = 40,
                cutoffPermutations = 20, threadsNumber = 8)
head(result2$RI)
plot(result2, type = "distances")
plot(result2, type = "cmatrix")


result3 <- mcfs(Pathogenicity~., data, projections=500, 
                projectionSize=75, splits=5, splitSetSize = 60,
                cutoffPermutations = 20, threadsNumber = 8)
head(result3$RI)
plot(result3, type = "distances")
plot(result3, type = "cmatrix")


#Test for 750
result4 <- mcfs(Pathogenicity~., data, projections=750, 
                projectionSize=25, splits=5, splitSetSize=20,
                cutoffPermutations =20, threadsNumber=8)
head(result4$RI)
plot(result4, type = "distances")
plot(result4, type = "cmatrix")
result4$RI
table(result4$RI$projections) #check the distribution to judge 
#4 chedck again 


result5 <- mcfs(Pathogenicity~., data, projections=750, 
                projectionSize=50, splits=5, splitSetSize = 40,
                cutoffPermutations = 20, threadsNumber = 8)
head(result5$RI)
plot(result5, type = "distances")
plot(result5, type = "cmatrix")


result6 <- mcfs(Pathogenicity~., data, projections=750, 
                projectionSize=75, splits=5, splitSetSize = 60,
                cutoffPermutations = 20, threadsNumber = 8)
head(result6$RI)
plot(result6, type = "distances")
plot(result6, type = "cmatrix")

# See results
cutoff_value<-19
r1 <- result1$RI[1:cutoff_value,c(2,6)]
r2 <- result2$RI[1:cutoff_value,c(2,6)]
r3 <- result3$RI[1:cutoff_value,c(2,6)]
r4 <- result4$RI[1:cutoff_value,c(2,6)]
r5 <- result5$RI[1:cutoff_value,c(2,6)]
r6 <- result6$RI[1:cutoff_value,c(2,6)]

combo_res <- cbind(r1,r2,r3,r4,r5,r6)
colnames(combo_res) <- c("500/25_attr","500/25_RI",
                         "500/50_attr","500/50_RI",
                         "500/75_attr","500/75_RI",
                         "750/25_attr","750/25_RI",
                         "750/50_attr","750/50_RI",
                         "750/75_attr","750/75_RI")

#Note: in all of the test case the same top 10 eatures are extracted
#Following the mcfs paper we select the 500 projections with 50 feautures size - result2
'''

