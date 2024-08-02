#Compare sigma2 values of other studies

#Panciroli et al. (2024) Mammal body mass (built upon previous datasets)

Panciroli_data <- read.csv('Panciroli2024Dataset.csv',header = T, row.names = NULL)
  head(Panciroli_data)

#Code adapted from Panciroli et al. (2024)

Panciroli_data <- Panciroli_data[ Panciroli_data$tree_names != "" , ]
Panciroli_data$Femur.Circumference <- as.numeric( Panciroli_data$Femur.Circumference )
Panciroli_data$Humerus.Circumference <- as.numeric( Panciroli_data$Humerus.Circumference )
Panciroli_data <- Panciroli_data[ complete.cases( Panciroli_data[ , c( "Femur.Circumference" , "Humerus.Circumference" ) ] ) , ]
Panciroli_data <- Panciroli_data[ Panciroli_data$Bi_Quad == "Q" , ]

##Read phylogenetic trees from Upham et al. 2019
Panciroli_trees <- read.nexus('Panciroli_trees.nex')
Panciroli_data[ , "FCHC" ] <- Panciroli_data$Femur.Circumference + Panciroli_data$Humerus.Circumference

##Phylogenetic regressions using new data

pgls.results.list.new.data <- list()
for( i in 1:100 ) {
  pgls.results.list.new.data[[ i ]] <- list()
  
  data.temp <- Panciroli_data[ Panciroli_data$Dataset == "A. New data (Panciroli et al 2023)" , ]
  data.temp$tree_names <- gsub(' ','_',data.temp$tree_names)
  tree.temp <- drop.tip( Panciroli_trees[[i]] , Panciroli_trees[[i]]$tip.label[ is.na( match( Panciroli_trees[[i]]$tip.label , data.temp$tree_names ) ) ] )	
  data.tree <- data.temp[ match( tree.temp$tip.label , data.temp$tree_names ) , ]
  rownames( data.tree ) <- data.tree$tree_names
  data.tree[ , "log_HC" ] <- log10( data.tree$Humerus.Circumference )
  data.tree[ , "log_FC" ] <- log10( data.tree$Femur.Circumference )
  data.tree[ , "log_HL" ] <- log10( data.tree$Humerus.Length )
  data.tree[ , "log_FL" ] <- log10( data.tree$Femur.Length )
  data.tree[ , "log_FCHC" ] <- log10( data.tree$Femur.Circumference + data.tree$Humerus.Circumference )
  data.tree[ , "log_BM" ] <- log10( data.tree$Body_mass_literature..g. )
  

  pgls.results.list.new.data[[ i ]][[1]] <- phylolm( log_BM ~ log_HC , data.tree , model='lambda',phy = Panciroli_trees[[i]]) 
  pgls.results.list.new.data[[ i ]][[2]] <- phylolm( log_BM ~ log_FC , data.tree , model='lambda',phy = Panciroli_trees[[i]]) 
  pgls.results.list.new.data[[ i ]][[3]] <- phylolm( log_BM ~ log_FCHC , data.tree , model='lambda',phy = Panciroli_trees[[i]]) 
  pgls.results.list.new.data[[ i ]][[4]] <- phylolm( log_BM ~ log_FL , data.tree , model='lambda',phy = Panciroli_trees[[i]]) 
  pgls.results.list.new.data[[ i ]][[5]] <- phylolm( log_BM ~ log_HL , data.tree , model='lambda',phy = Panciroli_trees[[i]]) 
  pgls.results.list.new.data[[ i ]][[6]] <- phylolm( log_FL ~ log_HL , data.tree , model='lambda',phy = Panciroli_trees[[i]]) 
  names( pgls.results.list.new.data[[ i ]] ) <- c("HC","FC","FCHC","FL","HL","FLHL")
}

#Sigma2 values
#Sigma2 individual traits
#Sigma2 for body mass
merror <- data.tree$log_BM * 0.01
names(merror) <- rownames(data.tree)

set.seed(123)
i= sample(1:length(Panciroli_trees),1,replace = F)
tree.temp <- drop.tip( Panciroli_trees[[i]] , Panciroli_trees[[i]]$tip.label[ is.na( match( Panciroli_trees[[i]]$tip.label , data.temp$tree_names ) ) ] )	

Panciroli_sigma_BM <- phylosig(tree.temp, 
                           x = data.tree[tree.temp$tip.label,'log_BM'] , 
                           method="lambda",test=T,
                           se=merror )

Panciroli_sigma_BM$sig2 #0.0142

#Sigma2 values for humerus circumference
merror <- data.tree$log_HC * 0.01
names(merror) <- rownames(data.tree)

set.seed(123)
i= sample(1:length(Panciroli_trees),1,replace = F)
tree.temp <- drop.tip( Panciroli_trees[[i]] , Panciroli_trees[[i]]$tip.label[ is.na( match( Panciroli_trees[[i]]$tip.label , data.temp$tree_names ) ) ] )	

Panciroli_sigma_HC <- phylosig(tree.temp, 
                               x = data.tree[tree.temp$tip.label,'log_HC'] , 
                               method="lambda",test=T,
                               se=merror )

Panciroli_sigma_HC$sig2 #0.0194

#Sigma2 values for femur circumference
merror <- data.tree$log_FC * 0.01
names(merror) <- rownames(data.tree)

set.seed(123)
i= sample(1:length(Panciroli_trees),1,replace = F)
tree.temp <- drop.tip( Panciroli_trees[[i]] , Panciroli_trees[[i]]$tip.label[ is.na( match( Panciroli_trees[[i]]$tip.label , data.temp$tree_names ) ) ] )	

Panciroli_sigma_FC <- phylosig(tree.temp, 
                               x = data.tree[tree.temp$tip.label,'log_FC'] , 
                               method="lambda",test=T,
                               se=merror )

Panciroli_sigma_FC$sig2 #0.0192

#Sigma2 values for humerus+femur circumference
merror <- data.tree$log_FCHC * 0.01
names(merror) <- rownames(data.tree)

set.seed(123)
i= sample(1:length(Panciroli_trees),1,replace = F)
tree.temp <- drop.tip( Panciroli_trees[[i]] , Panciroli_trees[[i]]$tip.label[ is.na( match( Panciroli_trees[[i]]$tip.label , data.temp$tree_names ) ) ] )	

Panciroli_sigma_FCHC <- phylosig(tree.temp, 
                               x = data.tree[tree.temp$tip.label,'log_FCHC'] , 
                               method="lambda",test=T,
                               se=merror )

Panciroli_sigma_FCHC$sig2 #0.002


#Sigma2 values for femur length
merror <- data.tree$log_FL * 0.01
names(merror) <- rownames(data.tree)

set.seed(123)
i= sample(1:length(Panciroli_trees),1,replace = F)
tree.temp <- drop.tip( Panciroli_trees[[i]] , Panciroli_trees[[i]]$tip.label[ is.na( match( Panciroli_trees[[i]]$tip.label , data.temp$tree_names ) ) ] )	

Panciroli_sigma_FL <- phylosig(tree.temp, 
                                 x = data.tree[tree.temp$tip.label,'log_FL'] , 
                                 method="lambda",test=T,
                                 se=merror )

Panciroli_sigma_FL$sig2 #0.00149

#Sigma2 values for humerus length
merror <- data.tree$log_HL * 0.01
names(merror) <- rownames(data.tree)

set.seed(123)
i= sample(1:length(Panciroli_trees),1,replace = F)
tree.temp <- drop.tip( Panciroli_trees[[i]] , Panciroli_trees[[i]]$tip.label[ is.na( match( Panciroli_trees[[i]]$tip.label , data.temp$tree_names ) ) ] )	

Panciroli_sigma_HL <- phylosig(tree.temp, 
                               x = data.tree[tree.temp$tip.label,'log_HL'] , 
                               method="lambda",test=T,
                               se=merror )

Panciroli_sigma_HL$sig2 #0.00149


#Sigma2 relationships
#Sigma2 values for body mass vs. humerus circumference
mean(sapply ( pgls.results.list.new.data, function(x) x$HC$sigma2))

#Sigma2 values for body mass vs. femur circumference
mean(sapply ( pgls.results.list.new.data, function(x) x$FC$sigma2))

#Sigma2 values for body mass vs. humerus+femur circumference
mean(sapply ( pgls.results.list.new.data, function(x) x$FCHC$sigma2))

#Sigma2 values for body mass vs. femur length
mean(sapply ( pgls.results.list.new.data, function(x) x$FL$sigma2))

#Sigma2 values for body mass vs. humerus length
mean(sapply ( pgls.results.list.new.data, function(x) x$HL$sigma2))

#Sigma2 values for femur length vs. humerus length
mean(sapply ( pgls.results.list.new.data, function(x) x$FLHL$sigma2))


#Other dataset, same code more or less
#Using Campione & Evans (2012) data (literature- and measured-based; only mammals, after Panciroli et al. 2024)

pgls.results.list.Campione_spec <- list()
pgls.results.list.Campione_lit <- list()

for( i in 1:100 ) {
  pgls.results.list.Campione_lit[[ i ]] <- list()
  pgls.results.list.Campione_spec[[ i ]] <- list()
  
  data.temp <- Panciroli_data[ Panciroli_data$Dataset == "Campione and Evans 2012" , ]
  data.temp$tree_names <- gsub(' ','_',data.temp$tree_names)
  tree.temp <- drop.tip( Panciroli_trees[[i]] , Panciroli_trees[[i]]$tip.label[ is.na( match( Panciroli_trees[[i]]$tip.label , data.temp$tree_names ) ) ] )	
  data.tree <- data.temp[ match( tree.temp$tip.label , data.temp$tree_names ) , ]
  rownames( data.tree ) <- data.tree$tree_names
  data.tree[ , "log_HC" ] <- log10( data.tree$Humerus.Circumference )
  data.tree[ , "log_FC" ] <- log10( data.tree$Femur.Circumference )
  data.tree[ , "log_FCHC" ] <- log10( data.tree$Femur.Circumference + data.tree$Humerus.Circumference )
  data.tree[ , "log_BMlit" ] <- log10( data.tree$Body_mass_literature..g. )
  data.tree[ , "log_BMspec" ] <- log10( data.tree$Body.Mass..g. )
  
  pgls.results.list.Campione_lit[[ i ]][[1]] <- phylolm( log_BMlit ~ log_HC , data.tree , model='lambda',phy = Panciroli_trees[[i]]) 
  pgls.results.list.Campione_lit[[ i ]][[2]] <- phylolm( log_BMlit ~ log_FC , data.tree , model='lambda',phy = Panciroli_trees[[i]]) 
  pgls.results.list.Campione_lit[[ i ]][[3]] <- phylolm( log_BMlit ~ log_FCHC , data.tree , model='lambda',phy = Panciroli_trees[[i]]) 
  names( pgls.results.list.Campione_lit[[ i ]] ) <- c("log_HC","log_FC","log_FCHC")
  
  pgls.results.list.Campione_spec[[ i ]][[1]] <- phylolm( log_BMspec ~ log_HC , data.tree , model='lambda',phy = Panciroli_trees[[i]]) 
  pgls.results.list.Campione_spec[[ i ]][[2]] <- phylolm( log_BMspec ~ log_FC , data.tree , model='lambda',phy = Panciroli_trees[[i]]) 
  pgls.results.list.Campione_spec[[ i ]][[3]] <- phylolm( log_BMspec ~ log_FCHC , data.tree , model='lambda',phy = Panciroli_trees[[i]]) 
  names( pgls.results.list.Campione_spec[[ i ]] ) <- c("log_HC","log_FC","log_FCHC")
  
}


#Sigma2 values
#Sigma2 individual traits
#Sigma2 for body mass-literature
merror <- data.tree$log_BMlit * 0.01
names(merror) <- rownames(data.tree)

set.seed(123)
i= sample(1:length(Panciroli_trees),1,replace = F)
tree.temp <- drop.tip( Panciroli_trees[[i]] , Panciroli_trees[[i]]$tip.label[ is.na( match( Panciroli_trees[[i]]$tip.label , data.temp$tree_names ) ) ] )	

Campione_sigma_BMlit <- phylosig(tree.temp, 
                               x = data.tree[tree.temp$tip.label,'log_BMlit'] , 
                               method="lambda",test=T,
                               se=merror )

Campione_sigma_BMlit$sig2 #0.01858

#Sigma2 for body mass-specimen
merror <- data.tree$log_BMspec * 0.01
names(merror) <- rownames(data.tree)

set.seed(123)
i= sample(1:length(Panciroli_trees),1,replace = F)
tree.temp <- drop.tip( Panciroli_trees[[i]] , Panciroli_trees[[i]]$tip.label[ is.na( match( Panciroli_trees[[i]]$tip.label , data.temp$tree_names ) ) ] )	

Campione_sigma_BMspec <- phylosig(tree.temp, 
                                 x = data.tree[tree.temp$tip.label,'log_BMspec'] , 
                                 method="lambda",test=T,
                                 se=merror )

Campione_sigma_BMspec$sig2 #0.01851


#Sigma2 for body humerus circumference
merror <- data.tree$log_HC * 0.01
names(merror) <- rownames(data.tree)

set.seed(123)
i= sample(1:length(Panciroli_trees),1,replace = F)
tree.temp <- drop.tip( Panciroli_trees[[i]] , Panciroli_trees[[i]]$tip.label[ is.na( match( Panciroli_trees[[i]]$tip.label , data.temp$tree_names ) ) ] )	

Campione_sigma_HC <- phylosig(tree.temp, 
                                  x = data.tree[tree.temp$tip.label,'log_HC'] , 
                                  method="lambda",test=T,
                                  se=merror )

Campione_sigma_HC$sig2 #0.00183


#Sigma2 for body femur circumference
merror <- data.tree$log_FC * 0.01
names(merror) <- rownames(data.tree)

set.seed(123)
i= sample(1:length(Panciroli_trees),1,replace = F)
tree.temp <- drop.tip( Panciroli_trees[[i]] , Panciroli_trees[[i]]$tip.label[ is.na( match( Panciroli_trees[[i]]$tip.label , data.temp$tree_names ) ) ] )	

Campione_sigma_FC <- phylosig(tree.temp, 
                              x = data.tree[tree.temp$tip.label,'log_FC'] , 
                              method="lambda",test=T,
                              se=merror )

Campione_sigma_FC$sig2 #0.0


#Sigma2 for body femur+humerus circumference
merror <- data.tree$log_FCHC * 0.01
names(merror) <- rownames(data.tree)

set.seed(123)
i= sample(1:length(Panciroli_trees),1,replace = F)
tree.temp <- drop.tip( Panciroli_trees[[i]] , Panciroli_trees[[i]]$tip.label[ is.na( match( Panciroli_trees[[i]]$tip.label , data.temp$tree_names ) ) ] )	

Campione_sigma_FCHC <- phylosig(tree.temp, 
                              x = data.tree[tree.temp$tip.label,'log_FCHC'] , 
                              method="lambda",test=T,
                              se=merror )

Campione_sigma_FCHC$sig2 #0.00164

#Sigma2 relationships
#Sigma2 values for body mass-literature vs. humerus circumference
mean(sapply ( pgls.results.list.Campione_lit, function(x) x$log_HC$sigma2))

#Sigma2 values for body mass-literature vs. femur circumference
mean(sapply ( pgls.results.list.Campione_lit, function(x) x$log_FC$sigma2))

#Sigma2 values for body mass-literature vs. humerus+femur circumference
mean(sapply ( pgls.results.list.Campione_lit, function(x) x$log_FCHC$sigma2))


#Sigma2 values for body mass-measured vs. humerus circumference
mean(sapply ( pgls.results.list.Campione_spec, function(x) x$log_HC$sigma2))

#Sigma2 values for body mass-measured vs. femur circumference
mean(sapply ( pgls.results.list.Campione_spec, function(x) x$log_FC$sigma2))

#Sigma2 values for body mass-measured vs. humerus+femur circumference
mean(sapply ( pgls.results.list.Campione_spec, function(x) x$log_FCHC$sigma2))




#Campione & Evans dataset (amniotes)

Campione_data <- read.csv('CampioneEvansData.csv',header = T, row.names = NULL)
  head(Campione_data)
Campione_tree <- read.tree('CampioneEvans_tree_timetree.nwk')

rownames(Campione_data) <- gsub(' ','_',Campione_data$Species)

#keep taxa on the tree

Campione_data <- Campione_data[Campione_tree$tip.label,]
Campione_data$logBM <- log10(Campione_data$Body.Mass..g.)
Campione_data$logHL <- log10(Campione_data$Humerus.Length)
Campione_data$logHC <- log10(Campione_data$Humerus.Circumference)
Campione_data$logFL <- log10(Campione_data$Femur.Length)
Campione_data$logFC <- log10(Campione_data$Femur.Circumference)
Campione_data$logFCHC <- log10(Campione_data$Femur.Circumference+Campione_data$Humerus.Circumference)

#PGLS
pgls.Campione.data <- list()

vars <- c('logHL','logHC','logFL','logFC','logFCHC')

for (i in 1:length(vars)){
  
  fm <- paste0('logBM~',vars[i])
  fm <- as.formula(fm)
  
  pgls.Campione.data[[i]] <- phylolm(fm, data=Campione_data,phy=Campione_tree,
                                     model='lambda',boot = 100)
  
}
names(pgls.Campione.data) <- vars

pgls.Campione.data$logFLHL <- phylolm(logFL~logHL, data=Campione_data,phy=Campione_tree,
                                      model='lambda',boot = 100)


#Sigma2 values
#Sigma2 individual traits
#Sigma2 for body mass
merror <- Campione_data$logBM * 0.01
names(merror) <- rownames(Campione_data)


Campione_sig2_BM <- phylosig(Campione_tree, 
                                 x = Campione_data[Campione_tree$tip.label,'logBM'] , 
                                 method="lambda",test=T,
                                 se=merror )

Campione_sig2_BM$sig2 #0.0125


#Sigma2 for humerus length
merror <- Campione_data$logHL * 0.01
names(merror) <- rownames(Campione_data)


Campione_sig2_HL <- phylosig(Campione_tree, 
                             x = Campione_data[Campione_tree$tip.label,'logHL'] , 
                             method="lambda",test=T,
                             se=merror )

Campione_sig2_HL$sig2 #0.00066


#Sigma2 for humerus circumference
merror <- Campione_data$logHC * 0.01
names(merror) <- rownames(Campione_data)


Campione_sig2_HC <- phylosig(Campione_tree, 
                             x = Campione_data[Campione_tree$tip.label,'logHC'] , 
                             method="lambda",test=T,
                             se=merror )

Campione_sig2_HC$sig2 #0.00184


#Sigma2 for femur length
merror <- Campione_data$logFL * 0.01
names(merror) <- rownames(Campione_data)


Campione_sig2_FL <- phylosig(Campione_tree, 
                             x = Campione_data[Campione_tree$tip.label,'logFL'] , 
                             method="lambda",test=T,
                             se=merror )

Campione_sig2_FL$sig2 #0.001016

#Sigma2 for femur circumference
merror <- Campione_data$logFC * 0.01
names(merror) <- rownames(Campione_data)


Campione_sig2_FC <- phylosig(Campione_tree, 
                             x = Campione_data[Campione_tree$tip.label,'logFC'] , 
                             method="lambda",test=T,
                             se=merror )

Campione_sig2_FC$sig2 #0.000486


#Sigma2 for femur+humerus circumference
merror <- Campione_data$logFCHC * 0.01
names(merror) <- rownames(Campione_data)


Campione_sig2_FCHC <- phylosig(Campione_tree, 
                             x = Campione_data[Campione_tree$tip.label,'logFCHC'] , 
                             method="lambda",test=T,
                             se=merror )

Campione_sig2_FCHC$sig2 #0.00167


#Sigma2 relationships
#Body mass vs. humerus length
pgls.Campione.data$logHL$sigma2

#Body mass vs. humerus circumference
pgls.Campione.data$logHC$sigma2

#Body mass vs. femur length
pgls.Campione.data$logFL$sigma2

#Body mass vs. femur circumference
pgls.Campione.data$logFC$sigma2

#Body mass vs. femur+humerus circumference
pgls.Campione.data$logFCHC$sigma2

#Femur length vs. humerus length
pgls.Campione.data$logFLHL$sigma2


#Benson et al. (2018) dataset #Dinosaurs

Benson_data <- read.csv('Benson2018Dataset.csv',header = T,row.names = NULL)
  head(Benson_data)

Benson_trees <- read.tree('Dinosaur trees cal3 Benson 2017.tre')  
    
#Prune data to include only taxa on the trees

Benson_data <- Benson_data[ Benson_data$Name_in_tree %in% Benson_trees[[1]]$tip.label , ]
  head(Benson_data)
rownames(Benson_data) <- Benson_data$Name_in_tree

Benson_data$logBM <- log10(Benson_data$Mass)
Benson_data$logHL <- log10(Benson_data$HL)
Benson_data$logHC <- log10(Benson_data$HCe)
Benson_data$logFL <- log10(Benson_data$FL)
Benson_data$logFC <- log10(Benson_data$FCe)
Benson_data$logFCHC <- log10(Benson_data$FCe+Benson_data$HCe)


pgls.Benson.data <- list()

for ( i in 1:length(Benson_trees)){
  
  Benson.tmp <- multi2di( Benson_trees[[i]])
  Benson.tmp$edge.length <- Benson.tmp$edge.length+0.001
  
pgls.tmp <- list()

  pgls.tmp$HL <- phylolm(logBM~logHL,data=Benson_data, 
                                    phy=Benson.tmp,model='lambda')  
  pgls.tmp$HC <- phylolm(logBM~logHC,data=Benson_data, 
                                    phy=Benson.tmp,model='lambda')  
  pgls.tmp$FL <- phylolm(logBM~logFL,data=Benson_data, 
                                    phy=Benson.tmp,model='lambda')  
  pgls.tmp$FC <-  phylolm(logBM~logFC,data=Benson_data, 
                                    phy=Benson.tmp,model='lambda')  
  pgls.tmp$FCHC <- phylolm(logBM~logFCHC,data=Benson_data, 
                                    phy=Benson.tmp,model='lambda')  
  pgls.tmp$FLHL <- phylolm(logFL~logHL,data=Benson_data, 
                           phy=Benson.tmp,model='lambda')  
  
  pgls.Benson.data[[i]] <- pgls.tmp

  
}


#Sigma2 values
#Sigma2 individual traits
#Sigma2 for body mass
merror <- Benson_data$logBM * 0.01
names(merror) <- rownames(Benson_data)

set.seed(123)
i= sample(1:length(Benson_trees),1,replace = F)
tree.temp <- Benson_trees[[i]]

Benson_sigma_BM <- phylosig(tree.temp, 
                               x = Benson_data[tree.temp$tip.label,'logBM'] , 
                               method="lambda",test=T,
                               se=merror )

Benson_sigma_BM$sig2 #0.0369

#Sigma2 for humerus length
merror <- Benson_data$logHL * 0.01
names(merror) <- rownames(Benson_data)

set.seed(123)
i= sample(1:length(Benson_trees),1,replace = F)
tree.temp <- Benson_trees[[i]]

Benson_sigma_HL <- phylosig(tree.temp, 
                            x = Benson_data[tree.temp$tip.label,'logHL'] , 
                            method="lambda",test=T,
                            se=merror )

Benson_sigma_HL$sig2 #0.00307


#Sigma2 for humerus circumference
merror <- Benson_data$logHC * 0.01
names(merror) <- rownames(Benson_data)

set.seed(123)
i= sample(1:length(Benson_trees),1,replace = F)
tree.temp <- Benson_trees[[i]]

Benson_sigma_HC <- phylosig(tree.temp, 
                            x = Benson_data[tree.temp$tip.label,'logHC'] , 
                            method="lambda",test=T,
                            se=merror )

Benson_sigma_HC$sig2 #0.00254


#Sigma2 for femur length
merror <- Benson_data$logFL * 0.01
names(merror) <- rownames(Benson_data)

set.seed(123)
i= sample(1:length(Benson_trees),1,replace = F)
tree.temp <- Benson_trees[[i]]

Benson_sigma_FL <- phylosig(tree.temp, 
                            x = Benson_data[tree.temp$tip.label,'logFL'] , 
                            method="lambda",test=T,
                            se=merror )

Benson_sigma_FL$sig2 #0.00413


#Sigma2 for femur circumference
merror <- Benson_data$logFC * 0.01
names(merror) <- rownames(Benson_data)

set.seed(123)
i= sample(1:length(Benson_trees),1,replace = F)
tree.temp <- Benson_trees[[i]]

Benson_sigma_FC <- phylosig(tree.temp, 
                            x = Benson_data[tree.temp$tip.label,'logFC'] , 
                            method="lambda",test=T,
                            se=merror )

Benson_sigma_FC$sig2 #0.00459


#Sigma2 for femur+humerus circumference
merror <- Benson_data$logFCHC * 0.01
names(merror) <- rownames(Benson_data)

set.seed(123)
i= sample(1:length(Benson_trees),1,replace = F)
tree.temp <- Benson_trees[[i]]

Benson_sigma_FCHC <- phylosig(tree.temp, 
                            x = Benson_data[tree.temp$tip.label,'logFCHC'] , 
                            method="lambda",test=T,
                            se=merror )

Benson_sigma_FCHC$sig2 #0.002248

#Sigma2 in the relationships
#Body mass vs. humerus length
mean(sapply ( pgls.Benson.data, function(x) x$HL$sigma2))

#Body mass vs. humerus circumference
mean(sapply ( pgls.Benson.data, function(x) x$HC$sigma2))

#Body mass vs. femur length
mean(sapply ( pgls.Benson.data, function(x) x$FL$sigma2))

#Body mass vs. femur circumference
mean(sapply ( pgls.Benson.data, function(x) x$FC$sigma2))

#Body mass vs. femur+humerus circumference
mean(sapply ( pgls.Benson.data, function(x) x$FCHC$sigma2))

#Femur length vs. humerus length
mean(sapply ( pgls.Benson.data, function(x) x$FLHL$sigma2))



#Iijima et al. (2018): extant crocs dataset
crocs <- read.csv('Iijima2018_crocs.csv',header = T,sep=';')
  head(crocs)
#get species means

crocs_dat <- data.frame('HL'=rep(NA,length(unique(crocs$Taxon))),
                    'FL'=rep(NA,length(unique(crocs$Taxon))),
                    'Trunk'=rep(NA,length(unique(crocs$Taxon))),
                    row.names = unique(crocs$Taxon))  
  
for ( i in 1:nrow(crocs_dat)){
  sp <- rownames(crocs_dat)[i]
  
  crocs_dat[i,] <- apply(crocs[crocs$Taxon == sp, c('Humerus.length','Femur.length','Trunk.length..C3.D15.')],
                         2,function(x) mean(x,na.rm=T))
  
  
}

crocs_dat <- log10(crocs_dat)
rownames(crocs_dat) <- gsub(' ','_',rownames(crocs_dat))

#Crocs tree from timetree.org
croc_tree <- read.tree('crocs_tree_timetree.nwk')

crocs_dat <- crocs_dat[croc_tree$tip.label,]


pgls.crocs <- list()

pgls.crocs$HL <- phylolm(Trunk~HL,data=crocs_dat,phy=croc_tree,model = 'lambda',boot=1000)
pgls.crocs$FL <- phylolm(Trunk~FL,data=crocs_dat,phy=croc_tree,model = 'lambda',boot=1000)
pgls.crocs$FLHL <- phylolm(FL~HL,data=crocs_dat,phy=croc_tree,model = 'lambda',boot=1000)


#Sigma2 values

#Sigma2 for Trunk length
merror <- crocs_dat$Trunk * 0.01
names(merror) <- rownames(crocs_dat)


crocs_sig2_TL <- phylosig(croc_tree, 
                               x = crocs_dat[croc_tree$tip.label,'Trunk'] , 
                               method="lambda",test=T,
                               se=merror )

crocs_sig2_TL$sig2 #0.00128


#Sigma2 for Humerus length
merror <- crocs_dat$HL * 0.01
names(merror) <- rownames(crocs_dat)


crocs_sig2_HL <- phylosig(croc_tree, 
                          x = crocs_dat[croc_tree$tip.label,'HL'] , 
                          method="lambda",test=T,
                          se=merror )

crocs_sig2_HL$sig2 #0.00105


#Sigma2 for Femur length
merror <- crocs_dat$FL * 0.01
names(merror) <- rownames(crocs_dat)


crocs_sig2_FL <- phylosig(croc_tree, 
                          x = crocs_dat[croc_tree$tip.label,'FL'] , 
                          method="lambda",test=T,
                          se=merror )

crocs_sig2_FL$sig2 #0.00105



#Sigma2 relationships
#Trunk length vs. humerus length
pgls.crocs$HL$sigma2

#Trunk length vs. femur length
pgls.crocs$FL$sigma2

#Femur length vs. humerus length
pgls.crocs$FLHL$sigma2



#Field et al. (2013): birds dataset

birds <- read.csv('Field2013_birds.csv',header = T, sep=';')
head(birds)

#get species means

birds_spp <- paste0(birds$Genus,"_",birds$species)
birds$sp_name <- birds_spp

  birds_dat <- data.frame('HL'=rep(NA,length(unique(birds_spp))),
                        'FL'=rep(NA,length(unique(birds_spp))),
                        'BM'=rep(NA,length(unique(birds_spp))),
                        row.names = unique(birds_spp))  

for ( i in 1:nrow(birds_dat)){
  sp <- rownames(birds_dat)[i]
  
  birds_dat[i,] <- apply(birds[birds$sp_name == sp, c('humerus_length..mm.','femur_length..mm.',
                                                    'mass..g.')],
                         2,function(x) mean(x,na.rm=T))
  
  
}

birds_dat <- log10(birds_dat)

#Birds tree from timetree.org
bird_tree <- read.tree('bird_tree_timetree.nwk')

birds_dat <- birds_dat[bird_tree$tip.label,]
birds_dat <- birds_dat[complete.cases(birds_dat),]

pgls.birds <- list()

pgls.birds$HL <- phylolm(BM~HL,data=birds_dat,phy=bird_tree,model = 'lambda',boot=1000)
pgls.birds$FL <- phylolm(BM~FL,data=birds_dat,phy=bird_tree,model = 'lambda',boot=1000)
pgls.birds$FLHL <- phylolm(FL~HL,data=birds_dat,phy=bird_tree,model = 'lambda',boot=1000)


#Sigma2 values

#Sigma2 for Body mass
merror <- birds_dat$BM * 0.01
names(merror) <- rownames(birds_dat)


birds_sig2_BM <- phylosig(bird_tree, 
                          x = birds_dat[bird_tree$tip.label,'BM'] , 
                          method="lambda",test=T,
                          se=merror )

birds_sig2_BM$sig2 #0.0059


#Sigma2 for Humerus length
merror <- birds_dat$HL * 0.01
names(merror) <- rownames(birds_dat)


birds_sig2_HL <- phylosig(bird_tree, 
                          x = birds_dat[bird_tree$tip.label,'HL'] , 
                          method="lambda",test=T,
                          se=merror )

birds_sig2_HL$sig2 #0.001114


#Sigma2 for Femur length
merror <- birds_dat$FL * 0.01
names(merror) <- rownames(birds_dat)


birds_sig2_FL <- phylosig(bird_tree, 
                          x = birds_dat[bird_tree$tip.label,'FL'] , 
                          method="lambda",test=T,
                          se=merror )

birds_sig2_FL$sig2 #0.000713


#Sigma2 relationships
#Body mass vs. humerus length
pgls.birds$HL$sigma2

#Body mass vs. femur length
pgls.birds$FL$sigma2

#Femur length vs. humerus length
pgls.birds$FLHL$sigma2


##Compare distributions of sigma2 values across different amniote groups

#Create object to store sigma2 values
#For FL vs. HL
FLHL_sig2 <- list(amniotes=pgls.Campione.data$logFLHL$bootstrap[,'sigma2'],  
                  turtles=LIMBS_pgls.models$general$lambda_ML$bootstrap[,'sigma2'],
                  crocs=pgls.crocs$FLHL$bootstrap[,'sigma2'],
                  dinosaurs=sapply ( pgls.Benson.data, function(x) x$FLHL$sigma2), 
                  birds=pgls.birds$FLHL$bootstrap[,'sigma2'],
                  mammals=sapply ( pgls.results.list.new.data, function(x) x$FLHL$sigma2))


#Comparisons
#with amniotes
t.test(FLHL_sig2$turtles,FLHL_sig2$amniotes)

#with crocs
t.test(FLHL_sig2$turtles,FLHL_sig2$crocs)

#with dinosaurs
t.test(FLHL_sig2$turtles,FLHL_sig2$dinosaurs)

#with birds
t.test(FLHL_sig2$turtles,FLHL_sig2$birds)

#with mammals
t.test(FLHL_sig2$turtles,FLHL_sig2$mammals)


#For Body size proxy vs. HL
HL_sig2 <- list(amniotes=pgls.Campione.data$logHL$bootstrap[,'sigma2'],  
                turtles=HUM_pgls.models$general$lambda_ML$bootstrap[,'sigma2'],
                crocs=pgls.crocs$HL$bootstrap[,'sigma2'],
                dinosaurs=sapply ( pgls.Benson.data, function(x) x$HL$sigma2), 
                birds=pgls.birds$HL$bootstrap[,'sigma2'],
                mammals=sapply ( pgls.results.list.new.data, function(x) x$HL$sigma2))


#Comparisons
#with amniotes
t.test(HL_sig2$turtles,HL_sig2$amniotes)

#with crocs
t.test(HL_sig2$turtles,HL_sig2$crocs)

#with dinosaurs
t.test(HL_sig2$turtles,HL_sig2$dinosaurs)

#with birds
t.test(HL_sig2$turtles,HL_sig2$birds)

#with mammals
t.test(HL_sig2$turtles,HL_sig2$mammals)


#For Body size proxy vs. FL
FL_sig2 <- list(amniotes=pgls.Campione.data$logFL$bootstrap[,'sigma2'],  
                turtles=FEM_pgls.models$general$lambda_ML$bootstrap[,'sigma2'],
                crocs=pgls.crocs$FL$bootstrap[,'sigma2'],
                dinosaurs=sapply ( pgls.Benson.data, function(x) x$FL$sigma2), 
                birds=pgls.birds$FL$bootstrap[,'sigma2'],
                mammals=sapply ( pgls.results.list.new.data, function(x) x$FL$sigma2))

#Comparisons
#with amniotes
t.test(FL_sig2$turtles,FL_sig2$amniotes)

#with crocs
t.test(FL_sig2$turtles,FL_sig2$crocs)

#with dinosaurs
t.test(FL_sig2$turtles,FL_sig2$dinosaurs)

#with birds
t.test(FL_sig2$turtles,FL_sig2$birds)

#with mammals
t.test(FL_sig2$turtles,FL_sig2$mammals)

