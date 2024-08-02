set.seed(123)

library(phytools)
library(paleotree)
library(phylolm)
library(rr2)
library(vioplot)
library(AICcmodavg)


#AICc function

AICc.phylolm<-function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){
  
  if(identical(nobs, NULL)) {n <- length(mod$residuals)} else {n <- nobs}
  LL <- logLik(mod)$logLik
  K <- logLik(mod)$df  #extract correct number of parameters included in model - this includes LM
  if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
  if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
  return(AICc)
}

##Read data

dat <- read.table('SupplementaryFileS1 (raw data).csv',header = T,row.names = NULL, sep=',')
head(dat)

cols <- c('darkorange','greenyellow','#6EB8CC','#D8A068','#8E8375','#7D87C9','gray80',
          'gold','#C68BC4')
names(cols) <- levels(as.factor(dat$Clade))

#To avoid kicking out fossils before the next step, add '1' to their Ratio columns
#dat$Ratio_F[is.na(dat$Ratio_F)] <- 1
#dat$Ratio_M[is.na(dat$Ratio_M)] <- 1

#Subset taxa with unknown body sizes, which will be predicted at the end
predict_subset <- subset(dat,is.na(SCL_mm))

#Subset taxa with known body sizes
size_subset <- subset(dat,!is.na(SCL_mm))

#Subset size data with only specimens >60% (female) AND >60% (male) of maximum SCL
#size_subset <-  subset(size_subset,Ratio_F >= 0.6 | Ratio_M >= 0.6 )
size_subset <-  subset(size_subset,To_use=='yes' )

#Data frame with species means and log-transformed variables
#For Humerus analysis
HUM_subset <- subset(size_subset,!is.na(HL_mm))

HUM_subset <- data.frame('SCL'=log10(aggregate(HUM_subset$SCL_mm~HUM_subset$Correct_taxonomy,FUN=mean)[,2]),
           'HL'=log10(aggregate(HUM_subset$HL_mm~HUM_subset$Correct_taxonomy,FUN=mean)[,2]),
           'Clade'=unlist(aggregate(HUM_subset$Clade~HUM_subset$Correct_taxonomy,FUN=unique)[,2]),
           'Terr_spec'=unlist(aggregate(HUM_subset$Terr_spec~HUM_subset$Correct_taxonomy,FUN=unique)[,2]),
           'Aq_spec'=unlist(aggregate(HUM_subset$Aq_spec~HUM_subset$Correct_taxonomy,FUN=unique)[,2]),
           'Trion'=unlist(aggregate(HUM_subset$Trion~HUM_subset$Correct_taxonomy,FUN=unique)[,2]),
           'Type'=aggregate(HUM_subset$Material~HUM_subset$Correct_taxonomy,FUN=unique)[,2],
           row.names = aggregate(HUM_subset$Material~HUM_subset$Correct_taxonomy,FUN=unique)[,1])



#List with different datasets for clade-specific regressions
HUM_list <- list('general'=HUM_subset,
                 'Chelidae'=subset(HUM_subset,Clade=='Chelidae'),
                 'Chelonioidea'=subset(HUM_subset,Clade=='Chelonioidea'),
                 'Chelydroidea'=subset(HUM_subset,Clade=='Chelydroidea'),
                 'Emysternia'=subset(HUM_subset,Clade=='Emysternia'),
                 'Geoemydidae'=subset(HUM_subset,Clade=='Geoemydidae'),
                 'Pelomedusoides'=subset(HUM_subset,Clade=='Pelomedusoides'),
                 'Testudinidae'=subset(HUM_subset,Clade=='Testudinidae'),
                 'Trionychia'=subset(HUM_subset,Clade=='Trionychia') )

#For Femur analysis
FEM_subset <- subset(size_subset,!is.na(FL_mm))

FEM_subset <- data.frame('SCL'=log10(aggregate(FEM_subset$SCL_mm~FEM_subset$Correct_taxonomy,FUN=mean)[,2]),
                         'FL'=log10(aggregate(FEM_subset$FL_mm~FEM_subset$Correct_taxonomy,FUN=mean)[,2]),
                         'Clade'=unlist(aggregate(FEM_subset$Clade~FEM_subset$Correct_taxonomy,FUN=unique)[,2]),
                         'Terr_spec'=unlist(aggregate(FEM_subset$Terr_spec~FEM_subset$Correct_taxonomy,FUN=unique)[,2]),
                         'Aq_spec'=unlist(aggregate(FEM_subset$Aq_spec~FEM_subset$Correct_taxonomy,FUN=unique)[,2]),
                         'Trion'=unlist(aggregate(FEM_subset$Trion~FEM_subset$Correct_taxonomy,FUN=unique)[,2]),
                         'Type'=aggregate(FEM_subset$Material~FEM_subset$Correct_taxonomy,FUN=unique)[,2],
                         row.names = aggregate(FEM_subset$Material~FEM_subset$Correct_taxonomy,FUN=unique)[,1])

#List with different datasets for clade-specific regressions
FEM_list <- list('general'=FEM_subset,
                 'Chelidae'=subset(FEM_subset,Clade=='Chelidae'),
                 'Chelonioidea'=subset(FEM_subset,Clade=='Chelonioidea'),
                 'Chelydroidea'=subset(FEM_subset,Clade=='Chelydroidea'),
                 'Emysternia'=subset(FEM_subset,Clade=='Emysternia'),
                 'Geoemydidae'=subset(FEM_subset,Clade=='Geoemydidae'),
                 'Pelomedusoides'=subset(FEM_subset,Clade=='Pelomedusoides'),
                 'Testudinidae'=subset(FEM_subset,Clade=='Testudinidae'),
                 'Trionychia'=subset(FEM_subset,Clade=='Trionychia') )


#For Stylopodia analysis
LIMB_subset <- subset(dat,!is.na(HL_mm) & !is.na(FL_mm) )
  LIMB_subset <-  subset(LIMB_subset,To_use=='yes' )

LIMB_subset <- data.frame(
  'HL'=log10(aggregate(LIMB_subset$HL_mm~LIMB_subset$Correct_taxonomy,FUN=mean)[,2]),
  'FL'=log10(aggregate(LIMB_subset$FL_mm~LIMB_subset$Correct_taxonomy,FUN=mean)[,2]),
  'Clade'=unlist(aggregate(LIMB_subset$Clade~LIMB_subset$Correct_taxonomy,FUN=unique)[,2]),
  'Terr_spec'=unlist(aggregate(LIMB_subset$Terr_spec~LIMB_subset$Correct_taxonomy,FUN=unique)[,2]),
  'Aq_spec'=unlist(aggregate(LIMB_subset$Aq_spec~LIMB_subset$Correct_taxonomy,FUN=unique)[,2]),
  'Trion'=unlist(aggregate(LIMB_subset$Trion~LIMB_subset$Correct_taxonomy,FUN=unique)[,2]),
  'Type'=aggregate(LIMB_subset$Material~LIMB_subset$Correct_taxonomy,FUN=unique)[,2],
  row.names = aggregate(LIMB_subset$Material~LIMB_subset$Correct_taxonomy,FUN=unique)[,1])



#List with different datasets for clade-specific regressions
LIMB_list <- list('general'=LIMB_subset,
                  'Chelidae'=subset(LIMB_subset,Clade=='Chelidae'),
                  'Chelonioidea'=subset(LIMB_subset,Clade=='Chelonioidea'),
                  'Chelydroidea'=subset(LIMB_subset,Clade=='Chelydroidea'),
                  'Emysternia'=subset(LIMB_subset,Clade=='Emysternia'),
                  'Geoemydidae'=subset(LIMB_subset,Clade=='Geoemydidae'),
                  'Pelomedusoides'=subset(LIMB_subset,Clade=='Pelomedusoides'),
                  'Testudinidae'=subset(LIMB_subset,Clade=='Testudinidae'),
                  'Trionychia'=subset(LIMB_subset,Clade=='Trionychia') )


#Read tree
#Sterli et al. 2018 topology
tree <- read.tree('Sterli_MCCT.tre')

#Add missing labels
tree_full <- bindPaleoTip(tree,tipLabel = 'Melanochelys_tricarinata',
                          nodeAttach = which(tree$tip.label=='Melanochelys_trijuga'),tipAge = 0,
                            positionBelow = 1)


tree_full <- bindPaleoTip(tree_full,tipLabel = 'Leiochelys_tokaryki',
                          nodeAttach = getMRCA(tree_full,c('Kinosternon_flavescens','Claudius_angustatus')),
                          tipAge = 66,
                          positionBelow = 1)


tree_full <- bindPaleoTip(tree_full,tipLabel = 'Oligochelone_rupelensis',
                          nodeAttach = getMRCA(tree_full,c('Osonachelus_decorata','Allopleuron_hofmanni')),
                          tipAge = 30,
                          positionBelow = 2)

#hypothetical position of 'Trionyx' singularis as stem-trionychid
tree_full <- bindPaleoTip(tree_full,tipLabel = 'Trionyx_singularis',
                          nodeAttach = getMRCA(tree_full,c('Trionyx_triunguis','Drazinderetes_tethyensis')),
                          tipAge = 63.3,
                          positionBelow = 20)


#hypothetical position of Ashleychelys palmeri as crownward-most stem-cheloniid
tree_full <- bindPaleoTip(tree_full,tipLabel = 'Ashleychelys_palmeri',
                          nodeAttach = getMRCA(tree_full,c('Chelonia_mydas','Lepidochelys_kempii')),
                          tipAge = 23.03,
                          positionBelow = 6)


#hypothetical position of Yuchelys nanyangensis as sister taxon to all nanhsiungchelyids
tree_full <- bindPaleoTip(tree_full,tipLabel = 'Yuchelys_nanyangensis',
                          nodeAttach = getMRCA(tree_full,c('Basilemys_variolosa','Nanhsiungchelys_wuchingensis')),
                          tipAge = 89,
                          positionBelow = 6)


#hypothetical position of Terlinguachelys fischbecki as sister taxon to all post-Santanachelys protostegids
tree_full <- bindPaleoTip(tree_full,tipLabel = 'Terlinguachelys_fischbecki',
                          nodeAttach = getMRCA(tree_full,c('Archelon_ischyros','Protostega_gigas')),
                          tipAge = 72,
                          positionBelow = 2)

#hypothetical position of Ptychogaster sp. as sister taxon to the Cuora-Mauremys group
tree_full <- bindPaleoTip(tree_full,tipLabel = 'Ptychogaster_sp',
                          nodeAttach = getMRCA(tree_full,c('Mauremys_leprosa','Echmatemys_sp')),
                          tipAge = 29,
                          positionBelow = 2)

tree_full <-  keep.tip(tree_full, 
                       unique(c(rownames(HUM_subset),
                                rownames(FEM_subset),
                                rownames(LIMB_subset)) ) )
  tree_full$root.time <- max(diag(vcv(tree_full)))

plot(tree_full,cex=0.4, no.margin = T)
dev.off()

#Rearrange tree to insert protostegids as stem-chelonioids
node_protost <- getMRCA(tree_full,tip = c('Archelon_ischyros','Santanachelys_gaffneyi','Rhinochelys_nammourensis','Desmatochelys_lowi'))
tree_protost <- extract.clade(tree_full, node = node_protost,root.edge = 1)
tree_protost$root.time <- max(diag(vcv(tree_protost)))

#drop protostegids from original tree
tree_full_tmp <- drop.tip(tree_full,tip=tree_protost$tip.label)

dateNodes(tree_full_tmp)[getMRCA(tree_full_tmp,c('Dermochelys_coriacea','Chelonia_mydas'))]

tree_full <-  bind.tree(x = tree_full_tmp , y = tree_protost, 
          where = getMRCA(tree_full_tmp,c('Dermochelys_coriacea','Chelonia_mydas')),
          position = 27)

tree_full$root.time <- max(diag(vcv(tree_full)))

plot(tree_full,cex=0.4)
axisPhylo(cex.axis=0.5)
abline(v=tree_full$root.time - c(23,66,145),lty=3)

dev.off()

#Prune tree to match datasets
  
#Humerus
HUM_trees <- list('general'=keep.tip(tree_full,rownames(HUM_list[['general']])),
                  'Chelidae'=keep.tip(tree_full,rownames(HUM_list[['Chelidae']])),
                  'Chelonioidea'=keep.tip(tree_full,rownames(HUM_list[['Chelonioidea']])),
                  'Chelydroidea'=keep.tip(tree_full,rownames(HUM_list[['Chelydroidea']])),
                  'Emysternia'=keep.tip(tree_full,rownames(HUM_list[['Emysternia']])),
                  'Geoemydidae'=keep.tip(tree_full,rownames(HUM_list[['Geoemydidae']])),
                  'Pelomedusoides'=keep.tip(tree_full,rownames(HUM_list[['Pelomedusoides']])),
                  'Testudinidae'=keep.tip(tree_full,rownames(HUM_list[['Testudinidae']])),
                  'Trionychia'=keep.tip(tree_full,rownames(HUM_list[['Trionychia']])) )

#Adjust root.time
for ( i in 1:length(HUM_trees)){HUM_trees[[i]]$root.time <- max(diag(vcv(HUM_trees[[i]]))) }


#Femur
FEM_trees <- list('general'=keep.tip(tree_full,rownames(FEM_list[['general']])),
                  'Chelidae'=keep.tip(tree_full,rownames(FEM_list[['Chelidae']])),
                  'Chelonioidea'=keep.tip(tree_full,rownames(FEM_list[['Chelonioidea']])),
                  'Chelydroidea'=keep.tip(tree_full,rownames(FEM_list[['Chelydroidea']])),
                  'Emysternia'=keep.tip(tree_full,rownames(FEM_list[['Emysternia']])),
                  'Geoemydidae'=keep.tip(tree_full,rownames(FEM_list[['Geoemydidae']])),
                  'Pelomedusoides'=keep.tip(tree_full,rownames(FEM_list[['Pelomedusoides']])),
                  'Testudinidae'=keep.tip(tree_full,rownames(FEM_list[['Testudinidae']])),
                  'Trionychia'=keep.tip(tree_full,rownames(FEM_list[['Trionychia']])) )

#Adjust root.time
for ( i in 1:length(FEM_trees)){FEM_trees[[i]]$root.time <- max(diag(vcv(FEM_trees[[i]]))) }


#Stylopodia
LIMB_trees <- list('general'=keep.tip(tree_full,rownames(LIMB_list[['general']])),
                   'Chelidae'=keep.tip(tree_full,rownames(LIMB_list[['Chelidae']])),
                   'Chelonioidea'=keep.tip(tree_full,rownames(LIMB_list[['Chelonioidea']])),
                   'Chelydroidea'=keep.tip(tree_full,rownames(LIMB_list[['Chelydroidea']])),
                   'Emysternia'=keep.tip(tree_full,rownames(LIMB_list[['Emysternia']])),
                   'Geoemydidae'=keep.tip(tree_full,rownames(LIMB_list[['Geoemydidae']])),
                   'Pelomedusoides'=keep.tip(tree_full,rownames(LIMB_list[['Pelomedusoides']])),
                   'Testudinidae'=keep.tip(tree_full,rownames(LIMB_list[['Testudinidae']])),
                   'Trionychia'=keep.tip(tree_full,rownames(LIMB_list[['Trionychia']])) )

#Adjust root.time
for ( i in 1:length(LIMB_trees)){LIMB_trees[[i]]$root.time <- max(diag(vcv(LIMB_trees[[i]]))) }



#Some counting stats

fossils <- subset(dat,Material=='fossil')

#Fossils with "SCL" (51)
unique(subset(fossils, !is.na(SCL_mm))$Correct_taxonomy)

#Fossils with "SCL" AND "HL" AND "FL" (29)
unique(subset(subset(fossils, !is.na(SCL_mm)) , !is.na(HL_mm) & !is.na(FL_mm) )$Correct_taxonomy)

#Fossils with "SCL" AND "HL" OR "SCL" AND "FL" (22)
unique(subset(subset(fossils, !is.na(SCL_mm)) , is.na(HL_mm) | is.na(FL_mm) )$Correct_taxonomy)

#Fossils with "SCL" AND "HL" BUT NOT "FL" ()
unique(subset(subset(fossils, !is.na(SCL_mm)) , !is.na(HL_mm) & is.na(FL_mm) )$Correct_taxonomy)

#Fossils with "SCL" AND "FL" BUT NOT "HL" ()
unique(subset(subset(fossils, !is.na(SCL_mm)) , is.na(HL_mm) & !is.na(FL_mm) )$Correct_taxonomy)


#Fossils without "SCL", but "HL" or "FL" (21)
unique(subset(subset(fossils, is.na(SCL_mm)) , !is.na(FL_mm) | !is.na(HL_mm))$Correct_taxonomy)

#Fossils without "SCL", but "HL" AND "FL" (15)
unique(subset(subset(fossils, is.na(SCL_mm)) , !is.na(FL_mm) & !is.na(HL_mm))$Correct_taxonomy)

#Fossils without "SCL", but "HL" AND NO "FL" ()
unique(subset(subset(fossils, is.na(SCL_mm)) , is.na(FL_mm) & !is.na(HL_mm))$Correct_taxonomy)

#Fossils without "SCL", but "FL" AND NO "HL" ()
unique(subset(subset(fossils, is.na(SCL_mm)) , !is.na(FL_mm) & is.na(HL_mm))$Correct_taxonomy)

fossils2 <- data.frame(SCL=rep(NA,length(unique(fossils$Correct_taxonomy))),
                      HL=rep(NA,length(unique(fossils$Correct_taxonomy))),
                      FL=rep(NA,length(unique(fossils$Correct_taxonomy))),
                      row.names = unique(fossils$Correct_taxonomy)
                      )

for ( i in 1:nrow(fossils2)){
  for ( j in 1:nrow(fossils)){
    
tax <-  fossils$Correct_taxonomy[j]

if ( tax == rownames(fossils2)[i])

fossils2[i,] <-   apply(fossils[fossils$Correct_taxonomy==tax, c('SCL_mm','HL_mm','FL_mm')],
                        2,function(x)mean(x,na.rm=T))

else next
      
  }
}

fossils2

#Fossils with SCL, FL and HL (33)
nrow(subset(fossils2, !is.na(SCL) & !is.na(HL) & !is.na(FL)))

#Fossils with SCL AND HL, but NO FL (15)
nrow(subset(fossils2, !is.na(SCL) & !is.na(HL) & is.na(FL)))

#Fossils with SCL AND FL, but NO HL (2*) Meiolania has HL too
nrow(subset(fossils2, !is.na(SCL) & is.na(HL) & !is.na(FL)))

#Fossils without SCL, but with HL AND FL (15)
nrow(subset(subset(fossils2,is.na(SCL)), !is.na(HL) & !is.na(FL)))

#Fossils without SCL, but with HL and NO FL (12*) *Meiolania has SCL
nrow(subset(subset(fossils2,is.na(SCL)), !is.na(HL) & is.na(FL)))

#Fossils without SCL, but with FL and NO HL (6*) *Meiolania has SCL
nrow(subset(subset(fossils2,is.na(SCL)), is.na(HL) & !is.na(FL)))
