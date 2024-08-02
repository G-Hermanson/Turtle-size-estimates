##STYLOPODIA PLOT##

{
  pdf('Stylopodia_plot.pdf',width=7, height =6.8, useDingbats = F, family = 'sans')
  
  layout ( matrix ( c(1,2,
                      3,4,
                      5,5),3,byrow = T))
  #layout.show(5)
  
  par(mar=c(2,4.5,2,3.5))
  
  #A
  plot (LIMB_list$general[,1:2],
        ylab='log10 Femur length (mm)',xlab='log10 Humerus length (mm)',
        cex=1,
        pch=NA,
        col= NA,
        xaxt='n',yaxt='n')
  axis(1,at=c(1.5,2,2.5))
  axis(2,at=c(1.5,2,2.5,3))
  
  # confidence interval
  LIMBS_CI_general <- summary(LIMBS_best_models$general)$coefficients[,c('lowerbootCI','upperbootCI')]
  LIMBS_CI_general <- t(LIMBS_CI_general)
  colnames(LIMBS_CI_general)[1] <- 'Int'
  LIMBS_CI_general[,1] <- coef(LIMBS_best_models$general)[1]
  
  conf.HL <- sort(LIMB_list$general$HL)
  preds = apply ( LIMBS_CI_general[,-3], 1, function(x) as.matrix(LIMBS_best_models$general$X[,-3] ) %*% x  )
  
  # fill in area between regression line and confidence interval
  polygon(x=c(rev(conf.HL), conf.HL), y=c(rev(sort(preds[,2])), sort(preds[,1])), 
          col = rgb(225/255, 225/255, 225/255, 0.5), border = NA)
  
  points(LIMB_list$general[,1:2],
         cex=1,
         pch=ifelse(LIMB_list$general$Type=='fossil',15,16),
         col= sapply(cols['stem'],make.transparent,alpha=0.95))
  
  abline(LIMBS_best_models$general,lty=2,lwd=1)
  
  points(LIMB_list$general[LIMB_list$general$Clade=='stem',1:2],
         cex=1,
         pch=ifelse(LIMB_list$general[LIMB_list$general$Clade=='stem','Type']=='fossil',22,21),
         bg=make.transparent('mediumpurple4' ,0.6),
         col=make.transparent('mediumpurple4',0.6))
  #colour terrestrial taxa
  points(LIMB_list$general[LIMB_list$general$Terr_spec=='1',1:2],
         cex=1,
         pch=ifelse(LIMB_list$general[LIMB_list$general$Terr_spec=='1','Type']=='fossil',22,21),
         col=make.transparent('black',0.7), bg=NA)
  
  legend('topleft',legend=c('PGLS fit','stem turtles'),bty='n',
         lty=c(2,NA),lwd=1.2,cex=0.9,col=c('black','mediumpurple4'),pch=c(NA,15),
         pt.cex=2)
  legend('bottomright',legend=c('extant','fossil'),
         pch=c(16,15),pt.cex=2,cex=0.9,
         col='gray95',
         bty='n',ncol=2)
  
  mtext('A',side = 3,cex = 1.2,adj = -0.17)
  
  #B
  plot (LIMB_list$general[,1:2],
        ylab='log10 Femur length (mm)',xlab='log10 Humerus length (mm)',
        cex=1,
        pch=NA,
        col= NA,
        xaxt='n',yaxt='n')
  axis(1,at=c(1.5,2,2.5))
  axis(2,at=c(1.5,2,2.5))
  
  # confidence interval
  conf.HL <- sort(LIMB_list$general$HL)
  preds = apply ( LIMBS_CI_general[,-3], 1, function(x) as.matrix(LIMBS_best_models$general$X[,-3] ) %*% x  )
  
  # fill in area between regression line and confidence interval
  polygon(x=c(rev(conf.HL), conf.HL), y=c(rev(sort(preds[,2])), sort(preds[,1])), 
          col = rgb(225/255, 225/255, 225/255, 0.5), border = NA)
  
  points(LIMB_list$general[,1:2],
         cex=1,
         pch=ifelse(LIMB_list$general$Type=='fossil',15,16),
         col= sapply(cols['stem'],make.transparent,alpha=0.95))
  abline(LIMBS_best_models$general,lty=2,lwd=1)
  
  
  points(LIMB_list$general[LIMB_list$general$Clade=='Pelomedusoides',1:2],
         cex=1,
         pch=ifelse(LIMB_list$general[LIMB_list$general$Clade=='Pelomedusoides','Type']=='fossil',22,21),
         bg=make.transparent( cols['Pelomedusoides'],0.6),
         col=make.transparent( cols['Pelomedusoides'],0.6))
  
  points( sort(predict(LIMBS_best_models$Pelomedusoides))~
            LIMB_list$Pelomedusoides[names(sort(predict(LIMBS_best_models$Pelomedusoides))),'HL'],type='l',
          lwd=2,
          col=cols['Pelomedusoides'])
  
  points(LIMB_list$general[LIMB_list$general$Clade=='Chelidae',1:2],
         cex=1,
         pch=ifelse(LIMB_list$general[LIMB_list$general$Clade=='Chelidae','Type']=='fossil',22,21),
         bg=make.transparent( cols['Chelidae'],0.6),
         col=make.transparent( cols['Chelidae'],0.6))
  
  points( sort(predict(LIMBS_best_models$Chelidae))~
            LIMB_list$Chelidae[names(sort(predict(LIMBS_best_models$Chelidae))),'HL'],type='l',lwd=2,
          col=cols['Chelidae'])
  
  legend('topleft',legend=c('Chelidae','Pelomedusoides'),
         pch=21,cex=0.9,pt.cex=2,
         pt.bg=cols[c('Chelidae','Pelomedusoides')],ncol=1,bty='n',
         col=sapply(cols[c('Chelidae','Pelomedusoides')],make.transparent,0.6))
  
  legend('bottomright',legend=c('extant','fossil'),
         pch=c(16,15),pt.cex=2,cex=0.9,
         col='gray95',
         bty='n',ncol=2)
  
  mtext('B',side = 3,cex = 1.2,adj = -0.17)
  
  #C
  plot (LIMB_list$general[,1:2],
        ylab='log10 Femur length (mm)',xlab='log10 Humerus length (mm)',
        cex=1,
        pch=NA,
        col= NA,
        xaxt='n',yaxt='n')
  axis(1,at=c(1.5,2,2.5))
  axis(2,at=c(1.5,2,2.5))
  
  # confidence interval
  conf.HL <- sort(LIMB_list$general$HL)
  preds = apply ( LIMBS_CI_general[,-3], 1, function(x) as.matrix(LIMBS_best_models$general$X[,-3] ) %*% x  )
  
  # fill in area between regression line and confidence interval
  polygon(x=c(rev(conf.HL), conf.HL), y=c(rev(sort(preds[,2])), sort(preds[,1])), 
          col = rgb(225/255, 225/255, 225/255, 0.5), border = NA)
  
  points(LIMB_list$general[,1:2],
         cex=1,
         pch=ifelse(LIMB_list$general$Type=='fossil',15,16),
         col= sapply(cols['stem'],make.transparent,alpha=0.95))
  abline(LIMBS_best_models$general,lty=2,lwd=1)
  
  points(LIMB_list$general[LIMB_list$general$Clade=='Chelydroidea',1:2],
         cex=1,
         pch=ifelse(LIMB_list$general[LIMB_list$general$Clade=='Chelydroidea','Type']=='fossil',22,21),
         bg=make.transparent( cols['Chelydroidea'],0.6),
         col=make.transparent( cols['Chelydroidea'],0.6))
  
  points( sort(predict(LIMBS_best_models$Chelydroidea))~
            LIMB_list$Chelydroidea[names(sort(predict(LIMBS_best_models$Chelydroidea))),'HL'],type='l',
          lwd=2,
          col=cols['Chelydroidea'])
  
  points(LIMB_list$general[LIMB_list$general$Clade=='Trionychia',1:2],
         cex=1,
         pch=ifelse(LIMB_list$general[LIMB_list$general$Clade=='Trionychia','Type']=='fossil',22,21),
         bg=make.transparent( cols['Trionychia'],0.6),
         col=make.transparent( cols['Trionychia'],0.6))
  
  tr <- sort(fitted(LIMBS_best_models$Trionychia))
  points( fitted(lm(tr ~ LIMB_list$Trionychia[names(tr),'HL'])) ~ LIMB_list$Trionychia[names(tr),'HL'],
          type='l', lwd=2, col=cols['Trionychia'])
  
  points(LIMB_list$general[LIMB_list$general$Clade=='Chelonioidea',1:2],
         cex=1,
         pch=ifelse(LIMB_list$general[LIMB_list$general$Clade=='Chelonioidea','Type']=='fossil',22,21),
         bg=make.transparent( cols['Chelonioidea'],0.6),
         col=make.transparent( cols['Chelonioidea'],0.6))
  
  points( sort(predict(LIMBS_best_models$Chelonioidea))~
            LIMB_list$Chelonioidea[names(sort(predict(LIMBS_best_models$Chelonioidea))),'HL'],type='l',
          lwd=2,
          col=cols['Chelonioidea'])
  
  legend('topleft',legend=c('Trionychia','Chelydroidea','Chelonioidea'),
         pch=21,cex=0.9,pt.cex=2,
         pt.bg=cols[c('Trionychia','Chelydroidea','Chelonioidea')],ncol=1,bty='n',
         col=sapply(cols[c('Trionychia','Chelydroidea','Chelonioidea')],make.transparent,0.6)  )
  
  legend('bottomright',legend=c('extant','fossil'),
         pch=c(16,15),pt.cex=2,cex=0.9,
         col='gray95',
         bty='n',ncol=2)
  
  mtext('C',side = 3,cex = 1.2,adj = -0.17)
  
  #D
  plot (LIMB_list$general[,1:2],
        ylab='log10 Femur length (mm)',xlab='log10 Humerus length (mm)',
        cex=1,
        pch=NA,
        col= NA,
        xaxt='n',yaxt='n')
  axis(1,at=c(1.5,2,2.5))
  axis(2,at=c(1.5,2,2.5))
  
  # confidence interval
  conf.HL <- sort(LIMB_list$general$HL)
  preds = apply ( LIMBS_CI_general[,-3], 1, function(x) as.matrix(LIMBS_best_models$general$X[,-3] ) %*% x  )
  
  # fill in area between regression line and confidence interval
  polygon(x=c(rev(conf.HL), conf.HL), y=c(rev(sort(preds[,2])), sort(preds[,1])), 
          col = rgb(225/255, 225/255, 225/255, 0.5), border = NA)
  
  points(LIMB_list$general[,1:2],
         cex=1,
         pch=ifelse(LIMB_list$general$Type=='fossil',15,16),
         col= sapply(cols['stem'],make.transparent,alpha=0.95))
  abline(LIMBS_best_models$general,lty=2,lwd=1)
  
  points(LIMB_list$general[LIMB_list$general$Clade=='Testudinidae',1:2],
         cex=1,
         pch=ifelse(LIMB_list$general[LIMB_list$general$Clade=='Testudinidae','Type']=='fossil',22,21),
         bg=make.transparent( cols['Testudinidae'],0.6),
         col=make.transparent( cols['Testudinidae'],0.6))
  
  points( sort(predict(LIMBS_best_models$Testudinidae))~
            LIMB_list$Testudinidae[names(sort(predict(LIMBS_best_models$Testudinidae))),'HL'],type='l',
          lwd=2,
          col=cols['Testudinidae'])
  
  points(LIMB_list$general[LIMB_list$general$Clade=='Geoemydidae',1:2],
         cex=1,
         pch=ifelse(LIMB_list$general[LIMB_list$general$Clade=='Geoemydidae','Type']=='fossil',22,21),
         bg=make.transparent( cols['Geoemydidae'],0.6),
         col=make.transparent( cols['Geoemydidae'],0.6))
  
  tr <- sort(fitted(LIMBS_best_models$Geoemydidae))
  points( fitted(lm(tr ~ LIMB_list$Geoemydidae[names(tr),'HL'])) ~ LIMB_list$Geoemydidae[names(tr),'HL'],
          type='l', lwd=2, col=cols['Geoemydidae'])
  
  points(LIMB_list$general[LIMB_list$general$Clade=='Emysternia',1:2],
         cex=1,
         pch=ifelse(LIMB_list$general[LIMB_list$general$Clade=='Emysternia','Type']=='fossil',22,21),
         bg=make.transparent( cols['Emysternia'],0.6),
         col=make.transparent( cols['Emysternia'],0.6))
  
  points( sort(predict(LIMBS_best_models$Emysternia))~
            LIMB_list$Emysternia[names(sort(predict(LIMBS_best_models$Emysternia))),'HL'],type='l',
          lwd=2,
          col=cols['Emysternia'])
  
  legend('topleft',legend=c('Emysternia','Geoemydidae','Testudinidae'),
         pch=21,cex=0.9,pt.cex=2,
         pt.bg=cols[c('Emysternia','Geoemydidae','Testudinidae')],ncol=1,bty='n',
         col=sapply(cols[c('Emysternia','Geoemydidae','Testudinidae')],make.transparent,0.6) )
  
  legend('bottomright',legend=c('extant','fossil'),
         pch=c(16,15),pt.cex=2,cex=0.9,
         col='gray95',
         bty='n',ncol=2)
  
  mtext('D',side = 3,cex = 1.2,adj = -0.17)
  
  #E
  par(mar=c(1.1,4.5,4.1,3.5))
  
  plot_order <- c('stem','Chelidae','Pelomedusoides','Trionychia','Chelonioidea','Chelydroidea',
                  'Emysternia','Geoemydidae','Testudinidae')
  
  cols.tmp <- cols
  cols.tmp['stem'] <- 'mediumpurple4'
  
  plot(as.numeric(as.factor(LIMB_list$general$Clade)) , LIMBS_resid[rownames(LIMB_list$general)] ,
       ylab='Residuals',yaxt=NULL,
       ylim=c(-0.3,0.3),xaxt='n',xlab='',pch=NA,
       xlim=c(-0.5,9.5))
  abline(h=0,col='black',lty=2,lwd=0.8)
  
  for ( i in 1:length(plot_order)){
    
    vioplot(  LIMBS_resid[rownames(LIMB_list$general)][LIMB_list$general$Clade==plot_order[i]] ,
              col='transparent' ,
              border='transparent',
              #border=sapply(cols[plot_order[i]],make.transparent,0.8),
              rectCol=sapply(cols.tmp[plot_order[i]],make.transparent,0.8),pchMed=21,
              colMed2=cols.tmp[plot_order[i]],colMed=sapply(cols.tmp[plot_order[i]],make.transparent,0.7),
              lineCol=cols.tmp[plot_order[i]],cex=1.5,
              horizontal = F,xlim=c(-0.3,0.3),add=T,at=i)  
    
  }
  
  for ( i in 2:length(plot_order)){
    
    quants <- quantile(sort(resid(LIMBS_best_models$general)[rownames(LIMB_list$general)[LIMB_list$general$Clade==plot_order[i]] ]),
                       probs = c(0.25,0.75))
    
    points( rep(i-0.35,
                length(resid(LIMBS_best_models$general)[rownames(LIMB_list$general)[LIMB_list$general$Clade==plot_order[i]] ] )),
            sort(resid(LIMBS_best_models$general)[rownames(LIMB_list$general)[LIMB_list$general$Clade==plot_order[i]] ]),
            type='l',lty=1,col='white',lwd=4 )
    
    points( rep(i-0.35,
                length(resid(LIMBS_best_models$general)[rownames(LIMB_list$general)[LIMB_list$general$Clade==plot_order[i]] ] )),
            sort(resid(LIMBS_best_models$general)[rownames(LIMB_list$general)[LIMB_list$general$Clade==plot_order[i]] ]),
            type='l',lty=1,col=cols.tmp[plot_order[i]] )
    
    segments(y0=quants[1],y1=quants[2],x0=i-0.35,x1=i-0.35,lwd=2.5, col=cols.tmp[plot_order[i]])
    
    points(i-0.35,median(sort(resid(LIMBS_best_models$general)[rownames(LIMB_list$general)[LIMB_list$general$Clade==plot_order[i]] ])),
           pch=21,cex=1.5,
           col=cols.tmp[plot_order[i]],bg='white')
    
  }
  
  legend('bottomleft', legend=c('general trend', 'clade-specific trend'),
         lty=c(1,1), pch=21, col='black', pt.bg = c('white','black'),bty='n', cex=0.9, pt.cex=2)
  
  mtext('E',side = 3,cex = 1.2,adj = -0.065,padj = -0.4)
  
  
  dev.off()
  
}



##HUMERUS PLOT##

{
pdf('Humerus_plot.pdf',width=7, height =6.8, useDingbats = F, family = 'sans')

layout ( matrix ( c(1,2,
                    3,4,
                    5,5),3,byrow = T))
#layout.show(5)

par(mar=c(2,4.5,2,3.5))

#A
plot (HUM_list$general[,2:1],
      ylab='log10 Carapace length (mm)',xlab='log10 Humerus length (mm)',
      cex=1,
      pch=NA,
      col= NA,
      xaxt='n',yaxt='n')
axis(1,at=c(1.5,2,2.5))
axis(2,at=c(2,2.5,3))

# confidence interval
conf.HL <- sort(HUM_list$general$HL)
preds = apply ( HUM_CI_general[,-3], 1, function(x) as.matrix(HUM_best_models$general$X[,-3] ) %*% x  )

# fill in area between regression line and confidence interval
polygon(x=c(rev(conf.HL), conf.HL), y=c(rev(sort(preds[,2])), sort(preds[,1])), 
        col = rgb(225/255, 225/255, 225/255, 0.5), border = NA)

points(HUM_list$general[,2:1],
       cex=1,
       pch=ifelse(HUM_list$general$Type=='fossil',15,16),
       col= sapply(cols['stem'],make.transparent,alpha=0.95))

abline(HUM_best_models$general,lty=2,lwd=1)

points(HUM_list$general[HUM_list$general$Clade=='stem',2:1],
       cex=1,
       pch=ifelse(HUM_list$general[HUM_list$general$Clade=='stem','Type']=='fossil',22,21),
       bg=make.transparent('mediumpurple4' ,0.6),
       col=make.transparent('mediumpurple4',0.6))
#colour terrestrial taxa
points(HUM_list$general[HUM_list$general$Terr_spec=='1',2:1],
       cex=1,
       pch=ifelse(HUM_list$general[HUM_list$general$Terr_spec=='1','Type']=='fossil',22,21),
       col=make.transparent('black',0.7), bg=NA)

legend('topleft',legend=c('PGLS fit','stem turtles'),bty='n',
       lty=c(2,NA),lwd=1.2,cex=0.9,col=c('black','mediumpurple4'),pch=c(NA,15),
       pt.cex=2)
legend('bottomright',legend=c('extant','fossil'),
       pch=c(16,15),pt.cex=2,cex=0.9,
       col='gray95',
       bty='n',ncol=2)

mtext('A',side = 3,cex = 1.2,adj = -0.17)

#B
plot (HUM_list$general[,2:1],
      ylab='log10 Carapace length (mm)',xlab='log10 Humerus length (mm)',
      cex=1,
      pch=NA,
      col= NA,
      xaxt='n',yaxt='n')
axis(1,at=c(1.5,2,2.5))
axis(2,at=c(2,2.5,3))

# confidence interval
conf.HL <- sort(HUM_list$general$HL)
preds = apply ( HUM_CI_general[,-3], 1, function(x) as.matrix(HUM_best_models$general$X[,-3] ) %*% x  )

# fill in area between regression line and confidence interval
polygon(x=c(rev(conf.HL), conf.HL), y=c(rev(sort(preds[,2])), sort(preds[,1])), 
        col = rgb(225/255, 225/255, 225/255, 0.5), border = NA)

points(HUM_list$general[,2:1],
       cex=1,
       pch=ifelse(HUM_list$general$Type=='fossil',15,16),
       col= sapply(cols['stem'],make.transparent,alpha=0.95))
abline(HUM_best_models$general,lty=2,lwd=1)


points(HUM_list$general[HUM_list$general$Clade=='Pelomedusoides',2:1],
       cex=1,
       pch=ifelse(HUM_list$general[HUM_list$general$Clade=='Pelomedusoides','Type']=='fossil',22,21),
       bg=make.transparent( cols['Pelomedusoides'],0.6),
       col=make.transparent( cols['Pelomedusoides'],0.6))

points( sort(predict(HUM_best_models$Pelomedusoides))~
          HUM_list$Pelomedusoides[names(sort(predict(HUM_best_models$Pelomedusoides))),'HL'],type='l',
        lwd=2,
        col=cols['Pelomedusoides'])

points(HUM_list$general[HUM_list$general$Clade=='Chelidae',2:1],
       cex=1,
       pch=ifelse(HUM_list$general[HUM_list$general$Clade=='Chelidae','Type']=='fossil',22,21),
       bg=make.transparent( cols['Chelidae'],0.6),
       col=make.transparent( cols['Chelidae'],0.6))

points( sort(predict(HUM_best_models$Chelidae))~
          HUM_list$Chelidae[names(sort(predict(HUM_best_models$Chelidae))),'HL'],type='l',lwd=2,
        col=cols['Chelidae'])

legend('topleft',legend=c('Chelidae','Pelomedusoides'),
       pch=21,cex=0.9,pt.cex=2,
       pt.bg=cols[c('Chelidae','Pelomedusoides')],ncol=1,bty='n',
       col=sapply(cols[c('Chelidae','Pelomedusoides')],make.transparent,0.6))

legend('bottomright',legend=c('extant','fossil'),
       pch=c(16,15),pt.cex=2,cex=0.9,
       col='gray95',
       bty='n',ncol=2)

mtext('B',side = 3,cex = 1.2,adj = -0.17)

#C
plot (HUM_list$general[,2:1],
      ylab='log10 Carapace length (mm)',xlab='log10 Humerus length (mm)',
      cex=1,
      pch=NA,
      col= NA,
      xaxt='n',yaxt='n')
axis(1,at=c(1.5,2,2.5))
axis(2,at=c(2,2.5,3))

# confidence interval
conf.HL <- sort(HUM_list$general$HL)
preds = apply ( HUM_CI_general[,-3], 1, function(x) as.matrix(HUM_best_models$general$X[,-3] ) %*% x  )

# fill in area between regression line and confidence interval
polygon(x=c(rev(conf.HL), conf.HL), y=c(rev(sort(preds[,2])), sort(preds[,1])), 
        col = rgb(225/255, 225/255, 225/255, 0.5), border = NA)

points(HUM_list$general[,2:1],
       cex=1,
       pch=ifelse(HUM_list$general$Type=='fossil',15,16),
       col= sapply(cols['stem'],make.transparent,alpha=0.95))
abline(HUM_best_models$general,lty=2,lwd=1)

points(HUM_list$general[HUM_list$general$Clade=='Chelydroidea',2:1],
       cex=1,
       pch=ifelse(HUM_list$general[HUM_list$general$Clade=='Chelydroidea','Type']=='fossil',22,21),
       bg=make.transparent( cols['Chelydroidea'],0.6),
       col=make.transparent( cols['Chelydroidea'],0.6))

points( sort(predict(HUM_best_models$Chelydroidea))~
          HUM_list$Chelydroidea[names(sort(predict(HUM_best_models$Chelydroidea))),'HL'],type='l',
        lwd=2,
        col=cols['Chelydroidea'])

points(HUM_list$general[HUM_list$general$Clade=='Trionychia',2:1],
       cex=1,
       pch=ifelse(HUM_list$general[HUM_list$general$Clade=='Trionychia','Type']=='fossil',22,21),
       bg=make.transparent( cols['Trionychia'],0.6),
       col=make.transparent( cols['Trionychia'],0.6))

tr <- sort(fitted(HUM_best_models$Trionychia))
points( fitted(lm(tr ~ HUM_list$Trionychia[names(tr),'HL'])) ~ HUM_list$Trionychia[names(tr),'HL'],
        type='l', lwd=2, col=cols['Trionychia'])

points(HUM_list$general[HUM_list$general$Clade=='Chelonioidea',2:1],
       cex=1,
       pch=ifelse(HUM_list$general[HUM_list$general$Clade=='Chelonioidea','Type']=='fossil',22,21),
       bg=make.transparent( cols['Chelonioidea'],0.6),
       col=make.transparent( cols['Chelonioidea'],0.6))

points( sort(predict(HUM_best_models$Chelonioidea))~
          HUM_list$Chelonioidea[names(sort(predict(HUM_best_models$Chelonioidea))),'HL'],type='l',
        lwd=2,
        col=cols['Chelonioidea'])

legend('topleft',legend=c('Trionychia','Chelydroidea','Chelonioidea'),
       pch=21,cex=0.9,pt.cex=2,
       pt.bg=cols[c('Trionychia','Chelydroidea','Chelonioidea')],ncol=1,bty='n',
       col=sapply(cols[c('Trionychia','Chelydroidea','Chelonioidea')],make.transparent,0.6)  )

legend('bottomright',legend=c('extant','fossil'),
       pch=c(16,15),pt.cex=2,cex=0.9,
       col='gray95',
       bty='n',ncol=2)

mtext('C',side = 3,cex = 1.2,adj = -0.17)

#D
plot (HUM_list$general[,2:1],
      ylab='log10 Carapace length (mm)',xlab='log10 Humerus length (mm)',
      cex=1,
      pch=NA,
      col= NA,
      xaxt='n',yaxt='n')
axis(1,at=c(1.5,2,2.5))
axis(2,at=c(2,2.5,3))

# confidence interval
conf.HL <- sort(HUM_list$general$HL)
preds = apply ( HUM_CI_general[,-3], 1, function(x) as.matrix(HUM_best_models$general$X[,-3] ) %*% x  )

# fill in area between regression line and confidence interval
polygon(x=c(rev(conf.HL), conf.HL), y=c(rev(sort(preds[,2])), sort(preds[,1])), 
        col = rgb(225/255, 225/255, 225/255, 0.5), border = NA)

points(HUM_list$general[,2:1],
       cex=1,
       pch=ifelse(HUM_list$general$Type=='fossil',15,16),
       col= sapply(cols['stem'],make.transparent,alpha=0.95))
abline(HUM_best_models$general,lty=2,lwd=1)

points(HUM_list$general[HUM_list$general$Clade=='Testudinidae',2:1],
       cex=1,
       pch=ifelse(HUM_list$general[HUM_list$general$Clade=='Testudinidae','Type']=='fossil',22,21),
       bg=make.transparent( cols['Testudinidae'],0.6),
       col=make.transparent( cols['Testudinidae'],0.6))

points( sort(predict(HUM_best_models$Testudinidae))~
          HUM_list$Testudinidae[names(sort(predict(HUM_best_models$Testudinidae))),'HL'],type='l',
        lwd=2,
        col=cols['Testudinidae'])

points(HUM_list$general[HUM_list$general$Clade=='Geoemydidae',2:1],
       cex=1,
       pch=ifelse(HUM_list$general[HUM_list$general$Clade=='Geoemydidae','Type']=='fossil',22,21),
       bg=make.transparent( cols['Geoemydidae'],0.6),
       col=make.transparent( cols['Geoemydidae'],0.6))

tr <- sort(fitted(HUM_best_models$Geoemydidae))
points( fitted(lm(tr ~ HUM_list$Geoemydidae[names(tr),'HL'])) ~ HUM_list$Geoemydidae[names(tr),'HL'],
        type='l', lwd=2, col=cols['Geoemydidae'])

points(HUM_list$general[HUM_list$general$Clade=='Emysternia',2:1],
       cex=1,
       pch=ifelse(HUM_list$general[HUM_list$general$Clade=='Emysternia','Type']=='fossil',22,21),
       bg=make.transparent( cols['Emysternia'],0.6),
       col=make.transparent( cols['Emysternia'],0.6))

points( sort(predict(HUM_best_models$Emysternia))~
          HUM_list$Emysternia[names(sort(predict(HUM_best_models$Emysternia))),'HL'],type='l',
        lwd=2,
        col=cols['Emysternia'])

legend('topleft',legend=c('Emysternia','Geoemydidae','Testudinidae'),
       pch=21,cex=0.9,pt.cex=2,
       pt.bg=cols[c('Emysternia','Geoemydidae','Testudinidae')],ncol=1,bty='n',
       col=sapply(cols[c('Emysternia','Geoemydidae','Testudinidae')],make.transparent,0.6) )

legend('bottomright',legend=c('extant','fossil'),
       pch=c(16,15),pt.cex=2,cex=0.9,
       col='gray95',
       bty='n',ncol=2)

mtext('D',side = 3,cex = 1.2,adj = -0.17)

#E
par(mar=c(1.1,4.5,4.1,3.5))

plot_order <- c('stem','Chelidae','Pelomedusoides','Trionychia','Chelonioidea','Chelydroidea',
                'Emysternia','Geoemydidae','Testudinidae')

plot(as.numeric(as.factor(HUM_list$general$Clade)) , HUM_resid[rownames(HUM_list$general)] ,
     ylab='Residuals',yaxt=NULL,
     ylim=c(-0.3,0.3),xaxt='n',xlab='',pch=NA,
     xlim=c(-0.5,9.5))
abline(h=0,col='black',lty=2,lwd=0.8)

for ( i in 1:length(plot_order)){
  
  vioplot(  HUM_resid[rownames(HUM_list$general)][HUM_list$general$Clade==plot_order[i]] ,
            col='transparent' ,
            border='transparent',
            #border=sapply(cols[plot_order[i]],make.transparent,0.8),
            rectCol=sapply(cols.tmp[plot_order[i]],make.transparent,0.8),pchMed=21,
            colMed2=cols.tmp[plot_order[i]],colMed=sapply(cols.tmp[plot_order[i]],make.transparent,0.7),
            lineCol=cols.tmp[plot_order[i]],cex=1.5,
            horizontal = F,xlim=c(-0.3,0.3),add=T,at=i)  
  
}

for ( i in 2:length(plot_order)){
  
  quants <- quantile(sort(resid(HUM_best_models$general)[rownames(HUM_list$general)[HUM_list$general$Clade==plot_order[i]] ]),
                     probs = c(0.25,0.75))
  
  points( rep(i-0.35,
              length(resid(HUM_best_models$general)[rownames(HUM_list$general)[HUM_list$general$Clade==plot_order[i]] ] )),
            sort(resid(HUM_best_models$general)[rownames(HUM_list$general)[HUM_list$general$Clade==plot_order[i]] ]),
          type='l',lty=1,col='white',lwd=4 )
  
  points( rep(i-0.35,
              length(resid(HUM_best_models$general)[rownames(HUM_list$general)[HUM_list$general$Clade==plot_order[i]] ] )),
            sort(resid(HUM_best_models$general)[rownames(HUM_list$general)[HUM_list$general$Clade==plot_order[i]] ]),
          type='l',lty=1,col=cols.tmp[plot_order[i]] )
  
  segments(y0=quants[1],y1=quants[2],x0=i-0.35,x1=i-0.35,lwd=2.5, col=cols.tmp[plot_order[i]])
  
  points(i-0.35,median(sort(resid(HUM_best_models$general)[rownames(HUM_list$general)[HUM_list$general$Clade==plot_order[i]] ])),
         pch=21,cex=1.5,
         col=cols.tmp[plot_order[i]],bg='white')
  
}

legend('bottomleft', legend=c('general trend', 'clade-specific trend'),
       lty=c(1,1), pch=21, col='black', pt.bg = c('white','black'),bty='n', cex=0.9, pt.cex=2)

mtext('E',side = 3,cex = 1.2,adj = -0.065,padj = -0.4)


dev.off()

}


##FEMUR PLOT##

{
pdf('Femur_plot.pdf',width=7, height =6.8, useDingbats = F, family = 'sans')
  
layout ( matrix ( c(1,2,
                    3,4,
                    5,5),3,byrow = T))

par(mar=c(2,4.5,2,3.5))


plot (FEM_list$general[,2:1],
      ylab='log10 Carapace length (mm)',xlab='log10 Femur length (mm)',
      cex=1,
      pch=NA,
      col= NA,
      xaxt='n',yaxt='n')
axis(1,at=c(1.5,2,2.5))
axis(2,at=c(2,2.5,3))

# confidence interval
conf.FL <- sort(FEM_list$general$FL)
preds = apply ( FEM_CI_general, 1, function(x) as.matrix(FEM_best_models$general$X ) %*% x  )

# fill in area between regression line and confidence interval
polygon(x=c(rev(conf.FL), conf.FL), y=c(rev(sort(preds[,2])), sort(preds[,1])), 
        col = rgb(225/255, 225/255, 225/255, 0.5), border = NA)

points(FEM_list$general[,2:1],
       cex=1,
       pch=ifelse(FEM_list$general$Type=='fossil',15,16),
       col= sapply(cols['stem'],make.transparent,alpha=0.95))

abline(FEM_best_models$general,lty=2,lwd=1)

points(FEM_list$general[FEM_list$general$Clade=='stem',2:1],
       cex=1,
       pch=ifelse(FEM_list$general[FEM_list$general$Clade=='stem','Type']=='fossil',22,21),
       bg=make.transparent('mediumpurple4' ,0.6),
       col=make.transparent('mediumpurple4',0.6))

legend('topleft',legend=c('PGLS fit','stem turtles'),bty='n',
       lty=c(2,NA),lwd=1.2,cex=0.9,col=c('black','mediumpurple4'),pch=c(NA,15),
       pt.cex=2)
legend('bottomright',legend=c('extant','fossil'),
       pch=c(16,15),pt.cex=2,cex=0.9,
       col=c('gray95','gray95'),
       bty='n',ncol=2)

mtext('A',side = 3,cex = 1.2,adj = -0.17)

#B Pleurodira
plot (FEM_list$general[,2:1],
      ylab='log10 Carapace length (mm)',xlab='log10 Femur length (mm)',
      cex=1,
      pch=NA,
      col= NA,
      xaxt='n',yaxt='n')
axis(1,at=c(1.5,2,2.5))
axis(2,at=c(2,2.5,3))

# confidence interval
conf.FL <- sort(FEM_list$general$FL)
preds = apply ( FEM_CI_general, 1, function(x) as.matrix(FEM_best_models$general$X ) %*% x  )

# fill in area between regression line and confidence interval
polygon(x=c(rev(conf.FL), conf.FL), y=c(rev(sort(preds[,2])), sort(preds[,1])), 
        col = rgb(225/255, 225/255, 225/255, 0.5), border = NA)

points(FEM_list$general[,2:1],
       cex=1,
       pch=ifelse(FEM_list$general$Type=='fossil',15,16),
       col= sapply(cols['stem'],make.transparent,alpha=0.95))

abline(FEM_best_models$general,lty=2,lwd=1)

points(FEM_list$general[FEM_list$general$Clade=='Pelomedusoides',2:1],
       cex=1,
       pch=ifelse(FEM_list$general[FEM_list$general$Clade=='Pelomedusoides','Type']=='fossil',22,21),
       bg=make.transparent( cols['Pelomedusoides'],0.6),
       col=make.transparent( cols['Pelomedusoides'],0.6))

points( sort(predict(FEM_best_models$Pelomedusoides))~
          FEM_list$Pelomedusoides[names(sort(predict(FEM_best_models$Pelomedusoides))),'FL'],type='l',
        lwd=2,
        col=cols['Pelomedusoides'])

points(FEM_list$general[FEM_list$general$Clade=='Chelidae',2:1],
       cex=1,
       pch=ifelse(FEM_list$general[FEM_list$general$Clade=='Chelidae','Type']=='fossil',22,21),
       bg=make.transparent( cols['Chelidae'],0.6),
       col=make.transparent( cols['Chelidae'],0.6))

points( sort(predict(FEM_best_models$Chelidae))~
          FEM_list$Chelidae[names(sort(predict(FEM_best_models$Chelidae))),'FL'],type='l',lwd=2,
        col=cols['Chelidae'])

legend('topleft',legend=c('Chelidae','Pelomedusoides'),
       pch=21,cex=0.9,pt.cex=2,
       pt.bg=cols[c('Chelidae','Pelomedusoides')],ncol=1,bty='n',
       col=sapply(cols[c('Chelidae','Pelomedusoides')],make.transparent,0.6))

legend('bottomright',legend=c('extant','fossil'),
       pch=c(16,15),pt.cex=2,cex=0.9,
       col='gray95',
       bty='n',ncol=2)

mtext('B',side = 3,cex = 1.2,adj = -0.17)

#C
plot (FEM_list$general[,2:1],
      ylab='log10 Carapace length (mm)',xlab='log10 Femur length (mm)',
      cex=1,
      pch=NA,
      col= NA,
      xaxt='n',yaxt='n')
axis(1,at=c(1.5,2,2.5))
axis(2,at=c(2,2.5,3))

# confidence interval
conf.FL <- sort(FEM_list$general$FL)
preds = apply ( FEM_CI_general, 1, function(x) as.matrix(FEM_best_models$general$X ) %*% x  )

# fill in area between regression line and confidence interval
polygon(x=c(rev(conf.FL), conf.FL), y=c(rev(sort(preds[,2])), sort(preds[,1])), 
        col = rgb(225/255, 225/255, 225/255, 0.5), border = NA)

points(FEM_list$general[,2:1],
       cex=1,
       pch=ifelse(FEM_list$general$Type=='fossil',15,16),
       col= sapply(cols['stem'],make.transparent,alpha=0.95))

abline(FEM_best_models$general,lty=2,lwd=1)

points(FEM_list$general[FEM_list$general$Clade=='Chelydroidea',2:1],
       cex=1,
       pch=ifelse(FEM_list$general[FEM_list$general$Clade=='Chelydroidea','Type']=='fossil',22,21),
       bg=make.transparent( cols['Chelydroidea'],0.6),
       col=make.transparent( cols['Chelydroidea'],0.6))

points( sort(predict(FEM_best_models$Chelydroidea))~
          FEM_list$Chelydroidea[names(sort(predict(FEM_best_models$Chelydroidea))),'FL'],type='l',
        lwd=2,
        col=cols['Chelydroidea'])

points(FEM_list$general[FEM_list$general$Clade=='Trionychia',2:1],
       cex=1,
       pch=ifelse(FEM_list$general[FEM_list$general$Clade=='Trionychia','Type']=='fossil',22,21),
       bg=make.transparent( cols['Trionychia'],0.6),
       col=make.transparent( cols['Trionychia'],0.6))
tr <- sort(fitted(FEM_best_models$Trionychia))
points( fitted(lm(tr ~ FEM_list$Trionychia[names(tr),'FL'])) ~ FEM_list$Trionychia[names(tr),'FL'],
        type='l', lwd=2, col=cols['Trionychia'])

points(FEM_list$general[FEM_list$general$Clade=='Chelonioidea',2:1],
       cex=1,
       pch=ifelse(FEM_list$general[FEM_list$general$Clade=='Chelonioidea','Type']=='fossil',22,21),
       bg=make.transparent( cols['Chelonioidea'],0.6),
       col=make.transparent( cols['Chelonioidea'],0.6))

points( sort(predict(FEM_best_models$Chelonioidea))~
          FEM_list$Chelonioidea[names(sort(predict(FEM_best_models$Chelonioidea))),'FL'],type='l',
        lwd=2,
        col=cols['Chelonioidea'])

legend('topleft',legend=c('Trionychia','Chelydroidea','Chelonioidea'),
       pch=21,cex=0.9,pt.cex=2,
       pt.bg=cols[c('Trionychia','Chelydroidea','Chelonioidea')],ncol=1,bty='n',
       col=sapply(cols[c('Trionychia','Chelydroidea','Chelonioidea')],make.transparent,0.6)  )

legend('bottomright',legend=c('extant','fossil'),
       pch=c(16,15),pt.cex=2,cex=0.9,
       col='gray95',
       bty='n',ncol=2)

mtext('C',side = 3,cex = 1.2,adj = -0.17)

#D
plot (FEM_list$general[,2:1],
      ylab='log10 Carapace length (mm)',xlab='log10 Femur length (mm)',
      cex=1,
      pch=NA,
      col= NA,
      xaxt='n',yaxt='n')
axis(1,at=c(1.5,2,2.5))
axis(2,at=c(2,2.5,3))

# confidence interval
conf.FL <- sort(FEM_list$general$FL)
preds = apply ( FEM_CI_general, 1, function(x) as.matrix(FEM_best_models$general$X ) %*% x  )

# fill in area between regression line and confidence interval
polygon(x=c(rev(conf.FL), conf.FL), y=c(rev(sort(preds[,2])), sort(preds[,1])), 
        col = rgb(225/255, 225/255, 225/255, 0.5), border = NA)

points(FEM_list$general[,2:1],
       cex=1,
       pch=ifelse(FEM_list$general$Type=='fossil',15,16),
       col= sapply(cols['stem'],make.transparent,alpha=0.95))

abline(FEM_best_models$general,lty=2,lwd=1)
points(FEM_list$general[FEM_list$general$Clade=='Testudinidae',2:1],
       cex=1,
       pch=ifelse(FEM_list$general[FEM_list$general$Clade=='Testudinidae','Type']=='fossil',22,21),
       bg=make.transparent( cols['Testudinidae'],0.6),
       col=make.transparent( cols['Testudinidae'],0.6))


points(FEM_list$general[FEM_list$general$Clade=='Geoemydidae',2:1],
       cex=1,
       pch=ifelse(FEM_list$general[FEM_list$general$Clade=='Geoemydidae','Type']=='fossil',22,21),
       bg=make.transparent( cols['Geoemydidae'],0.6),
       col=make.transparent( cols['Geoemydidae'],0.6))


points(FEM_list$general[FEM_list$general$Clade=='Emysternia',2:1],
       cex=1,
       pch=ifelse(FEM_list$general[FEM_list$general$Clade=='Emysternia','Type']=='fossil',22,21),
       bg=make.transparent( cols['Emysternia'],0.6),
       col=make.transparent( cols['Emysternia'],0.6))

points( sort(predict(FEM_best_models$Testudinidae))~
          FEM_list$Testudinidae[names(sort(predict(FEM_best_models$Testudinidae))),'FL'],type='l',
        lwd=2,
        col=cols['Testudinidae'])
points( sort(predict(FEM_best_models$Geoemydidae))~
          FEM_list$Geoemydidae[names(sort(predict(FEM_best_models$Geoemydidae))),'FL'],type='l',
        lwd=2,
        col=cols['Geoemydidae'])
points( sort(predict(FEM_best_models$Emysternia))~
          FEM_list$Emysternia[names(sort(predict(FEM_best_models$Emysternia))),'FL'],type='l',
        lwd=2,
        col=cols['Emysternia'])

legend('topleft',legend=c('Emysternia','Geoemydidae','Testudinidae'),
       pch=21,cex=0.9,pt.cex=2,
       pt.bg=cols[c('Emysternia','Geoemydidae','Testudinidae')],ncol=1,bty='n',
       col=sapply(cols[c('Emysternia','Geoemydidae','Testudinidae')],make.transparent,0.6) )

legend('bottomright',legend=c('extant','fossil'),
       pch=c(16,15),pt.cex=2,cex=0.9,
       col='gray95',
       bty='n',ncol=2)

mtext('D',side = 3,cex = 1.2,adj = -0.17)

#E
par(mar=c(1.1,4.5,4.1,3.5))

plot_order <- c('stem','Chelidae','Pelomedusoides','Trionychia','Chelonioidea','Chelydroidea',
                'Emysternia','Geoemydidae','Testudinidae')


plot(as.numeric(as.factor(FEM_list$general$Clade)) ~ FEM_resid[rownames(FEM_list$general)] ,
     ylab='Residuals',yaxt=NULL,
     ylim=c(-0.25,0.45),xaxt='n',xlab='',pch=NA,
     xlim=c(-0.5,9.5))
abline(h=0,col='black',lty=2,lwd=0.8)

for ( i in 1:length(plot_order)){
  
  vioplot(  FEM_resid[rownames(FEM_list$general)][FEM_list$general$Clade==plot_order[i]] ,
            col='transparent' ,
            border='transparent',rectCol=sapply(cols.tmp[plot_order[i]],make.transparent,0.8),pchMed=21,
            colMed2=cols.tmp[plot_order[i]],colMed=sapply(cols.tmp[plot_order[i]],make.transparent,0.7),
            lineCol=cols.tmp[plot_order[i]],cex=1.5,
            horizontal = F,xlim=c(-0.3,0.3),add=T,at=i)  
  
}

for ( i in 2:length(plot_order)){
  
  quants <- quantile(sort(resid(FEM_best_models$general)[rownames(FEM_list$general)[FEM_list$general$Clade==plot_order[i]] ]),
                     probs = c(0.25,0.75))
  
  points( rep(i-0.35,
              length(resid(FEM_best_models$general)[rownames(FEM_list$general)[FEM_list$general$Clade==plot_order[i]] ] )),
            sort(resid(FEM_best_models$general)[rownames(FEM_list$general)[FEM_list$general$Clade==plot_order[i]] ]),
          type='l',lty=1,col='white',lwd=4 )
  
  points( rep(i-0.35,
              length(resid(FEM_best_models$general)[rownames(FEM_list$general)[FEM_list$general$Clade==plot_order[i]] ] )),
            sort(resid(FEM_best_models$general)[rownames(FEM_list$general)[FEM_list$general$Clade==plot_order[i]] ]),
          type='l',lty=1,col=cols.tmp[plot_order[i]] )
  
  segments(y0=quants[1],y1=quants[2],x0=i-0.35,x1=i-0.35,lwd=2.5, col=cols.tmp[plot_order[i]])
  
  points(i-0.35,median(sort(resid(FEM_best_models$general)[rownames(FEM_list$general)[FEM_list$general$Clade==plot_order[i]] ])),
         pch=21,cex=1.5,
         col=cols.tmp[plot_order[i]],bg='white')
  
}

legend('bottomleft', legend=c('general trend', 'clade-specific trend'),
       lty=c(1,1), pch=21, col='black', pt.bg = c('white','black'),bty='n', cex=0.9, pt.cex=2)

mtext('E',side = 3,cex = 1.2,adj = -0.065,padj = -0.4)

dev.off()

}




#Violin plots showing distribution of sigma2 values across different amniote groups

groups_col <- c('darkslategray2','dimgray','lawngreen','gold1','orange','deeppink2')

#PLOTS

pdf('sigma2_comparisons_violin.pdf',width=7,height = 5,useDingbats = F)

layout(matrix(c(1,2,3,
                0,0,0),2,byrow = T),
       heights = c(0.5,0.35))


#For FL vs. HL

plot('n',xlim=c(0,6),ylim=range(unlist(FLHL_sig2)), xlab='', ylab='',xaxt='n',
     cex.axis=0.8)

axis(1,at=(1:6)-0.5,labels = names(FLHL_sig2),cex.axis=0.8)

rect(xleft=1,xright=2,ybottom=-0.000012,ytop=0.00035,border = NA, 
     col = adjustcolor('grey90',alpha.f = 0.5))

vioplot(FLHL_sig2,add=T,at=(1:6)-0.5,
        ylab=expression(sigma^2),main='Femur length vs. Humerus length' ,cex.main=1,
        col=adjustcolor( groups_col,alpha.f = 0.25) , rectCol='white',lineCol=groups_col,
        border=groups_col,pchMed=21, colMed=groups_col,colMed2=adjustcolor( groups_col,alpha.f = 0.5))


#For Body size proxy vs. HL

plot('n',xlim=c(0,6),ylim=range(unlist(HL_sig2)), xlab='', ylab=expression(sigma^2),xaxt='n',
     cex.axis=0.8)

axis(1,at=(1:6)-0.5,labels = names(HL_sig2),cex.axis=0.8)

rect(xleft=1,xright=2,ybottom=-0.00015,ytop=0.00418,border = NA, 
     col = adjustcolor('grey90',alpha.f = 0.5))

vioplot(HL_sig2,add=T,at=(1:6)-0.5,
        ylab=expression(sigma^2),main='Body size proxy vs. Humerus length' ,cex.main=1,
        col=adjustcolor( groups_col,alpha.f = 0.25) , rectCol='white',lineCol=groups_col,
        border=groups_col,pchMed=21, colMed=groups_col,colMed2=adjustcolor( groups_col,alpha.f = 0.5))


#For Body size proxy vs. FL

plot('n',xlim=c(0,6),ylim=range(unlist(FL_sig2)), xlab='', ylab=expression(sigma^2),xaxt='n',
     cex.axis=0.8)

axis(1,at=(1:6)-0.5,labels = names(FL_sig2),cex.axis=0.8)

rect(xleft=1,xright=2,ybottom=-0.0001,ytop=0.00325,border = NA, 
     col = adjustcolor('grey90',alpha.f = 0.5))

vioplot(FL_sig2,add=T,at=(1:6)-0.5,
        ylab=expression(sigma^2),main='Body size proxy vs. Femur length' ,cex.main=1,
        col=adjustcolor( groups_col,alpha.f = 0.25) , rectCol='white',lineCol=groups_col,
        border=groups_col,pchMed=21, colMed=groups_col,colMed2=adjustcolor( groups_col,alpha.f = 0.5))


dev.off()


############END



### Code to reproduce figures in the Supplementary Material ###

#LIMBS dataset
numbers <- setNames(as.numeric(as.factor(rownames(LIMB_list$general))),
                    rownames(LIMB_list$general))

#Create caption for supplementary figure, indicating which number corresponds to which species

write.table( (paste(paste0(numbers,'- ', gsub('_',' ',names(numbers)) ),collapse=', ')) , file='caption_limbs_pgls.txt',
             sep='')


pdf('LIMBS_suppl_plot.pdf', height = 7, width = 6.5, useDingbats = F)

par(mfrow=c(3,3))

#stem turtles
plot (LIMB_list$general[,1:2],
      ylab='log10 Femur length (mm)',xlab='log10 Humerus length (mm)',
      cex=1,
      pch=NA,
      col= NA,
      xaxt='n',yaxt='n', main='stem turtles')
axis(1,at=c(1.5,2,2.5))
axis(2,at=c(1.5,2,2.5,3))

#Add all points
points(LIMB_list$general[,1:2],
       cex=0.8,
       pch=ifelse(LIMB_list$general$Type=='fossil',15,16),
       col= sapply(cols['stem'],make.transparent,alpha=0.95))

#Add trendline
abline(LIMBS_best_models$general,lty=2,lwd=1)

#Add points
points(LIMB_list$general[LIMB_list$general$Clade=='stem',1:2],
       cex=1.8,
       pch=ifelse(LIMB_list$general[LIMB_list$general$Clade=='stem','Type']=='fossil',22,21),
       bg=make.transparent('mediumpurple4' ,0.4),
       col=make.transparent('mediumpurple4',0.4))

#Add names
text(LIMB_list$general[LIMB_list$general$Clade=='stem',1:2],
     cex=0.45 , lab = numbers[rownames(LIMB_list$general)[LIMB_list$general$Clade=='stem'] ] )


clades <- c('Chelidae','Pelomedusoides','Trionychia','Chelydroidea','Chelonioidea',
            'Emysternia','Geoemydidae','Testudinidae')

#Add other clades
for ( i in 1:length(clades)){
  
  
  plot (LIMB_list$general[,1:2],
        ylab='log10 Femur length (mm)',xlab='log10 Humerus length (mm)',
        cex=1,
        pch=NA,
        col= NA,
        xaxt='n',yaxt='n', main=clades[i])
  axis(1,at=c(1.5,2,2.5))
  axis(2,at=c(1.5,2,2.5,3))
  
  #Add all points
  points(LIMB_list$general[,1:2],
         cex=0.8,
         pch=ifelse(LIMB_list$general$Type=='fossil',15,16),
         col= sapply(cols['stem'],make.transparent,alpha=0.95))
  
  #Add trendline
  abline(LIMBS_best_models$general,lty=2,lwd=1)
  
  #Add points
  points(LIMB_list$general[LIMB_list$general$Clade==clades[i],1:2],
         cex=1.8,
         pch=ifelse(LIMB_list$general[LIMB_list$general$Clade==clades[i],'Type']=='fossil',22,21),
         bg=make.transparent(cols[clades[i]] ,0.6),
         col=make.transparent(cols[clades[i]],0.6))
  #Add names
  text(LIMB_list$general[LIMB_list$general$Clade==clades[i],1:2],
       cex=0.45 , lab = numbers[rownames(LIMB_list$general)[LIMB_list$general$Clade==clades[i]] ] )
  
}


dev.off()

####################################################

#HUMERUS dataset

numbers <- setNames(as.numeric(as.factor(rownames(HUM_list$general))),
                    rownames(HUM_list$general))

write.table( (paste(paste0(numbers,'- ', gsub('_',' ',names(numbers)) ),collapse=', ')) , file='caption_humerus_pgls.txt',
             sep='')


pdf('HUM_suppl_plot.pdf', height = 7, width = 6.5, useDingbats = F)

par(mfrow=c(3,3))

#stem turtles
plot (HUM_list$general[,2:1],
      ylab='log10 Carapace length (mm)',xlab='log10 Humerus length (mm)',
      cex=1,
      pch=NA,
      col= NA,
      xaxt='n',yaxt='n', main='stem turtles')
axis(1,at=c(1.5,2,2.5))
axis(2,at=c(1.5,2,2.5,3))

#Add all points
points(HUM_list$general[,2:1],
       cex=0.8,
       pch=ifelse(HUM_list$general$Type=='fossil',15,16),
       col= sapply(cols['stem'],make.transparent,alpha=0.95))

#Add trendline
abline(HUM_best_models$general,lty=2,lwd=1)

#Add points
points(HUM_list$general[HUM_list$general$Clade=='stem',2:1],
       cex=1.8,
       pch=ifelse(HUM_list$general[HUM_list$general$Clade=='stem','Type']=='fossil',22,21),
       bg=make.transparent('mediumpurple4' ,0.4),
       col=make.transparent('mediumpurple4',0.4))

#Add names
text(HUM_list$general[HUM_list$general$Clade=='stem',2:1],
     cex=0.45 , lab = numbers[rownames(HUM_list$general)[HUM_list$general$Clade=='stem'] ] )

#Add other clades
for ( i in 1:length(clades)){
  
  
  plot (HUM_list$general[,2:1],
        ylab='log10 Carapace length (mm)',xlab='log10 Humerus length (mm)',
        cex=1,
        pch=NA,
        col= NA,
        xaxt='n',yaxt='n', main=clades[i])
  axis(1,at=c(1.5,2,2.5))
  axis(2,at=c(1.5,2,2.5,3))
  
  #Add all points
  points(HUM_list$general[,2:1],
         cex=0.8,
         pch=ifelse(HUM_list$general$Type=='fossil',15,16),
         col= sapply(cols['stem'],make.transparent,alpha=0.95))
  
  #Add trendline
  abline(HUM_best_models$general,lty=2,lwd=1)
  
  #Add points
  points(HUM_list$general[HUM_list$general$Clade==clades[i],2:1],
         cex=1.8,
         pch=ifelse(HUM_list$general[HUM_list$general$Clade==clades[i],'Type']=='fossil',22,21),
         bg=make.transparent(cols[clades[i]] ,0.6),
         col=make.transparent(cols[clades[i]],0.6))
  #Add names
  text(HUM_list$general[HUM_list$general$Clade==clades[i],2:1],
       cex=0.45 , lab = numbers[rownames(HUM_list$general)[HUM_list$general$Clade==clades[i]] ] )
  
  }

dev.off()

####################################################

#FEMUR dataset

numbers <- setNames(as.numeric(as.factor(rownames(FEM_list$general))),
                    rownames(FEM_list$general))

write.table( (paste(paste0(numbers,'- ', gsub('_',' ',names(numbers)) ),collapse=', ')) , file='caption_femur_pgls.txt',
             sep='')

pdf('FEM_suppl_plot.pdf', height = 7, width = 6.5, useDingbats = F)


par(mfrow=c(3,3))

#stem turtles
plot (FEM_list$general[,2:1],
      ylab='log10 Carapace length (mm)',xlab='log10 Femur length (mm)',
      cex=1,
      pch=NA,
      col= NA,
      xaxt='n',yaxt='n', main='stem turtles')
axis(1,at=c(1.5,2,2.5))
axis(2,at=c(1.5,2,2.5,3))

#Add all points
points(FEM_list$general[,2:1],
       cex=0.8,
       pch=ifelse(FEM_list$general$Type=='fossil',15,16),
       col= sapply(cols['stem'],make.transparent,alpha=0.95))

#Add trendline
abline(FEM_best_models$general,lty=2,lwd=1)

#Add points
points(FEM_list$general[FEM_list$general$Clade=='stem',2:1],
       cex=1.8,
       pch=ifelse(FEM_list$general[FEM_list$general$Clade=='stem','Type']=='fossil',22,21),
       bg=make.transparent('mediumpurple4' ,0.4),
       col=make.transparent('mediumpurple4',0.4))

#Add names
text(FEM_list$general[FEM_list$general$Clade=='stem',2:1],
     cex=0.45 , lab = numbers[rownames(FEM_list$general)[FEM_list$general$Clade=='stem'] ] )

#Add other clades
for ( i in 1:length(clades)){
  
  
  plot (FEM_list$general[,2:1],
        ylab='log10 Carapace length (mm)',xlab='log10 Femur length (mm)',
        cex=1,
        pch=NA,
        col= NA,
        xaxt='n',yaxt='n', main=clades[i])
  axis(1,at=c(1.5,2,2.5))
  axis(2,at=c(1.5,2,2.5,3))
  
  #Add all points
  points(FEM_list$general[,2:1],
         cex=0.8,
         pch=ifelse(FEM_list$general$Type=='fossil',15,16),
         col= sapply(cols['stem'],make.transparent,alpha=0.95))
  
  #Add trendline
  abline(FEM_best_models$general,lty=2,lwd=1)
  
  #Add points
  points(FEM_list$general[FEM_list$general$Clade==clades[i],2:1],
         cex=1.8,
         pch=ifelse(FEM_list$general[FEM_list$general$Clade==clades[i],'Type']=='fossil',22,21),
         bg=make.transparent(cols[clades[i]] ,0.6),
         col=make.transparent(cols[clades[i]],0.6))
  #Add names
  text(FEM_list$general[FEM_list$general$Clade==clades[i],2:1],
       cex=0.45 , lab = numbers[rownames(FEM_list$general)[FEM_list$general$Clade==clades[i]] ] )
  
}

dev.off()

#END