##### STEP 1: PACKAGES AND WORKING DIRECTORY #####
library(dplyr)
library(ggbeeswarm)
library(ggplot2)
library(ggdist)
library(ggpubr)
library(ggExtra)
library(hexbin)
library(maveLLR)
library(patchwork)
library(PRROC)
library(pROC)
library(RColorBrewer)
library(rlang)
library(scales)
library(tidyverse)
library(tidyr)

# If you are on Andrew's computer:
workingDir='~/Dropbox/Andrew-Lorena/KCNQ1 Analysis/'

# If you are on Lorena's HOME computer:
workingDir='C:/Users/marga/Dropbox/Andrew-Lorena/KCNQ1 Analysis/'

# If you are on Lorena's LAB computer:
#workingDir='c:/Users/harveml1/Dropbox/Andrew-Lorena/KCNQ1 Analysis/'

# Otherwise, set the working directory as needed.

setwd(workingDir)

# If you have already run through the pipeline through STEP 3, load in tf4
tf4=read.csv("Current scores/KCNQ1Scores_9.11.25.csv",stringsAsFactors = FALSE)

# If you have already run through the pipeline through STEP 4, load in tf5
tf5=read.csv("Current scores/KCNQ1Scores2_9.11.25.csv",stringsAsFactors = FALSE)

##### STEP 2: FUNCTIONS #####
addKCNQ1Position=function(df,usePos=FALSE){
  start=c(104,120,150,186,196,224,245,258,298,312,322,363,389,509,537,587)
  end = c(115,143,181,194,216,241,257,285,311,321,361,384,395,532,567,620)
  type= c('VSD_helix','VSD_helix','VSD_helix','VSD_helix','VSD_helix','VSD_helix','Pore_helix','Pore_helix','Pore_helix','Pore_helix','Pore_helix','Cyto_helix','Cyto_helix','Cyto_helix','Cyto_helix','Cyto_helix')
  type2= c('S0','S1','S2','S2-3','S3','S4','S4-5','S5','P','Ploop','S6','HA','HA','HB','HC','HD')
  for(i in 1:nrow(df)){
    if(usePos){
      pos=df[i,'pos']
    } else{
      pos=df[i,'Position']
    }
    region='Other'
    region2='Other'
    for(j in 1:length(start)){
      s=start[j]
      e=end[j]
      if(pos>=s & pos<=e){
        region=type[j]
        region2=type2[j]
      }
    }
    df[i,'region2']=region
    df[i,'region3']=region2
  }
  return(df)
}
asinh_trans <- trans_new(
  #This is a pseudo-log transformation used for the plot(s) below.
  name = "asinh",
  transform = function(x) asinh(x),
  inverse = function(x) sinh(x),
  breaks = pretty_breaks()
)
brnichCounting=function(tf5_control,scoreName,lowerCutoff,upperCutoff,missenseOnly=TRUE,results,invertDirection=FALSE,strictClinvar=FALSE){
  if(is.na(results[1,1])){
    currRow=1
  } else{
    currRow = nrow(results)+1
  }
  if(missenseOnly){
    tf5_control=tf5_control[tf5_control$mut_type=='missense',]
  }
  if(strictClinvar){
    tf5_control$controlFinal=tf5_control$clinvar
  }
  tf5_control=tf5_control[!is.na(tf5_control[,scoreName]),]
  tf5_normal=tf5_control[tf5_control[,scoreName]>=lowerCutoff & tf5_control[,scoreName]<=upperCutoff,]
  tf5_abnormal=tf5_control[tf5_control[,scoreName]<lowerCutoff | tf5_control[,scoreName]>upperCutoff,]
  if(invertDirection){
    temp=tf5_normal
    tf5_normal=tf5_abnormal
    tf5_abnormal=temp
  }
  results[currRow,'scoreName']=scoreName
  results[currRow,'missenseOnly']=missenseOnly
  results[currRow,'strictClinvar']=strictClinvar
  results[currRow,'all_blb']=sum(tf5_control$controlFinal=='blb')
  results[currRow,'all_plp']=sum(tf5_control$controlFinal=='plp')
  results[currRow,'all_vus']=sum(tf5_control$controlFinal=='vus')
  results[currRow,'normal_blb']=sum(tf5_normal$controlFinal=='blb')
  results[currRow,'normal_plp']=sum(tf5_normal$controlFinal=='plp')
  results[currRow,'normal_vus']=sum(tf5_normal$controlFinal=='vus')
  results[currRow,'abnormal_blb']=sum(tf5_abnormal$controlFinal=='blb')
  results[currRow,'abnormal_plp']=sum(tf5_abnormal$controlFinal=='plp')
  results[currRow,'abnormal_vus']=sum(tf5_abnormal$controlFinal=='vus')
  brnichA=sum(tf5_normal$controlFinal=='blb')
  brnichB=sum(tf5_abnormal$controlFinal=='blb')
  brnichC=sum(tf5_normal$controlFinal=='plp')
  brnichD=sum(tf5_abnormal$controlFinal=='plp')
  results[currRow,'brnichA']=brnichA
  results[currRow,'brnichB']=brnichB
  results[currRow,'brnichC']=brnichC
  results[currRow,'brnichD']=brnichD
  
  # Brnich math
  P1=(brnichC+brnichD)/(brnichA+brnichB+brnichC+brnichD)
  P2path=brnichD/(brnichB+brnichD)
  if(brnichB==0){
    P2path=brnichD/(1+brnichB+brnichD)
  }
  P2benign=brnichC/(brnichA+brnichC)
  if(brnichC==0){
    P2path=(brnichC+1)/(1+brnichA+brnichC)
  }
  OddsPath_P = (P2path*(1-P1)) / ((1-P2path)*P1)
  OddsPath_B = (P2benign*(1-P1)) / ((1-P2benign)*P1) 
  OddsPath_P = round(OddsPath_P, 2)
  OddsPath_B = round(OddsPath_B, 4)
  if(OddsPath_P>350){
    PS3_level='veryStrong'
  }
  else if(OddsPath_P>18.7){
    PS3_level='strong'
  }
  else if(OddsPath_P>4.3){
    PS3_level='moderate'
  }
  else if(OddsPath_P>2.1){
    PS3_level='supporting'
  }
  else {
    PS3_level='indeterminate'
  }
  
  if(OddsPath_B<0.0029){
    BS3_level='veryStrong'
  }
  else if(OddsPath_B<0.053){
    BS3_level='strong'
  }
  else if(OddsPath_B<0.23){
    BS3_level='moderate'
  }
  else if(OddsPath_B<0.48){
    BS3_level='supporting'
  }
  else {
    BS3_level='indeterminate'
  }
  
  results[currRow,'P1']=P1
  results[currRow,'P2path']=P2path
  results[currRow,'P2benign']=P2benign
  results[currRow,'OddsPath_P']=OddsPath_P
  results[currRow,'OddsPath_B']=OddsPath_B
  results[currRow,'PS3_level']=PS3_level
  results[currRow,'BS3_level']=BS3_level
  
  return(results)
}
calcCutoffs=function(df,varName,sds,mut_type,hardCutoff=FALSE){
  allvals=na.omit(df[,varName])
  cat(paste0('Total # of variants:',length(allvals),'\n'))
  
  # Calculate normal cutoffs (lower+upper) based on distribution of synon variants
  df_syn=df[df$mut_type=='synon',]
  vals=na.omit(df_syn[,varName])
  syn_mean=mean(df_syn[,varName],na.rm=TRUE)
  syn_sd=sd(df_syn[,varName],na.rm=TRUE)
  lower_cutoff=round(syn_mean-sds*syn_sd,3)
  upper_cutoff=round(syn_mean+sds*syn_sd,3)
  if(hardCutoff){
    lower_cutoff=0.25
    upper_cutoff=2
  }
  hist(df_syn[,varName],breaks=50,main=varName)
  abline(v=lower_cutoff,col='green')
  abline(v=upper_cutoff,col='green')

  # Calculate LOF vs partial LOF cutoff based on distribution of nonsense variants
  df_nonsense=df[df[,mut_type]=='early_nonsense',]
  fullVsPartialLOF_cutoff=0.25
  hist(df_nonsense[,varName],breaks=50,main=varName)
  abline(v=fullVsPartialLOF_cutoff,col='red')
  
  # Calculate number of missense variants in each group (simple definition)
  lof=df[df$mut_type=='missense' & df[,varName]<fullVsPartialLOF_cutoff & !is.na(df[,varName]),]
  partial_lof=df[df$mut_type=='missense' & df[,varName]<lower_cutoff & df[,varName]>fullVsPartialLOF_cutoff & !is.na(df[,varName]),]
  normal=df[df$mut_type=='missense' & df[,varName]>lower_cutoff & df[,varName]<upper_cutoff & !is.na(df[,varName]),]
  gof=df[df$mut_type=='missense' & df[,varName]>upper_cutoff & !is.na(df[,varName]),]
  cat(paste0('Normal range:',lower_cutoff,'-',upper_cutoff,'\n'))
  cat(paste('Full vs partial loss of function cutoff:',fullVsPartialLOF_cutoff,'\n'))
  cat(paste(nrow(lof),'LOF missense variants\n'))
  cat(paste(nrow(partial_lof),'partial LOF missense variants\n'))
  cat(paste(nrow(normal),'normal missense variants\n'))
  cat(paste(nrow(gof),'GOF missense variants\n'))
  
  # Plotting missense variants
  missense=df[df$mut_type=='missense',]
  hist(missense[,varName],breaks=50,xlab=varName,main=paste('Missense',varName))
  abline(v=fullVsPartialLOF_cutoff,col='red')
  abline(v=lower_cutoff,col='green')
  abline(v=upper_cutoff,col='green')
}
calculate_summary_stats_with_ci <- function(data, group_var, summary_var) {
  data %>%
    group_by(!!sym(group_var)) %>%
    summarise(
      mean = mean(!!sym(summary_var)),
      sd = sd(!!sym(summary_var)),
      n = n(),
      ci_lower = mean - qt(0.975, df = n - 1) * (sd / sqrt(n)), # Lower bound of 95% CI
      ci_upper = mean + qt(0.975, df = n - 1) * (sd / sqrt(n))  # Upper bound of 95% CI
    )
}
convert_variant1to3 <- function(variant) {
  aa_map <- c(
    A = "Ala", R = "Arg", N = "Asn", D = "Asp", C = "Cys", 
    E = "Glu", Q = "Gln", G = "Gly", H = "His", I = "Ile", 
    L = "Leu", K = "Lys", M = "Met", F = "Phe", P = "Pro", 
    S = "Ser", T = "Thr", W = "Trp", Y = "Tyr", V = "Val",
    X = "Ter"
  )
  start_aa <- substr(variant, 1, 1)
  position <- as.numeric(gsub("[A-Za-z]", "", variant))
  end_aa <- substr(variant, nchar(variant), nchar(variant))
  start_aa_3 <- aa_map[[start_aa]]
  end_aa_3 <- aa_map[[end_aa]]
  hgvs_format <- paste0("p.", start_aa_3, position, end_aa_3)
  return(hgvs_format)
}
plotChargeSensitivity=function(df,pdfName,hisNeutral=TRUE){
  df <- df %>% 
    mutate(orig_charge = case_when(
      orig %in% c("K", "R")       ~  1,   # Lys, Arg
      orig %in% c("D", "E")            ~ -1,   # Asp, Glu
      TRUE                                ~  0    # all others → neutral
    ))
  df <- df %>% 
    mutate(new_charge = case_when(
      new %in% c("K", "R")            ~  1,   # Lys, Arg
      new %in% c("D", "E")            ~ -1,   # Asp, Glu
      TRUE                                ~  0    # all others → neutral
    ))
  df$deltaCharge=df$new_charge-df$orig_charge
  df$chargeClass='neutral'
  df[df$deltaCharge<=-1,"chargeClass"]='negative'
  df[df$deltaCharge>=1,"chargeClass"]='positive'
  df_neg=df[df$deltaCharge<=-1,]
  df_pos=df[df$deltaCharge>=1,]
  df_neutral=df[df$deltaCharge==0,]
  df_neg2 <- df_neg %>%
    group_by(pos) %>%
    summarise(trafficking_score = mean(trafficking_score, na.rm = TRUE)) %>% 
    ungroup()
  df_pos2 <- df_pos %>%
    group_by(pos) %>%
    summarise(trafficking_score = mean(trafficking_score, na.rm = TRUE)) %>% 
    ungroup()
  df_neutral2 <- df_neutral %>%
    group_by(pos) %>%
    summarise(trafficking_score = mean(trafficking_score, na.rm = TRUE)) %>% 
    ungroup()
  pdf(pdfName, width = 8, height = 6)
  plot(df_neutral2$pos,df_neutral2$trafficking_score,col='darkgray',pch=20,cex=2,ylim=c(0,1.5),xlab='Residue',ylab='Mean abundance score')
  points(df_pos2$pos,df_pos2$trafficking_score,col='orange',pch=20,cex=2)
  points(df_neg2$pos,df_neg2$trafficking_score,col='purple',pch=20,cex=2)
  legend("bottomright",legend = c("Negative", "Neutral","Positive"),col=c("purple", "darkgray", "orange"),
         pch=20, pt.cex=2, bty="n")
  df=df[,c('trafficking_score','chargeClass')]
  df <- na.omit(df)
  print(
    ggplot(df, aes(x = chargeClass, y = trafficking_score, fill = chargeClass)) +
      geom_violin(scale = 'width', width = 0.6) +  # Adjust violin width
      geom_quasirandom(dodge.width = 0.6, width = 0.3, size = 1, fill = 1, color = "black") +  # Smaller points
      labs(x = 'Charge class', y = 'Trafficking score', title = 'Charge sensitivity') +
      theme_classic() +
      geom_hline(yintercept = 0.25, linetype = "dashed") +
      geom_hline(yintercept = 0.61, linetype = "dashed") +
      geom_hline(yintercept = 1.29, linetype = "dashed")
  )
  dev.off()
  kruskal_test_result <- kruskal.test(trafficking_score ~ chargeClass, data = df)
  print(kruskal_test_result)
  print(kruskal_test_result$p.value)
}
plotPredictorVsMAVE=function(predictor,pdfName,tf4d_control,tf4c){
  pdf(pdfName)
  tf4d <- tf4c[!is.na(tf4c[, predictor]), ]
  plot(tf4d$function_score, tf4d[,predictor], pch='.',xlab="Function score",ylab=predictor,ylim=c(0,1),xlim=c(-0.5,2.5),main = paste("Function Score vs",predictor))
  ggplot(tf4d, aes(x = function_score, y = get(predictor))) +
    geom_bin2d(bins = 200) +  # Adjust `bins` to control resolution
    scale_fill_viridis_c() + # Use a visually appealing color scale
    labs(
      x = "Function score", 
      y = predictor, 
      title = paste("Scatter Density Plot: Function Score vs", predictor)
    ) +
    theme_minimal()
  
  library(scales)
  stretch_top <- function(x) x^2  # Squaring stretches values near 1 more than near 0
  stretch_top_inverse <- function(x) sqrt(x)
  stretch_transform <- trans_new(
    name = "stretch_top",
    transform = stretch_top,
    inverse = stretch_top_inverse
  )
  ggplot(tf4d, aes(x = function_score, y = get(predictor))) +
    geom_bin2d(bins = 100) +
    scale_fill_viridis_c() +
    scale_y_continuous(trans = stretch_transform, name = predictor) +  # Apply transformation to Y-axis
    scale_x_continuous(name = "Function score") +  # Leave X-axis unchanged
    labs(
      title = paste("Scatter Density Plot: Stretched Top End vs", predictor)
    ) +
    theme_minimal()
  
  hex_data <- hexbin(tf4d$function_score, tf4d[, predictor], xbins = 50)
  ggplot(data = as.data.frame(hex_data@cID), aes(x = hex_data@xcm, y = hex_data@ycm)) +
    geom_hex(stat = "identity", aes(fill = hex_data@count)) +
    scale_fill_viridis_c() +
    labs(
      x = "Function score", 
      y = predictor, 
      title = paste("Hexagonal Scatter Density Plot: Function Score vs", predictor)
    ) +
    theme_minimal()
  plot(tf4d$function_score, tf4d[,predictor], 
       pch = '.', 
       col = point_colors, 
       xlab = "Function Score", 
       ylab = predictor,
       ylim = c(0,1),
       xlim = c(-0.5,3),
       main = paste("Function Score vs",predictor))
  region_levels <- levels(tf4d_control$region2)   # Get the unique levels of region2
  for (region in unique(region_colors)) {
    subset_data <- tf4d[tf4d$region2 == region, ]
    plot(subset_data$function_score, subset_data[,predictor], 
         pch = '.', 
         col = 'black',  # You can adjust color if needed
         xlab = "Function Score", 
         ylab = predictor,
         ylim=c(0,1),
         xlim=c(-0.5,3),
         main = region)  # Title indicates the region
  }
  tf4d_control=tf4d_control[tf4d_control$controlFinal_missense=='blb'|tf4d_control$controlFinal_missense=='plp',]
  colors <- ifelse(tf4d_control$controlFinal_missense == "plp", "red", "blue")
  plot(tf4d_control$function_score, tf4d_control[,predictor], 
       pch = 20, 
       col = colors, 
       xlab = "Function Score", 
       ylab = predictor,
       ylim = c(0,1),
       xlim = c(-0.5,3),
       main = paste("Function Score vs",predictor))
  
  dev.off()
}
plot_mean_and_ci <- function(data, x_var, y_var, ci_lower_var, ci_upper_var, x_label, y_label, title, y_limits = NULL) {
  ggplot(data, aes(x = !!sym(x_var), y = !!sym(y_var))) +
    geom_bar(stat = "identity", fill = "skyblue", color = "black") + # Bar chart
    geom_errorbar(aes(ymin = !!sym(ci_lower_var), ymax = !!sym(ci_upper_var)),
                  width = 0.2, color = "black") + # Error bars for 95% CI
    geom_point(aes(y = !!sym(y_var)), color = "black", size = 3) + # Add mean dots
    labs(title = title,
         x = x_label,
         y = y_label) +
    theme_minimal() +
    (if (!is.null(y_limits)) ylim(y_limits[1], y_limits[2]) else NULL) # Conditional y-axis limits
}
summarize_clusters <- function(df, cluster_col = "cluster_n6_FINAL", ...) {
  crit_list <- rlang::list2(...)
  if (is.null(names(crit_list)) || any(names(crit_list) == "")) {
    names(crit_list) <- paste0("criteria", seq_along(crit_list))
  }
  
  lvls <- as.character(1:6)
  
  rows <- lapply(names(crit_list), function(nm) {
    idx <- rlang::eval_tidy(crit_list[[nm]], df)
    
    if (!is.logical(idx) || length(idx) != nrow(df)) {
      stop(sprintf("'%s' must be a logical vector of same length as nrow(df).", nm))
    }
    
    vals <- df[[cluster_col]][idx]
    
    tab <- as.data.frame(table(factor(vals, levels = as.numeric(lvls))), stringsAsFactors = FALSE)
    counts <- setNames(as.list(tab$Freq), lvls)
    na_count <- sum(is.na(vals))
    N <- length(vals)
    
    # derive description: text after "criteria_" if present
    desc <- sub("^criteria_?", "", nm)
    
    tibble(
      criteria = nm,
      description = desc,
      !!!counts,
      `NA` = na_count,
      N = N
    )
  })
  
  bind_rows(rows) %>%
    replace_na(list(`1`=0, `2`=0, `3`=0, `4`=0, `5`=0, `6`=0, `NA`=0))
}

##### STEP 3: MERGING SCORES AND DATASETS TOGETHER #####
# Load in datasets (Hom and Het)
homtraffic=read.csv('Current scores/KCNQ1_Hom_trafficking_3.19.25.csv')
homfunctional=read.csv('Current scores/KCNQ1_Hom_function_3.13.25.csv')
hettraffic=read.csv('Current scores/KCNQ1_Het_trafficking_3.19.25.csv')
hetfunctional=read.csv('Current scores/KCNQ1_Het_function_9.11.25.csv')

# Make histograms (Hom and Het)
hist(homtraffic$trafficking_score,breaks=30)
hist(homfunctional$function_score,breaks=30)
hist(hettraffic$het_trafficking_score,breaks=30)
hist(hetfunctional$het_function_score,breaks=30)

# Merge together datasets (Hom and Het)
temp1=merge(homtraffic,homfunctional,all=TRUE)
temp2=merge(temp1,hettraffic,all=TRUE)
tf1=merge(temp2,hetfunctional,all=TRUE) #13,392

summary(as.factor(tf1$mut_type))
#early_nonsense  late_nonsense   mid_nonsense       missense          synon 
#511             64              61                 12223             544 

tf1$mut_type2=tf1$mut_type
tf1[tf1$mut_type2=="mid_nonsense","mut_type2"]="early_nonsense"
summary(as.factor(tf1$mut_type2))
#early_nonsense  late_nonsense       missense          synon 
#572             64                  12223             544 
tf1=tf1[,c(1,2,14,3:13)]

# Merge in Clinvar
clinvar=read.csv('Validation/KCNQ1_clinvar_9.18.25.csv',stringsAsFactors = FALSE)
summary(as.factor(clinvar$clinvar))
#blb conflicting     missing         plp         vus 
#335         124          58         277         586 
clinvar=clinvar[,c(1,2,3,5,6)]
tf2=merge(tf1,clinvar,all.x=TRUE)

# Merge in our trafficking single variant validation data
a=read.csv('Validation/manual_KCNQ1_trafficking_validation_vumc_7.25.25.csv',stringsAsFactors = FALSE)
a2=a[,c('Mutation','AF647_norm')]
names(a2)=c('mutation','trafficking_1Var_thisStudy')
a2$trafficking_1Var_thisStudy=a2$trafficking_1Var_thisStudy/100
a2$het_trafficking_lit=NA
a2$paper='thisStudy'
tf3=merge(tf2,a2[,c('mutation','trafficking_1Var_thisStudy')],all.x=TRUE)

# Combine literature and our data trafficking data, average by study
curation=read.csv("Validation/KCNQ1 Trafficking and Patch Clamp Literature Curation 3.19.25 forR.csv",stringsAsFactors = FALSE)
curation2=curation[,c('mutation','trafficking','het_trafficking','paper')]
curation2=curation2[!is.na(curation2$trafficking),]
names(curation2)=c('mutation','trafficking_lit','het_trafficking_lit','paper')
curation2$trafficking_lit=curation2$trafficking_lit/100
curation2$het_trafficking_lit=curation2$het_trafficking_lit/100
names(a2)[2]='trafficking_lit'
manualTraffic=rbind(curation2,a2) #188 rows
manualTraffic=manualTraffic[manualTraffic$paper=='thisStudy' | manualTraffic$paper=='29532034' | manualTraffic$paper=='33600800' | manualTraffic$paper == '39969993',] #177
summary(as.factor(manualTraffic$paper))
#29532034  33600800  39969993 thisStudy 
#51        32        61        33 
51+32+61+33 #177
# So 177 measurements from 3 published papers and our study, 138 from literature.
write.csv(manualTraffic,'Validation/KCNQ1_trafficking_literature_full 9.11.25.csv',row.names=FALSE)
manualTraffic2 <- aggregate(trafficking_lit ~ mutation, data = manualTraffic, FUN = function(x) mean(x, na.rm = TRUE))
manualTraffic3 <- aggregate(het_trafficking_lit ~ mutation, data = manualTraffic, FUN = function(x) mean(x, na.rm = TRUE))
dim(manualTraffic2) #177 total variants with a manual trafficking measurement
tf3b=merge(tf3,manualTraffic2,all.x=TRUE)
tf3c=merge(tf3b,manualTraffic3,all.x=TRUE)
sum(!is.na(tf3c$trafficking_lit) & tf3c$mut_type=='missense') #158
#So 158 missense variants with a manual trafficking measurement and a MAVE score.

# Merge in our single variant patch clamp validation data
a2=read.csv('Validation/Syncropatch_KCNQ1_4.22.25.csv',stringsAsFactors = FALSE)
tf3d=merge(tf3c,a2,all.x=TRUE)
a2$paper='thisStudy'

# Combine literature patch clamp and our patch clamp data, average by study
curation=read.csv("Validation/KCNQ1 Trafficking and Patch Clamp Literature Curation 3.19.25 forR.csv",stringsAsFactors = FALSE)
curation2 <- curation %>%
  filter(!(is.na(peakCurrentNorm) & is.na(deltaV1.2act.Mut.WT.) & is.na(hetPeakCurrentNorm)))
curation2=curation2[,c('mutation','peakCurrentNorm','deltaV1.2act.Mut.WT.','hetPeakCurrentNorm','paper')]
names(curation2)=c('mutation','peakCurrent_lit','deltaV12act_lit','het_PeakCurrent_lit','paper')
curation2$peakCurrent_lit=curation2$peakCurrent_lit/100
curation2$het_PeakCurrent_lit=curation2$het_PeakCurrent_lit/100
a2=a2[,c('mutation','Tail.CD.mean','deltaV12act_thisStudy','het_PeakCurrent_thisStudy','paper')]
names(a2)=names(curation2)
manualFunction=rbind(curation2,a2) #318 rows
write.csv(manualFunction,'Validation/KCNQ1_function_literature_full_9.11.25.csv',row.names=FALSE)
manualFunction2 <- aggregate(peakCurrent_lit ~ mutation, data = manualFunction, FUN = function(x) mean(x, na.rm = TRUE))
manualFunction3 <- aggregate(deltaV12act_lit ~ mutation, data = manualFunction, FUN = function(x) mean(x, na.rm = TRUE))
manualFunction4 <- aggregate(het_PeakCurrent_lit ~ mutation, data = manualFunction, FUN = function(x) mean(x, na.rm = TRUE))
tf3e=merge(tf3d,manualFunction2,all.x=TRUE)
tf3f=merge(tf3e,manualFunction3,all.x=TRUE)
tf3g=merge(tf3f,manualFunction4,all.x=TRUE)

# Merge in gnomAD
gnomad=read.csv('Validation/KCNQ1_gnomAD_v4.1.0_forR.csv',stringsAsFactors = FALSE)
duplicated_rows <- gnomad %>% group_by(mutation) %>% filter(n() > 1) %>% ungroup() #26 duplicated rows
filtered_gnomad <- gnomad %>% group_by(mutation) %>% arrange(desc(AlleleFrequency)) %>%
  slice(1) %>% ungroup()
duplicated_rows2 <- filtered_gnomad %>% group_by(mutation) %>% filter(n() > 1) %>% ungroup() #0!
tf4=merge(tf3g,filtered_gnomad,all.x=TRUE)
tf4$AlleleFrequency[is.na(tf4$AlleleFrequency)]=0
tf4$GroupMaxFAFfrequency[is.na(tf4$GroupMaxFAFfrequency)]=0

# Add in region columns
tf4$region=tf4$mut_type
tf4[tf4$mut_type=='missense','region']='missense_cytosolic'
tf4[tf4$pos>=104 & tf4$pos<=(361) & tf4$mut_type=='missense','region']='missense_transmembrane'
tf4=addKCNQ1Position(tf4,usePos=TRUE)
tf4[!tf4$mut_type=='missense','region2']=tf4[!tf4$mut_type=='missense','mut_type']
tf4[!tf4$mut_type=='missense','region3']=tf4[!tf4$mut_type=='missense','mut_type']
tf4[tf4$region2=='Other','region2']='other_missense'
tf4[tf4$region3=='Other','region3']='other_missense'
summary(as.factor(tf4$mut_type))
#early_nonsense  late_nonsense   mid_nonsense       missense          synon 
#511             64             61                  12223             544
summary(as.factor(tf4$region))
#early_nonsense          late_nonsense           mid_nonsense     missense_cytosolic missense_transmembrane  synon
#511                    64                      61                7524              4699                     544
summary(as.factor(tf4$region2))
#Cyto_helix early_nonsense  late_nonsense   mid_nonsense other_missense     Pore_helix          synon      VSD_helix 
#2147            511             64             61           6077           1960            544           2039 
summary(as.factor(tf4$region3))
#early_nonsense             HA             HB             HC             HD  late_nonsense   mid_nonsense other_missense 
#511            492            449            577            629             64             61           6077 
#P          Ploop             S0             S1             S2           S2-3             S3             S4 
#262            187            195            416            603            171            367            287 
#S4-5             S5             S6          synon 
#241            520            750            544 

#Making extra column: whether it is in a drug-binding position (ML277 and/or R-L3)
# We identified 6 key drug-binding residues: W248, L251, V255, Y267, F335, and F339
drugBindingPos <- c(248, 251, 255, 267, 335, 339)
drugBindingTF = tf4$mut_type=='missense' & tf4$pos %in% drugBindingPos
tf4$drugBindingTF= drugBindingTF
sum(tf4$drugBindingTF) # 112

# How many variants at drug binding positions have reduced function?
tf4_drugBindingT=tf4[drugBindingTF,]
tf4_drugBindingT$mutation
hist(tf4_drugBindingT$function_score, breaks=30)
toKeep=tf4_drugBindingT$function_score<.58
tf4_drugBindingT_low=tf4_drugBindingT[toKeep,]
dim(tf4_drugBindingT) #112
dim(tf4_drugBindingT_low) #96
# So 96 of 112 (86% of drug binding variants had scores below 0.58 (reduced function)

# Turning missense function scores at drug-binding positions to NA
tf4[tf4$drugBindingTF,'function_score']=NA
tf4[tf4$drugBindingTF,'function_score_sem']=NA
tf4[tf4$drugBindingTF,'het_function_score']=NA
tf4[tf4$drugBindingTF,'het_function_score_sem']=NA
# Abundance scores at these positions can remain because drugs aren't used for that assay

# Defining strong DN variants as meeting 3 or 4 of LOF in the 4 assays:
criterion1=tf4$trafficking_score<0.25 & !is.na(tf4$trafficking_score)
criterion2=tf4$function_score<0.25 & !is.na(tf4$function_score)
criterion3=tf4$het_trafficking_score<0.25 & !is.na(tf4$het_trafficking_score)
criterion4=tf4$het_function_score<0 & !is.na(tf4$het_function_score)
sum(criterion1) #2901
sum(criterion2) #1748
sum(criterion3) #924
sum(criterion4) #1023
tf4$criterionTotal=criterion1+criterion2+criterion3+criterion4
summary(as.factor(tf4$criterionTotal))    
#0    1    2    3    4 
#10152  1248   951   762   290 
tf4$criterionAll4=tf4$criterionTotal==4
tf4$criterion3or4=tf4$criterionTotal==3|tf4$criterionTotal==4
sum(tf4$criterion3or4) #1052
sum(tf4$criterion3or4 & tf4$mut_type=='missense') #838
sum(tf4$criterionAll4) #290
sum(tf4$criterionAll4 & tf4$mut_type=='missense') #290
#838-290=548, so 548 variants met 3/4 criteria for DN

# Add in mut_type3 for k-means clustering plot
tf4$mut_type3=tf4$mut_type
tf4[tf4$mut_type2=="missense" & tf4$criterionAll4,"mut_type3"]="severe_DN"
tf4[tf4$mut_type2=="missense" & tf4$pos>=659 & tf4$pos<=665,"mut_type3"]="PY_missense"
sum(tf4$mut_type3=='severe_DN') #290
sum(tf4$mut_type3=='PY_missense') #132
summary(as.factor(tf4$mut_type3))
#early_nonsense  late_nonsense   mid_nonsense       missense    PY_missense severe_DN  synon 
#511             64              61                  11801            132      290     544
tf4=tf4[,c(1:3,39,4:38)]

# Add in SpliceAI scores
spliceai=read.csv('Validation/kcnq1_spliceai_output3.csv',header=FALSE)
spliceai=spliceai[,c(1,2,4,5,9,10,11,12,17)]
names(spliceai)=c('chr','genomicPos','ref','alt','spliceai_accLoss','spliceai_donorLoss','spliceai_accGain','spliceai_donorGain','spliceai_merged')
head(spliceai)
spliceai <- spliceai %>%
  mutate(hgvs_g = paste0("chr", chr, ":", genomicPos, ref, ">", alt))
spliceai=spliceai[,c(10,5:9)]

hist(spliceai$spliceai_merged,breaks=1000,ylim=c(0,100))
allpossible=read.csv("Current scores/KCNQ1_Hom_function_3.13.25_allPossibleSNVs.csv")
allpossible2=merge(allpossible,spliceai,all.x=TRUE)
write.csv(allpossible2,"Current scores/KCNQ1_Hom_function_3.13.25_allPossibleSNVs2.csv",row.names = F)
dim(allpossible) #6084
dim(spliceai) #1212354
dim(allpossible2) #6084
sum(is.na(allpossible2$spliceai_merged)) #0
sum(!is.na(allpossible2$spliceai_merged)) #6084
# So we got spliceai scores for every possible KCNQ1 SNV.
hist(allpossible2$spliceai_merged,breaks=1000)
hist(allpossible2$spliceai_merged,breaks=1000,ylim=c(0,50))
# Spot check of A344A
# c.1032G>A, chr11:2583545G>A
spliceai[spliceai$hgvs_g=='chr11:2583545G>A',] #spliceAI merged 0.68
allpossible2[allpossible2$hgvs_g=='chr11:2583545G>A',]
spliceai_SNV=allpossible2[,c(2,19:23)]
spliceai_max <- spliceai_SNV %>%
  group_by(mutation) %>%
  slice_max(order_by = spliceai_merged, n = 1, with_ties = FALSE) %>% 
  ungroup()
spliceai_max=as.data.frame(spliceai_max)
tf4=merge(tf4,spliceai_max,all.x=TRUE)

# Write out merged file
write.csv(tf4,"Current scores/KCNQ1Scores_9.11.25.csv",row.names=FALSE)

# Make a smaller version for DN heatmap
tf4_dn=tf4[,c("mutation","mut_type","mut_type2","orig","pos","new","criterionTotal")]
tf4_dn[tf4_dn$criterionTotal<=2,'criterionTotal2']=0
tf4_dn[tf4_dn$criterionTotal==3,'criterionTotal2']=1
tf4_dn[tf4_dn$criterionTotal==4,'criterionTotal2']=2
tf4_dn=tf4_dn[,c("mutation","mut_type","mut_type2","orig","pos","new","criterionTotal2")]
names(tf4_dn)[7]='criterionTotal'
write.csv(tf4_dn,"Current scores/KCNQ1Scores_noDrugBindingPos_forDNheatmap_9.11.25.csv",row.names=FALSE)

#Making violin plots of MAVE score versus mut_type
pdf('Validation/ViolinPlotsforMAVEscores_9.11.25.pdf')
#Hom Trafficking Score vs Mutation Type - Violin plots
ggplot(tf4, aes(x = mut_type2, y = trafficking_score, fill = mut_type2)) +
  geom_violin(trim = TRUE, alpha = 0.7) + 
  labs(title = "Violin Plot of Hom Trafficking Score by Mutation Type", x = "Mutation Type", y = "Hom Trafficking Score") +
  scale_fill_brewer(palette = "Set3") + 
  scale_color_brewer(palette = "Set3") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10), 
        legend.position = "none") +
  geom_hline(yintercept = 0.25, linetype = "dashed") +
  geom_hline(yintercept = 0.61, linetype = "dashed") +
  geom_hline(yintercept = 1.29, linetype = "dashed")

#Hom Function Score vs Mutation Type - Violin plots
ggplot(tf4, aes(x = mut_type, y = function_score, fill = mut_type)) +
  geom_violin(trim = TRUE, alpha = 0.7) + 
  labs(
    title = "Violin Plot of Hom Function Score by Mutation Type",
    x = "Mutation Type",
    y = "Hom Function Score"
  ) +
  scale_fill_brewer(palette = "Set3") + 
  scale_color_brewer(palette = "Set3") + 
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 14), 
    axis.title = element_text(size = 12), 
    axis.text = element_text(size = 10), 
    legend.position = "none") +
  geom_hline(yintercept = 0.25, linetype = "dashed") +
  geom_hline(yintercept = 0.58, linetype = "dashed") +
  geom_hline(yintercept = 1.391, linetype = "dashed")
  
  
#Het Trafficking Score vs Mutation Type - Violin plots
ggplot(tf4, aes(x = mut_type2, y = het_trafficking_score, fill = mut_type2)) +
  geom_violin(trim = TRUE, alpha = 0.7) + 
  labs(
    title = "Violin Plot of Het Trafficking Score by Mutation Type",
    x = "Mutation Type",
    y = "Het Trafficking Score"
  ) +
  scale_fill_brewer(palette = "Set3") + 
  scale_color_brewer(palette = "Set3") + 
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 14), 
    axis.title = element_text(size = 12), 
    axis.text = element_text(size = 10), 
    legend.position = "none") +
  geom_hline(yintercept = 0.25, linetype = "dashed") +
  geom_hline(yintercept = 0.676, linetype = "dashed") +
  geom_hline(yintercept = 1.325, linetype = "dashed")

#Het Function Score vs Mutation Type - Violin plots
ggplot(tf4, aes(x = mut_type, y = het_function_score, fill = mut_type)) +
  geom_violin(trim = TRUE, alpha = 0.7) + 
  labs(
    title = "Violin Plot of Het Function Score by Mutation Type",
    x = "Mutation Type",
    y = "Het Function Score"
  ) +
  scale_fill_brewer(palette = "Set3") + 
  scale_color_brewer(palette = "Set3") + 
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 14), 
    axis.title = element_text(size = 12), 
    axis.text = element_text(size = 10), 
    legend.position = "none")+
    geom_hline(yintercept = 0.25, linetype = "dashed") +
    geom_hline(yintercept = 2.0, linetype = "dashed")
dev.off()

#Counting the number of mutation types (mut_types) in each MAVE experiment
#Hom Trafficking Mutation Type Summary
sum(tf4$mut_type2=='early_nonsense' & tf4$trafficking_score,na.rm=TRUE) #558
sum(tf4$mut_type2=='late_nonsense' & tf4$trafficking_score,na.rm=TRUE) #63
sum(tf4$mut_type2=='synon' & tf4$trafficking_score,na.rm=TRUE) #541
sum(tf4$mut_type2=='missense' & tf4$trafficking_score,na.rm=TRUE) #12089
#Hom Function Mutation Type Summary
sum(tf4$mut_type=='early_nonsense' & tf4$function_score,na.rm=TRUE) #440
sum(tf4$mut_type=='mid_nonsense' & tf4$function_score,na.rm=TRUE) #50
sum(tf4$mut_type=='late_nonsense' & tf4$function_score,na.rm=TRUE) #58
sum(tf4$mut_type=='synon' & tf4$function_score,na.rm=TRUE) #494
sum(tf4$mut_type=='missense' & tf4$function_score,na.rm=TRUE) #10970
#Het Trafficking Mutation Type Summary
sum(tf4$mut_type2=='early_nonsense' & tf4$het_trafficking_score,na.rm=TRUE) #463
sum(tf4$mut_type2=='late_nonsense' & tf4$het_trafficking_score,na.rm=TRUE) #45
sum(tf4$mut_type2=='synon' & tf4$het_trafficking_score,na.rm=TRUE) #470
sum(tf4$mut_type2=='missense' & tf4$het_trafficking_score,na.rm=TRUE) #10260
#Het Function Mutation Type Summary
sum(tf4$mut_type=='early_nonsense' & tf4$het_function_score,na.rm=TRUE) #480
sum(tf4$mut_type=='mid_nonsense' & tf4$het_function_score,na.rm=TRUE) #50
sum(tf4$mut_type=='late_nonsense' & tf4$het_function_score,na.rm=TRUE) #58
sum(tf4$mut_type=='synon' & tf4$het_function_score,na.rm=TRUE) #523
sum(tf4$mut_type=='missense' & tf4$het_function_score,na.rm=TRUE) #11483

##### STEP 4: ADDING IN MERGED SCORES AND LLR tf4-->tf5 #####
# Adding in 3 merged scores
# These are calculated in a separate python Jupyter notebook.
tf4=read.csv("Current scores/KCNQ1Scores_9.11.25.csv",stringsAsFactors = FALSE)
mergedScores=read.csv("Current scores/combination-mave-analysis_Classifiers_All_10-10-25.csv",stringsAsFactors = FALSE)
mergedScores2=mergedScores[,c('mutation','LR_prob','RF_prob','GNB_prob')]
tf5=merge(tf4,mergedScores2,all.x=TRUE)

# LLR calculation (Hom function), all variants
tf5b = tf5[!is.na(tf5$controlFinal) & (tf5$controlFinal=='blb' | tf5$controlFinal=='plp'), ]
tf5b$controlFinal <- factor(tf5b$controlFinal)
bw=0.2 #bandwidth of 0.2
posScores_func=tf5b[tf5b$controlFinal=='plp','function_score']
posScores_func=posScores_func[!is.na(posScores_func)] #192
hist(posScores_func,breaks=100)
abline(v=-0.3)
abline(v=1.5)
# Outlier filtering
posScores_func=posScores_func[posScores_func > -0.3 & posScores_func < 1.5]
negScores_func=tf5b[tf5b$controlFinal=='blb','function_score']
hist(negScores_func,breaks=100)
negScores_func=negScores_func[!is.na(negScores_func)] #29
data_func <- data.frame(
  score = c(posScores_func, negScores_func),
  type = factor(rep(c("Positive", "Negative"), c(length(posScores_func), length(negScores_func))))
)
ggplot(data_func, aes(x = score, fill = type, color = type)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  labs(title = "Histogram of Hom Func Positive and Negative Scores", x = "Scores", y = "Frequency") +
  theme_minimal()
llr_func <- buildLLR.kernel(posScores_func,negScores_func,bw=bw,kernel="gaussian")
# Manually edited buildLLR.kernel code to allow benign moderate at a symmetric level to path moderate
pdf('Validation/LLR_HomFunction_allVariants_9.11.25.pdf')
drawDensityLLR(c(posScores_func,negScores_func),llr_func$llr,llr_func$posDens,llr_func$negDens,posScores_func,negScores_func)
dev.off()
# Merging LLR scores back in with tf5
tf5$llr <- llr_func$llr(tf5$function_score)
hist(tf5$llr)
tf5 <- tf5 %>%
  mutate(llr_category = case_when(
    is.na(llr) ~ NA_character_,  # Assign NA if llr is missing
    llr > 2.54 ~ "Path_veryStrong",
    llr > 1.27 ~ "Path_strong",
    llr > 0.63 ~ "Path_moderate",
    llr > 0.31 ~ "Path_supporting",
    llr >= -0.31 & llr <= 0.31 ~ "No_evidence",  # Values near 0
    llr < -2.54 ~ "Benign_veryStrong",
    llr < -1.27 ~ "Benign_strong",
    llr < -0.63 ~ "Benign_moderate",
    llr < -0.31 ~ "Benign_supporting"
  ))
summary(as.factor(tf5$llr_category))
#Benign_moderate Benign_supporting  No_evidence   Path_moderate   Path_strong   Path_supporting   NA's 
#       6660              1477        1079             707          1757               333        1390 
# Either Strong: 1757
# Either Moderate: 707+6660 #7367
# Either Supporting: 1477+333 #1810

# LLR calculation (Hom function), missense variants
tf5b = tf5[tf5$mut_type=='missense' & !is.na(tf5$controlFinal) & (tf5$controlFinal=='blb' | tf5$controlFinal=='plp'), ]
tf5b$controlFinal <- factor(tf5b$controlFinal)
bw=0.2 #bandwidth of 0.2
posScores_func=tf5b[tf5b$controlFinal=='plp','function_score']
posScores_func=posScores_func[!is.na(posScores_func)] #192
hist(posScores_func,breaks=100)
abline(v=-0.3)
abline(v=1.5)
# Outlier filtering
posScores_func=posScores_func[posScores_func > -0.3 & posScores_func < 1.5]
negScores_func=tf5b[tf5b$controlFinal=='blb','function_score']
hist(negScores_func,breaks=100)
negScores_func=negScores_func[!is.na(negScores_func)] #29
data_func <- data.frame(
  score = c(posScores_func, negScores_func),
  type = factor(rep(c("Positive", "Negative"), c(length(posScores_func), length(negScores_func))))
)
ggplot(data_func, aes(x = score, fill = type, color = type)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  labs(title = "Histogram of Hom Func Positive and Negative Scores", x = "Scores", y = "Frequency") +
  theme_minimal()
llr_func <- buildLLR.kernel(posScores_func,negScores_func,bw=bw,kernel="gaussian")
# Manually edited buildLLR.kernel code to allow benign moderate at a symmetric level to path moderate
pdf('Validation/LLR_HomFunction_missenseVariants_9.11.25.pdf')
drawDensityLLR(c(posScores_func,negScores_func),llr_func$llr,llr_func$posDens,llr_func$negDens,posScores_func,negScores_func)
dev.off()
# Merging LLR scores back in with tf5
tf5$llr_missenseOnly <- llr_func$llr(tf5$function_score)
hist(tf5$llr_missenseOnly)

# Add K-means clusters to tf5
clusters=read.csv('./Current scores/combination-mave-analysis_PCA-Kmeans_10-10-25.csv')
for(cluster in 0:5){
  print(cluster)
  currCluster=clusters[clusters$cluster_n6==cluster,]
  #print(head(currCluster))
  print(dim(currCluster))
  print(median(currCluster$trafficking_score),na.rm=TRUE)
  print(median(currCluster$function_score))
  print(median(currCluster$het_trafficking_score))
  print(median(currCluster$het_function_score))
  print('')
}
# Manual inspection to rename and reorder clusters:
# Cluster 0: Normal2 -->Cluster 5
# Cluster 1: DN -->Cluster 1
# Cluster 2: Normal1 -->Cluster 4
# Cluster 3: PartialLOF -->Cluster 3
# Cluster 4: GOF -->Cluster 6
# Cluster 5: Haploinsufficient -->Cluster 2

clusters$cluster_n6_FINAL = clusters$cluster_n6
clusters[clusters$cluster_n6==0,'cluster_n6_FINAL']=5
clusters[clusters$cluster_n6==1,'cluster_n6_FINAL']=1
clusters[clusters$cluster_n6==2,'cluster_n6_FINAL']=4
clusters[clusters$cluster_n6==3,'cluster_n6_FINAL']=3
clusters[clusters$cluster_n6==4,'cluster_n6_FINAL']=6
clusters[clusters$cluster_n6==5,'cluster_n6_FINAL']=2
summary(as.factor(clusters$cluster_n6_FINAL))
#   1    2    3    4    5    6 
#1411  847 1440 3547 2945  126 
map <- c(
  "1" = "DN",
  "2" = "Haploinsufficient",
  "3" = "PartialLOF",
  "4" = "Normal1",
  "5" = "Normal2",
  "6" = "GOF"
)
clusters$cluster_n6_name_FINAL <- map[as.character(clusters$cluster_n6_FINAL)]
write.csv(clusters,'./Current scores/combination-mave-analysis_PCA-Kmeans_10-10-25_withFinalClusterNames.csv',row.names=FALSE)

clusters_forMerge=clusters[,c('mutation','cluster_n6_FINAL','cluster_n6_name_FINAL')]
tf5=merge(tf5,clusters_forMerge,all.x=TRUE)
summary(as.factor(tf5$cluster_n6_name_FINAL))
#DN     GOF Haploinsufficient    Normal1    Normal2    PartialLOF     NA's 
#1411   126             847       3547        2945        1440       3087
write.csv(tf5,"Current scores/KCNQ1Scores2_9.11.25.csv",row.names=FALSE)


##### STEP 5: VARIANT COUNTS #####
#Our analysis assayed 13,403 single residue variants, representing 94% of all possible variants.
dim(tf5) #13403
13403/(676*21) #94%
summary(as.factor(tf5$mut_type))
#early_nonsense  late_nonsense   mid_nonsense       missense          synon 
#511             64             61          12223            544 
511+64+61 #636 nonsense variants

#A small number (21/541 or 3.8%) of synonymous variants had MAVE scores below 0.5.
sum(!is.na(tf5$trafficking_score) & tf5$mut_type=='synon' & tf5$trafficking_score<0.5) #21
sum(!is.na(tf5$trafficking_score) & tf5$mut_type=='synon') #541

# How many variants had scores in all 4 assays?
sum(!is.na(tf5$trafficking_score) & !is.na(tf5$function_score) & !is.na(tf5$het_trafficking_score) & !is.na(tf5$het_function_score)) 
#10316

##### STEP 6: HOW MANY LOF AND GOF VARIANTS DO WE HAVE? #####
calcCutoffs(tf5,'trafficking_score',1.96,'mut_type')
#Total # of variants:13252
#Normal range:0.61-1.29
#Full vs partial loss of function cutoff: 0.25 
#2321 LOF missense variants
#1333 partial LOF missense variants
#8130 normal missense variants
#305 GOF missense variants
#2321+1333 = 3654 
# In homozygous experiments, we identified 3,654 variants with decreased surface abundance 
# and 305 variants with gain-of-abundance.

calcCutoffs(tf5,'function_score',1.96,'mut_type')
#Total # of variants:12013
#Normal range:0.58-1.391
#Full vs partial loss of function cutoff: 0.25 
#1362 LOF missense variants
#1259 partial LOF missense variants
#8010 normal missense variants
#339 GOF missense variants
#1362+1259 = 2621 LOF or partial LOF missense variants
#Additionally, we found 2,621 variants with decreased function and 339 variants gain-of-function variants.

calcCutoffs(tf5,'het_trafficking_score',1.96,'mut_type')
#Normal range:0.676-1.325
#Full vs partial loss of function cutoff: 0.25 
#921 LOF missense variants
#1856 partial LOF missense variants
#6964 normal missense variants
#529 GOF missense variants

calcCutoffs(tf5,'het_function_score',1.96,'mut_type')
#Normal range:0.363-1.717
#Full vs partial loss of function cutoff: 0.25 
#1418 LOF missense variants
#354 partial LOF missense variants
#9169 normal missense variants
#542 GOF missense variants

##### STEP 7: TABLE FOR MANUAL TRAFFICKING DATA #####
a <- read.csv('Validation/manual_KCNQ1_trafficking_validation_vumc_7.25.25.csv', stringsAsFactors = FALSE)
subset_a <- a[, c("Mutation", "DMS_Trafficking", "SD_trafficking", "Reason", "mutType")]
write.csv(subset_a, "Validation/KCNQ1_SingleVariantSurfaceAbundanceMeasurements.csv", row.names = FALSE)

##### STEP 8: CIRCLE PLOT OF CLINVAR COUNTS #####
clinvar=read.csv('Validation/KCNQ1_clinvar_9.18.25.csv',stringsAsFactors = FALSE)

#all variants (missense, nonsense, synonymous)
summary(as.factor(clinvar$clinvar))
#blb conflicting     missing         plp         vus 
#335         124          58         277         586 

#missense variants only
summary(as.factor(clinvar$clinvar[clinvar$class == "missense"]))
#blb conflicting     missing         plp         vus 
#13         105          58         221         576 

data <- data.frame(
  category = c("blb","plp","vus", "conflicting", "missing"),
  value = c(13,221,576,105,58)
)
data$fraction <- data$value / sum(data$value)
data$ymax <- cumsum(data$fraction)
data$ymin <- c(0, head(data$ymax, n = -1))
data[,c('category','fraction')]
#category   fraction
#1         blb 0.01336074
#2         plp 0.22713258
#3         vus 0.59198356
#4 conflicting 0.10791367
#5     missing 0.05960946
# 59% of missense variants in Clinvar are VUS, 11% are conflicting
# Total number of missense variants in Clinvar: 13+221+576+105+58 = 973. 973/12,844 = 7.5%

pdf('Validation/clinvar_donut_raw_10.7.25.pdf')
ggplot(data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2, fill = category)) +
  geom_rect() +
  coord_polar(theta = "y") +
  xlim(c(0, 4)) +
  theme_void() +
  theme(legend.position = "right") +
  labs(fill = "Category", title = "Donut Plot")
dev.off()

##### STEP 9: NUMBER OF CLINVAR VARIANTS WITH AVAILABLE MAVE SCORES #####
summary(as.factor(tf5$clinvar)) # all variants, true clinvar
#blb conflicting     missing         plp         vus        NA's 
#299         116          56         272         574       12086 
summary(as.factor(tf5$clinvar[tf5$mut_type == "missense"])) #missense variants, true clinvar
#blb conflicting     missing         plp         vus        NA's 
#12          99          56         218         566       11272 
summary(as.factor(tf5$controlFinal)) # all variants, clinvar+
#blb conflicting     missing         plp         vus        NA's 
#320         102          56         272         567       12086 
summary(as.factor(tf5$controlFinal_missense)) # missense variants, clinvar+
#blb conflicting     missing         plp         vus        NA's 
#28          90          56         218         559       12452 

##### STEP 10: VIOLIN PLOTS OF MAVE SCORES VS CLINVAR #####
tf5_control=tf5[tf5$controlFinal_missense=='blb'|tf5$controlFinal_missense=='plp'|tf5$controlFinal_missense=='vus'|tf5$controlFinal_missense=='conflicting',]
tf5_control=tf5_control[!is.na(tf5_control$controlFinal_missense),]
tf5_control[tf5_control$controlFinal_missense=='conflicting','controlFinal_missense']='vus'
tf5_control[tf5_control$controlFinal_missense=='vus','controlFinal_missense']='vus_confl'
write.csv(tf5_control,'Validation/KCNQ1_scores_forLLR_9.11.25.csv',row.names=FALSE)
summary(as.factor(tf5_control$controlFinal_missense)) #28 BLB, 218 PLP, 649 VUS_confl

# Counts for each MAVE:
tf5_control_traf=tf5_control[!is.na(tf5_control$trafficking_score),]
summary(as.factor(tf5_control_traf$controlFinal_missense))
#blb       plp vus_confl 
#28       217       645 
tf5_control_func=tf5_control[!is.na(tf5_control$function_score),]
summary(as.factor(tf5_control_func$controlFinal_missense))
#blb       plp vus_confl 
#28       204       602 
tf5_control_hettraf=tf5_control[!is.na(tf5_control$het_trafficking_score),]
summary(as.factor(tf5_control_hettraf$controlFinal_missense))
#blb       plp vus_confl 
#26       205       565 
tf5_control_hetfunc=tf5_control[!is.na(tf5_control$het_function_score),]
summary(as.factor(tf5_control_hetfunc$controlFinal_missense))
#blb       plp vus_confl 
#28       204       625 

pdf('Validation/Clinvar_violin_9.11.25.pdf')
ggplot(tf5_control, aes(x = controlFinal_missense, y = trafficking_score, fill = controlFinal_missense)) +
  geom_violin(scale='width') +
  geom_quasirandom(fill = 1, color = "black") +
  labs(x = 'Clinvar', y = 'Trafficking score', title = 'Clinvar comparison (trafficking scores)') +
  theme_classic()+
  geom_hline(yintercept = 0.25, linetype = "dashed") +
  geom_hline(yintercept = 0.61, linetype = "dashed") +
  geom_hline(yintercept = 1.29, linetype = "dashed")
ggplot(tf5_control, aes(x = controlFinal_missense, y = function_score, fill = controlFinal_missense)) +
  geom_violin(scale='width') +
  geom_quasirandom(fill = 1, color = "black") +
  labs(x = 'Clinvar', y = 'Function score', title = 'Clinvar comparison (function scores)') +
  theme_classic()+
  geom_hline(yintercept = 0.25, linetype = "dashed") +
  geom_hline(yintercept = 0.58, linetype = "dashed") +
  geom_hline(yintercept = 1.391, linetype = "dashed")
ggplot(tf5_control, aes(x = controlFinal_missense, y = het_trafficking_score, fill = controlFinal_missense)) +
  geom_violin(scale='width') +
  geom_quasirandom(fill = 1, color = "black") +
  labs(x = 'Clinvar', y = 'het_trafficking score', title = 'Clinvar comparison (het_trafficking scores)') +
  theme_classic()+
  geom_hline(yintercept = 0.25, linetype = "dashed") +
  geom_hline(yintercept = 0.676, linetype = "dashed") +
  geom_hline(yintercept = 1.325, linetype = "dashed")
ggplot(tf5_control, aes(x = controlFinal_missense, y = het_function_score, fill = controlFinal_missense)) +
  geom_violin(scale='width') +
  geom_quasirandom(fill = 1, color = "black") +
  labs(x = 'Clinvar', y = 'het_function score', title = 'Clinvar comparison (het_function scores)') +
  theme_classic()+
  geom_hline(yintercept = 0.25, linetype = "dashed") +
  geom_hline(yintercept = 2.0, linetype = "dashed") 
  
# Plot Clinvar violin plots with all classes of variants including missing, conflicting, NA (all four assays)
ggplot(tf5, aes(x = controlFinal_missense, y = trafficking_score, fill = controlFinal_missense)) +
  geom_violin(scale='width') +
  geom_quasirandom(fill = 1, color = "black") +
  labs(x = 'Clinvar', y = 'Trafficking score', title = 'Clinvar comparison (trafficking scores)') +
  theme_classic()
ggplot(tf5, aes(x = controlFinal_missense, y = function_score, fill = controlFinal_missense)) +
  geom_violin(scale='width') +
  geom_quasirandom(fill = 1, color = "black") +
  labs(x = 'Clinvar', y = 'Function score', title = 'Clinvar comparison (function scores)') +
  theme_classic()
ggplot(tf5, aes(x = controlFinal_missense, y = het_trafficking_score, fill = controlFinal_missense)) +
  geom_violin(scale='width') +
  geom_quasirandom(fill = 1, color = "black") +
  labs(x = 'Clinvar', y = 'het_trafficking score', title = 'Clinvar comparison (het_trafficking scores)') +
  theme_classic()
ggplot(tf5, aes(x = controlFinal_missense, y = het_function_score, fill = controlFinal_missense)) +
  geom_violin(scale='width') +
  geom_quasirandom(fill = 1, color = "black") +
  labs(x = 'Clinvar', y = 'het_function score', title = 'Clinvar comparison (het_function scores)') +
  theme_classic()
dev.off()

pdf('Validation/Clinvar_violin_4truthSets_9.11.25.pdf')
# Missense Clinvar
tf5_control=tf5[!is.na(tf5$clinvar) & !tf5$clinvar=='missing' & tf5$mut_type=='missense',]
tf5_control[tf5_control$clinvar=='conflicting' | tf5_control$clinvar=='vus','clinvar']='vus_confl'
ggplot(tf5_control, aes(x = clinvar, y = function_score, fill = clinvar)) +
  geom_violin(scale='width') +
  geom_quasirandom(fill = 1, color = "black") +
  labs(x = 'Clinvar', y = 'Function score', title = 'Missense, Clinvar function scores') +
  theme_classic()+
  geom_hline(yintercept = 0.25, linetype = "dashed") +
  geom_hline(yintercept = 0.58, linetype = "dashed") +
  geom_hline(yintercept = 1.391, linetype = "dashed")
# Missense Clinvar+
tf5_control=tf5[tf5$mut_type=='missense' & !is.na(tf5$controlFinal) & !tf5$controlFinal=='missing',]
tf5_control[tf5_control$controlFinal=='conflicting' | tf5_control$controlFinal=='vus','controlFinal']='vus_confl'
ggplot(tf5_control, aes(x = controlFinal, y = function_score, fill = controlFinal)) +
  geom_violin(scale='width') +
  geom_quasirandom(fill = 1, color = "black") +
  labs(x = 'Clinvar+', y = 'Function score', title = 'Missense Clinvar+ function scores') +
  theme_classic()+
  geom_hline(yintercept = 0.25, linetype = "dashed") +
  geom_hline(yintercept = 0.58, linetype = "dashed") +
  geom_hline(yintercept = 1.391, linetype = "dashed")
# All Clinvar
tf5_control=tf5[!is.na(tf5$clinvar) & !tf5$clinvar=='missing',]
tf5_control[tf5_control$clinvar=='conflicting' | tf5_control$clinvar=='vus','clinvar']='vus_confl'
ggplot(tf5_control, aes(x = clinvar, y = function_score, fill = clinvar)) +
  geom_violin(scale='width') +
  geom_quasirandom(fill = 1, color = "black") +
  labs(x = 'Clinvar', y = 'Function score', title = 'All, Clinvar function scores') +
  theme_classic()+
  geom_hline(yintercept = 0.25, linetype = "dashed") +
  geom_hline(yintercept = 0.58, linetype = "dashed") +
  geom_hline(yintercept = 1.391, linetype = "dashed")
# All Clinvar+
tf5_control=tf5[!is.na(tf5$controlFinal) & !tf5$controlFinal=='missing',]
tf5_control[tf5_control$controlFinal=='conflicting' | tf5_control$controlFinal=='vus','controlFinal']='vus_confl'
ggplot(tf5_control, aes(x = controlFinal, y = function_score, fill = controlFinal)) +
  geom_violin(scale='width') +
  geom_quasirandom(fill = 1, color = "black") +
  labs(x = 'Clinvar+', y = 'Function score', title = 'All, Clinvar+ function scores') +
  theme_classic()+
  geom_hline(yintercept = 0.25, linetype = "dashed") +
  geom_hline(yintercept = 0.58, linetype = "dashed") +
  geom_hline(yintercept = 1.391, linetype = "dashed")
dev.off()

##### STEP 11: STRONG DN VS CLINVAR #####
# Scatterplot with Faceting | Jitter Points | BLB/PLP/VUS | Het Trafficking & Function Scores | Strong-DN & Severe-DN
pdf('Validation/StrongDN_ClinVar_vs_HetMAVEScores_JitterPlot_9.11.25.pdf')
tf5_filtered <- tf5_control %>%
  filter(!is.na(controlFinal_missense)) %>%
  mutate(criterion_category = case_when(
    criterionTotal %in% c(0, 1, 2) ~ "0-2",
    criterionTotal == 3 ~ "3",
    criterionTotal == 4 ~ "4",
    TRUE ~ "Other"
  ))
plot1 <- ggplot(tf5_filtered, aes(x = criterion_category, y = het_trafficking_score, color = criterion_category)) +
  geom_jitter(size = 2.0, alpha = 1, width = 0.2, stroke = 0) +
  scale_color_manual(values = c("0-2" = "grey", "3" = "#5D3A9B", "4" = "#E6399B")) +
  facet_wrap(~ controlFinal_missense) +
  labs(x = 'Criterion Total Category', y = 'Het Trafficking score', title = 'Strong-DN in ClinVar vs Het Trafficking Scores') +
  theme_classic()
plot2 <- ggplot(tf5_filtered, aes(x = criterion_category, y = het_function_score, color = criterion_category)) +
  geom_jitter(size = 2.0, alpha = 1, width = 0.2, stroke = 0) +
  scale_color_manual(values = c("0-2" = "grey", "3" = "#4F4FA8", "4" = "#FFA500")) +
  facet_wrap(~ controlFinal_missense) +
  labs(x = 'Criterion Total Category', y = 'Het Function score', title = 'Strong-DN in ClinVar vs Het Function Scores') +
  ylim(NA, 3.1) + 
  theme_classic()
combined_plot <- plot1 / plot2
print(combined_plot)
dev.off()

variant_counts <- tf5_filtered %>%
  group_by(controlFinal_missense, criterion_category) %>%
  summarise(count = n(), .groups = 'drop')
print(variant_counts)
variant_counts_trafficking <- tf5_filtered %>%
  group_by(controlFinal_missense, criterion_category) %>%
  summarise(count = sum(!is.na(het_trafficking_score)), .groups = 'drop')
variant_counts_function <- tf5_filtered %>%
  group_by(controlFinal_missense, criterion_category) %>%
  summarise(count = sum(!is.na(het_function_score)), .groups = 'drop')
print(variant_counts_trafficking)
#Counting how many variants are in each MAVE, each classification, and each criterion category.
#controlFinal_missense criterion_category count
#1 blb          0-2                   26
#2 plp          0-2                  124
#3 plp          3                     40
#4 plp          4                     26
#5 vus_confl    0-2                  466
#6 vus_confl    3                     24
#7 vus_confl    4                     11

print(variant_counts_function)
#controlFinal_missense criterion_category count
#1 blb          0-2                   29
#2 plp          0-2                  125
#3 plp          3                     42
#4 plp          4                     26
#5 vus_confl    0-2                  521
#6 vus_confl    3                     24
#7 vus_confl    4                     11
#42+26 = 68, so 68 plp variants had strong/severe DN effects
#24+11 = 35, so 35 vus_confl variants had strong/severe DN effects

# Selectivity filter DN scores GYGD 314-317
selfilter=tf5[tf5$mut_type=='missense' & tf5$pos>=314 & tf5$pos<=317,]
selfilter_dn=selfilter[selfilter$criterion3or4,]
nrow(selfilter_dn) #44
# 44/4=11, so there is a mean of 11 strong/severe DN missense variants across the 4 SF residues
sum(selfilter_dn$controlFinal_missense=='plp',na.rm=TRUE) #10

##### STEP 12: VIOLIN PLOTS JLNS+NEONATAL #####
jlns=read.csv('Validation/jlns_variants_forR.csv',stringsAsFactors = F)
for(var in jlns$mutation){
  if(any(tf5$mutation==var)){
    jlns[jlns$mutation==var,'function_score']=tf5[tf5$mutation==var,'function_score']
    jlns[jlns$mutation==var,'criterion3or4']=tf5[tf5$mutation==var,'criterion3or4']
    jlns[jlns$mutation==var,'clinvar_actual']=tf5[tf5$mutation==var,'clinvar_actual']
    jlns[jlns$mutation==var,'spliceai_merged']=tf5[tf5$mutation==var,'spliceai_merged']
  }
}
jlns2=jlns[jlns$spliceai_merged<0.15,] # removing variants with moderate/high spliceAI
# D202N in particular had a spliceAI of 0.18 and function score of 1.37
sum(!is.na(jlns2$function_score),na.rm=TRUE) #41 JLNS variants with a function score 
sum(jlns2$criterion3or4,na.rm=TRUE) #17 of these had strong/severe DN effects (17/41 41%)
sum(!is.na(jlns2$function_score) & jlns2$function_score<0.58) #34/41 (83%) of JLNS variants had LOF DMS scores
jlns2$jlns=TRUE
pdf('Validation/jlns_violin2.pdf')
ggplot(jlns2, aes(x = jlns, y = function_score, fill = jlns)) +
  geom_violin(scale='width') +
  geom_quasirandom(fill = 1, color = "black") +
  labs(x = 'JLNS', y = 'Function score', title = 'JLNS') +
  theme_classic()+
  geom_hline(yintercept = 0.25, linetype = "dashed") +
  geom_hline(yintercept = 0.58, linetype = "dashed") +
  geom_hline(yintercept = 1.391, linetype = "dashed")
dev.off()

##### STEP 13: BRNICH MATH -COUNTING THE NUMBER OF Q1 VARIANTS IN EACH CLINVAR CLASSIFICATION - B/LB, CONFLICTING, MISSING, P/LP, VUS, AND NA.#####
tf5_control=tf5[!is.na(tf5$controlFinal) & (tf5$controlFinal=='plp'|tf5$controlFinal=='blb'|tf5$controlFinal=='vus'),]
results=data.frame(scoreName=NA,missenseOnly=NA,strictClinvar=NA,all_blb=NA,all_plp=NA,all_vus=NA,normal_blb=NA,normal_plp=NA,normal_vus=NA,abnormal_plp=NA,abnormal_blb=NA,abnormal_vus=NA,brnichA=NA,brnichB=NA,brnichC=NA,brnichD=NA)
# Missense only; controlFinal
results=brnichCounting(tf5_control,'trafficking_score',0.61,1.29,missenseOnly=TRUE,results)
results=brnichCounting(tf5_control,'function_score',0.58,1.391,missenseOnly=TRUE,results)
results=brnichCounting(tf5_control,'het_trafficking_score',0.676,1.325,missenseOnly=TRUE,results)
results=brnichCounting(tf5_control,'het_function_score',0.363,1.717,missenseOnly=TRUE,results)
results=brnichCounting(tf5_control,'GNB_prob',0.5,999,missenseOnly=TRUE,invertDirection=TRUE,results)
results=brnichCounting(tf5_control,'LR_prob',0.75,999,missenseOnly=TRUE,invertDirection=TRUE,results)
results=brnichCounting(tf5_control,'RF_prob',0.8,999,missenseOnly=TRUE,invertDirection=TRUE,results)
# All variants including synon and nonsense; controlFinal
results=brnichCounting(tf5_control,'trafficking_score',0.61,1.29,missenseOnly=FALSE,results)
results=brnichCounting(tf5_control,'function_score',0.58,1.391,missenseOnly=FALSE,results)
results=brnichCounting(tf5_control,'het_trafficking_score',0.676,1.325,missenseOnly=FALSE,results)
results=brnichCounting(tf5_control,'het_function_score',0.363,1.717,missenseOnly=FALSE,results)
results=brnichCounting(tf5_control,'GNB_prob',0.5,999,missenseOnly=FALSE,invertDirection=TRUE,results)
results=brnichCounting(tf5_control,'LR_prob',0.75,999,missenseOnly=FALSE,invertDirection=TRUE,results)
results=brnichCounting(tf5_control,'RF_prob',0.8,999,missenseOnly=FALSE,invertDirection=TRUE,results)
# Missense only; strictClinvar
results=brnichCounting(tf5_control,'trafficking_score',0.61,1.29,missenseOnly=TRUE,results,strictClinvar=TRUE)
results=brnichCounting(tf5_control,'function_score',0.58,1.391,missenseOnly=TRUE,results,strictClinvar=TRUE)
results=brnichCounting(tf5_control,'het_trafficking_score',0.676,1.325,missenseOnly=TRUE,results,strictClinvar=TRUE)
results=brnichCounting(tf5_control,'het_function_score',0.363,1.717,missenseOnly=TRUE,results,strictClinvar=TRUE)
results=brnichCounting(tf5_control,'GNB_prob',0.5,999,missenseOnly=TRUE,invertDirection=TRUE,results,strictClinvar=TRUE)
results=brnichCounting(tf5_control,'LR_prob',0.75,999,missenseOnly=TRUE,invertDirection=TRUE,results,strictClinvar=TRUE)
results=brnichCounting(tf5_control,'RF_prob',0.8,999,missenseOnly=TRUE,invertDirection=TRUE,results,strictClinvar=TRUE)
# All variants including synon and nonsense; strictClinvar
results=brnichCounting(tf5_control,'trafficking_score',0.61,1.29,missenseOnly=FALSE,results,strictClinvar=TRUE)
results=brnichCounting(tf5_control,'function_score',0.58,1.391,missenseOnly=FALSE,results,strictClinvar=TRUE)
results=brnichCounting(tf5_control,'het_trafficking_score',0.676,1.325,missenseOnly=FALSE,results,strictClinvar=TRUE)
results=brnichCounting(tf5_control,'het_function_score',0.363,1.717,missenseOnly=FALSE,results,strictClinvar=TRUE)
results=brnichCounting(tf5_control,'GNB_prob',0.5,999,missenseOnly=FALSE,invertDirection=TRUE,results,strictClinvar=TRUE)
results=brnichCounting(tf5_control,'LR_prob',0.75,999,missenseOnly=FALSE,invertDirection=TRUE,results,strictClinvar=TRUE)
results=brnichCounting(tf5_control,'RF_prob',0.8,999,missenseOnly=FALSE,invertDirection=TRUE,results,strictClinvar=TRUE)
results
write.csv(results,'Validation/Brnich_summary_forExcel.csv',row.names=FALSE)

# Brnich Summary #
results_preferred=results[results$missenseOnly & !results$strictClinvar,]
results_preferred[,c(1,13:23)]
#scoreName brnichA brnichB brnichC brnichD        P1    P2path  P2benign OddsPath_P OddsPath_B PS3_level     BS3_level
#1     trafficking_score      26       2      59     158 0.8857143 0.9875000 0.6941176      10.19     0.2928  moderate    supporting
#2        function_score      28       0      36     168 0.8793103 0.9940828 0.5625000      23.06     0.1765    strong      moderate
#3 het_trafficking_score      25       1      76     129 0.8874459 0.9923077 0.7524752      16.36     0.3856  moderate    supporting
#4    het_function_score      27       1      97     107 0.8793103 0.9907407 0.7822581      14.69     0.4931  moderate indeterminate
#5              GNB_prob      22       3      15     151 0.8691099 0.9805195 0.4054054       7.58     0.1027  moderate      moderate
#6               LR_prob      20       5      17     149 0.8691099 0.9675325 0.4594595       4.49     0.1280  moderate      moderate
#7               RF_prob      20       5      19     147 0.8691099 0.9671053 0.4871795       4.43     0.1431  moderate      moderate

# MAVE Trafficking: Comparing scores with clinical classifications showed that 26/28 (93%) of benign/likely benign variants had normal abundance, while 158/217 (73%) of pathogenic/likely pathogenic variants showed reduced abundance (Figure 2D). 
# MAVE Function: Function scores were more closely correlated to clinical classifications than surface abundance scores: 28/28 (100%) B/LB variants had normal function scores, while 168/204 (82%) P/LP variants showed partial loss-of-function or loss-of-function (Figure 3I).

# MERGED Random Forest
ggplot(tf5_control, aes(x = RF_prob, fill = controlFinal_missense)) +
  geom_histogram(bins = 100, position = "identity", alpha = 0.6) +
  scale_fill_manual(values = c("blb" = "blue", "plp" = "red")) +
  labs(title = "Histogram of RF_prob", x = "LR_prob", y = "Count") +
  theme_minimal()
sum(tf5_control$controlFinal_missense=='blb' & tf5_control$RF_prob<=0.8,na.rm=TRUE) #20
sum(tf5_control$controlFinal_missense=='blb' & tf5_control$RF_prob>0.8,na.rm=TRUE) #5
sum(tf5_control$controlFinal_missense=='plp' & tf5_control$RF_prob<=0.8,na.rm=TRUE) #19
sum(tf5_control$controlFinal_missense=='plp' & tf5_control$RF_prob>0.8,na.rm=TRUE) #147
# so Random forest is capturing 148/(20+148)=88% of path clinvar variants.

##### STEP 14: BARPLOT LOF/GOF PROPORTIONS BY POSITION FOR HEATMAPS #####
# Barplot - Proportion of LOF/GOF by Position in Hom Trafficking Map.
scores=tf5[tf5$mut_type=='missense',c("pos",'trafficking_score')]
scores$lowScore=scores$trafficking_score<0.61
scores$highScore=scores$trafficking_score>1.29
lowScoreResult <- scores %>%
  group_by(pos) %>%
  summarize(true_count = sum(lowScore, na.rm = TRUE))
lowScoreResult$true_count=-1*lowScoreResult$true_count
highScoreResult <- scores %>%
  group_by(pos) %>%
  summarize(true_count = sum(highScore, na.rm = TRUE))
allResult=cbind(lowScoreResult,highScoreResult[,c(2)])
names(allResult)=c('pos','lowScore','highScore')
allResult_long <- reshape2::melt(allResult, id.vars = "pos", variable.name = "ScoreType", value.name = "Score")
pdf('Validation/Proportion_LOFandGOF_Barplot_9.11.25.pdf')
ggplot(allResult_long, aes(x = pos, y = Score, fill = ScoreType)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.6) +
  scale_fill_manual(values = c("lowScore" = "red", "highScore" = "blue")) +
  theme_minimal() +
  labs(title = "Bar Graph of Hom Trafficking Scores", x = "Position", y = "Score")

# Barplot - Proportion of LOF/GOF by Position in Hom Function Map.
scores=tf5[tf5$mut_type=='missense',c("pos",'function_score')]
scores$lowScore=scores$function_score<0.58
scores$highScore=scores$function_score>1.391
lowScoreResult = scores %>%
  group_by(pos) %>%
  summarize(true_count = sum(lowScore, na.rm = TRUE))
lowScoreResult$true_count=-1*lowScoreResult$true_count
highScoreResult <- scores %>%
  group_by(pos) %>%
  summarize(true_count = sum(highScore, na.rm = TRUE))
allResult=cbind(lowScoreResult,highScoreResult[,c(2)])
names(allResult)=c('pos','lowScore','highScore')
allResult_long <- reshape2::melt(allResult, id.vars = "pos", variable.name = "ScoreType", value.name = "Score")
ggplot(allResult_long, aes(x = pos, y = Score, fill = ScoreType)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.6) +
  scale_fill_manual(values = c("lowScore" = "red", "highScore" = "blue")) +
  theme_minimal() +
  labs(title = "Bar Graph of Hom Function Scores", x = "Position", y = "Score")

# Barplot - Proportion of LOF/GOF by Position in Het Trafficking Map.
scores=tf5[tf5$mut_type=='missense',c("pos",'het_trafficking_score')]
scores$lowScore=scores$het_trafficking_score<0.676
scores$highScore=scores$het_trafficking_score>1.325
lowScoreResult = scores %>%
  group_by(pos) %>%
  summarize(true_count = sum(lowScore, na.rm = TRUE))
lowScoreResult$true_count=-1*lowScoreResult$true_count
highScoreResult <- scores %>%
  group_by(pos) %>%
  summarize(true_count = sum(highScore, na.rm = TRUE))
allResult=cbind(lowScoreResult,highScoreResult[,c(2)])
names(allResult)=c('pos','lowScore','highScore')
allResult_long <- reshape2::melt(allResult, id.vars = "pos", variable.name = "ScoreType", value.name = "Score")
ggplot(allResult_long, aes(x = pos, y = Score, fill = ScoreType)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.6) +
  scale_fill_manual(values = c("lowScore" = "red", "highScore" = "blue")) +
  theme_minimal() +
  labs(title = "Bar Graph of Het Trafficking Scores Scores", x = "Position", y = "Score")

# Barplot - Proportion of LOF/GOF by Position in Het Function Map.
scores=tf5[tf5$mut_type=='missense',c("pos",'het_function_score')]
scores$lowScore=scores$het_function_score<0.363
scores$highScore=scores$het_function_score>1.717
lowScoreResult = scores %>%
  group_by(pos) %>%
  summarize(true_count = sum(lowScore, na.rm = TRUE))
lowScoreResult$true_count=-1*lowScoreResult$true_count
highScoreResult <- scores %>%
  group_by(pos) %>%
  summarize(true_count = sum(highScore, na.rm = TRUE))
allResult=cbind(lowScoreResult,highScoreResult[,c(2)])
names(allResult)=c('pos','lowScore','highScore')
allResult_long <- reshape2::melt(allResult, id.vars = "pos", variable.name = "ScoreType", value.name = "Score")
ggplot(allResult_long, aes(x = pos, y = Score, fill = ScoreType)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.6) +
  scale_fill_manual(values = c("lowScore" = "red", "highScore" = "blue")) +
  theme_minimal() +
  labs(title = "Bar Graph of Het Function Scores Scores", x = "Position", y = "Score")
dev.off()

# Barplot - Proportion of Strong/Severe-DN by Position.
tf5_dn=tf5[,c("mutation","mut_type","mut_type2","orig","pos","new","criterionTotal")]
tf5_dn[tf5_dn$criterionTotal<=2,'criterionTotal2']=0
tf5_dn[tf5_dn$criterionTotal==3,'criterionTotal2']=1
tf5_dn[tf5_dn$criterionTotal==4,'criterionTotal2']=2
tf5_dn=tf5_dn[,c("mutation","mut_type","mut_type2","orig","pos","new","criterionTotal2")]
names(tf5_dn)[7]='criterionTotal'
all_positions <- unique(tf5_dn$pos)
scoreCounts <- tf5_dn %>%
  filter(criterionTotal %in% c(1, 2)) %>%
  group_by(pos, criterionTotal) %>%
  summarize(count = n(), .groups = 'drop')
scoreCounts <- scoreCounts %>%
  right_join(data.frame(pos = all_positions), by = "pos") 
scoreCounts$criterionTotal <- as.factor(scoreCounts$criterionTotal)
pdf('Validation/Proportion_StrongSevereDN_Barplot_10.9.25.pdf')
ggplot(scoreCounts, aes(x = as.factor(pos), y = count, fill = criterionTotal)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.7, na.rm = TRUE) +
  scale_fill_manual(values = c("1" = "orange", "2" = "red"), na.translate = FALSE) +
  geom_segment(aes(x = as.factor(pos), xend = as.factor(pos), y = 0, yend = 0), 
               color = "black", size = 1) +
  scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 10)]) +
  theme_minimal() +
  labs(title = "Stacked Bar Graph of CriterionTotal Counts per Position",
       x = "Position",
       y = "Count",
       fill = "Criterion Total") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

##### STEP 15: CHARGE SENSITIVITY #####
# OPTION 1 --don't include initial methionine.
nterm=tf5[(tf5$mut_type=='missense'|tf5$mut_type=='synon') & !is.na(tf5$trafficking_score) & tf5$pos<=100  & tf5$pos>1,]
plotChargeSensitivity(nterm,'Validation/Nterm_chargeSensitivity.pdf')

# OPTION 2--include initial methionine.
nterm=tf5[(tf5$mut_type=='missense'|tf5$mut_type=='synon') & !is.na(tf5$trafficking_score) & tf5$pos<=100,]
plotChargeSensitivity(nterm,'Validation/Nterm_chargeSensitivity_withM.pdf')

postHA=tf5[(tf5$mut_type=='missense'|tf5$mut_type=='synon') & !is.na(tf5$trafficking_score) & tf5$pos<=508  & tf5$pos>396,]
plotChargeSensitivity(postHA,'Validation/postHA_chargeSensitivity.pdf')
#p=1.0e-204, Kruskal-Wallis test

plotChargeSensitivity(tf5,'wholeGene_chargeSensitivity.pdf')
#p=2.5e-105, Kruskal-Wallis test

# Violin plot of residues 2-37
nterm37=tf5[(tf5$mut_type=='missense'|tf5$mut_type=='synon') & !is.na(tf5$trafficking_score) & tf5$pos<=37  & tf5$pos>1,]
plotChargeSensitivity(nterm37,'nTerm37_chargeSensitivity.pdf')
#p=1.5e-85, Kruskal-Wallis test

##### STEP 16: PRECISION RECALL PLOTS (AUC) #####
tf5b = tf5[!is.na(tf5$controlFinal_missense) & (tf5$controlFinal_missense=='blb' | tf5$controlFinal_missense=='plp'), ]
tf5b$controlFinal_missense <- factor(tf5b$controlFinal_missense)
tf5b = tf5[!is.na(tf5$controlFinal_missense) & (tf5$controlFinal_missense=='blb' | tf5$controlFinal_missense=='plp'), ]
tf5b$controlFinal_missense <- factor(tf5b$controlFinal_missense)

#AUC Scores
#Hom Traff AUC 0.88 , CI 0.84-0.93
#Hom Func AUC 0.94, CI 0.91-0.97 
#Het Traff AUC 0.80, CI 0.74-0.87
#Het Func AUC 0.81, CI 0.74-0.87

# Hom traf - ROC
tf5b_traf=tf5b[!is.na(tf5b$trafficking_score),]
scores <- tf5b_traf$trafficking_score
scores = abs(0.95-scores)
labels <- tf5b_traf$controlFinal_missense
roc_obj_hom_traff <- roc(labels, scores)

# Hom func - ROC
tf5b_func=tf5b[!is.na(tf5b$function_score),]
scores <- tf5b_func$function_score
scores = abs(0.9855-scores)
labels <- tf5b_func$controlFinal_missense
roc_obj_hom_func <- roc(labels, scores)

# Het traf - ROC
tf5b_traf=tf5b[!is.na(tf5b$het_trafficking_score),]
scores <- tf5b_traf$het_trafficking_score
scores = abs(1.0005-scores)
labels <- tf5b_traf$controlFinal_missense
roc_obj_het_traff <- roc(labels, scores)

# Het func - ROC
tf5b_func=tf5b[!is.na(tf5b$het_function_score),]
scores <- tf5b_func$het_function_score
scores = abs(1-scores)
labels <- tf5b_func$controlFinal_missense
roc_obj_het_func <- roc(labels, scores)

# Gaussian NB - ROC
tf5b_GNB=tf5b[!is.na(tf5b$GNB_prob),]
scores <- tf5b_GNB$GNB_prob
labels <- tf5b_GNB$controlFinal_missense
roc_obj_GNB <- roc(labels, scores)

# Logistic Regression - ROC
tf5b_LR=tf5b[!is.na(tf5b$LR_prob),]
scores <- tf5b_LR$LR_prob
labels <- tf5b_LR$controlFinal_missense
roc_obj_LR <- roc(labels, scores)

# Random Forest - ROC
tf5b_RF=tf5b[!is.na(tf5b$RF_prob),]
scores <- tf5b_RF$RF_prob
labels <- tf5b_RF$controlFinal_missense
roc_obj_RF <- roc(labels, scores)

# Plotting ROC curves 
pdf('Validation/ROC_Curves_9.11.25.pdf')
# Each of the 4 main assays+Gaussian NV
plot(roc_obj_hom_traff, main = "ROC Curve", col = "red", lwd = 2, xlim = c(1.1, 0.1))
auc_hom_traff <- auc(roc_obj_hom_traff)
ci_hom_traff <- ci(roc_obj_hom_traff)
lines(roc_obj_hom_func, col = "#0072B2", lwd = 2)
auc_hom_func <- auc(roc_obj_hom_func)
ci_hom_func <- ci(roc_obj_hom_func)
lines(roc_obj_het_traff, col = "#6B8E23", lwd = 2)
auc_het_traff <- auc(roc_obj_het_traff)
ci_het_traff <- ci(roc_obj_het_traff)
lines(roc_obj_het_func, col = "orange", lwd = 2)
auc_het_func <- auc(roc_obj_het_func)
ci_het_func <- ci(roc_obj_het_func)
legend("bottomright", legend = c(paste("Hom Traff (AUC =", round(auc_hom_traff, 2),')'),
                                 paste("Hom Func (AUC =", round(auc_hom_func, 2),')'),
                                 paste("Het Traff (AUC =", round(auc_het_traff, 2),')'),
                                 paste("Het Func (AUC =", round(auc_het_func, 2),')')),
       col = c("red", "#0072B2", "#6B8E23", "orange"), lwd = 2)

# Plotting 3 Merged predictors
plot(roc_obj_GNB, main = "ROC Curve", col = "purple", lwd = 2, xlim = c(1.1, 0.1))
auc_GNB <- auc(roc_obj_GNB)
ci_GNB <- ci(roc_obj_GNB)
lines(roc_obj_LR, col = "green", lwd = 2)
auc_LR <- auc(roc_obj_LR)
ci_LR <- ci(roc_obj_LR)
lines(roc_obj_RF, col = "blue", lwd = 2)
auc_RF <- auc(roc_obj_RF)
ci_RF <- ci(roc_obj_RF)
legend("bottomright", legend = c(paste("Gaussian NB (AUC =", round(auc_GNB, 2),')'),
                                 paste("Logistic Regression (AUC =", round(auc_LR, 2),')'),
                                 paste("Random Forest (AUC =", round(auc_RF, 2),')')),
       col = c("purple", "green", "blue"), lwd = 2)
dev.off()

##### STEP 17: PLOTTING 4 DMS SCORES AGAINST EACH OTHER #####
tf5_missense <- tf5[tf5$mut_type == 'missense', ]

# Plot Hom trafficking vs Hom function scores
cor.test(tf5$trafficking_score, tf5$function_score, method="spearman") #all variants: rho=0.60
cor.test(tf5_missense$trafficking_score, tf5_missense$function_score, method="spearman") #missense variants: rho=0.57
pdf('Validation/HomTrafficking_vs_HomFunction_raw_9.11.25.pdf')
plot(tf5$trafficking_score,tf5$function_score,pch='.',xlab='Hom trafficking score',ylab='Hom function score',main='All variants',xlim=c(-0.5,3),ylim=c(-0.5,3))
text(0.25,2.75,"rho=0.60")
plot(tf5_missense$trafficking_score, tf5_missense$function_score, pch='.', xlab='Hom trafficking score', ylab='Hom function score', main='Missense variants', xlim=c(-0.1, 1.55), ylim=c(-0.5, 3), xaxt = 'n')
axis(1, at = seq(-0.5, 1.5, by = 0.5)) # Custom x-axis labels
text(0.25, 2.75, "rho=0.57")
abline(v=0.25)
abline(v=0.61)
abline(v=1.29)
abline(h=0.25)
abline(h=0.58)
abline(h=1.391)
abline(-.5,1)
dev.off()

# Plot Het trafficking vs Het function scores
cor.test(tf5$het_trafficking_score, tf5$het_function_score, method="spearman", use="complete.obs")#all variants: rho=0.37
cor.test(tf5_missense$het_trafficking_score, tf5_missense$het_function_score, method="spearman") #missense variants: rho=0.46
pdf('Validation/HetTrafficking_vs_HetFunction_raw_9.11.25.pdf')
plot(tf5$het_trafficking_score,tf5$het_function_score,pch='.',xlab='Het Trafficking score',ylab='Het Function score',main='Het All variants',xlim=c(0,2),ylim=c(-2,5))
text(0.25,2.75,"rho=0.37")
plot(tf5_missense$het_trafficking_score,tf5_missense$het_function_score,pch='.',xlab='Trafficking score',ylab='Function score',main='Het Missense variants',xlim=c(0,2),ylim=c(-2,5))
text(0.25,2.75,"rho=0.46")
dev.off()

# Plot Hom vs Het trafficking scores
cor.test(tf5$trafficking_score, tf5$het_trafficking_score, method="spearman", use="complete.obs")#all variants: rho=0.56
cor.test(tf5_missense$trafficking_score, tf5_missense$het_trafficking_score, method="spearman", use="complete.obs")#missense: rho=0.69
pdf('Validation/HomTrafficking_vs_HetTrafficking_raw_9.11.25.pdf')
plot(tf5$trafficking_score,tf5$het_trafficking_score,pch='.',xlab='Hom trafficking score',ylab='Het trafficking score',main='Hom vs Het All variants')
text(0.25,2.75,"rho=0.69")
points(tf5$trafficking_score[tf5$trafficking_score < 0.25 & tf5$het_trafficking_score > 0.676], 
       tf5$het_trafficking_score[tf5$trafficking_score < 0.25 & tf5$het_trafficking_score > 0.676], 
       col='purple', pch='.')
points(tf5$trafficking_score[tf5$trafficking_score < 0.25 & tf5$het_trafficking_score >= 0.25 & tf5$het_trafficking_score <= 0.676], 
       tf5$het_trafficking_score[tf5$trafficking_score < 0.25 & tf5$het_trafficking_score >= 0.25 & tf5$het_trafficking_score <= 0.676], 
       col='blue', pch='.')
points(tf5$trafficking_score[tf5$trafficking_score < 0.25 & tf5$het_trafficking_score < 0.25], 
       tf5$het_trafficking_score[tf5$trafficking_score < 0.25 & tf5$het_trafficking_score < 0.25], 
       col='green', pch='.')
abline(v=0.25)
abline(h=0.25)
abline(h=0.676)
#plot(tf5_missense$trafficking_score,tf5_missense$het_trafficking_score,pch='.',xlab='Hom trafficking score',ylab='Het trafficking score',main='Hom vs Het Missense variants')
#text(0.25,1.75,"rho=0.69")
dev.off()

# Plot Hom vs Het function scores
cor.test(tf5$function_score, tf5$het_function_score, method="spearman", use="complete.obs")#all variants: rho=0.58
cor.test(tf5_missense$function_score, tf5_missense$het_function_score, method="spearman", use="complete.obs")#missense: rho=0.54
pdf('Validation/HomFunction_vs_HetFunction_10.14.25.pdf')
plot(tf5$function_score,tf5$het_function_score,pch='.',xlab='Hom function score',ylab='Het function score',main='Hom vs Het All variants')
text(0.25,2.75,"rho=0.58")
plot(tf5_missense$function_score,tf5_missense$het_function_score,pch='.',xlab='Hom function score',ylab='Het function score',main='Hom vs Het Missense variants')
text(0.25,1.75,"rho=0.54")
dev.off()

# Strong DN plots
# Hom vs Het trafficking: calculate deviation from best fit line:
tf5_LOT=tf5[tf5$trafficking_score<0.25 & !is.na(tf5$het_trafficking_score) & !is.na(tf5$trafficking_score) & (tf5$mut_type2=='early_nonsense' | tf5$mut_type2=='missense'),]
tf5_LOT$region='missense_cytosolic'
tf5_LOT[tf5_LOT$pos>=104 & tf5_LOT$pos<=(361),'region']='missense_transmembrane'
tf5_LOT[tf5_LOT$mut_type2=='early_nonsense','region']='early_nonsense'
tf5_LOT$region2b=tf5_LOT$region2
tf5_LOT[tf5_LOT$region2b=='mid_nonsense','region2b']='early_nonsense'
tf5_LOT[tf5_LOT$region2b=='early_nonsense','region2b']='A_early_nonsense'
tf5_LOT[tf5_LOT$region2b=='VSD_helix','region2b']='B_VSD_helix'
tf5_LOT[tf5_LOT$region2b=='Pore_helix','region2b']='C_Pore_helix'
tf5_LOT[tf5_LOT$region2b=='Cyto_helix','region2b']='D_Cyto_helix'
tf5_LOT[tf5_LOT$region2b=='other_missense','region2b']='E_other_missense'

summary(as.factor(tf5_LOT$region2b))
#A_early_nonsense      B_VSD_helix     C_Pore_helix     D_Cyto_helix E_other_missense 
#463                   613              757              514              204 

summary(as.factor(tf5_LOT$region3))
#early_nonsense             HA             HB             HC             HD   mid_nonsense other_missense              P 
#428            171            102            120            121             35            204            145 
#Ploop             S0             S1             S2           S2-3             S3             S4           S4-5 
#143             75             83            235             98             99             23             26 
#S5             S6 
#194            249 

# Dominant Negative Variants | Het Trafficking Score | Early nonsense, Missense Cytosolic, Missense Transmembrane
pdf("validation/DN_HetTraffickingScore_ByRegion_Cytosolic_vs_Transmembrane_9.11.25.pdf")
ggplot(tf5_LOT, aes(x = region2b, y = het_trafficking_score, fill = region2b)) +
  geom_violin(trim = TRUE, alpha = 0.7) + # Violin plot with color mapped to region
  labs(
    title = "Violin Plot of Het Trafficking Score by Region",
    x = "Region",
    y = "Het Trafficking Score"
  ) +
  scale_fill_brewer(palette = "Set2") + # Color palette for violins
  theme_classic() + # Classic theme for white background
  geom_hline(yintercept = 0.676) +
  geom_hline(yintercept = 0.25) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "none"
  )
dev.off()
dim(tf5_LOT) #2551
sum(tf5_LOT$het_trafficking_score<(0.25) & !is.na(tf5_LOT$het_trafficking_score)) #848 #dominant negative
sum(tf5_LOT$het_trafficking_score>=0.25 & tf5_LOT$het_trafficking_score<(0.676)) #876 #intermediate
sum(tf5_LOT$het_trafficking_score>=0.676) #827 haploinsufficient
#Of 2551 variants with strong loss of homozygous abundance, 848 had DN abundance, 827 had haploins abund. and remaining 876 were intermediate

# Breaking it down by location 
#(Strong/Severe-DN variants with <0.25 Het Trafficking)
sum(tf5_LOT$region2b=='A_early_nonsense' & tf5_LOT$het_trafficking_score<(0.25)) #0
sum(tf5_LOT$region2b=='B_VSD_helix' & tf5_LOT$het_trafficking_score<(0.25)) #375
sum(tf5_LOT$region2b=='C_Pore_helix' & tf5_LOT$het_trafficking_score<(0.25)) #335
sum(tf5_LOT$region2b=='D_Cyto_helix' & tf5_LOT$het_trafficking_score<(0.25)) #84
sum(tf5_LOT$region2b=='E_other_missense' & tf5_LOT$het_trafficking_score<(0.25)) #54

#(Strong/Severe-DN variants with >0.676 Het Trafficking)
sum(tf5_LOT$region2b=='A_early_nonsense' & tf5_LOT$het_trafficking_score>(0.676)) #448
sum(tf5_LOT$region2b=='B_VSD_helix' & tf5_LOT$het_trafficking_score>(0.676)) #8
sum(tf5_LOT$region2b=='C_Pore_helix' & tf5_LOT$het_trafficking_score>(0.676)) #120
sum(tf5_LOT$region2b=='D_Cyto_helix' & tf5_LOT$het_trafficking_score>(0.676)) #189
sum(tf5_LOT$region2b=='E_other_missense' & tf5_LOT$het_trafficking_score>(0.676)) #62

# Rank plot Hom Trafficking (Log Scale)
veryDN=tf5[tf5$criterion3or4 & tf5$mut_type=='missense',] #838
summary(as.factor(veryDN$region))
#missense_cytosolic missense_transmembrane 
#212                    626

summary(as.factor(veryDN$region2))
#Cyto_helix other_missense     Pore_helix     VSD_helix 
#191             52            273            322 
52/838 #6.2% other
100-(100*sum(veryDN$region2=='other_missense')/nrow(veryDN)) 
#94% of strong/severe DN variants were helical
(273+322)/838 #71% strong/severe DN in transmembrane helices
191/838 #23% strong/severe DN in cytoplasmic helices

#CountsForPaper
summary(as.factor(veryDN$region2)) #322 VSD_helix
summary(as.factor(tf5_LOT$region2)) #613 VSD_helix
#so 322/613=52.5% of VSD_helix LOF variants had DN effects.
# Of the 613 loss of abundance variants in the VSD domain, 322 (52.5%) 
#    had a dominant negative effect on abundance. The VSD had the highest fraction of
#    DN effects compared to the pore domain or the cytoplasmic domains.

# DN variants were spread throughout the protein, but were enriched in certain domains,
# especially helical variants in the VSD (Figure 4)

summary(as.factor(veryDN$region2))
#VSD_helix, Pore_helix, and Cyto_helix variants represented 322/838, 273/838, and 191/838, respectively
#(322+273+191)/838=93.8%
#93.8% of the 838 strong/severe DN missense variants were in helical regions, 595 (71.0%) in transmembrane helices and 191 (22.8%) in cytoplasmic helices.
sum(veryDN$criterion3or4) #838
sum(veryDN$criterionAll4) #290
#838-290=548 exactly 3
veryDNpos=veryDN$pos
veryDNFreq <- as.data.frame(table(veryDNpos))
colnames(veryDNFreq) <- c("Position", "veryDNMutants")
veryDNFreq$Position <- as.numeric(as.character(veryDNFreq$Position))
all_numbers <- data.frame(Position = 1:676)
merged_df <- merge(all_numbers, veryDNFreq, by = "Position", all.x = TRUE)
merged_df[is.na(merged_df)] <- 0
veryDNFreqSorted=merged_df[order(merged_df$veryDNMutants,decreasing=TRUE),]
veryDNFreqSorted$rank=1:nrow(veryDNFreqSorted)

# Counting number of PLP with strong+severe DN
veryDNplp=tf5[tf5$criterion3or4 & tf5$mut_type=='missense' & tf5$controlFinal_missense=='plp' & !is.na(tf5$controlFinal_missense),] #65
dim(veryDNplp) #69
count_criteria_3_plp <- sum(veryDNplp$criterionTotal == 3, na.rm = TRUE)
count_criteria_4_plp <- sum(veryDNplp$criterionTotal == 4, na.rm = TRUE)
print(paste("Count of criterionTotal == 3:", count_criteria_3_plp))  #43
print(paste("Count of criterionTotal == 4:", count_criteria_4_plp))  #26
#Of the 202 PLPs, 69 were strong/severe DN.
#43 PLP variants had Strong DN effects.
#26 PLP variants had Severe DN effects.

# Counting number of VUS with strong+severe DN
veryDNvus <- tf5[tf5$criterion3or4 & tf5$mut_type == 'missense' & !is.na(tf5$controlFinal_missense) & (tf5$controlFinal_missense == 'vus' | tf5$controlFinal_missense == 'conflicting'), ]
dim(veryDNvus) #35
count_criteria_3_vus <- sum(veryDNvus$criterionTotal == 3, na.rm = TRUE)
count_criteria_4_vus <- sum(veryDNvus$criterionTotal == 4, na.rm = TRUE)
print(paste("Count of criterionTotal == 3:", count_criteria_3_vus)) #27
print(paste("Count of criterionTotal == 4:", count_criteria_4_vus)) #11
#Of the 479 VUS, 35 were strong/severe DN.
#27 VUS/Conflicting variants had Strong DN effects. 
#11 VUS/Conflicting variants had Severe DN effects. 

# Make a table with pos, mutation, clinvar classification, DN effect. 
all_DN=tf5[tf5$criterion3or4 & (tf5$controlFinal_missense=='plp' | tf5$controlFinal_missense=='vus' | tf5$controlFinal_missense=='conflicting') & !is.na(tf5$controlFinal_missense),]
all_DN=all_DN[,c('pos','mutation','region3','controlFinal_missense','criterionTotal')]
all_DN$DN_effect <- ifelse(all_DN$criterionTotal == 3, 'strong-DN', 
                           ifelse(all_DN$criterionTotal == 4, 'severe-DN', NA))
all_DN <- all_DN %>%
  select(-criterionTotal) %>%
  arrange(pos)
print(all_DN)
write.csv(all_DN,'Validation/Strong and Severe DN Variants_9.11.25.csv',row.names=FALSE,)

# Plotting DN mutants on structure (strong+severe DN)
veryDNFreqSorted2 = merged_df[order(merged_df$Position, decreasing = FALSE),]
write.table(veryDNFreqSorted2,'dnFreq_9.11.25.txt',row.names=FALSE,col.names=FALSE) #move to PymolHeatmap folder

# Plotting DN mutants on structure (P/LP only+strong+severe DN)
veryDNplppos=veryDNplp$pos
veryDNplpFreq <- as.data.frame(table(veryDNplppos))
colnames(veryDNplpFreq) <- c("Position", "veryDNMutants")
veryDNplpFreq$Position <- as.numeric(as.character(veryDNplpFreq$Position))
all_numbers <- data.frame(Position = 1:676)
merged_df2 <- merge(all_numbers, veryDNplpFreq, by = "Position", all.x = TRUE)
merged_df2[is.na(merged_df2)] <- 0
veryDNFreqSorted3=merged_df2[order(merged_df2$veryDNMutants,decreasing=TRUE),]
veryDNFreqSorted3b=merged_df2[order(merged_df2$Position,decreasing=FALSE),]
write.table(veryDNFreqSorted3b,'dnplpFreq_9.11.25.txt',row.names=FALSE,col.names=FALSE) #move to PymolHeatmap folder

# Figure out exact helix locations of DN mutants for table
dnHigh=veryDNFreqSorted[veryDNFreqSorted$veryDNMutants>=10,]
dnHigh=addKCNQ1Position(dnHigh)
dnHigh2=dnHigh[order(dnHigh$Position),]
#Position veryDNMutants rank    region2 region3
#111      111            14    5  VSD_helix      S0
#114      114            14    6  VSD_helix      S0
#160      160            12   13  VSD_helix      S2
#171      171            11   15  VSD_helix      S2
#174      174            13    9  VSD_helix      S2
#176      176            13   10  VSD_helix      S2
#177      177            13   11  VSD_helix      S2
#178      178            13   12  VSD_helix      S2
#179      179            16    2  VSD_helix      S2
#184      184            15    4      Other   Other
#189      189            14    7  VSD_helix    S2-3
#193      193            11   16  VSD_helix    S2-3
#284      284            10   19 Pore_helix      S5
#314      314            17    1 Pore_helix   Ploop
#316      316            11   17 Pore_helix   Ploop
#317      317            10   20 Pore_helix   Ploop
#320      320            10   21 Pore_helix   Ploop
#346      346            10   22 Pore_helix      S6
#375      375            16    3 Cyto_helix      HA
#376      376            10   23 Cyto_helix      HA
#380      380            14    8 Cyto_helix      HA
#567      567            10   24 Cyto_helix      HC
#568      568            12   14      Other   Other
#588      588            11   18 Cyto_helix      HD

##### STEP 18: GATING MUTANTS  #####
# Plotting gating mutants on structure (high trafficking+low function)
tf5$diff_score=tf5$trafficking_score-tf5$function_score
gatingMutantDF=tf5[!is.na(tf5$diff_score) & tf5$diff_score>0.5 & tf5$trafficking_score>0.61 & tf5$function_score<0.58,] #339
#CountsForPaper
nrow(gatingMutantDF) #339 gating mutants (all missense)
"Notably, we identified 339 variants that exhibited loss-of-function without loss-of-abundance."
summary(as.factor(gatingMutantDF$controlFinal_missense)) #22 plp
#conflicting     missing         plp         vus        NA's 
#2           4          22           9         302 
gatingMutantDF[gatingMutantDF$mutation=='A341V'|gatingMutantDF$mutation=='G269S'|gatingMutantDF$mutation=='R366W',]
"This included 21 known pathogenic gating variants, including A341V"
plot(gatingMutantDF$trafficking_score,gatingMutantDF$function_score,pch='.',xlab='Trafficking score',ylab='Function score',main='Gating variants',xlim=c(-0.5,3),ylim=c(-0.5,3))
gatingMutants=gatingMutantDF$pos
gatingFreq <- as.data.frame(table(gatingMutants))
colnames(gatingFreq) <- c("Position", "GatingMutants")
all_numbers <- data.frame(Position = 1:676)
merged_df <- merge(all_numbers, gatingFreq, by = "Position", all.x = TRUE)
merged_df$GatingMutants[is.na(merged_df$GatingMutants)] <- 0
gatingFreqSorted=merged_df[order(merged_df$GatingMutants,decreasing=TRUE),]
merged_df[merged_df$GatingMutants>=10,"GatingMutants"]=10
write.table(merged_df,'gatingFreq_9.11.25.txt',row.names=FALSE,col.names=FALSE) #move to Structure/Gating

#Plotting trafficking vs function for control variants with line for gating variants
tf5_control=tf5[!is.na(tf5$controlFinal_missense) & (tf5$controlFinal_missense=='plp'|tf5$controlFinal_missense=='blb'),]
tf_B = tf5[tf5$controlFinal_missense=='blb' & !is.na(tf5$controlFinal_missense),]
tf_P = tf5[tf5$controlFinal_missense=='plp' & !is.na(tf5$controlFinal_missense),]
pdf('Validation/trafficking_vs_function_PB_9.11.25.pdf')
plot(tf5_control$trafficking_score, tf5_control$function_score, type='n', xlab = "Trafficking Score", ylab = "Function Score",
     xlim=c(-0.1,1.5),ylim=c(-0.5,1.5),main = "Trafficking vs Function Score")
points(tf_P$trafficking_score, tf_P$function_score, col = "red",  pch = 20)
points(tf_B$trafficking_score, tf_B$function_score, col = "blue", pch = 17)
abline(-.5,1)
abline(v=0.61)
abline(h=0.58)
gatingMutantDF[gatingMutantDF$mutation=='A341V'|gatingMutantDF$mutation=='G269S'|gatingMutantDF$mutation=='R366W',c('mutation','trafficking_score','function_score')]
#mutation trafficking_score function_score
#A341V         1.17      0.16
#G269S         0.88      0.34
#R366W         1.01      0.31
#points(c(1.17,0.88,1.01),c(0.16,0.34,0.31),col='purple')
points(c(1.17,1.01),c(0.16,0.31),col='purple') #A341V and R366W
dev.off()

# Mutants that disrupt calmodulin binding (16556865)
tf5[tf5$mutation=='R366W'|tf5$mutation=='S373P'|tf5$mutation=='W392R',c('mutation','trafficking_score','function_score')]
tf5[tf5$mutation=='R366P',c('mutation','trafficking_score','function_score')]
#Traf 0.97, function 0.24
#also gating 
tf5[tf5$mutation=='R366P',] #R366P and R366W disrupt calmodulin binding+have very abnormal current properties.

##### STEP 19: PATCH CLAMP VS DMS FUNCTION #####
# Hom Patch Clamp vs Hom DMS Function
cor.test(tf5$peakCurrent_lit, tf5$function_score, method = "spearman") #rho=0.64
sum(!is.na(tf5$peakCurrent_lit)) #223 variants
highlight_variants = c('R14D','R14T','R25E','L38N','E39F','G119V','S143F','A223L','R228A','L282P','L353P','K467A','K467R','E473I','P477L','I514T','P570L','L659I','P660A','T666Y')
tf5$color = ifelse(tf5$mutation %in% highlight_variants, '#C54B8C', 'black')
tf5$shape = ifelse(tf5$mutation %in% highlight_variants, 17, 16) # 17 is the code for triangle, 16 for circle
pdf('Validation/DMS_vs_patchClamp_9.11.25.pdf')
ggplot(tf5, aes(x = peakCurrent_lit, y = function_score)) +
  scale_x_continuous(trans = asinh_trans) +
  scale_y_continuous(trans = asinh_trans) +
  geom_point(aes(color = color, shape = shape)) +
  scale_color_identity() +
  scale_shape_identity() +
  stat_smooth(method = "loess", se = TRUE) +
  xlab("Manual Current Density (% of WT)") +
  ylab("DMS Function Score") +
  theme_minimal() +
  geom_hline(yintercept = 0, color = "darkgrey") +
  geom_hline(yintercept = 1, color = "darkgrey") +
  geom_vline(xintercept = 0, color = "darkgrey") +
  geom_vline(xintercept = 1, color = "darkgrey") +
  ggtitle('Lit manual vs DMS Function scores (n=223 variants)') +
  annotate("text", x = 2, y = -0.2, label = "Spearman Rho=0.64", color = "black", size = 5)
dev.off()

# Het Literature vs DMS Het Trafficking
cor.test(tf5$het_trafficking_lit, tf5$het_trafficking_score, method = "spearman") #rho=0.87
sum(!is.na(tf5$het_trafficking_lit)) #24 variants
pdf('Validation/HetDMSTrafficking_vs_HetTraffickingLit_10.15.25.pdf')
ggplot(tf4, aes(x = het_trafficking_lit, y = het_trafficking_score)) +
  scale_x_continuous(trans = asinh_trans()) +
  scale_y_continuous(trans = asinh_trans()) +
  geom_point() +
  stat_smooth(method = "loess", se = TRUE) +
  xlab("Manual Het Literature Score (% of WT)") +
  ylab("DMS Het Trafficking Score") +
  theme_minimal() +
  geom_hline(yintercept = 0, color = "darkgrey") +
  geom_hline(yintercept = 1, color = "darkgrey") +
  geom_vline(xintercept = 0, color = "darkgrey") +
  geom_vline(xintercept = 1, color = "darkgrey") +
  ggtitle('Literature Het manual trafficking vs DMS Het Trafficking scores (n=24 variants)')+
  annotate("text", x = 2, y = -0.2, label = "Spearman Rho=0.87", color = "black", size = 5)
dev.off()


# Het Patch Clamp vs Het DMS Function
cor.test(tf5$het_PeakCurrent_lit, tf5$het_function_score, method = "spearman") #rho=0.44
sum(!is.na(tf5$het_PeakCurrent_lit)) #64 variants
pdf('Validation/HetDMS_vs_HetpatchClampLit_9.11.25.pdf')
ggplot(tf5, aes(x = het_PeakCurrent_lit, y = het_function_score)) +
  scale_x_continuous(trans = asinh_trans()) +
  scale_y_continuous(trans = asinh_trans()) +
  geom_point() +
  stat_smooth(method = "loess", se = TRUE) +
  xlab("Manual Het Peak current (% of WT)") +
  ylab("DMS Het Function Score") +
  theme_minimal() +
  geom_hline(yintercept = 0, color = "darkgrey") +
  geom_hline(yintercept = 1, color = "darkgrey") +
  geom_vline(xintercept = 0, color = "darkgrey") +
  geom_vline(xintercept = 1, color = "darkgrey") +
  ggtitle('Literature Het manual vs DMS Het Trafficking scores (n=64 variants)')+
  annotate("text", x = 2, y = -0.2, label = "Spearman Rho=0.44", color = "black", size = 5)+
  ylim(-1.1, 3)+
  annotate("point", x = 0.33, y = 0.27, color = "red", size = 3)  # Red point at (0.3360182, 0.195)
dev.off()

G314S_TF=(tf5$mutation=="G314S")
tf5[G314S_TF,] #het_function_score 0.2444856, het_peakCurrent_lit 0.195

##### STEP 20: MANUAL TRAFFICKING VS DMS TRAFFICKING #####
# Pearson correlation calculation
# Hom trafficking - Pearson correlation calculation
tf5_missense=tf5[tf5$mut_type=='missense',]
cor.test(tf5_missense$trafficking_score, tf5_missense$trafficking_lit, method="spearman") #rho=0.87 (spearman)
sum(!is.na(tf5_missense$trafficking_lit)) #158 variants

#Het trafficking - Pearson correlation calculation (all points black)
tf5_missense=tf5[tf5$mut_type=='missense',]
cor.test(tf5_missense$het_trafficking_score, tf5_missense$het_trafficking_lit, method="spearman") #rho=0.87 (spearman)
sum(!is.na(tf5_missense$het_trafficking_lit)) #24 variants

#Het function - Pearson correlation calculation
tf5_missense=tf5[tf5$mut_type=='missense',]
cor.test(tf5_missense$het_function_score, tf5_missense$het_PeakCurrent_lit, method="spearman") #rho=0.47 (spearman)
sum(!is.na(tf5_missense$het_PeakCurrent_lit)) #63 variants

#Hom Trafficking - Plot MAVE scores vs literature
tf5_missense_under4=tf5_missense[!is.na(tf5_missense$trafficking_lit) & tf5_missense$trafficking_lit<4,] #removing 1 outlier point, missense only
pdf('Validation/Lit_vs_DMS_trafficking_9.11.25.pdf')
ggplot(tf5_missense_under4, aes(x = trafficking_lit, y = trafficking_score)) +
  geom_point() +
  scale_x_continuous(trans = asinh_trans()) +
  scale_y_continuous(trans = asinh_trans())+
  geom_hline(yintercept = 0, color = "darkgrey") +
  geom_hline(yintercept = 1, color = "darkgrey") +
  geom_vline(xintercept = 0, color = "darkgrey") +
  geom_vline(xintercept = 1, color = "darkgrey") +
  stat_smooth(method = "loess", se = TRUE) +
  ggtitle('Manual vs DMS trafficking scores (n=157 variants)') +
  annotate("text", x = 2, y = 0.3, label = "Spearman Rho=0.87", color = "black", size = 5) +
  theme_minimal()
dev.off()

#Hom Trafficking - Plot MAVE scores vs literature (rosy's data points are in red, literature data in black)
highlight_variants = c("R14D",'R14T','R25E','L38N','E39F','G119V','S143F','A223L','Q260K','L282P','A352T','L353P','S409R','K467A','K467R','E473I','P477L','E487K','L506V','I514T','N551Y','P570L','L659I','Q664W','T666Y')
tf5_missense_under4$highlight = ifelse(tf5_missense_under4$mutation %in% highlight_variants, '#C54B8C', 'black')
pdf('Validation/Lit_vs_DMS_trafficking_RED_triangles_for_Rosy_data_9.11.25.pdf')
ggplot(tf5_missense_under4, aes(x = trafficking_lit, y = trafficking_score)) +
  geom_point(aes(color = highlight, shape = highlight), size = 3) + # Apply color and shape to points
  scale_x_continuous(trans = asinh_trans()) +
  scale_y_continuous(trans = asinh_trans()) +
  geom_hline(yintercept = 0, color = "darkgrey") +
  geom_hline(yintercept = 1, color = "darkgrey") +
  geom_vline(xintercept = 0, color = "darkgrey") +
  geom_vline(xintercept = 1, color = "darkgrey") +
  stat_smooth(method = "loess", se = TRUE, color = "black") + # Force single smooth line
  ggtitle('Manual vs DMS trafficking scores (n=157 variants)') +
  annotate("text", x = 2, y = 0.3, label = "Spearman Rho=0.87", color = "black", size = 5) +
  scale_color_manual(values = c("#C54B8C" = "#C54B8C", "black" = "black")) + # Define colors
  scale_shape_manual(values = c("#C54B8C" = 17, "black" = 16)) + # Define shapes (triangle for red, circle for black)
  theme_minimal() +
  theme(legend.position = "none") # Remove legend if not needed
dev.off()

# Making a bargraph of single variant abundance validation data
a=read.csv('Validation/manual_KCNQ1_trafficking_validation_vumc_7.25.25.csv',stringsAsFactors = FALSE)
a2=a[(!a$Mutation=='WT' & !a$Mutation=='unstained' & !a$Mutation=='empty'),c('Mutation','AF647_norm','SE_trafficking','mutType')]
names(a2)=c('mutation','trafficking_1Var_thisStudy','trafficking_1Var_thisStudy_se','mutType')
a2_missense=a2[a2$mutType=='missense',] #we measured 26 missense variants and 7 synonymous variants
a3=merge(a2_missense,tf5[,c('mutation','trafficking_score','trafficking_score_sem','pos')],all.x=TRUE)
a3$trafficking_1Var_thisStudy=a3$trafficking_1Var_thisStudy/100
a3$trafficking_1Var_thisStudy_se=a3$trafficking_1Var_thisStudy_se/100
a3=a3[!is.na(a3$trafficking_score),]
a3=a3[order(a3$pos),]
a3$mutation <- factor(a3$mutation, levels =a3$mutation)
# How many of the missense variants we measured by flow had final DMS trafficking scores?
sum(!is.na(tf5_missense$trafficking_1Var_thisStudy) & !is.na(tf5_missense$trafficking_score)) #24 of 26

# How many synonymous variants had scores outside the normal range?
tf5_synon=tf5[tf5$mut_type=='synon' & !is.na(tf5$trafficking_score),] #541 synon scores total
hist(tf5_synon$trafficking_score,breaks=100)
sum(tf5_synon$trafficking_score<0.5) #21 below 0.5
21/541 #3.8%

# Barplot for Figure 2
pdf('Validation/trafficking_rosy_barplot_missense_9.11.25.pdf')
ggplot(a3, aes(x = mutation, y = trafficking_1Var_thisStudy)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_errorbar(aes(ymin = trafficking_1Var_thisStudy - trafficking_1Var_thisStudy_se, 
                    ymax = trafficking_1Var_thisStudy + trafficking_1Var_thisStudy_se), width = 0.2) +
  labs(title = "Manual Trafficking Scores", x = "Mutation", y = "Manual Trafficking Score") +
  theme_minimal()+
  scale_y_continuous(trans = asinh_trans)
ggplot(a3, aes(x = mutation, y = trafficking_score)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_errorbar(aes(ymin = trafficking_score - trafficking_score_sem, 
                    ymax = trafficking_score + trafficking_score_sem), width = 0.2) +
  labs(title = "DMS Trafficking Scores", x = "Mutation", y = "DMS Trafficking Score") +
  theme_minimal()+
  scale_y_continuous(trans = asinh_trans)
dev.off()

# Synonymous variant barplot
pdf('Validation/trafficking_rosy_barplot_synon_9.11.25.pdf')
a2_synon=a2[a2$mutType=='synon',] #we measured 26 missense variants and 7 synonymous variants
a2_synon$mutation <- factor(a2_synon$mutation, levels =a2_synon$mutation)
ggplot(a2_synon, aes(x = mutation, y = trafficking_1Var_thisStudy)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_errorbar(aes(ymin = trafficking_1Var_thisStudy - trafficking_1Var_thisStudy_se, 
                    ymax = trafficking_1Var_thisStudy + trafficking_1Var_thisStudy_se), width = 0.2) +
  labs(title = "Mean Trafficking Score by Mutation", x = "Mutation", y = "Mean Trafficking Score") +
  theme_minimal()
dev.off()

# Merge in real gnomAD
pdf('Validation/KCNQ1_scores_vs_gnomad_9.11.25.pdf')
plot(tf5$function_score,
     log10(tf5$AlleleFrequency),
     col = "black",
     pch = 1,
     xlab = "MAVE Function Score",
     ylab = expression(log[10]*"(Allele Frequency)"))
abline(v = 0.25,  lty = 2)
abline(v = 0.58,  lty = 2)
abline(v = 1.391, lty = 2)
dev.off()

highAF=tf5[tf5$AlleleFrequency>3*10^-5,]
sum(highAF$function_score<=1.391 & highAF$function_score>=0.58 & !is.na(highAF$function_score))
#35/36

# Choose which AF to use; switch to GroupMaxFAFfrequency if you prefer
tf <- subset(tf5, !is.na(function_score) & !is.na(AlleleFrequency))
tf$log10_AF <- log10(tf$AlleleFrequency)
x_lim <- c(-0.5, 3.0)
y_lim <- c(-6.5, -2.0)
v_guides <- c(0.5, 1.0, 1.25)
p <- ggplot(tf, aes(x = function_score, y = log10_AF)) +
  geom_point(size = 2, alpha = 0.7, color = "#E7B8AE") +   # peach points
  geom_vline(xintercept = v_guides, linetype = "dashed",
             linewidth = 0.4, color = "grey55") +
  labs(x = "MAVE Function Score",
       y = expression(log[10]*"(Allele Frequency)")) +
  coord_cartesian(xlim = x_lim, ylim = y_lim, expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(color = "grey20"),
    axis.text  = element_text(color = "grey20")
  ) +
  annotate("text", x = x_lim[1] + 0.02*diff(x_lim), 
           y = y_lim[2] + 0.05*diff(y_lim),
           label = "B", fontface = "bold", size = 6, hjust = 0, vjust = 1)
# size controls the relative height/width of the margins; tweak 'alpha' for fill
p_marg <- ggMarginal(
  p,
  type   = "density",
  margins = "both",
  size   = 6,
  alpha  = 0.25,
  fill   = "#E7B8AE",
  color  = "#E7B8AE"
)

p_marg


##### STEP 21: BARPLOT SYNCROPATCH VALIDATION DATA #####
a=read.csv('Validation/Syncropatch_KCNQ1_4.22.25.csv',stringsAsFactors = FALSE)
a2=a[,c('mutation','Tail.CD.mean','Tail.CD.SE')]
a3=merge(a2,tf5[,c('mutation','function_score','function_score_sem','pos')],all.x=TRUE)
a4=a3[!is.na(a3$function_score),]
a4=a4[order(a4$pos),]
a4$mutation <- factor(a4$mutation, levels =a4$mutation)

# Barplot for Figure 2_RosyData
pdf('Validation/function_SP_barplot_9.11.25.pdf')
ggplot(a4, aes(x = mutation, y = Tail.CD.mean)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_errorbar(aes(ymin = Tail.CD.mean - Tail.CD.SE, 
                    ymax = Tail.CD.mean + Tail.CD.SE), width = 0.2) +
  labs(title = "Manual Function Scores", x = "Mutation", y = "Manual Function Score") +
  theme_minimal()+
  scale_y_continuous(trans = asinh_trans)
ggplot(a4, aes(x = mutation, y = function_score)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_errorbar(aes(ymin = function_score - function_score_sem, 
                    ymax = function_score + function_score_sem), width = 0.2) +
  labs(title = "DMS Function Scores", x = "Mutation", y = "DMS Function Score") +
  theme_minimal()+
  scale_y_continuous(trans = asinh_trans)
dev.off()

##### STEP 22: AVERAGE SCORES FOR STRUCTURE HEATMAPS #####
# Calculating average trafficking scores for structure heatmap
tf5_missense=tf5[tf5$mut_type=='missense',]
tf5_missense$trans_score2=tf5_missense$trafficking_score
hist(tf5_missense$trans_score2)
df_avg <- tf5_missense %>%
  group_by(pos) %>%
  summarize(avg_trans_score2 = mean(trans_score2, na.rm = TRUE))
all_positions <- data.frame(pos = 1:676)
df_complete <- all_positions %>%
  left_join(df_avg, by = "pos") %>%
  mutate(avg_trans_score2 = ifelse(is.na(avg_trans_score2), 1, avg_trans_score2))
hist(df_complete$avg_trans_score2,breaks=100)
df_complete$avg_trans_score2[df_complete $avg_trans_score2>2]=2
df_complete $avg_trans_score2[df_complete $avg_trans_score2 <0]=0
df_complete[df_complete$pos==676,"avg_trans_score2"]=2
hist(df_complete$avg_trans_score2,breaks=100)
write.table(df_complete,'KCNQ1traf_scores.txt',col.names=FALSE,row.names=FALSE) #move to Structure/PymolHeatmap/

# Calculating average function scores
tf5_missense=tf5[tf5$mut_type=='missense',]
tf5_missense$trans_score2=tf5_missense$function_score
hist(tf5_missense$trans_score2)
df_avg <- tf5_missense %>%
  group_by(pos) %>%
  summarize(avg_trans_score2 = mean(trans_score2, na.rm = TRUE))
all_positions <- data.frame(pos = 1:676)
df_complete <- all_positions %>%
  left_join(df_avg, by = "pos") %>%
  mutate(avg_trans_score2 = ifelse(is.na(avg_trans_score2), 1, avg_trans_score2))
hist(df_complete$avg_trans_score2,breaks=100)
df_complete$avg_trans_score2[df_complete $avg_trans_score2>2]=2
df_complete $avg_trans_score2[df_complete $avg_trans_score2 <0]=0
hist(df_complete$avg_trans_score2,breaks=100)
write.table(df_complete,'KCNQ1func_scores.txt',col.names=FALSE,row.names=FALSE) #move to Structure/PymolHeatmap/

##### STEP 23: PORE DISTANCE #####
pore=read.csv('Validation/pore_distances.csv',stringsAsFactors = FALSE)
pore=pore[pore$chain=='A',]
pore=pore[,c('resi','distance_A')]
names(pore)=c('pos','pore_distance')
pore2=pore[pore$pos>=258 & pore$pos<=361,]
tf_pore=tf5[,c('pos','mutation','trafficking_score','function_score','mut_type')]
tf_pore2=merge(tf_pore,pore2)
tf_pore2=tf_pore2[tf_pore2$mut_type=='missense',]
plot(tf_pore2$pore_distance,tf_pore2$function_score,pch='.')

# Put each observation into a 0‑5, 5‑10 … 25‑30 Å bin
tf_pore2 <- tf_pore2 %>% 
  mutate(pore_bin = cut(
    pore_distance,
    breaks = c(0,5,10,30),              # 0,5,10,…,30
    labels = c("0–5", "5–10", "10+"),
    include.lowest = TRUE, right = FALSE))

# Violin plot (with a thin boxplot inside for the median/IQR)
pdf('Validation/pore_distance.pdf')
ggplot(tf_pore2, aes(x = pore_bin, y = function_score)) +
  geom_violin(trim = TRUE, fill = "steelblue", alpha = 0.6) +
  geom_boxplot(width = 0.10, outlier.shape = NA, alpha = 0.8) +
  labs(x = "Pore–axis distance (Å bins)",
       y = "Function score") +
  theme_classic(base_size = 12)+
  geom_hline(yintercept = 0.58, linetype = "dashed", linewidth = 0.4, color = "gray")
dev.off()  

sum(!is.na(tf_pore2$function_score) & tf_pore2$function_score<0.58 & tf_pore2$pore_distance<5) #118
sum(!is.na(tf_pore2$function_score) & tf_pore2$pore_distance<5) #129
# So 118/129 (91%) of missense variants in pore domain within 5A of pore had reduced function


##### STEP 24: AVERAGE MISSENSE SCORE FOR HC/HD COILED COIL DOMAINS #####
tf5_chargedP=tf5[tf5$new=='E'|tf5$new=='D'|tf5$new=='H'|tf5$new=='R'|tf5$new=='K'|tf5$new=='P',]
avg_trafficking_chargedP <- tf5_chargedP %>%
  group_by(pos) %>%
  summarize(avg_trafficking_chargedP = mean(trafficking_score, na.rm = TRUE))
avg_trafficking_chargedP=avg_trafficking_chargedP[avg_trafficking_chargedP$pos<=616 & avg_trafficking_chargedP$pos>=598,]
pdf('Validation/HD_coiled_coil_charged_9.11.25.pdf')
plot(avg_trafficking_chargedP$pos,avg_trafficking_chargedP$avg_trafficking_chargedP,ylim=c(0,1.1),xlab='Residue',ylab='Mean Score (charged mutations)')
text(599,0.05,"V599")
text(602,0.05,"L602")
text(606,0.05,"L606")
text(609,0.05,"I609")
text(613,0.05,"L613")
dev.off()

##### STEP 25: NONSENSE SCORE VS POSITION IN PROTEIN #####
traffic_nonsense=tf5[tf5$mut_type=='early_nonsense'| tf5$mut_type=='mid_nonsense'|tf5$mut_type=='late_nonsense',]
traffic_nonsense=traffic_nonsense[!is.na(traffic_nonsense$trafficking_score),]
functional_nonsense=tf5[tf5$mut_type=='early_nonsense'| tf5$mut_type=='mid_nonsense'|tf5$mut_type=='late_nonsense',]
functional_nonsense=functional_nonsense[!is.na(functional_nonsense$function_score),]
het_traffic_nonsense=tf5[tf5$mut_type=='early_nonsense'| tf5$mut_type=='mid_nonsense'|tf5$mut_type=='late_nonsense',]
het_traffic_nonsense=traffic_nonsense[!is.na(traffic_nonsense$het_trafficking_score),]
het_functional_nonsense=tf5[tf5$mut_type=='early_nonsense'| tf5$mut_type=='mid_nonsense'|tf5$mut_type=='late_nonsense',]
het_functional_nonsense=functional_nonsense[!is.na(functional_nonsense$het_function_score),]
#239-307 is mid_nonsense
pdf('Validation/KCNQ1_nonsense_vs_pos_pretty_9.11.25.pdf')
ggplot(traffic_nonsense, aes(x = pos, y = trafficking_score, color = mut_type)) +
  geom_point() + # Add points
  labs(title="Hom Trafficking MAVE",x = "Position in Protein", y = "Variant Score") + 
  theme_minimal()
ggplot(functional_nonsense, aes(x = pos, y = function_score, color = mut_type)) +
  geom_point() + # Add points
  labs(title="Hom Functional MAVE",x = "Position in Protein", y = "Variant Score") + 
  theme_minimal()
ggplot(het_traffic_nonsense, aes(x = pos, y = het_trafficking_score, color = mut_type)) +
  geom_point() + # Add points
  labs(title="Het Trafficking MAVE",x = "Position in Protein", y = "Variant Score") + 
  theme_minimal()
ggplot(het_functional_nonsense, aes(x = pos, y = het_function_score, color = mut_type)) +
  geom_point() + # Add points
  labs(title="Het Functional MAVE",x = "Position in Protein", y = "Variant Score") + 
  theme_minimal()
dev.off()

##### STEP 26: PLOT DMS VS PREDICTORS #####
comp=read.csv('Validation/KCNQ1_compPred.csv')
tf5b=merge(tf5,comp,all.x=TRUE)
tf5c=tf5b[tf5b$mut_type=='missense',]
pdf('./Validation/KCNQ1_AlphaMissense_vs_MAVE_withColors.pdf')
plot(tf5c$function_score,tf5c$alpha_missense,pch='.')
point_colors <- ifelse(tf5c$function_score < 0.58 & tf5c$alpha_missense > 0.906, "red", 
                       ifelse(tf5c$function_score >= 0.58 & tf5c$function_score <= 1.391 & 
                                tf5c$alpha_missense > 0.906, "orange",
                              ifelse(tf5c$function_score >= 0.58 & tf5c$function_score <= 1.391 & 
                                       tf5c$alpha_missense < 0.169, "purple", "black")))
plot(tf5c$function_score,tf5c$alpha_missense,pch='.',col=point_colors)
abline(h=0.792, lty=2,col='grey')
abline(h=0.906, lty=2,col='grey')
abline(h=0.99, lty=2,col='grey')
abline(h=0.169, lty=2,col='grey')
abline(h=0.099, lty=2,col='grey')
abline(v=0.25,lty=2,col='grey')
abline(v=0.58,lty=2,col='grey')
abline(v=1.391,lty=2,col='grey')
dev.off()

pdf('Validation/compPrediction_vs_functionScores.pdf')
plot(tf5c$function_score,tf5c$alpha_missense,pch='.')
plot(tf5c$function_score,tf5c$CPT1,pch='.')
plot(tf5c$function_score,tf5c$VARITY_R_LOO,pch='.')
plot(tf5c$function_score,tf5c$REVEL,pch='.')
cor.test(tf5c$function_score, tf5c$alpha_missense, method="spearman") #rho=-0.49
cor.test(tf5c$function_score, tf5c$CPT1, method="spearman") #rho=-0.51
cor.test(tf5c$function_score, tf5c$VARITY_R_LOO, method="spearman") #rho=-0.48
cor.test(tf5c$function_score, tf5c$REVEL, method="spearman") #rho=-0.47
dev.off()

tf5c_control=tf5c[!is.na(tf5c$controlFinal_missense) & (tf5c$controlFinal_missense=='blb'|tf5c$controlFinal_missense=='plp'|tf5c$controlFinal_missense=='vus'|tf5c$controlFinal_missense=='conflicting'),]
tf5c_control[tf5c_control$controlFinal_missense=='conflicting','controlFinal_missense']='vus'
tf5c_control[tf5c_control$controlFinal_missense=='vus','controlFinal_missense']='vus_confl'

plotPredictorVsMAVE('alpha_missense','./Validation/KCNQ1_AlphaMissense_vs_MAVE.pdf',tf5c,tf5c)
cor.test(tf5c$function_score, tf5c$alpha_missense, method="spearman") #all variants: rho=-0.49
plotPredictorVsMAVE('CPT1','./Validation/KCNQ1_CPT1_vs_MAVE.pdf',tf5c,tf5c)
plotPredictorVsMAVE('VARITY_R','./Validation/KCNQ1_VARITY_R_LOO_vs_MAVE.pdf',tf5c,tf5c)
pdf('./Validation/KCNQ1_computationPredictor_violin.pdf')
ggplot(tf5c_control, aes(x = controlFinal_missense, y = alpha_missense, fill = controlFinal_missense)) +
  geom_violin(scale='width') +
  geom_quasirandom(fill = 1, color = "black") +
  labs(x = 'Clinvar', y = 'Score', title = 'Clinvar comparison alpha_missense') +
  theme_classic()
ggplot(tf5c_control, aes(x = controlFinal_missense, y = CPT1, fill = controlFinal_missense)) +
  geom_violin(scale='width') +
  geom_quasirandom(fill = 1, color = "black") +
  labs(x = 'Clinvar', y = 'Score', title = 'Clinvar comparison CPT1') +
  theme_classic()
ggplot(tf5c_control, aes(x = controlFinal_missense, y = VARITY_R_LOO, fill = controlFinal_missense)) +
  geom_violin(scale='width') +
  geom_quasirandom(fill = 1, color = "black") +
  labs(x = 'Clinvar', y = 'Score', title = 'Clinvar comparison VARITY_R_LOO') +
  theme_classic()
dev.off()
cor.test(tf5c$function_score, tf5c$alpha_missense, method="spearman") #all variants: rho=-0.49
cor.test(tf5c$function_score, tf5c$CPT1, method="spearman") #all variants: rho=-0.51
cor.test(tf5c$function_score, tf5c$VARITY_R_LOO, method="spearman") #all variants: rho=-0.48

##### STEP 27: CLUSTERS BY CLASS/LOCATION #####
criteria_missense=tf5$mut_type=='missense'
criteria_synon=tf5$mut_type=='synon'
criteria_earlyNonsense=tf5$mut_type=='early_nonsense'
criteria_midNonsense=tf5$mut_type=='mid_nonsense'
criteria_lateNonsense=tf5$mut_type=='late_nonsense'
criteria_clinvarplus_missense_plp=tf5$controlFinal_missense=='plp'
criteria_clinvarplus_missense_blb=tf5$controlFinal_missense=='blb'
criteria_clinvarplus_all_plp=tf5$controlFinal=='plp'
criteria_clinvarplus_all_blb=tf5$controlFinal=='blb'
criteria_clinvar_missense_plp=tf5$clinvar=='plp' & tf5$mut_type=='missense'
criteria_clinvar_missense_blb=tf5$clinvar=='blb' & tf5$mut_type=='missense'
criteria_clinvar_all_plp=tf5$clinvar=='plp'
criteria_clinvar_all_blb=tf5$clinvar=='blb'
criteria_PY=tf5$pos>=659 & tf5$pos<=665 & tf5$mut_type=='missense' #92
criteria_severeDN = tf5$criterionAll4
criteria_strongDN = tf5$criterion3or4 & !(tf5$criterionAll4)

out_counts <- summarize_clusters(
  tf5,
  cluster_col = "cluster_n6_FINAL",
  criteria_missense = criteria_missense,
  criteria_synon = criteria_synon,
  criteria_earlyNonsense = criteria_earlyNonsense,
  criteria_midNonsense = criteria_midNonsense,
  criteria_lateNonsense = criteria_lateNonsense,
  criteria_clinvarplus_missense_plp = criteria_clinvarplus_missense_plp,
  criteria_clinvarplus_missense_blb = criteria_clinvarplus_missense_blb,
  criteria_clinvarplus_all_plp = criteria_clinvarplus_all_plp,
  criteria_clinvarplus_all_blb = criteria_clinvarplus_all_blb,
  criteria_clinvar_missense_plp = criteria_clinvar_missense_plp,
  criteria_clinvar_missense_blb = criteria_clinvar_missense_blb,
  criteria_clinvar_all_plp = criteria_clinvar_all_plp,
  criteria_clinvar_all_blb = criteria_clinvar_all_blb,
  criteria_PY=criteria_PY,
  criteria_severeDN=criteria_severeDN,
  criteria_strongDN=criteria_strongDN
)

out_counts=as.data.frame(out_counts)
write.csv(out_counts,'./Validation/Cluster_counts_by_class_10.20.25.csv',row.names=FALSE)

##### STEP 28: BIOBANK ANALYSES #####
# Read in datasets and calculate summary statistics for QTc
lqts_biobank=data.frame(function_class=c('Loss','PartialLoss','Normal','Gain'),lqt_fraction=c(0.081,0.019,0.005,0))
lqts_biobank$function_class <- factor(lqts_biobank$function_class,
                                      levels = c('Loss', 'PartialLoss', 'Normal', 'Gain'))
qtc_biobank=read.csv('./Validation/KCNQ1_Biobank_QTc_forR.csv')
summary(as.factor(qtc_biobank$function_score))
qtc_biobank$function_score <- factor(qtc_biobank$function_score,
                                      levels = c('LOF', 'partialLOF', 'normal', 'GOF'))
summary_stats_raw <- calculate_summary_stats_with_ci(qtc_biobank, "function_score", "QTc_raw")
as.data.frame(summary_stats_raw)
#function_score     mean       sd    n ci_lower ci_upper
#1            LOF 456.6667 26.34266   63 450.0324 463.3010
#2     partialLOF 438.1284 24.27737   74 432.5038 443.7530
#3         normal 430.7942 29.37519 2760 429.6978 431.8906
#4            GOF 432.4143 25.34958   70 426.3699 438.4587
456.7-430.8 #25.9 ms prolongation of LOF vs normal

# Plots
pdf('./Validation/Biobank_LQT_barplot.pdf')
ggplot(lqts_biobank, aes(x = function_class, y = lqt_fraction, fill = function_class)) +
  geom_bar(stat = "identity") +
  labs(title = "Barplot of LQT Fractions by Function Class",
       x = "Function Class",
       y = "LQT Fraction") +
  ylim(0, 0.1) + 
  theme_minimal() + 
  theme(legend.position = "none")
plot_mean_and_ci(
  data = summary_stats_raw,
  x_var = "function_score",
  y_var = "mean",
  ci_lower_var = "ci_lower",
  ci_upper_var = "ci_upper",
  x_label = "Function Score",
  y_label = "Mean QTc Raw",
  title = "QTc by Function Score (95%CI)",
  y_limits = c(420, 470) # Optional limits
)
dev.off()

##### STEP 29: WALSH VUS #####
walshVUS=tf5[tf5$mutation=='P570L'|tf5$mutation=='L282P'|tf5$mutation=='G119V'|tf5$mutation=='L353P'|tf5$mutation=='S143F'|tf5$mutation=='I514T',c('mutation','trafficking_score','function_score','peakCurrent_thisStudy','deltaV12act_thisStudy','llr_category')]
walshVUS=walshVUS[order(walshVUS$llr,decreasing=TRUE),]
write.csv(walshVUS,'KCNQ1_walshVUS.csv',row.names=FALSE)

##### STEP 30: MAVE DB FILES #####
# Trafficking
mavedb1=tf5[,c('mutation','trafficking_score','trafficking_score_sem')]
mavedb1$hgvs_pro=sapply(mavedb1$mutation,convert_variant1to3)
mavedb1=mavedb1[,c(4,2,3)]
names(mavedb1)=c('hgvs_pro','score','SE')
write.csv(mavedb1,'Current scores/MaveDB/KCNQ1_mavedb1_traf.csv',row.names=FALSE)

# Function
mavedb2=tf5[,c('mutation','function_score','function_score_sem')]
mavedb2$hgvs_pro=sapply(mavedb2$mutation,convert_variant1to3)
mavedb2=mavedb2[,c(4,2,3)]
names(mavedb2)=c('hgvs_pro','score','SE')
write.csv(mavedb2,'Current scores/MaveDB/KCNQ1_mavedb2_func.csv',row.names=FALSE)

# Het Trafficking
mavedb3=tf5[,c('mutation','het_trafficking_score','het_trafficking_sem')]
mavedb3$hgvs_pro=sapply(mavedb3$mutation,convert_variant1to3)
mavedb3=mavedb3[,c(4,2,3)]
names(mavedb3)=c('hgvs_pro','score','SE')
write.csv(mavedb3,'Current scores/MaveDB/KCNQ1_mavedb3_hettraf.csv',row.names=FALSE)

# Het Function
mavedb4=tf5[,c('mutation','het_function_score','het_function_score_sem')]
mavedb4$hgvs_pro=sapply(mavedb4$mutation,convert_variant1to3)
mavedb4=mavedb4[,c(4,2,3)]
names(mavedb4)=c('hgvs_pro','score','SE')
write.csv(mavedb4,'Current scores/MaveDB/KCNQ1_mavedb4_hetfunc.csv',row.names=FALSE)

##### STEP 31: SpliceAI variant counts #####
tf5_ms=tf5[tf5$mut_type=='missense'|tf5$mut_type=='synon',]
spliceai_high=tf5_ms[tf5_ms$spliceai_merged>=0.5 & !is.na(tf5_ms$spliceai_merged),] #44
spliceai_med=tf5_ms[tf5_ms$spliceai_merged>=0.2 & tf5_ms$spliceai_merged<0.5 & !is.na(tf5_ms$spliceai_merged),] #75
# So there are 44 missense/synon variants that have a spliceAI score above 0.5 
#    and 75 that have a score between 0.2 and 0.5.
spliceai_medHigh=tf5_ms[tf5_ms$spliceai_merged>=0.2 & !is.na(tf5_ms$spliceai_merged),] #119
spliceai_medHigh=spliceai_medHigh[order(spliceai_medHigh$pos),c(1,2,6,40:44)]
spliceai_medHigh=spliceai_medHigh[,-3]
write.csv(spliceai_medHigh,'Validation/kcnq1_spliceai_medHigh.csv',row.names=FALSE)
summary(as.factor(spliceai_medHigh$mut_type))
#missense    synon 
#82           37

