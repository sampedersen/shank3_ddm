library(writexl)

for (i in 1:2) {
  file_path <- file.choose()
  file_name <- basename(file_path)
  
  if (grepl("^HE", file_name)) {
    load(file_path)
    gg <- "HE"
    message("HE group loaded")
    
  } else if (grepl("^WT", file_name)) {
    load(file_path)
    gg <- "WT"
    message("WT group loaded")
    
  } else {
    next  # skip file if neither HE nor WT
  }
  
  Summary_df <- as.data.frame(Summary)
  Summary_df <- cbind(Metric = rownames(Summary), Summary_df)
  output_path <- file.path(dirname(file_path), paste0(file_name, "_Summary.xlsx"))
  
  write_xlsx(Summary_df, output_path)
}



path = "C:\\Users\\Feede\\Documents\\Sam_Temp\\Shank3\\Outputs\\April-23\\16-27\\"
file_name = "HE_April-23_16-27.RData"
load(paste0(path,"HE_April-23_16-27.RData"))
output = file.path(dirname(path), paste0(file_name, "_Summary.xlsx"))
write_xlsx(Summary_df,output)



n_chains = 6
params=NULL
for (i in 1:2) {
  file_path <- file.choose()
  file_name <- basename(file_path)
  if (grepl("^HE", file_name)) {
    load(file_path)
    gg = "HE"
    message("HE group loaded")
  } else if (grepl("^WT", file_name)) {
    load(file_path)
    gg = "WT"
    message("WT group loaded")
  }
  Summary %>% as.data.frame() %>%
    filter(psrf > 1.1)
  for(j in 1:n_chains){
    chain=rbind(Results$mcmc[j])
  }
  subj=unique(Data$idxP)
  for (s in subj){
    subj_idx = unique(Data$subj_idx[Data$idxP ==s])
    tmp_res = data.frame(
      idxP = s,
      subj_idx = subj_idx,
      group = gg,
      wOffer = mean(chain[,c( paste( c("b1.p[",toString(s),"]"), collapse = ""))]),
      boundary = mean(chain[,c( paste( c("alpha.p[",toString(s),"]"), collapse = ""))]),
      nDT = mean(chain[,c( paste( c("theta.p[",toString(s),"]"), collapse = ""))]),
      bias = mean(chain[,c( paste( c("bias.p[",toString(s),"]"), collapse = ""))])
    )
    params=rbind(params,tmp_res)
    
    params %>%
      pivot_longer(cols = c(wOffer,boundary,nDT,bias)) %>%
      group_by(group,name) %>%
      summarise(tvalue = t.test(value ~ session, paired = TRUE)$statistic,
                pvalue = t.test(value ~ session, paired = TRUE)$p.value) %>%
      filter(pvalue < .05)
    
    params %>%
      pivot_longer(cols = c(wOffer,boundary,nDT,bias)) %>%
      mutate(group = factor(group, levels = c("HC","BD"))) %>%
      group_by(session,name) %>%
      summarise(tvalue = t.test(value ~ group)$statistic,
                pvalue = t.test(value ~ group)$p.value) %>%
      filter(pvalue < .05)
    
    summary(aov(data = params,
                formula = wOffer ~ group))
    TukeyHSD(aov(data = params,
                 formula = wOffer ~ group),which = "group")
    
    
    ggplot(params,aes(x = group, y = wOffer, color = group, group = group)) +
      theme_pubr(base_size = 18) +
      geom_point(position = position_dodge2(width=.25), alpha = .5,size=3) +
      geom_hline(yintercept =0,linetype="dashed",color="grey",linewidth=1) +
      stat_summary(geom = "line",position = position_dodge2(width=.25),linewidth=1.5) +
      stat_summary(position = position_dodge2(width=.25),size=1.5,linewidth=1.5) +
      scale_color_brewer(type="qual",palette = 4) +
      labs(x = "Dx",
           y = "Influence of Offer on Drift Rate",
           color = "Dx")
  }
  
  }

# Load the .RData file and capture the loaded object names
loaded_objects <- load(file_path)

# Now you have:
# - file_path: the path of the file you loaded
# - loaded_objects: a character vector of object names loaded into the environment

# Example output
print(paste("Loaded from:", file_path))
print("Objects loaded:")
print(loaded_objects)


print(file_name)

