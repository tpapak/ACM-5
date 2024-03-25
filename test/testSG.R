rm(list=ls())
source("studygraph.R")
library(readxl)
ACM_5_with_networks <- 
  read_excel("data/mortality/ACM-5%-networks-overlap-correction-TP.xlsx", 
                                                   sheet = "data-ACM")

treatment_order <- read_excel("treatment_order.xlsx", 
                              skip = 1)

studyTable <- ACM_5_with_networks %>%
  filter(!is.na(effect)) %>%
  mutate(sef = ifelse(is.na(se_Adj),1, se_Adj))  %>% 
  mutate(var = (se*sef)^2)

studyIds <- unique(studyTable$id)
errorstudies <- c()
warnmess <- c()
studyInconsistency <- data.frame()
studyGraphs <- lapply(studyIds, function(stid){
  res <- {}
  tryCatch({
    stgr <- fullGraph(stid, studyTable ,rand=F)
    res <- stgr
  },error = function(cond){
    # print(c("error in study",stid))
    errorstudies <<- c(errorstudies,stid)
    # message(conditionMessage(cond))
  })
  return(res)
})
SGs <- Filter(function(g){!is.null(g)},studyGraphs)
names(SGs)<-lapply(SGs,function(x){return(x$study)}) %>% unlist()
netReport(studyTable, SGs, 2, outcome="mortality")

# pairwiseReport(studyTable, SGs, 7, outcome="mortality")

# networks <- netReport(studyTable, SGs, 1, outcome="mortality", subgroup="sex")
# # st <- studyTable %>% filter(id=="115") %>% filter(`treat 1`!="MUFA")
# st <- studyTable %>% filter(`treat 2`!="PRO")
# sg <- fullGraph("44",st)
# plotSG(sg)

# netid <- 2
# 
# sg <- fullGraph("191",studyTable)
# n2s <- prepareStudyForNetwork(sg,netid,studyTable)
# n2sr <- rowsFromGraph(n2s)
# n2sr
