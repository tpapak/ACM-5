library(matlib)
library(dplyr)
library(igraph)
library(MASS)
getArms <- function(study,tr1name="treat 1",tr2name="treat 2"){
  res <- unique(unlist(c(study[tr1name],study[tr2name])))
  res
}

getVertices <- function(strows) {
  verticesnames <- c(strows$`treat 1`,strows$`treat 2`) %>% 
    unique() %>% sort()
  vertices <- data.frame(id=1:length(verticesnames), label=verticesnames)
  return(vertices)
}

getVertexId <- function(name,vertices){
  subset(vertices, label==name)$id 
}

getVertexLabel <- function(vid,vertices){
  subset(vertices, id==vid)$label
}

edgeFromRow <- function(row,vertices) {
  fromV <- getVertexId(row$`treat 1`,vertices)
  toV <- getVertexId(row$`treat 2`,vertices)
  return(c(fromV,toV))
}

studyGraph <- function(studyid, data) {
  strows <- subset(data, id==studyid) %>%
    dplyr::select(id, `treat 1`, `treat 2`, effect, var)
  vertices <- getVertices(strows)
  edges <- sapply(1:nrow(strows),function(r){
    e <- edgeFromRow(strows[r,],vertices)
    return(e)
  }) %>% c()
  stgr <- igraph::graph(edges, directed=F)
  vertex_attr(stgr, name="label")<-vertices$label
  # plot(stgr)
  edgeList <- ends(stgr,E(stgr))
  edgeList <- cbind(edgeList,1:nrow(edgeList)) #edge index
  vars <- mapply(function(ei){
    labelfrom <- getVertexLabel(edgeList[ei,1],vertices)
    labelto <- getVertexLabel(edgeList[ei,2],vertices)
    vi <- strows[strows$`treat 1`==labelfrom & strows$`treat 2`==labelto | 
                 strows$`treat 1`==labelto & strows$`treat 2`==labelfrom , ]
    # print(c("vi",vi))
    return(vi$var)
  },1:nrow(edgeList))  %>% unlist()
  effs <- mapply(function(ei){
    labelfrom <- getVertexLabel(edgeList[ei,1],vertices)
    labelto <- getVertexLabel(edgeList[ei,2],vertices)
    vi <- strows[strows$`treat 1`==labelfrom & strows$`treat 2`==labelto | 
                 strows$`treat 1`==labelto & strows$`treat 2`==labelfrom , ]
    # print(c("vi",vi))
    return(vi$effect)
  },1:nrow(edgeList))  %>% unlist()
  treat1 <- mapply(function(ei){
    labelfrom <- getVertexLabel(edgeList[ei,1],vertices)
    labelto <- getVertexLabel(edgeList[ei,2],vertices)
    vi <- strows[strows$`treat 1`==labelfrom & strows$`treat 2`==labelto | 
                 strows$`treat 1`==labelto & strows$`treat 2`==labelfrom , ]
    # print(c("vi",vi))
    return(vi$`treat 1`)
  },1:nrow(edgeList))  %>% unlist()
  graph_attr(stgr,name="study")<-studyid
  edge_attr(stgr,name="variance",index=E(stgr))<-vars 
  edge_attr(stgr,name="effect",index=E(stgr))<-effs
  edge_attr(stgr,name="treat1",index=E(stgr))<-treat1
  return(stgr)
}

#Get minimum spanning tree by weighting the comparisons variance
studySpanningTree <- function(studygraph, rand=F){
  #study weights
  if(rand){
    spt <- subgraph.edges(studygraph,sample_spanning_tree(studygraph))
  }else{
    ws <- E(studygraph)$variance
    spt <- mst(studygraph, weights = ws)
  }
  return(spt)
}

#form A matrix where AxVi=Vij and thus Vi=A^-1 x Vij
varianceGraph <- function(studygraph, spt){
  #Get the edge with the smallest variance so as to avoid negative variances
  loopedges <- difference(studygraph,spt)
  if(ecount(loopedges)>0){
    loopedgeId <- which.min(E(loopedges)$variance)
    loopedge <- E(loopedges)[loopedgeId]
    vg <- add_edges(spt, ends(loopedges,loopedge), variance=loopedge$variance)
  }else{
    vss <- data.frame(vid=V(spt)
                      ,label=V(spt)$label
                      ,degree=degree(spt,V(spt))
                      ,strength=strength(spt,V(spt),weights=E(spt)$variance))%>%
              filter(degree>1) %>% 
              arrange(strength, degree)
    cornerV <- vss[1,]
    eij <- neighbors(spt, cornerV$vid)
    vi <- eij[1]
    vj <- eij[2]
    eijs <- get.edge.ids(spt, c(cornerV$vid,vi,cornerV$vid,vj))
    varij <- E(spt)[eijs]$variance %>% max()
    vg <- add_edges(spt, c(vi,vj), variance=varij)
  }
  return(vg)
}

varianceMatrix <- function(variancegraph){
  vg <- variancegraph
  vm <- matrix(nrow=ecount(vg),ncol=vcount(vg),0)
  vmm <- Reduce(function(acc,ei){
    evs <- ends(vg,E(vg))[ei,]
    acc[ei,evs]<-1
    return(acc)
  },1:ecount(vg),vm)
 return(vmm) 
}

treatVars <- function(studygraph, spt){
  if(ecount(spt)>1){
    vg <- varianceGraph(studygraph, spt)
    mv <- varianceMatrix(vg)
    vij <- matrix(ncol=1, E(vg)$variance)
    ##If there higher than triangle order polygons the design matrix
    #is not invertible so we use the pseudoinverse of MASS::ginv
    vars <- MASS::ginv(mv) %*% vij
    #Remove all loops so no inconsistencies happen 
    if(any(vars<0)){
        vg <- varianceGraph(spt, spt)
        loopedges <- difference(studygraph, spt)
        edgelbels <- paste(E(loopedges)$comparison,collapse=",")
        wmes <- paste( "found negative variances in study "
                     , studygraph$study
                     , " the following comparisons will be ignored in order to get consistent variances "
                     , edgelbels
                     , sep="")
        warning(wmes)
        mv <- varianceMatrix(vg)
        vij <- matrix(ncol=1, E(vg)$variance)
        vars <- MASS::ginv(mv) %*% vij
      }
  }else{
    vars <- matrix(ncol=1,rep(E(spt)$variance/2,2))
  }
  if(any(vars<0)){stop("found negative variances")}
  return(vars)
}

#the input is a spanning tree
treatEffects <- function(sp){
  rootVertex <- 1 
  sdfs <- dfs(sp,rootVertex,father=T)
  efm <- matrix(nrow=ecount(sp),ncol=vcount(sp),0)
  efmm <- Reduce(function(acc,ei){
    evs <- ends(sp,E(sp))[ei,]
    if(evs[2] != rootVertex){
      if(sdfs$father[evs[2]]==evs[1]){
        vi <- evs[1]
        vj <- evs[2]
      }else{
        vi <- evs[2]
        vj <- evs[1]
      }
      if(E(sp)$treat1[ei] == V(sp)[vi]$label){
        acc[ei,vi] <- -1
        acc[ei,vj] <- 1
      }else{
        acc[ei,vi] <- 1
        acc[ei,vj] <- -1
      }
    }
    return(acc)
  },1:ecount(sp),efm)[,-1]
  eij <- matrix(ncol=1, E(sp)$effect)
  if(ecount(sp)>1){
    effs <- inv(efmm) %*% eij
  }else{
    effs<-efmm*E(sp)$effect
  }
  return(effs)
}

labelEdges <- function (gr){
  getEdgeLabel <- function(gr, eid){
    vi <- ends(gr, E(gr))[eid,1]
    vj <- ends(gr, E(gr))[eid,2]
    paste(V(gr)$label[vi],V(gr)$label[vj],sep=":")
  }
  edgelabels <- lapply(1:ecount(gr), function(ei){
    return(getEdgeLabel(gr,ei))
  }) %>% unlist()
  edge_attr(gr, name="comparison") <- edgelabels
  return(gr)
}

fullGraph <- function(studyid, data, rand=F){
  gr <- studyGraph(studyid, data)
  edge_attr(gr,name="oldvariance")<-E(gr)$variance 
  gr <- labelEdges(gr)
  spt <- studySpanningTree(gr, rand)
  vars <- treatVars(gr, spt)
  vertex_attr(gr, name="vi")<-vars
  effs <- treatEffects(spt)
  vertex_attr(gr, name="effi")<-c(0,effs)
  directed_edges <- lapply(1:ecount(gr), function(ei){
    eij <- ends(gr,E(gr))[ei,]
    if(E(gr)$treat1[ei] == V(gr)[eij[1]]$label){
      rese <- eij
    }else{
      rese <- c(eij[2], eij[1])
    }
    return(rese)
  }) 
  dgr <- graph(unlist(directed_edges),directed=T)
  graph_attr(dgr,name="study")<-studyid
  edge.attributes(dgr)<-edge.attributes(gr)
  vertex.attributes(dgr)<-vertex.attributes(gr)
  dgr <- inconsistency(dgr, spt)
  dgr <- labelEdges(dgr)
  return(dgr)
}

#get the difference of the reported 
inconsistency <- function(gr, spt){
  loopedges <- difference(as.undirected(gr), spt)
  diffs <- lapply(1:ecount(gr), function(ei){
    eij <- ends(gr,E(gr))[ei,]
    vi <- eij[1]
    vj <- eij[2]
    if(get.edge.ids(loopedges,c(vi,vj))>0){
      #found effect
      teff <- V(gr)$effi[vj]-V(gr)$effi[vi]
      #reported effect
      reff <- E(gr)$effect[ei]
      var <- E(gr)$variance[ei]
      diff <- reff - teff
      res <- abs(diff)
    }else{
      res <- NA
    }
    return(res)
  }) %>% unlist()
  varRess <- lapply(1:ecount(gr), function(ei){
    eij <- ends(gr,E(gr))[ei,]
    vi <- eij[1]
    vj <- eij[2]
    if(get.edge.ids(loopedges,c(vi,vj))>0){
      #found effect
      teff <- V(gr)$vi[vj]+V(gr)$vi[vi]
      #reported effect
      var <- E(gr)$oldvariance[ei]
      diff <- teff - var
      res <- abs(diff)
    }else{
      res <- NA
    }
    return(res)
  }) %>% unlist()
  qis <- lapply(1:ecount(gr), function(ei){
    eij <- ends(gr,E(gr))[ei,]
    vi <- eij[1]
    vj <- eij[2]
    if(get.edge.ids(loopedges,c(vi,vj))>0){
      #found effect
      teff <- V(gr)$effi[vj]-V(gr)$effi[vi]
      #reported effect
      reff <- E(gr)$effect[ei]
      var <- E(gr)$variance[ei]
      diff <- reff - teff
      res <- diff^2/var
    }else{
      res <- NA
    }
    return(res)
  }) %>% unlist()
  vars <- lapply(1:ecount(gr), function(ei){
    var <- E(gr)$variance[ei]
    eij <- ends(gr,E(gr))[ei,]
    vi <- eij[1]
    vj <- eij[2]
    if(get.edge.ids(loopedges,c(vi,vj))>0){
      res <- 1/var
    }else{
      res <- NA
    }
    return(res)
  }) %>% unlist()
  edge_attr(gr, name="diff")<-diffs
  edge_attr(gr, name="varRes")<-varRess
  nw <- sum (vars,na.rm=T)
  Q <- sum(qis, na.rm=T)
  dofs <- ecount(loopedges)
  pvalue <- pchisq(Q, dofs, lower.tail = F)
  graph_attr(gr, name="Q")<-Q
  graph_attr(gr, name="dofs")<-dofs
  graph_attr(gr, name="pvalue")<-pvalue
  return(gr)
}

plotVarRes <- function(gr,circ=F){
  if(circ==F){
  plot( gr, edge.label=round(E(gr)$varRes,4)
      , main=paste(gr$study,"Variance Residuals")
      )
  }else{
    plot( gr, edge.label=round(E(gr)$varRes,4)
        , main=paste(gr$study,"Variance Residuals")
        , layout=layout_in_circle(sg)
        )
  }
}

plotDiffs <- function(gr,circ=F){
  if(circ==F){
  plot( gr, edge.label=round(E(gr)$diff,4)
      , main=paste(gr$study,"Residuals Effects")
      )
  }else{
    plot( gr, edge.label=round(E(gr)$diff,4)
        , main=paste(gr$study,"Residuals Effects")
        , layout=layout_in_circle(sg)
        )
  }
}

plotVars <- function(gr,circ=F){
  if(circ==F){
    plot(gr, edge.label=round(E(sg)$variance,4)
        , main=paste(gr$study,"variances")
    )
  }else{
    plot(gr, edge.label=round(E(sg)$variance,4)
       , layout=layout_in_circle(gr)
        , main=paste(gr$study,"variances")
    )
         
    
  }
}

plotEffects <- function(gr,circ=F){
  if(circ==F){
    plot(gr
        , edge.label=round(E(sg)$effect,4)
        , main=paste(gr$study,"effects")
    )
  }else{
    plot(gr
         , edge.label=round(E(sg)$effect,4)
          , main=paste(gr$study,"effects")
         , layout=layout_in_circle(gr)
    )
  }
}

prepareStudyForNetwork <- function(gr, netid, studyTable){
  sid <- gr$study
  nt1col <- paste("N",netid,"treat1",sep="")
  nt2col <- paste("N",netid,"treat2",sep="")
  strows <- studyTable %>% filter(id == sid)
  strows <- strows[!is.na(strows[nt1col]),]
  
  
  netnodes <- c(strows[nt1col],strows[nt2col]) %>% unlist () %>% 
    unique() %>% sort()
  newvertices <- data.frame(vid=1:length(netnodes),label=netnodes)
  cg <- make_full_graph(nrow(newvertices),directed=F)
  newEdges <- ends(cg,E(cg)) %>% t() %>% as.vector()
  newSG <- graph(newEdges, directed=T)
  graph_attr(newSG, name="study")<-sid
  vertex_attr(newSG, name="label")<-newvertices$label
  fullSG <- Reduce(function(acc, i){
    nv <- V(newSG)$label[i]
    vrs1 <- strows[strows[nt1col]==nv,"treat 1"] %>% as.vector()
    vrs2 <- strows[strows[nt2col]==nv,"treat 2"] %>% as.vector()
    mergedNodes <- c(vrs1,vrs2) %>% unlist()
    vids <- which(V(gr)$label %in% mergedNodes)
    varis <- V(gr)$vi[vids]
    effis <- V(gr)$effi[vids]
    nw <- sum(sapply(X=varis,FUN=function(v){return(1/v)}))
    newEffi <- sum(sapply(X=1:length(effis)
                  ,FUN=function(vid){return(effis[vid]/varis[vid])})) / nw
    vertex_attr(acc,index=i,name="effi")<-newEffi
    vertex_attr(acc,index=i,name="vari")<-1/nw
    return(acc)
  }, 1:nrow(newvertices),newSG)
  newEffs <- mapply(function(ei){
    eij <- ends(fullSG,ei)
    vi <- eij[1]
    vj <- eij[2]
    effij <- V(fullSG)$effi[vj] -
             V(fullSG)$effi[vi]
    return(effij) 
  }, 1:ecount(newSG))
  newVariances <- mapply(function(ei){
    eij <- ends(fullSG,ei)
    vi <- eij[1]
    vj <- eij[2]
    varij <- V(fullSG)$vari[vj] +
             V(fullSG)$vari[vi]
    return(varij) 
  }, 1:ecount(newSG))
  edge_attr(fullSG,name="effect")<-newEffs
  edge_attr(fullSG,name="variance")<-newVariances 
 return(fullSG) 
}

#Get rows of a multiarm study for netmeta consumption
rowsFromGraph <- function(fg){
  rows <- mapply(function(ei){
    eij <- ends(fg,ei)
    t1 <- V(fg)$label[eij[1]]
    t2 <- V(fg)$label[eij[2]]
    ri <- c( studlab=fg$study
           , TE=E(fg)$effect[ei]
           , seTE=sqrt(E(fg)$variance[ei])
           , treat1=t1
           , treat2=t2)
    return(ri)
  },1:ecount(fg))
  res <- data.frame(studlab=rows[1,]
                    ,TE=as.numeric(rows[2,])
                    ,seTE=as.numeric(rows[3,])
                    ,treat1=rows[4,]
                    ,treat2=rows[5,])
  return(res)
}

# library(readxl)
# ACM_5_with_networks <- read_excel("data/ACM-5%-with-networks-without-Argos.xlsx")
# #
# #Insert variance
# studyTable <- ACM_5_with_networks %>%
#   filter(!is.na(effect)) %>%
#   mutate(var = se^2)
# 
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
