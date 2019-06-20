

coexp.gdp = paste0(gdp.coexp(),"/data/")

#If mode == "FULL" then all networks are expected to exist, including
#those in supplementary. If not, only those in /data/rdsnet/ will
#be used

coexp.initDbGTEX = function(){
  coexp.nets$gtexv6 <<- NULL
  coexp.data$gtexv6 <<- NULL
  nets = coexp.getGTExTissues()
  for(net in nets){
    #cat(paste0("Loading GTEx ",net,"\n"))
    coexp.nets$gtexv6[[net]] <<- coexp.getGTExNet(net)
    coexp.data$gtexv6[[net]] <<- paste0(gdp.coexp(),"/supplementary/rdsnets/gtexv6/",net,".resids.rds")
    coexp.go$gtexv6[[net]] <<- paste0(coexp.nets$gtexv6[[net]],"_gprof.csv")
    coexp.ctype$gtexv6[[net]] <<- paste0(coexp.nets$gtexv6[[net]],".celltype.csv")
    stopifnot(file.exists(coexp.nets$gtexv6[[net]]))
    stopifnot(file.exists(coexp.data$gtexv6[[net]]))
    stopifnot(file.exists(coexp.go$gtexv6[[net]]))
    stopifnot(file.exists(coexp.ctype$gtexv6[[net]]))
  }
}


coexp.initDb = function(coexp.mode = "FULL"){
  coexp.nets <<- NULL
  coexp.data <<- NULL
  coexp.ctype <<- NULL
  coexp.go <<- NULL
  
  if(coexp.mode == "GTEX")
    coexp.initDbGTEX()    
  
  #Init gtexv6 nets
  if(coexp.mode == "FULL"){
      
    coexp.initDbGTEX()    
    #Init exonic RNA-seq nets
    coexp.nets$exonic <<- NULL
    coexp.nets$exonic$SNIG <<- paste0(gdp.coexp(),"/supplementary/rdsnets/exonic/netSNIG.7.12.it.30.rds")
    coexp.nets$exonic$PUTM <<- paste0(gdp.coexp(),"/supplementary/rdsnets/exonic/netPUTM.8.13.it.30.rds")
    coexp.data$exonic <<- NULL
    coexp.data$exonic$SNIG <<- paste0(gdp.coexp(),"/supplementary/rdsnets/exonic/resids.SNIG.7.rds")
    coexp.data$exonic$PUTM <<- paste0(gdp.coexp(),"/supplementary/rdsnets/exonic/resids.PUTM.8.rds")
    coexp.go$exonic <<- NULL
    coexp.go$exonic$SNIG <<- paste0(gdp.coexp(),"/supplementary/rdsnets/exonic/netSNIG.7.12.it.30.rds_gprofiler.csv")
    coexp.go$exonic$PUTM <<- paste0(gdp.coexp(),"/supplementary/rdsnets/exonic/netPUTM.8.13.it.30.rds_gprofiler.csv")
    coexp.ctype$exonic <<- NULL
    coexp.ctype$exonic$SNIG <<- paste0(gdp.coexp(),"/supplementary/rdsnets/exonic/netSNIG.7.12.it.30.rds.cell.type.csv")
    coexp.ctype$exonic$PUTM <<- paste0(gdp.coexp(),"/supplementary/rdsnets/exonic/netPUTM.8.13.it.30.rds.cell.type.csv")
    
    #Init microarray nets
    nets = coexp.getMicTissues()
    coexp.nets[["micro19K"]] <<- NULL
    coexp.data[["micro19K"]] <<- NULL
    
    for(net in nets){
      #cat(paste0("Loading UKBEC microarray ",net,"\n"))
      coexp.nets$micro19K[[net]] <<- paste0(gdp.coexp(),"/supplementary/rdsnets/micro19K/net",net,
                                          ".12.signed.it.20.rds_cor_pca.rds")
      stopifnot(file.exists(coexp.nets$micro19K[[net]]))
      coexp.data$micro19K[[net]] <<- paste0(gdp.coexp(),"/supplementary/rdsnets/micro19K/",net,
                                          ".mic.expr.data.19K.rds")
      stopifnot(file.exists(coexp.data$micro19K[[net]]))
      coexp.go$micro19K[[net]] <<- paste0(gdp.coexp(),"/supplementary/rdsnets/micro19K/net",net,
                                        ".12.signed.it.20.rds_cor_pca.rds_gprofiler.csv")
      stopifnot(file.exists(coexp.go$micro19K[[net]]))
      coexp.ctype$micro19K[[net]] <<- paste0(gdp.coexp(),"/supplementary/rdsnets/micro19K/net",net,
                                           ".12.signed.it.20.rds",".USER_terms.csv")
      stopifnot(file.exists(coexp.ctype$micro19K[[net]]))
    }
  }
  #Init ROS/MAP nets
  coexp.nets$rosmap$cogdxad <<- paste0(gdp.coexp(),"/data/rdsnets/rosmap/netad.8.it.50.rds")
  coexp.nets$rosmap$cogdxprobad <<- paste0(gdp.coexp(),"/data/rdsnets/rosmap/netprobad.11.it.50.rds")
  coexp.nets$rosmap$cogdnotxad <<- paste0(gdp.coexp(),"/data/rdsnets/rosmap/netnotad.8.it.50.rds")
  
  coexp.data$rosmap$cogdxad <<- paste0(gdp.coexp(),"/data/rdsnets/rosmap/fpkm.casectrl.qc.qn.combat.covs.svas2.res.cogdx.ad.rds")
  coexp.data$rosmap$cogdxprobad <<- paste0(gdp.coexp(),"/data/rdsnets/rosmap/fpkm.casectrl.qc.qn.combat.covs.svas2.res.cogdx.probad.rds")
  coexp.data$rosmap$cogdnotxad <<- paste0(gdp.coexp(),"/data/rdsnets/rosmap/fpkm.casectrl.qc.qn.combat.covs.svas2.res.cogdx.notad.rds")
  
  coexp.go$rosmap$cogdxad <<- paste0(gdp.coexp(),"/data/rdsnets/rosmap/netad.8.it.50.rds_gprofiler.csv")
  coexp.go$rosmap$cogdxprobad <<- paste0(gdp.coexp(),"/data/rdsnets/rosmap/netprobad.11.it.50.rds_gprofiler.csv")
  coexp.go$rosmap$cogdnotxad <<- paste0(gdp.coexp(),"/data/rdsnets/rosmap/netnotad.8.it.50.rds_gprofiler.csv")
  coexp.ctype$rosmap$cogdxad <<- paste0(gdp.coexp(),"/data/rdsnets/rosmap/netad.8.it.50.rds.celltype.csv")
  coexp.ctype$rosmap$cogdxprobad <<- paste0(gdp.coexp(),"/data/rdsnets/rosmap/netprobad.11.it.50.rds.celltype.csv")
  coexp.ctype$rosmap$cogdnotxad <<- paste0(gdp.coexp(),"/data/rdsnets/rosmap/netnotad.8.it.50.rds.celltype.csv")
  
  
  coexp.nets$test$test <<- paste0(gdp.coexp(),"/data/rdsnets/test/smallNetSNIG.rds")
  coexp.nets$rosmap$all <<- paste0(gdp.coexp(),"/supplementary/rdsnets/rosmap/netROSMAPSingle.6.it.50.rds")
  coexp.data$rosmap$all <<- paste0(gdp.coexp(),"/supplementary/rdsnets/rosmap/fpkm.casectrl.qc.qn.combat.covs.svas2.res.rds")
  coexp.go$rosmap$all <<- paste0(gdp.coexp(),"/supplementary/rdsnets/rosmap/netROSMAPSingle.6.it.50.rds_gprofiler.csv")
  coexp.ctype$rosmap$all <<- paste0(gdp.coexp(),"/supplementary/rdsnets/rosmap/netROSMAPSingle.6.it.50.rds.celltype.csv")
  
  cat(paste0("Networks from ",paste0(names(coexp.nets),collapse=",")," loaded\n"))
  cat(paste0("Expression data from ",paste0(names(coexp.data),collapse=",")," loaded\n"))
}

coexp.getCovariates = function(tissue,which.one,cov){
  if(which.one == "rosmap"){
    expr.data = coexp.getExprDataFromTissue(tissue=tissue,which.one=which.one)
    key = read.csv(paste0(gdp.coexp(),"/data/rdsnets/rosmap/ROSMAP_IDkey.csv"))
    covs = read.csv(paste0(gdp.coexp(),"/data/rdsnets/rosmap/ROSMAP_clinical.csv"))
    ids = rownames(expr.data)
    mask = key$projid[match(ids,key$mrna_id)]
    nonmatchingids = ids[is.na(mask)]
    goodids = NULL
    for(id in nonmatchingids){
      subids = str_split(id,"_")
      recid = paste0(subids[[1]][1],"_",subids[[1]][2])
      goodids = c(goodids,recid)
    }
    mask[is.na(mask)] = key$projid[match(goodids,key$mrna_id)]
    samples = mask
    #samples = rosmap.fromRNAseqID2ProjectID(rownames(expr.data))
    gender = as.factor(covs$msex[match(samples,covs$projid)])
    pmi = as.numeric(covs$pmi[match(samples,covs$projid)])
    braaksc = as.factor(covs$braaksc[match(samples,covs$projid)])
    cogdx = as.factor(covs$cogdx[match(samples,covs$projid)])
    educ = as.numeric(covs$educ[match(samples,covs$projid)])
    ceradsc = as.factor(covs$ceradsc[match(samples,covs$projid)])
    age = as.character(covs$age_death[match(samples,covs$projid)])
    
    #Impute
    pmi[is.na(pmi)] = mean(pmi[!is.na(pmi)])
    age[grep("90\\+",age)] = "90"
    age = as.numeric(age)
    race = as.factor(covs$race[match(samples,covs$projid)])
    
    batch = str_split(rownames(expr.data),"_")
    batch = as.factor(unlist(lapply(batch,function(x){return(x[[3]])})))
    
    toreturn = data.frame(batch,gender,pmi,age,race,braaksc,cogdx,educ,ceradsc)
    
    rownames(toreturn) = rownames(expr.data)
    return(toreturn)
  }else if(which.one == "gtexv6"){
    expr.data = coexp.getExprDataFromTissue(tissue=tissue,which.one=which.one)
    covs = read.table(paste0(gdp.coexp(),"/supplementary/rdsnets/gtexv6/gtexv6covariates.txt"),
                      header=T,
                      sep="\t")
    ids = rownames(expr.data)
    covs$ID = gsub("-","\\.",covs$ID)
    return(covs[match(ids,covs$ID),])
  }else if(which.one == "exonic"){
    expr.data = coexp.getExprDataFromTissue(tissue=tissue,which.one=which.one)
    ids = rownames(expr.data)
    covs = read.csv(paste0(gdp.coexp(),"/supplementary/rdsnets/exonic/all_samples_data.csv"),stringsAsFactors=F)
    covs = covs[match(ids,covs$A.CEL_file),c(2:6,8,9)]
    colnames(covs) = c("Age","PMI","RIN","Gender","tissue","causeofdeath","batch")
    return(covs)
  }else if(which.one == "micro19K"){
    expr.data = coexp.getExprDataFromTissue(tissue=tissue,which.one=which.one)
    ids = rownames(expr.data)
    covs = read.csv(paste0(gdp.coexp(),"/supplementary/rdsnets/exonic/all_samples_data.csv"),stringsAsFactors=F)
    covs$U.SD_No = gsub("/","_",covs$U.SD_No)
    covs = covs[match(ids,covs$U.SD_No),c(2:6,8,9)]
    colnames(covs) = c("Age","PMI","RIN","Gender","tissue","causeofdeath","batch")
    return(covs)
  }
  return(NULL)
}


coexp.getMicTissues = function(){
#	load("/Users/juanbot/Dropbox/kcl/WGCNA/data_no_outliers.RData")
#	tissues = names(dat.clean)
#	rm(dat.clean)
	return(c("CRBL", "FCTX", "HIPP", "MEDU", "OCTX", "PUTM", "SNIG", "TCTX", 
					"THAL", "WHMT"))
}

coexp.getGTExNet = function(tissue){
	
	the.dir = paste0(gdp.coexp(),"/supplementary/rdsnets/gtexv6/")
	files = list.files(the.dir)
	
	net.file = files[grep(paste0("net",tissue,"\\.\\d+\\.it\\.50\\.rds$"),files)]
	if(length(net.file) == 0)
		return(NULL)
	
	return(paste0(the.dir,net.file))	
}

coexp.getGTExTissues = function(){
	
	the.dir = paste0(gdp.coexp(),"/supplementary/rdsnets/gtexv6/")
	files = list.files(the.dir)
	
	net.files = files[grep("net\\w+\\.\\d+\\.it\\.50\\.rds$",files)]
	net.files = gsub("net","",net.files)
	net.files = gsub(".\\d+\\.it\\.50\\.rds$","",net.files)
	
	return(net.files)	
}


coexp.getInfoFromTissue = function(which.one,tissue,what,...){
  if(what == "net")
    return(coexp.getNetworkFromTissue(which.one=which.one,tissue=tissue))
  if(what == "expr")
    return(coexp.getExprDataFromTissue(which.one=which.one,tissue=tissue))
  if(what == "go")
    return(coexp.getGOFromTissue(which.one=which.one,tissue=tissue))
  if(what == "ctype")
    return(coexp.getCellTypeFromTissue(which.one=which.one,tissue=tissue))
  
}

coexp.getNetworkFromTissue = function(tissue="SNIG",which.one="exonic",only.file=F,genes=NULL,modules){

  if(which.one == "new"){
    cat("New network, reading from",tissue,"\n")
    return(readRDS(tissue))
    
  }
    
	if(which.one == "onthefly"){
		stopifnot(!is.null(genes))
		net = NULL
		net$moduleColors = modules
		names(net$moduleColors) = genes
		net$MEs = NULL
		return(net)
	}
	
	file = coexp.nets[[which.one]][[tissue]]
	if(only.file)
		return(file)
	else{
		if(file.exists(file))
			return(readRDS(file))
	}
	return(NULL)
				
	cat(paste0("When getting the network for tissue ",tissue," we don't know ",which.one,"\nReturning NULL\n"))
	return(NULL)
}

coexp.getNetworkEigengenes = function(tissue,which.one){
  mes = coexp.getNetworkFromTissue(tissue=tissue,which.one=which.one)$MEs
  class(mes) = c("data.frame","meigengenes") 
  return(mes)
}

coexp.getModulesMostRelevantGenes = function(tissues="SNIG",
                                            modules="red",which.ones,
                                            n=10,cutoff=-1){
  toreturn = NULL
  for(tissue in tissues){
    i = which(tissues == tissue)
    expr.data.file = coexp.getExprDataFromTissue(tissue=tissue,
                                                 which.one=which.ones[i],only.file=T)
    mm = coexp.getMM(which.one=which.ones[i],tissue=tissue,
                     table.format = T,genes=NULL,
                     expr.data.file=expr.data.file)
    
    mm = mm[mm$module == modules[i],]
    mm = mm[order(mm$mm,decreasing=TRUE),]
    ncount = 0
    if(n > 0){
      mask = 1:n
      ncount = n
    }
    else{
      mask = mm$mm > cutoff
      ncount = sum(mask)
    } 
    toreturn = rbind(toreturn,c(rep(tissue,ncount),rep(modules[i],ncount),rep(which.ones[i],ncount),
                                names(mm$mm[mask]),mm$mm[mask]))
      
  }
  toreturn
}

coexp.getModuleMostRelevantGenes = function(tissue="SNIG",
		module="red",which.one,n=10,cutoff=-1,expr.data.file=NULL){
  
	mm = coexp.getMM(which.one=which.one,tissue=tissue,table.format = T,genes=NULL,expr.data.file=expr.data.file)
	
	mm = mm[mm$module == module,]
	mm = mm[order(mm$mm,decreasing=TRUE),]
	if(n > 0)
		return(mm[1:n,])
	return(mm[mm$mm > cutoff,])
}

coexp.getExprDataFromTissue = function(tissue="SNIG",which.one="rnaseq",only.file=F){
	
  if(which.one == "new"){
    cat("Reading new expression data from",tissue,"\n")
    if(only.file)
      return(tissue)
    return(readRDS(tissue))
  }
    
	file = coexp.data[[which.one]][[tissue]]
	if(only.file)
		return(file)
	else{
		if(file.exists(file)){
		  if(which.one == "rosmap"){
		    expr.data = data.frame(readRDS(file))
		    colnames(expr.data) = gsub("\\.[0-9]+","",colnames(expr.data))
		    return(expr.data)
		  }
		  return(readRDS(file))
		}
			
	}
	return(NULL)
	
	cat(paste0("When getting the expression data for tissue ",tissue," we don't know ",which.one,"\nReturning NULL\n"))
	return(NULL)
}

coexp.getModulesFromTissue = function(tissue="SNIG",which.one="rnaseq",in.cluster=F){
	return(unique(coexp.getNetworkFromTissue(tissue,which.one)$moduleColors))
}

coexp.getGenesFromModule = function(tissue="SNIG",which.one="rnaseq",module="black"){
	net = coexp.getNetworkFromTissue(tissue=tissue,which.one=which.one)
	if(is.null(module)) return(names(net$moduleColors))
	return(names(net$moduleColors)[net$moduleColors == module])
}

coexp.getModuleFromGenes = function(tissue="SNIG",which.one="exonic",genes){
	net = coexp.getNetworkFromTissue(tissue=tissue,which.one=which.one)
	return(net$moduleColors[match(genes,names(net$moduleColors))])
}

coexp.getGOFromTissues = function(tissues,which.ones,modules){
  n = length(modules)
  toreturn = NULL
  for(i in 1:n){
    cat("Working with",which.ones[i],", ",tissues[i]," and ",modules[i],"\n")
    localtr = coexp.getGOFromTissue(tissue=tissues[i],
                                    which.one=which.ones[i],
                                    module=modules[i])
    key = paste0(which.ones[i],"_",tissues[i],"_",modules[i])
    toreturn = rbind(toreturn,cbind(rep(key,nrow(localtr)),localtr))
  }
  colnames(toreturn)[1] = "key"
  return(toreturn)
}


coexp.getGOFromTissue = function(tissue="SNIG",which.one="rnaseq",module=NULL){
  if(file.exists(coexp.go[[which.one]][[tissue]])){
    go = read.csv(coexp.go[[which.one]][[tissue]],stringsAsFactors=F)
    if(!is.null(module))
      return(go[go$query.number == module,])
    return(go)
  }
  cat(paste0("When getting the network GO for tissue ",tissue," we don't know ",
             which.one,"\nReturning NULL\n"))
  return(NULL)
}

coexp.getCellTypeFromTissue = function(tissue="SNIG",which.one="rnaseq",module=NULL){
  if(file.exists(coexp.ctype[[which.one]][[tissue]])){
    ct = data.frame(read.csv(coexp.ctype[[which.one]][[tissue]],stringsAsFactors=F))
    if(!is.null(module)){
      if(which.one == "micro19K"){
        ct = ct[ct$InputCategories == module,]
        ctvec = ct[,4]
        names(ctvec) = ct[,2]
        
      }else{
        ctvec = ct[,module]
        names(ctvec) = ct[,1]
        
      }
      return(ctvec)
    }
      
    else return(ct)
  }
    
  cat(paste0("When getting the network cell type signals for tissue ",tissue,
             " we don't know ",which.one,"\nReturning NULL\n"))
  return(NULL)
}

