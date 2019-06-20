# TODO: Add comment
# 
# Author: juanbot
###############################################################################

library(preprocessCore)
require(flashClust)
library(limma)
library(stringr)
library(swamp)
library(sva)
library(outliers)

source("coexpinit.R")

coexp.trasposeDataFrame = function(file.in,first.c.is.name=T){
	if(typeof(file.in) == "character")
		data.in = readRDS(file.in)
	else{
		data.in = file.in
		rm(file.in)
	}
	
	if(first.c.is.name){
		data.t = as.data.frame(cbind(apply(data.in[,-1],MARGIN=1,function(x){ return(as.numeric(x))})))
		colnames(data.t) = data.in[,1]
		rownames(data.t) = colnames(data.in[-1])
	}else{
		data.t = as.data.frame(cbind(apply(data.in,MARGIN=1,function(x){ return(as.numeric(x))})))
		colnames(data.t) = rownames(data.in)
		rownames(data.t) = colnames(data.in)	
	}
	return(data.t)
}

coexp.rosmapRepairSamples = function(samples,reload=F){
	
	if(reload){
		rpkms <<- read.table(paste0(gdp.rosmap(),"/ROSMAP_RNAseq_FPKM_gene.tsv"),
				header=T,row.names=1)
		colnames(rpkms) <<- gsub("X","",colnames(rpkms))
		rpkms <<- rpkms[,-1]
		cat("Reading raw expression data\n")
		print(rpkms[1:5,1:5])
	}
	
	allids = colnames(rpkms)
	allprjids = coexp.rosmapFromRNAseqID2ProjectID(allids)
	#Now get covariates for cases and ctrls
	covs = read.csv(paste0(gdp.rosmap(),"/ROSMAP_clinical.csv"))
	plates_78 = allids[is.na(match(allids,colnames(rpkms)))]
	plates_78[paste0(plates_78,"_8") %in% colnames(rpkms)] = 
			paste0(plates_78[paste0(plates_78,"_8") %in% colnames(rpkms)],"_8")
	plates_78[paste0(plates_78,"_7") %in% colnames(rpkms)] = 
			paste0(plates_78[paste0(plates_78,"_7") %in% colnames(rpkms)],"_7")
	allids[is.na(match(allids,colnames(rpkms)))] = plates_78
	return(allids)
}

coexp.rosmapFromRNAseqID2ProjectID = function(ids){
	key = read.csv(paste0(gdp.rosmap(),"/ROSMAP_IDkey.csv"))
	covs = read.csv(paste0(gdp.rosmap(),"/ROSMAP_clinical.csv"))
	mask = key$projid[match(ids,key$mrna_id)]
	nonmatchingids = ids[is.na(mask)]
	goodids = NULL
	for(id in nonmatchingids){
		subids = str_split(id,"_")
		recid = paste0(subids[[1]][1],"_",subids[[1]][2])
		goodids = c(goodids,recid)
	}
	mask[is.na(mask)] = key$projid[match(goodids,key$mrna_id)]
	return(mask)
	
}

coexp.rosmapPlotMDS = function(rpkms.net,path,covvars,label,n.mds=-1){
	covs = coexp.rosmapGetCovariates(ids=colnames(rpkms.net))
	intersect.s = intersect(rownames(covs),colnames(rpkms.net))
	covs = covs[intersect.s,]
	rpkms.net = rpkms.net[,intersect.s]
	stopifnot(identical(colnames(rpkms.net),rownames(covs)))
	
	for(covvar in covvars){
		cat(paste0("Working on MDS plot for ",covvar,"\n"))
		pdf(paste0(path,covvar,".pdf"),height=8,width=30)
		if(n.mds > 0)
			mask = sample(ncol(rpkms.net),n.mds)
		else
			mask = c(1:ncol(rpkms.net))
		colors = rainbow(length(levels(covs[,covvar])))
		plotMDS(rpkms.net[,mask],col=colors[as.numeric(covs[,covvar])],
				main=paste0("ROS/MAP MDS using ",covvar," ",label))
		legend("topright",fill=colors,
				legend=levels(covs[,covvar]))
		dev.off()
	}
	
}

coexp.rosmapGetCovariates = function(ids=NULL){
	
	if(is.null(ids)){
		expr.data = readRDS(paste0(gdp.rosmap(),"/fpkm.qc.qn.combat.covs.svas2.res.rds"))
		samples = coexp.rosmapFromRNAseqID2ProjectID(rownames(expr.data))
		
	}else
		samples = coexp.rosmapFromRNAseqID2ProjectID(ids)
	#samples = unique(samples)
	covs = read.csv(paste0(gdp.rosmap(),"/ROSMAP_clinical.csv"),stringsAsFactors=F)
	gender = as.factor(covs$msex[match(samples,covs$projid)])
	pmi = as.numeric(covs$pmi[match(samples,covs$projid)])
	braaksc = as.factor(covs$braaksc[match(samples,covs$projid)])
	cogdx = as.factor(covs$cogdx[match(samples,covs$projid)])
	educ = as.numeric(covs$educ[match(samples,covs$projid)])
	ceradsc = as.factor(covs$ceradsc[match(samples,covs$projid)])
	age = covs$age_death[match(samples,covs$projid)]
	
	#Impute
	pmi[is.na(pmi)] = mean(pmi[!is.na(pmi)])
	age[grep("90\\+",age)] = "90"
	age = as.numeric(age)
	race = as.factor(covs$race[match(samples,covs$projid)])
	
	if(is.null(ids))
		batch = str_split(rownames(rownames(expr.data)),"_")
	else
		batch = str_split(ids,"_")
	batch = as.factor(unlist(lapply(batch,function(x){return(x[[3]])})))
	
	toreturn = data.frame(batch,gender,pmi,age,race,braaksc,cogdx,educ,ceradsc)
	
	if(is.null(ids))
		rownames(toreturn) = rownames(expr.data)
	else
		rownames(toreturn) = ids
	
	return(toreturn)
}

coexp.rosmapFromFPKM2Residuals = function(reload=F,
		max.outliers=1,
		drop.sva.axes=NULL,
		covnames = c("batch","braaksc","cogdx","gender","race","ceradsc"),
		n.mds=100,
		do.plot=T){
	
	if(reload){
		rpkms <<- read.table(paste0(gdp.rosmap(),"/ROSMAP_RNAseq_FPKM_gene.tsv"),
				header=T,row.names=1)
		colnames(rpkms) <<- gsub("X","",colnames(rpkms))
		rpkms <<- rpkms[,-1]
		cat("Reading raw expression data\n")
		print(rpkms[1:5,1:5])
	}
	
	#Filter as in DE, getting same gene set for both datasets
	expressed.genes = rowSums(rpkms > 0.1)  > (0.8 * ncol(rpkms))
	rpkms.net = rpkms[expressed.genes,]
	cat("After filtering we keep",nrow(rpkms.net),"genes\n")
	
	colnames(rpkms.net) = coexp.rosmapRepairSamples(reload=F)
	#Get the data right!
	cat("Using",ncol(rpkms.net),"samples for QC\n")
	
	#DEBUG
	#rpkms.net = rpkms.net[1:1000,]
	
	covs = coexp.rosmapGetCovariates(ids=colnames(rpkms.net))
	intersect.s = intersect(rownames(covs),colnames(rpkms.net))
	covs = covs[intersect.s,]
	rpkms.net = rpkms.net[,intersect.s]
	stopifnot(identical(colnames(rpkms.net),rownames(covs)))
	
	if(do.plot){
		coexp.rosmapPlotMDS(rpkms.net=rpkms.net,
				path=paste0(gdp.rosmap(),"/fpkm.mds."),
				covvars=covnames,label="unprocessed",
				n.mds=n.mds)
		
		
		pdf(paste0(gdp.rosmap(),"/fpkm.outliers.pdf"),height=8,width=30)
		plot(flashClust(dist(coexp.trasposeDataFrame(rpkms.net,F)), method = "average"),
				main=paste0("Possible outliers for ROS/MAP cases + ctrls"))
		dev.off()
		
	}
	
	#Lets study outliers
	outlist = NULL
	if(max.outliers > 0){
		sampdist = as.matrix(dist(coexp.trasposeDataFrame(rpkms.net,F)))
		for(i in 1:max.outliers){
			result = grubbs.test(apply(sampdist,1,sum))
			if(result$p.value < 0.05){
				sampdist = sampdist[rownames(sampdist) !=  names(result$p.value),
						colnames(sampdist) !=  names(result$p.value)]
				cat("We should drop out",names(result$p.value),"\n")
				outlist = c(outlist,names(result$p.value))
				pdf(paste0(gdp.rosmap(),"/fpkm.outliers",i,".pdf"),height=8,width=30)
				plot(flashClust(dist(coexp.trasposeDataFrame(rpkms.net[,colnames(rpkms.net) != names(result$p.value)],F)), method = "average"),
						main=paste0("Removed outlier ",names(result$p.value)," for ROS/MAP cases + ctrls"))
				dev.off()
				
			}else
				break
		}
		
		if(do.plot){
			coexp.rosmapPlotMDS(rpkms.net=rpkms.net,
					path=paste0(gdp.rosmap(),"/fpkm.qc.mds."),
					covvars=covnames,label="QCed",
					n.mds=n.mds)
		}	
	}
	
	
	if(!is.null(outlist))
		write.csv(as.data.frame(outlist),paste0(gdp.rosmap(),"/fpkm.outliers.csv"))
	print("We get rid of the following samples")
	print(outlist)
	
	rpkms.net = rpkms.net[,!(colnames(rpkms.net) %in% outlist)]
	print("We're left with")
	print(dim(rpkms.net))
	#samples = samples[!(samples %in% outlist)]
	#Work the covariates
	covs = coexp.rosmapGetCovariates(ids=colnames(rpkms.net))
	intersect.s = intersect(rownames(covs),colnames(rpkms.net))
	covs = covs[intersect.s,]
	rpkms.net = rpkms.net[,intersect.s]
	stopifnot(identical(colnames(rpkms.net),rownames(covs)))
	the.genes = rownames(rpkms.net)
	the.samples = colnames(rpkms.net)
	
	#Normalize
	rpkms.net = normalize.quantiles(as.matrix(rpkms.net))
	#Ready to roll!!		
	colnames(rpkms.net) = the.samples
	rownames(rpkms.net) = the.genes
	
	saveRDS(coexp.trasposeDataFrame(rpkms.net,F),paste0(gdp.rosmap(),"/fpkm.qc.qn.rds"))
	
	debug.mask = rep(T,nrow(rpkms.net))
	#debug.mask[1:10000] = F
	#Lets correct for batch effects
	bch.corrected = ComBat(dat=rpkms.net,batch=covs$batch,
			mod=model.matrix(~1,data=covs),par.prior=T,prior.plots=T)
	#Is this a right approach???
	bch.corrected = bch.corrected - min(bch.corrected)
	out.file = paste0(gdp.rosmap(),"/fpkm.qc.qn.combat.rds")
	cat("Saving batch corrected expression at",out.file,"\n")
	saveRDS(bch.corrected,out.file)
	
	coexp.rosmapPlotMDS(rpkms.net=bch.corrected,
			path=paste0(gdp.rosmap(),"/fpkm.qc.qn.combat.mds."),
			covvars=covnames,label="QCed+Combat",
			n.mds=n.mds)
	
	if(do.plot){
		pdf(paste0(gdp.rosmap(),"/fpkm.qc.qn.pdf"))
		pcres = prince(as.matrix(rpkms.net[debug.mask,]),covs[,colnames(covs) != "diseasestatus"],top=20)
		coexp.princePlot(prince=pcres,main="ROS/MAP All samples, batch uncorrected FPKMQN")
		dev.off()
		pdf(paste0(gdp.rosmap(),"/fpkm.qc.qn.combat.pdf"))
		pcres = prince(as.matrix(bch.corrected),covs[,colnames(covs) != "diseasestatus"],top=20)
		coexp.princePlot(prince=pcres,main="ROS/MAP All samples, batch corrected FPKMQN")
		dev.off()
		
	}
	
	#Estimate the number of surrogate variables
	#Given that samples = c(rosmap.ctrlsID(),rosmap.casesID())
	#cc.design = as.factor(c(rep(1,length(rosmap.ctrlsID())),rep(0,length(rosmap.casesID()))))
	mm = model.matrix(~ gender + pmi + age + race,data=covs)
	nullmm = model.matrix(~ 1,data=covs)
	#print(Sys.time())
	#nsvas = num.sv(rpkms.net,mm,method="leek")
	cat("Launching svaseq\n")
	print(Sys.time())
	svas = svaseq(dat=as.matrix(bch.corrected),mod=mm,mod0=nullmm)
	print(Sys.time())
	out.file = paste0(gdp.rosmap(),"/fpkm.qc.qn.combat.sva.rds")
	cat("Saving surrogate variables at",out.file,"\n")
	saveRDS(svas,out.file)
	
	#Now we control whether the surrogate variables
	#have some relation with the covariates used in the mm model
	#through a heat map of correlation p-values
	numeric.covs = covs[,colnames(covs) != "diseasestatus"]
	linp = matrix(ncol=svas$n.sv,nrow=ncol(numeric.covs))
	rownames(linp) = colnames(numeric.covs)
	colnames(linp) = paste0("SV",1:svas$n.sv)
	linp[] = 0
	for(cov in 1:ncol(numeric.covs)){	
		for(sva in 1:svas$n.sv){
			if(svas$n.sv == 1)
				axis = svas$sv
			else 
				axis = svas$sv[,sva]
			
			linp[cov,sva] = cor.test(as.numeric(numeric.covs[,cov]),axis)$p.value
		}
	}
	smallest = -10
	linp10 <- log10(linp)
	linp10 <- replace(linp10, linp10 <= smallest, smallest)
	tonote <- signif(linp, 1)
	pdf(paste0(gdp.rosmap(),"/fpkm.qc.qn.combat.sva.pdf"))
	
	heatmap.2(linp10, Colv = F, Rowv = F, dendrogram = "none", 
			trace = "none", symbreaks = F, symkey = F, breaks = seq(-20, 0, length.out = 100), 
			key = T, colsep = NULL, rowsep = NULL, sepcolor = "black", 
			sepwidth = c(0.05, 0.05), main = "ROS/MAP, correlation of SVs and covariates, FPKM QN", 
			labCol = colnames(linp10), 
			labRow = colnames(covs), 
			xlab = "Surrogate variables")
	dev.off()
	
	covs.rs = covs[,match(c("gender","pmi","age","race"),colnames(covs))]
	
	##Apply data correction
	print("Getting residuals")
	resids <- apply(bch.corrected, 1, function(y){
				lm( y ~ . , data=cbind(covs.rs,svas$sv))$residuals 
			})
	rownames(resids) = colnames(bch.corrected)
	print("Saving residuals")
	saveRDS(resids,paste0(gdp.rosmap(),"/fpkm.casectrl.qc.qn.combat.sva.res.rds"))
	
	pdf(paste0(gdp.rosmap(),"/fpkm.qc.qn.combat.sva.res.pdf"))
	pcres = prince(as.matrix(coexp.trasposeDataFrame(resids,F)),covs.rs,top=20)
	coexp.princePlot(prince=pcres,main="ROS/MAP All samples, residuals (gender, pmi, age, race)")
	dev.off()
	
	
	coexp.rosmapPlotMDS(rpkms.net=coexp.trasposeDataFrame(resids,F),
			path=paste0(gdp.rosmap(),"/fpkm.qc.qn.combat.sva.res.mds."),
			covvars=covnames,label="residuals",
			n.mds=n.mds)
	print("Done!")
}

coexp.princePlot = function (prince, label = colnames(prince$o), smallest = -20, 
		note = F, notecol = "black", notecex = 1, breaks = seq(-20, 
				0, length.out = 100), col = heat.colors(99), margins = c(5, 
				7), key = T, cexRow = 1, cexCol = 1, xlab = "Principal Components (Variation)", 
		colsep = NULL, rowsep = NULL, sepcolor = "black", sepwidth = c(0.05, 
				0.05), Rsquared = F, breaksRsquared = seq(0, 1, length.out = 100),main) 
{
	if (class(prince) != "prince") {
		stop("prince is not an object generated by the prince function")
	}
	if (smallest > 0) {
		stop("smallest has to be less than 0")
	}
	require(gplots)
	linp10 <- log10(prince$linp)
	linp10 <- replace(linp10, linp10 <= smallest, smallest)
	tonote <- signif(prince$linp, 1)
	prop <- round(prince$prop, 0)
	if (Rsquared == T) {
		linp10 <- prince$rsquared
		breaks <- breaksRsquared
		tonote <- round(prince$rsquared, 2)
		col <- col[length(col):1]
	}
	heatmap.2(linp10, Colv = F, Rowv = F, dendrogram = "none", 
			trace = "none", symbreaks = F, symkey = F, breaks = breaks, 
			key = key, col = col, cexRow = cexRow, cexCol = cexCol, 
			colsep = colsep, rowsep = rowsep, sepcolor = sepcolor, 
			sepwidth = sepwidth, main = main, labCol = paste(1:ncol(linp10), 
					" (", prop, ")", sep = ""), margins = margins, labRow = label, 
			xlab = xlab, cellnote = if (note == T) {
						tonote
					}
					else {
						matrix(ncol = ncol(prince$linp), nrow = nrow(prince$linp))
					}, notecol = notecol, notecex = notecex)
}

