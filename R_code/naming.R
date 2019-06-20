

coexp.biomart.conv.table = paste0(gdp.coexp(),"/data/EnsemblNamesGRCh38.txt")

coexp.fromAny2Ensembl = function(genes){
	if(substr(genes[1],1,1) == "X"){
		return(coexp.fromGeneName2Ensembl(coexp.fromXtIDToGeneSymbols19K(genes)))
	}else if(substr(genes[1],1,4) != "ENSG"){
		return(coexp.fromGeneName2Ensembl(genes))
	}
	return(genes)
}


coexp.fromEntrez2Ensembl <- function(genes,
		table.file=biotype.all.file,
		ignore.unknown=FALSE,
		which.are.unknown=FALSE){
	
	#table.file=biotype.file,ignore.unknown=FALSE){
	#table.file=biotype.all.file,ignore.unknown=FALSE){
	
	function.table <- read.table(table.file,header=TRUE,stringsAsFactors=FALSE,sep=",")
	gene.names <- function.table$ensembl_gene_id[match(genes,function.table$external_gene_id)]
	
	if(which.are.unknown)
		return(is.na(gene.names))
	
	if(ignore.unknown){
		gene.names <- na.omit(gene.names)
	}else{
		unnamed.genes.count <- sum(is.na(gene.names))
		if(unnamed.genes.count > 0){
			print(paste0("Genes without name ",unnamed.genes.count))
			if(unnamed.genes.count > 5)
				print(paste0("Genes ",paste0(genes[is.na(gene.names)][1:5],
										collapse=", "), " and more... don't have a name"))
			else
				print(paste0("Genes ",paste0(genes[is.na(gene.names)],collapse=", "), " don't have a name"))
			gene.names[is.na(gene.names)] <- genes[is.na(gene.names)]	
		}
		
	}
	
	return(gene.names)
}

coexp.fromGeneName2Ensembl <- function(genes,ignore.unknown=FALSE,which.are.unknown=FALSE){
	if(!exists("coexp.utils.gene.names.conv.table"))
		coexp.loadConvTable()
	
	gene.names <- coexp.utils.gene.names.conv.table$Ensembl[match(genes,coexp.utils.gene.names.conv.table$Gene)]
	
	if(which.are.unknown)
		return(is.na(gene.names))
	
	if(ignore.unknown){
		gene.names <- na.omit(gene.names)
	}else{
		unnamed.genes.count <- sum(is.na(gene.names))
		if(unnamed.genes.count > 0){
			print(paste0("Genes without name ",unnamed.genes.count))
			print(paste0("Genes ",paste0(genes[is.na(gene.names)][1:unnamed.genes.count],collapse=", "), 
							" don't have a name"))
			gene.names[is.na(gene.names)] <- genes[is.na(gene.names)]	
		}
	}
	return(gene.names)
}

coexp.fromAny2GeneName = function(genes){
	if(substr(genes[1],1,1) == "X"){
		return(coexp.fromXtIDToGeneSymbols19K(genes))
	}else if(substr(genes[1],1,4) == "ENSG"){
		return(coexp.fromEnsembl2GeneName(genes))
	}
	return(genes)
}

coexp.fromXtIDToGeneSymbols19K = function(xids){
	
	trans.table <- read.csv(paste0(gdp.coexp(),"/supplementary/rdsnets/micro19K/annot_19K.csv"),stringsAsFactors=F)
	gene.symbols = trans.table$Gene_Symbol[match(xids,trans.table$XtID)]
	gene.symbols[is.na(gene.symbols)] = xids[is.na(gene.symbols)]
	return(gene.symbols)
}

coexp.loadConvTable = function(in.table=coexp.biomart.conv.table){
	the.table = read.table(in.table,header=F,stringsAsFactors=F)
	colnames(the.table) = c("Ensembl","Gene")
	coexp.utils.gene.names.conv.table <<- the.table
}

coexp.fromEnsembl2GeneName <- function(genes,ignore.unknown=FALSE,which.are.unknown=FALSE){
	if(!exists("coexp.utils.gene.names.conv.table"))
		coexp.loadConvTable()
	
	gene.names <- coexp.utils.gene.names.conv.table$Gene[match(genes,coexp.utils.gene.names.conv.table$Ensembl)]
	
	if(which.are.unknown)
		return(is.na(gene.names))
	
	if(ignore.unknown){
		gene.names <- na.omit(gene.names)
	}else{
		unnamed.genes.count <- sum(is.na(gene.names))
		if(unnamed.genes.count > 0){
			cat("Genes without name ",unnamed.genes.count,"\n")
			if(unnamed.genes.count > 3)
				unnamed.genes.count = 10
			cat("Genes ",paste0(genes[is.na(gene.names)][1:unnamed.genes.count],collapse=", "), " and more... don't have a name\n")
			gene.names[is.na(gene.names)] <- genes[is.na(gene.names)]	
		}
	}
	return(gene.names)
}

