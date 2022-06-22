#' Creates a features pair file from a GFF file
#'
#' This function creates a features pair file from a GFF file
#'
#' @param gffFile Path to the input GFF file
#' @param outFile Path to the output features pair file, if not provided the output will be gffFile.pairs
#' @param forceOverwrite If output file exists, overwrite the existing file. (default FALSE)
#' @return Path to the generated features pair file  
#' @keywords internal

get_pairs_from_gff3<-function(gffFile, outFile, forceOverwrite=FALSE){
	if (!base::file.exists(gffFile)) { 
		stop(paste0("Input gff3 file ",gffFile," not found."))
	}

	if(missing(outFile)){
		foutput=paste0(gffFile,".pairs")
	} else {
		foutput=outFile
	}

	if (base::file.exists(foutput) && forceOverwrite==FALSE) { 
		## ToDo: Remove message if outFile is missing
		if(!missing(outFile)){
			message ("Pairs file for this GFF already exists. Use forceOverwrite=TRUE to overwrite the existing file."); 
		}
		return(foutput) 
	}


	cnames <- c("seqname","source","feature","start","end","score","strand","frame","attribute")
	myGFFdata <- utils::read.delim(gffFile, header=F,  sep="\t", comment.char = "#", quote="", blank.lines.skip=T) %>%    `colnames<-`(cnames)  

	myIDPattern<-"^(.*;)?( *ID *= *([^;]+))(.*)$"
	selected_id_lines<-grep(myIDPattern,myGFFdata$attribute,perl=T,ignore.case=T)

	resID <- data.frame(FeatureID=myGFFdata[selected_id_lines,"feature"], ID = gsub(myIDPattern,"\\3",myGFFdata[selected_id_lines,"attribute"],perl=T,ignore.case=T))


	myParentPattern<-"^(.*;)?( *Parent *= *([^;]+))(.*)$"
	selected_parent_lines<-grep(myParentPattern,myGFFdata$attribute,perl=T,ignore.case=T)

	resParent <- data.frame(Feature=myGFFdata[selected_parent_lines,"feature"], ID = gsub(myParentPattern,"\\3",myGFFdata[selected_parent_lines,"attribute"],perl=T,ignore.case=T))

	## INCLUDE MULTIPLE PARENT OCCURRENCES
	withMulti<-grep(",",resParent$ID)
	if(length(withMulti)>0){
		newItems<-do.call("rbind",lapply(withMulti,function(x){t(sapply(strsplit(resParent$ID[x],",")[[1]],function(y){c(resParent$Feature[x],y)}))}))
		colnames(newItems) <- c("Feature", "ID")
		resParent<-rbind(resParent[-withMulti,],newItems)
	}

	mergedDF<-merge(x = resParent, y = resID, by = "ID", all.x = TRUE)

	countPairs<-mergedDF %>% dplyr::count(.data$Feature,.data$FeatureID)

	## Robust to missing parent IDS 
	if(any(is.na(countPairs$FeatureID))){
		message("Warning: there are missing Parent IDs and will be ignored")
		# Option 1: (active) For features with missing parents, parent feature changed from NA to empty string in pairs
		countPairs[is.na(countPairs$FeatureID),]$FeatureID <- ""
		# Option 2: (disabled) Features with missing parents are removed from pairs file
		#countPairs<-countPairs[!is.na(countPairs$FeatureID),]
	}
	countElems<-data.frame(ELEMENT=myGFFdata[-selected_parent_lines,"feature"], PARENT=rep("",nrow(myGFFdata)-length(selected_parent_lines))) %>% dplyr::count(.data$ELEMENT, .data$PARENT)

	colnames(countPairs)<-c("ELEMENT","PARENT", "N")
	colnames(countElems)<-c("ELEMENT","PARENT", "N")

	DF<-rbind(countPairs,countElems)

	utils::write.table(DF,file=foutput,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t", append=FALSE)

	rm(list=setdiff(ls(), "foutput"))	
	gc(reset=T)

	return(foutput)
}

#' Creates a features path file from a GFF file
#'
#' This function creates a features path file from a GFF file
#'
#' @param gffFile Path to the input GFF file
#' @param outFile Path to the output features path file, if not provided the output will be gffFile.paths
#' @param forceOverwrite If output file exists, overwrite the existing file. (default FALSE)
#' @return Path to the generated features path file
#'
#' @importFrom rlang .data
#' @keywords internal
gff3_to_paths<-function(gffFile, outFile, forceOverwrite=FALSE){ 
	if (!base::file.exists(gffFile)) { 
		stop(paste0("Input gff3 file ",gffFile," not found."))
	}

	if(missing(outFile)){
		foutput=paste0(gffFile,".paths")
	} else {
		foutput=outFile
	}

	if (base::file.exists(foutput) && forceOverwrite==FALSE) {
		if(!missing(outFile)){
			message ("Paths file for this GFF already exists. Use forceOverwrite=TRUE to overwrite the existing file."); 
		}
	}
	myGFFdata<-readLines(gffFile)
	message(paste0("Creating ",foutput," file for ",gffFile))
	CommentPattern<-"(^#(.*)$)|(^#([\t ]*)$)|(^[\t ]*$)"
	comment_lines<-grep(CommentPattern,myGFFdata,perl=T)
	if(length(comment_lines)>0){
		myGFFdata<-myGFFdata[-comment_lines]
	}

	ParentAndIDPattern<-"^([^\t]+)\t[^\t]+\t(([^\t]+\t){3})[^\t]+\t([^\t])+\t(.*)[\t;]ID *= *([^;]+)(.*);Parent *= *([^;]+)(.*)$"
	parent_and_id_lines<-grep(ParentAndIDPattern,myGFFdata,perl=T,ignore.case=T)
	parentIdPairs<-gsub(ParentAndIDPattern,"\\1\t\\2\\4\t\\6\t\\8",myGFFdata[parent_and_id_lines],perl=T,ignore.case=T)

	IDAndParentPattern<-"^([^\t]+)\t[^\t]+\t(([^\t]+\t){3})[^\t]+\t([^\t])+\t(.*)[\t;]Parent *= *([^;]+)(.*);ID *= *([^;]+)(.*)$"
	id_and_parent_lines<-grep(IDAndParentPattern,myGFFdata,perl=T,ignore.case=T)
	if(length(id_and_parent_lines>0)){
		IdparentPairs<-gsub(IDAndParentPattern,"\\1\t\\2\\4\t\\8\t\\6",myGFFdata[id_and_parent_lines],perl=T,ignore.case=T)
		parentIdPairs<-c(parentIdPairs,IdparentPairs)
		parent_and_id_lines<-c(parent_and_id_lines,id_and_parent_lines)
	}

	if(length(parentIdPairs)==0){
		resParentID<-data.frame(Chr=character(),Feature=character(),Start=character(),End=character(),Strand=character(), ID=character(),ParentID=character())
	} else {
		vParentIdPairs<-unlist(strsplit(parentIdPairs,"\t"))
		resParentID <- cbind.data.frame(split(vParentIdPairs, rep(1:7, times=length(vParentIdPairs)/7)), stringsAsFactors=F)
		colnames(resParentID) <- c("Chr","Feature","Start","End","Strand", "ID","ParentID")
		## CHECK FOR MULTIPLE PARENT OCCURRENCES
		withMulti<-grep(",",resParentID$ParentID)
		if(length(withMulti)>0){
			newItems<- resParentID[withMulti,] %>% dplyr::mutate(ParentID = strsplit(as.character(.data$ParentID), ",")) %>% tidyr::unnest(.data$ParentID)
			resParentID<-rbind(resParentID[-withMulti,],newItems)
		}
		resParentID<-resParentID[order(resParentID$ParentID),]
	}

	if(length(parent_and_id_lines)>0){
		GFFdata<-myGFFdata[-parent_and_id_lines]
	} else {
		GFFdata<-myGFFdata
	}
	ParentPattern<-"^([^\t]+)\t[^\t]+\t(([^\t]+\t){3})[^\t]+\t([^\t])+\t(.*)[\t;]Parent *= *([^;]+)(.*)$"

	parent_lines<-grep(ParentPattern,GFFdata,perl=T,ignore.case=T)
	parentPairs<-gsub(ParentPattern,"\\1\t\\2\\4\t\\6",GFFdata[parent_lines],perl=T,ignore.case=T)
	if(length(parentPairs)==0){
		resParent<-data.frame(Chr=character(),Feature=character(),Start=character(),End=character(),Strand=character(), ParentID=character(),ID=character())
	} else {
		vParentPairs<-unlist(strsplit(parentPairs,"\t"))
		resParent <- cbind.data.frame(split(vParentPairs, rep(1:6, times=length(vParentPairs)/6)), stringsAsFactors=F)
		colnames(resParent) <- c("Chr","Feature","Start","End","Strand", "ParentID")
		## CHECK FOR MULTIPLE PARENT OCCURRENCES
		withMulti<-grep(",",resParent$ParentID)
		if(length(withMulti)>0){
			newItems<-do.call("rbind",lapply(withMulti,function(x){as.data.frame(t(sapply(strsplit(as.character(resParent$ParentID[x]),",")[[1]],function(y){as.character(c(resParent[x,-6],y))})))}))
			colnames(newItems) <- c("Chr","Feature","Start","End","Strand", "ParentID")
			resParent<-rbind(resParent[-withMulti,],newItems)
		}
		resParent<-resParent[order(resParent$ParentID),]

		resParent$ID<-paste0(resParent$Feature,":",rownames(resParent))
	}
	if(length(parent_lines)>0){
		GFFdata<-GFFdata[-parent_lines]
	}

	idPattern<-"^([^\t]+)\t[^\t]+\t(([^\t]+\t){3})[^\t]+\t([^\t])+\t(.*)[\t;]ID *= *([^;]+)(.*)$"
	selected_lines<-grep(idPattern,GFFdata,perl=T,ignore.case=T)
	idPairs<-gsub(idPattern,"\\1\t\\2\\4\t\\6",GFFdata[selected_lines],perl=T,ignore.case=T)
	if(length(idPairs)==0){
		resID<-data.frame(Chr=character(),Feature=character(),Start=character(),End=character(),Strand=character(), ID=character(),ParentID=character())
	} else {
		vIdPairs<-unlist(strsplit(idPairs,"\t"))
		resID <- cbind.data.frame(split(vIdPairs, rep(1:6, times=length(vIdPairs)/6)), stringsAsFactors=F)
		colnames(resID) <- c("Chr","Feature","Start","End","Strand", "ID")

		resID$ParentID<-rep("",nrow(resID))
	}
	if(length(selected_lines)>0){
		if(length(GFFdata[-selected_lines])>0){
			message(paste0("There are ",length(GFFdata[-selected_lines])," remaining elements in the GFF file"))
		}
	} else {
		if(length(GFFdata)>0){
			message(paste0("There are ",length(GFFdata)," remaining elements in the GFF file"))
		}
	}
	finalColumns<-c("Chr","Feature","Start","End","Strand", "ID","ParentID")


	fullData<-rbind(resParentID[,c("Feature","ID","ParentID")],resParent[,c("Feature","ID","ParentID")],resID[,c("Feature","ID","ParentID")])

	mergedDF<-rbind(resParentID[,finalColumns],resParent[,finalColumns],resID[,finalColumns])
	colnames(mergedDF)<-c(finalColumns)

	counter<-""
	while(any(!is.na(mergedDF[,paste0("ParentID",counter)]) & (mergedDF[,paste0("ParentID",counter)]!=""))){
		if(counter==""){
			nextCounter<-1
		} else{
			nextCounter<-counter+1
		}

		# print(paste0("Checking level ",nextCounter))
		mergedDF<-merge(x = mergedDF, y = fullData, by.x = paste0("ParentID",counter), by.y = "ID", all.x = TRUE,suffixes=c("",as.character(nextCounter)))
		counter<-nextCounter
		repeatedIDs<-which(!is.na(mergedDF[,paste0("ParentID",counter)]) & (mergedDF[,paste0("ParentID",counter)]!="") & (mergedDF[,paste0("ParentID",counter)]==mergedDF[,"ParentID"]))
		if(length(repeatedIDs) >0){
			warning(paste0(length(repeatedIDs), " duplicated ID found in the input GFF file. e.g. ",mergedDF[repeatedIDs[1],"ParentID"]))
			break;
		}
	}
	colOrder<-c(which(colnames(mergedDF)=="Chr"),which(colnames(mergedDF)=="Start"),which(colnames(mergedDF)=="End"),which(colnames(mergedDF)=="Strand"),which(colnames(mergedDF)=="Feature"),which(colnames(mergedDF)=="ID"))
	for (i in c(1:(counter))){
		if(i==1){
			p<-""
		} else {
			p<-i-1
		}
		colOrder<-c(colOrder,which(colnames(mergedDF)==paste0("Feature",i)))
		colOrder<-c(colOrder,which(colnames(mergedDF)==paste0("ParentID",p)))
	}
	mergedDF<-mergedDF[,colOrder]

	# SORT FEATURES 
	testTree<-get_features(gffFile, outFormat="tree",fileType = "GFF3")

	featureList<-rev(unique(rev((sapply(data.tree::Traverse(testTree,traversal="level"), function(node){{node$name}})[-1]))))

	featureOrder<-factor(featureList, levels=featureList)
	
	mergedDF<- mergedDF %>% dplyr::arrange(.data$Chr, as.numeric(.data$Start), dplyr::desc(as.numeric(.data$End)), factor(.data$Feature, levels = featureList))
	
	
	pathList<-as.data.frame(t(apply(mergedDF,1,function(x){c(x["Chr"],x["Start"],x["End"],x["Strand"],paste(x[5:length(x)][!is.na(x[5:length(x)]) & x[5:length(x)]!=""],collapse=";"))})))
	message("Copying paths file") 
	
	utils::write.table(pathList,file=foutput,col.names=FALSE,quote=FALSE,sep="\t", row.names=FALSE,append=FALSE)

	rm(list=setdiff(ls(), "foutput"))	
	gc(reset=T)

	return(foutput)
}




#' Generates SAF formatted data from a features path file for the given blocks and features
#'
#' This function creates a SAF file from a features path file for the given blocks 
#' and features
#'
#' @param pathsFile Path to the input features path file 
#' @param groupBy Vector of features to group by in feature count
#' @param block Vector of features to be used as block in feature count
#' @return SAF formatted data frame for the given blocks and features
#' @keywords internal
saf_from_paths<-function(pathsFile,groupBy=c("mRNA","gene"),block=c("exon","CDS")) {

	message("Loading File")
	thelines<-as.vector(readLines(pathsFile))
	DF<-data.frame("GeneID"=character(0),"Chr"=character(0),"Start"=numeric(0),"End"=numeric(0),"Strand"=character(0),"Notes"=character(0), stringsAsFactors=FALSE)
	message("Getting SAF")
	if(length(block)==0){
		searchLine=paste(paste0("\t",groupBy,";([^;]+)$"),collapse="|")
		selected_lines<-thelines[which(grepl(searchLine,thelines,perl=TRUE))]
		THE_RESULT<-sapply(selected_lines,function(x) {if(is.null(x)) { "FEATURENOTFOUND" } else { x }})
		for (g in seq(from=1,to=length(groupBy),by=1)) {
			search_this_group=paste0("^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(",groupBy[g],");([^;]+)")
			mz<-gsub(search_this_group,"\\6\t\\1\t\\2\t\\3\t\\4\t\\5",THE_RESULT,perl=TRUE)
			DF<-rbind(DF,do.call(rbind.data.frame,strsplit(mz,'\t')),stringsAsFactors=FALSE)
			colnames(DF)<-c("GeneID","Chr","Start","End","Strand","Notes")
		}

	} else{

		searchBlock=paste(paste0("\t",block,";"),collapse="|")
		searchGroup=paste(paste0("^([^;]+).*[;\t]",groupBy,";([^;]+).+"),collapse="|")
		# message("The search: ",searchBlock)
		# message("The search2: ",searchGroup)
		selected_lines<-thelines[which(grepl(searchBlock,thelines,perl=TRUE))]
		selected_lines<-selected_lines[which(grepl(searchGroup,selected_lines,perl=TRUE))]
		THE_RESULT<-lapply(selected_lines,function(x) {if(is.null(x)) { "FEATURENOTFOUND" } else { x }})
		THE_RESULT<-grep("FEATURENOTFOUND",THE_RESULT,fixed=TRUE,value=TRUE,invert=TRUE)
		for (g in seq(from=1,to=length(groupBy),by=1)) {
			search_this_group=paste0("^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t;]+);.*;(",groupBy[g],");([^\t;]+).*")
			mz<-gsub(search_this_group,"\\7\t\\1\t\\2\t\\3\t\\4\t\\5->\\6",THE_RESULT,perl=TRUE)
			mz<-mz[grep("[^\\t]+\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+\\t[^\\t]+",mz,perl=TRUE)]
			m<-paste(mz,collapse="\n")
			m2<-utils::read.table(text=m,sep="\t",quote="")
			if (!is.null(m2[1,6])) {
				names(m2)<-names(DF)
				DF<-rbind(DF,m2,stringsAsFactors=FALSE)
			}
			else {
				message(groupBy[g]," not found")
			}
		}
	}
	rm(list=setdiff(ls(), "DF"))	
	gc(reset=T)
	
	return (DF[order(DF$Chr,DF$GeneID,DF$Start),])
}


#' Summarizes the number of elements of each type in each chromosome of a gff file
#'
#' This function summarizes the number of features of each type in each chromosome  
#' of a gff file and returns the statistics
#'
#' @param inFile Path to the input gff file
#' @return A tibble with the summary data
#'
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' test_gff3<-system.file("extdata", "AthSmall.gff3", package="Rgff")
#' gff_stats_by_chr(test_gff3)

gff_stats_by_chr <-function(inFile) {
	nFields<-max(utils::count.fields(inFile,sep="\t", comment.char = "#",quote=""))	

	if(nFields > 9){
		message("Warning: input file contains lines with more than 9 fields")
	}
	
	gffTable<-utils::read.table(inFile,sep="\t", comment.char = "#", colClasses=c("character","NULL","character","numeric","numeric","NULL", "NULL","NULL","NULL"), stringsAsFactors=FALSE,quote="")
	colnames(gffTable)<-c("Chr", "FeatureType","Start","End")

	summaryData<-gffTable %>%
	dplyr::group_by(.data$Chr,.data$FeatureType) %>%
	dplyr::summarise(AvgLen = mean(.data$End-.data$Start), MaxLen = max(.data$End-.data$Start), MinLen = min(.data$End-.data$Start), n=dplyr::n())	

	rm(list=setdiff(ls(), "summaryData"))	
	gc(reset=T)

	return (summaryData)
}


#' Summarizes the number of features of each type in a gff file
#'
#' This function summarizes the number of features of each type in  
#' a gff file and returns the statistics
#'
#' @param inFile Path to the input gff file
#' @return A tibble with the summary data
#'
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' test_gff3<-system.file("extdata", "AthSmall.gff3", package="Rgff")
#' gff_stats(test_gff3)
gff_stats <-function(inFile) {
	nFields<-max(utils::count.fields(inFile,sep="\t", comment.char = "#",quote=""))	

	if(nFields > 9){
		message("Warning: input file contains lines with more than 9 fields")
	}
	
	gffTable<-utils::read.table(inFile,sep="\t", comment.char = "#", colClasses=c("character","NULL","character","numeric","numeric","NULL", "NULL","NULL","NULL"), stringsAsFactors=FALSE,quote="")
	colnames(gffTable)<-c("Chr", "FeatureType","Start","End")

	summaryData<-gffTable %>%
	dplyr::group_by(.data$FeatureType) %>%
	dplyr::summarise(AvgLen = mean(.data$End-.data$Start), MaxLen = max(.data$End-.data$Start), MinLen = min(.data$End-.data$Start), n=dplyr::n())	

	if (!("chromosome" %in% summaryData$FeatureType)){
		summaryData<-summaryData %>% tibble::add_row(FeatureType="chromosome",n=length(unique(gffTable$Chr)),.before=1)
	}

	rm(list=setdiff(ls(), "summaryData"))	
	gc(reset=T)
	
	return (summaryData)
}



#' Creates a SAF file from a GFF3 for the given pairs of blocks and features,
#' note that groupBy and block vectors must have the same length
#'
#' This function creates a SAF file from a GFF3 for the given blocks 
#' and features
#'
#' @param gffFile Path to the input GFF file
#' @param outFile Path to the output SAF file,  if not provided the output path will be the input path (without extension) with the suffix ".groupBy1-block1.groupBy2-block2(...).saf"
#' @param forceOverwrite If output file exists, overwrite the existing file. (default FALSE)
#' @param features Vector of pairs of features, separated by '>' (see sep argument), to be used respectively as "group by" and "block"
#' @param sep Separator of each "group by" and "block" provided in the feature argument (default '>')
#' @return Path to the generated SAF file
#' @keywords internal
saf_from_gff3<-function(gffFile, outFile, forceOverwrite=FALSE, features=c("gene > exon"), sep=">") {

	escapedSep<-gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", sep)

	# Read features adnd check 
	featureDF<-get_features(gffFile, outFormat="data.frame",fileType = "GFF3")
	featureList<-c()
	blockList<-c()
	for (i in c(1:length(features))){
		groupByAndBlock<-trimws(strsplit(features[i],escapedSep)[[1]])
		groupBy=groupByAndBlock[1]
		block=""
		
		if(length(groupByAndBlock)>2 ){
			stop(paste0("Invalid feature: ",features[i],". Each feature element must be composed by at most two elements with the syntax groupBy ",sep," block."))
		}

		if(!(groupBy %in% (featureDF$FEATURES))){
			stop(paste0("Feature: '",groupBy,"' is not defined in ",gffFile,". \nRun get_features(\"",gffFile,"\") to see the feature tree"))
		}
		
		if(length(groupByAndBlock)>1 ){
			block=groupByAndBlock[2]

			
			if(featureDF[featureDF$FEATURES == groupBy,]$BLOCKS == "" || !(block %in% (strsplit(as.character(featureDF[featureDF$FEATURES == groupBy,]$BLOCKS)," ")[[1]]))){
				stop(paste0("Block '",block, "' is not part of the feature '", groupBy,"' in ",gffFile,". \nRun get_features(\"",gffFile,"\") to see the feature tree"))
			}

			
		} 
		
		featureList<-c(featureList,groupBy)
		blockList<-c(blockList,block)		
	}

	if(missing(outFile)){
		fSuffix<-""
		for (i in c(1:length(featureList))){
				fSuffix<-paste0(fSuffix,".",featureList[i])
				if(blockList[i] != ""){
					fSuffix<-paste0(fSuffix,"-",blockList[i])
				}
		}
		foutput=paste0(sub(pattern = "(.*)\\..{0,4}$", replacement = "\\1", gffFile),fSuffix,".saf")
	} else {
		foutput=outFile
	}


	if(base::file.exists(foutput) && forceOverwrite==FALSE){
		message("Output SAF file exists. Use forceOverwrite=TRUE to overwrite the existing file.")
		return(foutput)
	}
	
	PATHSfile=paste0(gffFile,".paths")
	if (!base::file.exists(PATHSfile)) {
		message ("Creating PATHS file... please be patient...")
		gff3_to_paths(gffFile)
	}

	

	for (i in c(1:length(featureList))){
		if(blockList[i] == ""){
			blockParam=c()
		} else {
			blockParam=c(blockList[i])
		}
		if(i==1){
			safDF<-unique(saf_from_paths(PATHSfile,groupBy=c(featureList[i]),block=blockParam))
		} else {

			safDF<-rbind(safDF,unique(saf_from_paths(PATHSfile,groupBy=c(featureList[i]),block=blockParam)))
		}
	}	
	safDF<-unique(safDF[order(safDF[,2], safDF[,3], safDF[,4], safDF[,1] ),])

	utils::write.table(safDF,foutput,sep="\t",quote=FALSE,row.names = FALSE)

	rm(list=setdiff(ls(), "foutput"))	
	gc(reset=T)
	
	return(foutput)
}	


#' Plots or exports an image of the feature tree from a gff file
#'
#' This function plots the feature tree from a gff file or, if an output file name is provided, 
#' exports an image of in the desired format ("png", "pdf" or "svg"). 
#' Packages "DiagrammeR", "DiagrammeRsvg" and "rsvg" must be installed to use this function.
#'
#' @param inFile Path to the input gff file 
#' @param outFile Path to the output features image file, if not provided the tree will be plotted
#' @param includeCounts Include number of occurrences of each subfeature 
#' @param fileType Version of the input file (GTF/GFF3). If not provided it will be determined from the file name.
#' @param exportFormat Output image format when it is not possible to deduce it from the extension of outFile ("png", "pdf" or "svg"). Default, "png"
#' @return Path of the output features image file
#'
#' @export
#'
#' @examples
#' test_gff3<-system.file("extdata", "AthSmall.gff3", package="Rgff")
#' plot_features(test_gff3)

plot_features<- function(inFile, outFile, includeCounts=FALSE, fileType=c("AUTO","GFF3","GTF"), exportFormat=c("png","pdf","svg")){
	if (!requireNamespace("DiagrammeR", quietly = TRUE)) {    
		stop("Package \"DiagrammeR\" is required for plotting a `data.tree::Node` object. Please install it.")
	}

	fileType <- match.arg(fileType)
	exportFormat <- match.arg(fileType)

	if (!base::file.exists(inFile)) { 
		stop(paste0("Input gff file ",inFile," not found."))
	}
	detectedExt<-tools::file_ext(inFile)

	if(fileType == "AUTO"){
		if(tolower(detectedExt) == "gtf"){
			fileExt="gtf"
			fileType = "GTF"
		} else if( (tolower(detectedExt) == "gff") || (tolower(detectedExt) == "gff3")) {
			fileExt=tolower(detectedExt)
			fileType = "GFF3"			
		} else {
			message(paste0("Unknown file type", detectedExt,": Defaulting to gff3"))
			fileExt="gff3"		
			fileType = "GFF3"				
		}
	}

	
	pairsFile=paste0(inFile,".pairs")
	if (!base::file.exists(pairsFile)) {
		message ("Creating pairs file...")
		get_pairs_from_gff(inFile,fileType=fileType)
	}

	myTree<-get_features(inFile,outFormat="tree", includeCounts=includeCounts, fileType=fileType ) 

	if(includeCounts) {	
		data.tree::SetNodeStyle(myTree, label = function(x) {paste0(x$name,":", x$N)} )
	}

	if(missing(outFile)){		
		return(plot(myTree))
	} else {

		acceptedFormats<-c("png","pdf","svg")

		if (!requireNamespace("DiagrammeRsvg", quietly = TRUE) || !requireNamespace("rsvg", quietly = TRUE)) {    
			stop("Packages \"DiagrammeRsvg\" and \"rsvg\" are required for exporting a `data.tree::Node` object. Please install them.")
		}

		detectedOutExt<-tools::file_ext(outFile)
		
		if(tolower(detectedOutExt) %in% acceptedFormats){
			foutput=outFile			
			if(tolower(exportFormat) != tolower(detectedOutExt)){
				message(paste0("fileType automatically set to ",tolower(detectedOutExt)," based on output file name"))
				exportFormat<-tolower(detectedOutExt)
			}

		} else {
			if(!(tolower(exportFormat) %in% acceptedFormats)){
				message(paste0("Invalid fileType (",exportFormat,"). Allowed formats are ",paste(acceptedFormats,collapse=", ")))
				return()
			}

			foutput=paste0(inFile,".",exportFormat)
		}
		
		data.tree::ToDiagrammeRGraph(myTree) %>% DiagrammeR::export_graph(file_name=foutput, file_type=exportFormat)

		rm(list=setdiff(ls(), "foutput"))	
		gc(reset=T)
	
		return(foutput) 
	}

}



#' Analyses the feature type hierarchy of a gff file 
#'
#' Based on the feature type hierarchy a gff file, this function creates and returns 
#' a feature tree or a feature dependency table. 
#'
#' @param inFile Path to the input GTF/GFF3 features  file 
#' @param includeCounts Include number of occurrences of each feature and subfeature 
#' @param outFormat Output format of the function. Available formats are: tree (DEFAULT), data.frame and JSON.
#' @param fileType Version of the input file (GTF/GFF3). Default AUTO: determined from the file name.
#' @return Depending on the outFormat selected returns a feature tree (tree), a feature dependency table as data.frame (data.frame) or a feature dependency table as JSON object (JSON)
#'
#' @export
#'
#' @examples
#' test_gff3<-system.file("extdata", "AthSmall.gff3", package="Rgff")
#' get_features(test_gff3)
get_features<- function(inFile, includeCounts=FALSE, outFormat=c("tree", "data.frame", "JSON"),fileType = c("AUTO", "GFF3","GTF")){

	fileType <- match.arg(fileType)
	outFormat <- match.arg(outFormat)
	
	pairsFile=paste0(inFile,".pairs")
	if (!base::file.exists(pairsFile)) {
		message ("Creating pairs file...")
		get_pairs_from_gff(inFile,fileType=fileType)
	}

	if(is.data.frame(pairsFile)){
		myPairData<-pairsFile
	} else if (is.vector(pairsFile) && is.character(pairsFile) && length(pairsFile)==1) {
		if(file.exists(pairsFile)){
			myPairData<-utils::read.table(pairsFile,header=T,sep="\t")
		} else {
			stop("Only Data frame of pairs or a valid path of pairs file are allowed")
		}
	}
	if(!all(c("ELEMENT","PARENT") %in% names(myPairData))){
		stop("Missing ELEMENT or PARENT columns in pair data")
	} 
	
	myTree<-data.tree::as.Node(myPairData,mode="network") 
	if(outFormat == "tree"){
		return(myTree)
	} else {
		if(includeCounts) {
			Descendants <- function(node) {
				descData<-data.tree::ToDataFrameTree(node,"name","N")[-1,c("name","N")]		
				aggrData<-descData %>% dplyr::group_by(.data$name) %>% dplyr::summarise(N =sum(.data$N))
				return(paste(paste(aggrData$name, aggrData$N, sep=":"), collapse=" "))
			}
		} else {
			Descendants <- function(node) {
				descData<-unique(data.tree::ToDataFrameTree(node,"name")[-1,c("name")])
				return(paste(descData, collapse=" "))
			}
		}

		myTree$Do(function(node) node$descendants <- Descendants(node))

		qFeatures<-as.data.frame((myTree$Get("descendants"))[-1])
		colnames(qFeatures)<-c("BLOCKS")

		colnames(qFeatures)<-c("BLOCKS")
		if(includeCounts){
			qNames<-dplyr::left_join(data.frame(name=names((myTree$Get("descendants"))[-1])),data.tree::ToDataFrameTree(myTree,"name","N")[-1,c("name","N")] %>% dplyr::group_by(.data$name) %>% dplyr::summarise(N =sum(.data$N)))
			qFeatures$FEATURES<-with(qNames, paste0(name, ":", N))
			
		} else {
			qFeatures$FEATURES<-names((myTree$Get("descendants"))[-1])
		}

		qFeatures<-unique(qFeatures)


		if(outFormat == "data.frame"){
			return (qFeatures)
		} else if(outFormat == "JSON"){
			qFeaturesList <- as.list(as.data.frame(t(qFeatures)))
			outText<-RJSONIO::toJSON(list(features=qFeaturesList),.escapeEscapes = FALSE)
	
			return(outText)
		} 
	}
}

#' Sorts a GFF3 file
#'
#' This function produces a sorted GFF3 file from an unsorted GFF3 file.
#' The output file will be sorted by Chromosome, Start, End (reverse) and feature (based on the precedency in feature tree)
#'
#' @param gffFile Path to the input GFF3 file
#' @param outFile Path to the output GFF3 file, if not provided the output will be  the input path (without extension) with the suffix sorted.gff3
#' @param forceOverwrite If output file exists, overwrite the existing file. (default FALSE)
#' @return Path to the sorted GFF3 file
#'
#' @importFrom rlang .data
#' @keywords internal
sort_gff3<-function(gffFile, outFile, forceOverwrite=FALSE){
	if (!base::file.exists(gffFile)) { 
		stop(paste0("Input gff3 file ",gffFile," not found."))
	}

	if(missing(outFile)){
		foutput=paste0(sub(pattern = "(.*)\\..{0,4}$", replacement = "\\1", gffFile),".sorted.gff3")
	} else {
		foutput=outFile
	}

	if (base::file.exists(foutput) && forceOverwrite==FALSE) { 
		message ("Sorted file for this GFF already exists. Use forceOverwrite=TRUE to overwrite the existing file."); 
		return(foutput) 
	}



	fileConn<-file(gffFile,"r")
	gffHeader<-c()
	readingHeader<-TRUE
	CommentPattern<-"(^#(.*)$)|(^#([\t ]*)$)"
	BlankPattern<-"(^[\t ]*$)"
	while(readingHeader){
		hLine<-readLines(fileConn,n=1)
		if(identical(hLine, character(0))){
			readingHeader<-FALSE	
		} else if(grepl(CommentPattern,hLine,perl=T)){
			gffHeader<-c(gffHeader,hLine)
		} else if(grepl(BlankPattern,hLine,perl=T)){
		} else {
			readingHeader<-FALSE
		}
	}
	close(fileConn)

	testTree<-get_features(gffFile,outFormat="tree", fileType="GFF3" ) 

	featureList<-rev(unique(rev((sapply(data.tree::Traverse(testTree,traversal="level"), function(node){{node$name}})[-1]))))

	featureOrder<-factor(featureList, levels=featureList)
	
	gffTable<-utils::read.delim(gffFile, comment.char="#",  header = FALSE, sep = "\t", blank.lines.skip = TRUE)

	orderedGff<- gffTable %>% dplyr::arrange(.data$V1, .data$V4, dplyr::desc(.data$V5), factor(.data$V3, levels = featureList))


	fileConn<-file(foutput,"w+")

	if(length(gffHeader)> 0){
		writeLines(gffHeader, fileConn)
	}
	utils::write.table(orderedGff, fileConn, append=TRUE, sep="\t", quote=F, row.names = FALSE, col.names=FALSE)

	close(fileConn)	

	rm(list=setdiff(ls(), "foutput"))	
	gc(reset=T)
	
	return(foutput)

}

#' Creates a features pair file from a GTF file
#'
#' This function creates a features pair file from a GTF file
#'
#' @param gtfFile Path to the input GTF file
#' @param outFile Path to the output features pair file, if not provided the output will be gtfFile.pairs
#' @param forceOverwrite If output file exists, overwrite the existing file. (default FALSE)
#' @return Path to the generated features pair file  
#'
#' @importFrom rlang .data
#' @keywords internal
get_pairs_from_gtf<-function(gtfFile, outFile, forceOverwrite=FALSE){
	if (!base::file.exists(gtfFile)) { 
		stop(paste0("Input GTF file ",gtfFile," not found."))
	}

	if(missing(outFile)){
		foutput=paste0(gtfFile,".pairs")
	} else {
		foutput=outFile
	}

	if (base::file.exists(foutput) && forceOverwrite==FALSE) { 
		## ToDo: Remove message if outFile is missing
		if(!missing(outFile)){
			message ("Pairs file for this GTF already exists. Use forceOverwrite=TRUE to overwrite the existing file."); 
		}
		return(foutput) 
	}

	cnames <- c("seqname","source","feature","start","end","score","strand","frame","attribute")
	myGTFdata <- utils::read.delim(gtfFile, header=F,  sep="\t", comment.char = "#", quote="", blank.lines.skip=T) %>%    `colnames<-`(cnames)  

	featureList<- myGTFdata %>% dplyr::group_by(.data$feature) %>% dplyr::summarise(N = dplyr::n())
	parentElem<-sapply(featureList$feature,function(x){ if(x == "transcript") {"gene"} else { if ( x == "gene") {""} else {"transcript"}}})

	DF<-data.frame(ELEMENT=featureList$feature, PARENT= parentElem, N=featureList$N)


	# CHECK FOR MISSING TRANSCRIPT OR GENE PARENTS
	transcriptPattern<-" *transcript_id +\"?([^\"]+)\"? *(;|$)"
	transcriptPatternFull<-".* *transcript_id +\"?([^\"]+)\"? *(;|$).*"

	

	addTranscripts<-!any(grepl(transcriptPattern,myGTFdata[myGTFdata$feature == "transcript","attribute"],perl=T,ignore.case=T))
	if(addTranscripts){
			selected_transcriptid_lines<-grep(transcriptPattern,myGTFdata[,"attribute"],perl=T,ignore.case=T)
		transcriptIDs<-unique(gsub(transcriptPatternFull,"\\1",myGTFdata[selected_geneid_lines,"attribute"],perl=T,ignore.case=T))
		DF[nrow(DF)+1,]<-c("transcript","gene",length(transcriptIDs))
	}

	genePattern<-" *gene_id +\"?([^\"]+)\"? *(;|$)"
	genePatternFull<-".* *gene_id +\"?([^\"]+)\"? *(;|$).*"

	addGenes<-!any(grepl(genePattern,myGTFdata[myGTFdata$feature == "gene","attribute"],perl=T,ignore.case=T))
	if(addGenes){
		selected_geneid_lines<-grep(genePattern,myGTFdata[,"attribute"],perl=T,ignore.case=T)		
		geneIDs<-unique(gsub(genePatternFull,"\\1",myGTFdata[selected_geneid_lines,"attribute"],perl=T,ignore.case=T))
		DF[nrow(DF)+1,]<-c("gene","",length(geneIDs))		
	}

	# END CHECK FOR MISSING TRANSCRIPT OR GENE PARENTS

	
	utils::write.table(DF,file=foutput,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t", append=FALSE)

	rm(list=setdiff(ls(), "foutput"))	
	gc(reset=T)

	return(foutput)
}
	

#' Creates a features path file from a GTF file
#'
#' This function creates a features path file from a GTF file
#'
#' @param gtfFile Path to the input GTF file
#' @param outFile Path to the output features path file, if not provided the output will be gtfFile.paths
#' @param forceOverwrite If output file exists, overwrite the existing file. (default FALSE)
#' @return Path to the generated features path file
#'
#' @importFrom rlang .data
#' @keywords internal
gtf_to_paths<-function(gtfFile, outFile, forceOverwrite=FALSE){ 
	if (!base::file.exists(gtfFile)) { 
		stop(paste0("Input gtf file ",gtfFile," not found."))
	}

	if(missing(outFile)){
		foutput=paste0(gtfFile,".paths")
	} else {
		foutput=outFile
	}

	if (base::file.exists(foutput) && forceOverwrite==FALSE) {
		if(!missing(outFile)){
			message ("Paths file for this GTF already exists. Use forceOverwrite=TRUE to overwrite the existing file."); 
		}
	}


	cnames <- c("seqname","source","feature","start","end","score","strand","frame","attribute")
	myGTFdata <- utils::read.delim(gtfFile, header=F,  sep="\t", comment.char = "#", quote="", blank.lines.skip=T) %>%    `colnames<-`(cnames)  

	fLevels<-unique(myGTFdata$feature)
	tIndex<-which(fLevels=="transcript")
	if(length(tIndex)>0){
		fLevels<-c(fLevels[tIndex],fLevels[-tIndex])
	}
	
	gIndex<-which(fLevels=="gene")
	if(length(gIndex)>0){
		fLevels<-c(fLevels[gIndex],fLevels[-gIndex])
	}

	
	myGTFdata$feature<-factor(myGTFdata$feature, ordered=T,levels=fLevels)

	myGTFdata<- myGTFdata %>% dplyr::arrange(.data$seqname, .data$start, dplyr::desc(.data$end), factor(.data$feature, levels = levels(myGTFdata$feature)))
	
	transcriptPattern<-" *transcript_id +\"?([^\"]+)\"? *(;|$)"
	transcriptPatternFull<-".* *transcript_id +\"?([^\"]+)\"? *(;|$).*"

	genePattern<-" *gene_id +\"?([^\"]+)\"? *(;|$)"
	genePatternFull<-".* *gene_id +\"?([^\"]+)\"? *(;|$).*"


	selected_geneid_lines<-grep(genePattern,myGTFdata[,"attribute"],perl=T,ignore.case=T)
	selected_transcriptid_lines<-grep(transcriptPattern,myGTFdata[,"attribute"],perl=T,ignore.case=T)

	geneIDs<-rep("",nrow(myGTFdata))
	transcriptIDs<-rep("",nrow(myGTFdata))

	geneIDs[selected_geneid_lines]<-gsub(genePatternFull,"\\1",myGTFdata[selected_geneid_lines,"attribute"],perl=T,ignore.case=T)
	transcriptIDs[selected_transcriptid_lines]<-gsub(transcriptPatternFull,"\\1",myGTFdata[selected_transcriptid_lines,"attribute"],perl=T,ignore.case=T)


	featureCounter<-rep("",nrow(myGTFdata))
	cFeatures<-fLevels[fLevels != "transcript" & fLevels != "gene"]

	for (cf in cFeatures){
		fLines<-which(myGTFdata$feature == cf)
		featureCounter[fLines]<-c(1:length(fLines))
	}

	myGTFdata <- myGTFdata %>% dplyr::mutate(gene_id=geneIDs, transcript_id=transcriptIDs, fcounter=featureCounter)

	# CHECK FOR MISSING TRANSCRIPT OR GENE PARENTS
	addGenes<-!any(grepl(genePattern,myGTFdata[myGTFdata$feature == "gene","attribute"],perl=T,ignore.case=T))
	if(addGenes){
		geneList<-data.frame(geneID=gsub(genePatternFull,"\\1",myGTFdata[selected_geneid_lines,"attribute"],perl=T,ignore.case=T),start=myGTFdata[selected_geneid_lines,"start"],end=myGTFdata[selected_geneid_lines,"end"] )

		geneDF <- geneList %>% dplyr::group_by(.data$geneID) %>% dplyr::summarise(start=min(.data$start), end=max(.data$end)) %>% dplyr::mutate(ASSIGNED=0)

	}

	addTranscripts<-!any(grepl(transcriptPattern,myGTFdata[myGTFdata$feature == "transcript","attribute"],perl=T,ignore.case=T))
	if(addTranscripts){
		transcriptList<-data.frame(transcriptID=gsub(transcriptPatternFull,"\\1",myGTFdata[selected_geneid_lines,"attribute"],perl=T,ignore.case=T),start=myGTFdata[selected_geneid_lines,"start"],end=myGTFdata[selected_geneid_lines,"end"] )

		transcriptDF <- transcriptList %>% dplyr::group_by(.data$transcriptID) %>% dplyr::summarise(start=min(.data$start), end=max(.data$end)) %>% dplyr::mutate(ASSIGNED=0)
	}

	
	fileConn<-base::file(foutput, 'w')
	apply(myGTFdata,1, function(x){
		fType<-x["feature"]
		if (fType == "gene") {
			cat(paste0(x["seqname"],"\t",as.numeric(x["start"]),"\t",as.numeric(x["end"]),"\t",x["strand"],"\tgene;",x["gene_id"],"\n"),file=fileConn,append=TRUE, sep="")					
		} else {
			if (addGenes && x["gene_id"] != ""){
				if(geneDF[geneDF$geneID==x["gene_id"],]$ASSIGNED == 0){
					cat(paste0(x["seqname"],"\t",geneDF[geneDF$geneID==x["gene_id"],]$start,"\t",geneDF[geneDF$geneID==x["gene_id"],]$end,"\t",x["strand"],"\tgene;",x["gene_id"],"\n"),file=fileConn,append=TRUE, sep="")
					geneDF[geneDF$geneID==x["gene_id"],]$ASSIGNED == 1
				}
			}
		
			if (fType == "transcript") {
				cat(paste0(x["seqname"],"\t",as.numeric(x["start"]),"\t",as.numeric(x["end"]),"\t",x["strand"],"\ttranscript;",x["transcript_id"],";gene;",x["gene_id"],"\n"),file=fileConn,append=TRUE, sep="")
			} else {
				if (addTranscripts && x["transcript_id"] != ""){			
					if(transcriptDF[transcriptDF$transcriptID==x["transcript_id"],]$ASSIGNED == 0){
						cat(paste0(x["seqname"],"\t",transcriptDF[transcriptDF$transcriptID==x["transcript_id"],]$start,"\t",transcriptDF[transcriptDF$transcriptID==x["transcript_id"],]$end,"\t",x["strand"],"\ttranscript;",x["transcript_id"],";gene;",x["gene_id"],"\n"),file=fileConn,append=TRUE, sep="")
						transcriptDF[transcriptDF$transcriptID==x["transcript_id"],]$ASSIGNED == 1
					}
				}
				cat(paste0(x["seqname"],"\t",as.numeric(as.numeric(x["start"])),"\t",as.numeric(as.numeric(x["end"])),"\t",x["strand"],"\t",fType,";",fType,"_",x["fcounter"],";transcript;",x["transcript_id"],";gene;",x["gene_id"],"\n"),file=fileConn,append=TRUE, sep="")
			}
		}
	})	
	close(fileConn)

	rm(list=setdiff(ls(), "foutput"))	
	gc(reset=T)

	return(foutput)
	
}



#' Sorts a GTF file
#'
#' This function produces a sorted GTF file from an unsorted GTF file.
#' The default order is by Chromosome, Start, End (reverse) and feature (gene -> transcript -> OTHERS)
#'
#' @param gtfFile Path to the input GTF file
#' @param outFile Path to the output GTF file, if not provided the output will be  the input path (without extension) with the suffix sorted.gtf
#' @param forceOverwrite If output file exists, overwrite the existing file. (default FALSE)
#' @return Path to the sorted GTF file
#'
#' @importFrom rlang .data
#' @keywords internal
sort_gtf<-function(gtfFile, outFile, forceOverwrite=FALSE){
	if (!base::file.exists(gtfFile)) { 
		stop(paste0("Input gtf file ",gtfFile," not found."))
	}

	if(missing(outFile)){
		foutput=paste0(sub(pattern = "(.*)\\..{0,4}$", replacement = "\\1", gtfFile),".sorted.gtf")
	} else {
		foutput=outFile
	}

	if (base::file.exists(foutput) && forceOverwrite==FALSE) { 
		message ("Sorted file for this GTF already exists. Use forceOverwrite=TRUE to overwrite the existing file."); 
		return(foutput) 
	}



	fileConn<-file(gtfFile,"r")
	gtfHeader<-c()
	readingHeader<-TRUE
	CommentPattern<-"(^#(.*)$)|(^#([\t ]*)$)"
	BlankPattern<-"(^[\t ]*$)"
	while(readingHeader){
		hLine<-readLines(fileConn,n=1)
		if(identical(hLine, character(0))){
			readingHeader<-FALSE	
		} else if(grepl(CommentPattern,hLine,perl=T)){
			gtfHeader<-c(gtfHeader,hLine)
		} else if(grepl(BlankPattern,hLine,perl=T)){
		} else {
			readingHeader<-FALSE
		}
	}
	close(fileConn)

	cnames <- c("seqname","source","feature","start","end","score","strand","frame","attribute")
	myGTFdata <- utils::read.delim(gtfFile, header=F,  sep="\t", comment.char = "#", quote="", blank.lines.skip=T) %>%    `colnames<-`(cnames)  

	fLevels<-unique(myGTFdata$feature)
	tIndex<-which(fLevels=="transcript")
	if(length(tIndex)>0){
		fLevels<-c(fLevels[tIndex],fLevels[-tIndex])
	}
	
	gIndex<-which(fLevels=="gene")
	if(length(gIndex)>0){
		fLevels<-c(fLevels[gIndex],fLevels[-gIndex])
	}

	
	myGTFdata$feature<-factor(myGTFdata$feature, ordered=T,levels=fLevels)

	orderedGTFdata<- myGTFdata %>% dplyr::arrange(.data$seqname, .data$start, dplyr::desc(.data$end), factor(.data$feature, levels = levels(myGTFdata$feature)))
	
	fileConn<-file(foutput,"w+")

	if(length(gtfHeader)> 0){
		writeLines(gtfHeader, fileConn)
	}
	utils::write.table(orderedGTFdata, fileConn, append=TRUE, sep="\t", quote=F, row.names = FALSE, col.names=FALSE)

	close(fileConn)	

	rm(list=setdiff(ls(), "foutput"))	
	gc(reset=T)
	
	return(foutput)

}

#####################################################
##  GTF/GFF WRAPPERS
#####################################################

#' Sorts a GTF/GFF3 file
#'
#' This function produces a sorted gff file from an unsorted gff file.
#' The default order is by Chromosome, Start, End (reverse) and feature (based on the precedency in feature tree)
#'
#' @param inFile Path to the input gff file
#' @param outFile Path to the output sorted file, if not provided the output will be the input path (without extension) with the suffix sorted.gtf/gff3
#' @param fileType Version of the input file (GTF/GFF3). Default AUTO: determined from the file name.
#' @param forceOverwrite If output file exists, overwrite the existing file. (default FALSE)
#' @return Path to the sorted feature file
#'
#' @export
#'
#' @examples
#' test_gff3<-system.file("extdata", "eden.gff3", package="Rgff")
#' sort_gff(test_gff3)


sort_gff<-function(inFile, outFile, fileType=c("AUTO","GFF3","GTF"), forceOverwrite=FALSE){
	fileType <- match.arg(fileType)

	if (!base::file.exists(inFile)) { 
		stop(paste0("Input file ",inFile," not found."))
	}
	detectedExt<-tools::file_ext(inFile)

	if(fileType == "AUTO"){
		if(tolower(detectedExt) == "gtf"){
			fileExt="gtf"
			fileType = "GTF"
		} else if( (tolower(detectedExt) == "gff") || (tolower(detectedExt) == "gff3")) {
			fileExt=tolower(detectedExt)
			fileType = "GFF3"			
		} else {
			message(paste0("Unknown file type", detectedExt,": Defaulting to gff3"))
			fileExt="gff3"		
			fileType = "GFF3"				
		}
	}

	if(missing(outFile)){
		foutput=paste0(sub(pattern = "(.*)\\..{0,4}$", replacement = "\\1", inFile),".sorted.",fileExt)
	} else {
		foutput=outFile
	}

	if (base::file.exists(foutput) && forceOverwrite==FALSE) { 
		message ("Sorted file for this annotation file already exists. Use forceOverwrite=TRUE to overwrite the existing file."); 
		return(foutput) 
	}
	message(foutput)

	if(fileType == "GTF"){
		return(sort_gtf(inFile, foutput, forceOverwrite))
	} else {
		return(sort_gff3(inFile, foutput, forceOverwrite))
	}
}




#' Creates a SAF file from a GTF/GFF3 features for the given pairs of blocks and features
#'
#' This function creates a SAF file from a GTF/GFF3 features for the given blocks 
#' and features
#'
#' @param inFile Path to the input gff file
#' @param outFile Path to the output SAF file, if not provided the output path will be the input path with the suffix ".feature1-block1.feature2-block2(...).saf"
#' @param fileType Version of the input file (GTF/GFF3). Default AUTO: determined from the file name.
#' @param forceOverwrite If output file exists, overwrite the existing file. (default FALSE)
#' @param features Vector of pairs of features/blocks, separated by '>' (see sep argument). In the case of features without defined blocks, only the feature is needed (see example)
#' @param sep Separator of each "feature" and "block" provided in the feature argument (default '>')
#' @return Path to the generated SAF file
#'
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' test_gff3<-system.file("extdata", "AthSmall.gff3", package="Rgff")
#' ## Default usage, extract gene features by exon blocks 
#' saf_from_gff(test_gff3)
#' ## Define only feature without block to count reads within the whole genomic locus 
#' saf_from_gff(test_gff3, features=c("gene"))
#' ## Define multiple features for counting readsoverlapping only in exonic regions 
#' saf_from_gff(test_gff3, features=c("gene > exon", "ncRNA_gene > exon"))
saf_from_gff<-function(inFile, outFile, fileType=c("AUTO","GFF3","GTF"), forceOverwrite=FALSE, features=c("gene > exon"), sep=">") {
	fileType <- match.arg(fileType)

	if (!base::file.exists(inFile)) { 
		stop(paste0("Input annotation file ",inFile," not found."))
	}
	detectedExt<-tools::file_ext(inFile)

	if(fileType == "AUTO"){
		if(tolower(detectedExt) == "gtf"){
			fileExt="gtf"
			fileType = "GTF"
		} else if( (tolower(detectedExt) == "gff") || (tolower(detectedExt) == "gff3")) {
			fileExt=tolower(detectedExt)
			fileType = "GFF3"			
		} else {
			message(paste0("Unknown file type", detectedExt,": Defaulting to gff3"))
			fileExt="gff3"		
			fileType = "GFF3"				
		}
	}


	escapedSep<-gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", sep)

	# Read features and check 
	featureDF<-get_features(inFile, outFormat="data.frame",fileType = fileType)

	featureList<-c()
	blockList<-c()
	for (i in c(1:length(features))){
		groupByAndBlock<-trimws(strsplit(features[i],escapedSep)[[1]])
		groupBy=groupByAndBlock[1]
		block=""
		
		if(length(groupByAndBlock)>2 ){
			stop(paste0("Invalid feature: ",features[i],". Each feature element must be composed by at most two elements with the syntax feature ",sep," block."))
		}

		if(!(groupBy %in% (featureDF$FEATURES))){
			stop(paste0("Feature: '",groupBy,"' is not defined in ",inFile,". \nRun get_features(\"",inFile,"\") to see the feature tree"))
		}
		
		if(length(groupByAndBlock)>1 ){
			block=groupByAndBlock[2]

			
			if(featureDF[featureDF$FEATURES == groupBy,]$BLOCKS == "" || !(block %in% (strsplit(as.character(featureDF[featureDF$FEATURES == groupBy,]$BLOCKS)," ")[[1]]))){
				stop(paste0("Block '",block, "' is not part of the feature '", groupBy,"' in ",inFile,". \nRun get_features(\"",inFile,"\") to see the feature tree"))
			}

			
		} 
		
		featureList<-c(featureList,groupBy)
		blockList<-c(blockList,block)		
	}

	if(missing(outFile)){
		fSuffix<-""
		for (i in c(1:length(featureList))){
				fSuffix<-paste0(fSuffix,".",featureList[i])
				if(blockList[i] != ""){
					fSuffix<-paste0(fSuffix,"-",blockList[i])
				}
		}
		foutput=paste0(inFile,fSuffix,".saf")
	} else {
		foutput=outFile
	}


	if(base::file.exists(foutput) && forceOverwrite==FALSE){
		message("Output SAF file exists. Use forceOverwrite=TRUE to overwrite the existing file.")
		return(foutput)
	}
	
	PATHSfile=paste0(inFile,".paths")
	if (!base::file.exists(PATHSfile)) {
		message ("Creating PATHS file... please be patient...")
		if(fileType == "GFF3"){
			gff3_to_paths(inFile)
		} else {
			gtf_to_paths(inFile)			
		}
	}

	

	for (i in c(1:length(featureList))){
		if(blockList[i] == ""){
			blockParam=c()
		} else {
			blockParam=c(blockList[i])
		}
		if(i==1){
			safDF<-unique(saf_from_paths(PATHSfile,groupBy=c(featureList[i]),block=blockParam))
		} else {

			safDF<-rbind(safDF,unique(saf_from_paths(PATHSfile,groupBy=c(featureList[i]),block=blockParam)))
		}
	}	
	safDF<-unique(safDF[order(safDF[,2], safDF[,3], safDF[,4], safDF[,1] ),])

	utils::write.table(safDF,foutput,sep="\t",quote=FALSE,row.names = FALSE)

	rm(list=setdiff(ls(), "foutput"))	
	gc(reset=T)
	
	return(foutput)
}	


#' Creates a features pair file from a gff file
#'
#' This function creates a features pair file from a gff file
#'
#' @param inFile Path to the input gff file
#' @param outFile Path to the output features pair file, if not provided the output will be inFile.pairs
#' @param fileType Version of the input file (GTF/GFF3). Default AUTO: determined from the file name.
#' @param forceOverwrite If output file exists, overwrite the existing file. (default FALSE)
#' @return Path to the generated features pair file  
#'
#' @keywords internal

get_pairs_from_gff<-function(inFile, outFile, fileType=c("AUTO","GFF3","GTF"), forceOverwrite=FALSE){
	fileType <- match.arg(fileType)

	if (!base::file.exists(inFile)) { 
		stop(paste0("Input file ",inFile," not found."))
	}
	detectedExt<-tools::file_ext(inFile)

	if(fileType == "AUTO"){
		if(tolower(detectedExt) == "gtf"){
			fileExt="gtf"
			fileType = "GTF"
		} else if( (tolower(detectedExt) == "gff") || (tolower(detectedExt) == "gff3")) {
			fileExt=tolower(detectedExt)
			fileType = "GFF3"			
		} else {
			message(paste0("Unknown file type", detectedExt,": Defaulting to gff3"))
			fileExt="gff3"		
			fileType = "GFF3"				
		}
	}

	if(missing(outFile)){
		foutput=paste0(inFile,".pairs")
	} else {
		foutput=outFile
	}

	if (base::file.exists(foutput) && forceOverwrite==FALSE) { 
		## ToDo: Remove message if outFile is missing
		if(!missing(outFile)){
			message ("Pairs file for this feature file already exists. Use forceOverwrite=TRUE to overwrite the existing file."); 
		}
		return(foutput) 
	}

	if(fileType == "GFF3"){
		return(get_pairs_from_gff3(inFile, outFile, forceOverwrite))
	} else {
		return(get_pairs_from_gtf(inFile, outFile, forceOverwrite))
	}
}


#' Converts a GTF file into a GFF3 file
#'
#' This function converts a GTF file into a GFF3 file mantaining the feature hierarchy 
#' defined by the gene_id and transcript_id attributes. The remaining attributes of each feature 
#' will be kept with the same name and value. 
#'
#' @param gtfFile Path to the input GTF file
#' @param outFile Path to the output GFF3 file, inf not provided the output will be gtfFile.gff3
#' @param forceOverwrite If output file exists, overwrite the existing file. (default FALSE)
#' @return Path to the generated GFF3 file
#'
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' test_gtf<-system.file("extdata", "AthSmall.gtf", package="Rgff")
#' gtf_to_gff3(test_gtf)
#' }
gtf_to_gff3<-function(gtfFile, outFile, forceOverwrite=FALSE){

	if (!requireNamespace("rtracklayer", quietly = TRUE)) {    
		stop("Package \"rtracklayer\" is required to convert from GTF to GFF3. Please install it.")
	}

	if(missing(outFile)){
		foutput=paste0(gtfFile,".gff3")
	} else {
		foutput=outFile
	}

	if(base::file.exists(foutput) && forceOverwrite==FALSE){
		message("Output GFF3 file exists. Use forceOverwrite=TRUE to overwrite the existing file")
		return(foutput)
	}


	# startTime <- as.numeric(Sys.time())

	myGTFdata<-rtracklayer::import.gff2(gtfFile)

	GroupedGenes <- (as.data.frame(myGTFdata) %>% dplyr::group_by(.data$gene_id) %>% dplyr::summarise(gtypes = paste0(unique(.data$type), collapse = " ")) )

	GroupedTranscripts <- (as.data.frame(myGTFdata) %>% dplyr::group_by(.data$transcript_id) %>% dplyr::summarise(ttypes = paste0(unique(.data$type), collapse = " ")) )


	############# CHECK MISSING TRANSCRIPTS
	addTranscript<-!any(grepl("transcript ",unique(GroupedTranscripts$ttypes)))
	if(addTranscript){

			TranscriptsCoords <- as.data.frame(myGTFdata) %>% dplyr::group_by(.data$transcript_id,.data$seqnames, .data$strand, .data$source) %>% dplyr::summarise(start = min(.data$start), end = max(.data$end))
			TranscriptsCoords$type<-"transcript"
			## Todo ADD UNIQUE FIELDS?
			presentTranscripts<-grep("transcript ",GroupedTranscripts$ttypes)
			if(length(presentTranscripts) > 0){
				TranscriptsCoords<-TranscriptsCoords[TranscriptsCoords$transcript_id %in% GroupedTranscripts[-presentTranscripts,]$transcript_id]
			} 
			

			

			TranscriptsCoordsGR<-GenomicRanges::makeGRangesFromDataFrame(TranscriptsCoords,
							 keep.extra.columns=TRUE,
							 ignore.strand=FALSE,
							 seqinfo=NULL,
							 seqnames.field=c("seqnames"),
							 start.field="start",
							 end.field="end",
							 strand.field="strand")
			
			for( i in 1:nrow(TranscriptsCoords)){
				TranscriptPos<-match(TranscriptsCoords[i,"transcript_id"], myGTFdata$transcript_id)-1
				myGTFdata<-S4Vectors::append(myGTFdata,TranscriptsCoordsGR[i,],TranscriptPos)			
			}
	}

	############# CHECK MISSING GENES
	addGenes<-!any(grepl("gene ",unique(GroupedGenes$gtypes)))
	if(addGenes){

			GenesCoords <- as.data.frame(myGTFdata) %>% dplyr::group_by(.data$gene_id,.data$seqnames, .data$strand, .data$source) %>% dplyr::summarise(start = min(.data$start), end = max(.data$end))
			GenesCoords$type<-"gene"
			## Todo ADD UNIQUE FIELDS?
			presentGenes<-grep("gene ",GroupedGenes$gtypes)
			if(length(presentGenes) > 0){
				GenesCoords<-GenesCoords[GenesCoords$gene_id %in% GroupedGenes[-presentGenes,]$gene_id,]
			} 
			


			GenesCoordsGR<-GenomicRanges::makeGRangesFromDataFrame(GenesCoords,
							 keep.extra.columns=TRUE,
							 ignore.strand=FALSE,
							 seqinfo=NULL,
							 seqnames.field=c("seqnames"),
							 start.field="start",
							 end.field="end",
							 strand.field="strand")
			
			myGTFdata<-S4Vectors::append(myGTFdata,GenesCoordsGR)
	}

	
	tIndex<-which(levels(myGTFdata$type)=="transcript")
	newLevels<-c("transcript",levels(myGTFdata$type)[-tIndex])
	gIndex<-which(newLevels=="gene")
	newLevels<-c("gene",newLevels[-gIndex])

	# sort(myGTFdata, by = ~ seqnames + start + rev(end))
	
	# sort(myGTFdata)
	orderedGff<- as.data.frame(myGTFdata) %>% dplyr::arrange(.data$seqnames, .data$start, dplyr::desc(.data$end), factor(.data$type, levels = newLevels))


	myGTFdata<-GenomicRanges::makeGRangesFromDataFrame(orderedGff,
					 keep.extra.columns=TRUE,
					 ignore.strand=FALSE,
					 seqnames.field=c("seqnames"),
					 start.field="start",
					 end.field="end",
					 strand.field="strand")
			
			
	#		for( i in 1:nrow(GenesCoords)){
	#			GenePos<-match(GenesCoords[i,"gene_id"], myGTFdata$gene_id)-1
	#			myGTFdata<-S4Vectors::append(myGTFdata,GenesCoordsGR[i,],GenePos)			
	#		}

	fLevels<-unique(myGTFdata$type)
	featureCounter<-rep("",length(myGTFdata))
	cFeatures<-fLevels[fLevels != "transcript" & fLevels != "gene"]

	for (cf in cFeatures){
		fLines<-which(myGTFdata$type == cf)
		featureCounter[fLines]<-c(1:length(fLines))
	}
	S4Vectors::mcols(myGTFdata)$featureCounter<-featureCounter



	#apply(myGTFdata, MARGIN = 1, function(x) {if(x$type=="gene"){x$gene_id} else if(x$type=="transcript"){x$transcript_id} else {NA}})

	#IDs<-apply(as.data.frame(S4Vectors::mcols(myGTFdata)), MARGIN = 1, function(x) { if (as.character(x["type"])=="gene"){as.character(x["gene_id"])} else if(as.character(x["type"]) =="transcript"){as.character(x["transcript_id"])} else {NA}})
	IDs<-apply(as.data.frame(S4Vectors::mcols(myGTFdata)), MARGIN = 1, function(x) { 
		if (as.character(x["type"])=="gene"){
			as.character(x["gene_id"])
		} else if(as.character(x["type"]) =="transcript"){
			as.character(x["transcript_id"])
		} else {
			paste0(as.character(x["type"]),"_",x["featureCounter"])
		}
		
	})

	S4Vectors::mcols(myGTFdata)$featureCounter<-NULL


	# Find repeated IDs
	n_occur <- data.frame(table(IDs))
	if(nrow(n_occur[n_occur$Freq > 1,])>0){
			geneRows<-which(S4Vectors::mcols(myGTFdata)$type=="gene")
			IDs[geneRows]<-paste0("gene:",IDs[geneRows])

			transcriptRows<-which(S4Vectors::mcols(myGTFdata)$type=="transcript")
			IDs[transcriptRows]<-paste0("transcript:",IDs[transcriptRows])
			
			geneIdRows<-which(!is.na(S4Vectors::mcols(myGTFdata)$gene_id))			
			S4Vectors::mcols(myGTFdata)$gene_id[geneIdRows] = paste0("gene:",S4Vectors::mcols(myGTFdata)$gene_id[geneIdRows])
			
			transcriptIdRows<-which(!is.na(S4Vectors::mcols(myGTFdata)$transcript_id))
			S4Vectors::mcols(myGTFdata)$transcript_id[transcriptIdRows] = paste0("transcript:",S4Vectors::mcols(myGTFdata)$transcript_id[transcriptIdRows])
	} 
	
	n_occur <- data.frame(table(IDs))
	if(nrow(n_occur[n_occur$Freq > 1,])){
		warning(paste0("There are repeated IDs in the conversion:",n_occur[n_occur$Freq > 1,]))
	}
	

	S4Vectors::mcols(myGTFdata)$ID<-IDs
	ParentIDs<-apply(as.data.frame(S4Vectors::mcols(myGTFdata)), MARGIN = 1, function(x) { if (as.character(x["type"])=="gene"){NA} else if(as.character(x["type"]) =="transcript"){as.character(x["gene_id"])} else {as.character(x["transcript_id"])}})
	S4Vectors::mcols(myGTFdata)$Parent<-ParentIDs

	S4Vectors::mcols(myGTFdata)$gene_id<-NULL
	S4Vectors::mcols(myGTFdata)$transcript_id<-NULL
	
	# Replace special characters
	S4Vectors::mcols(myGTFdata)<-as.data.frame(S4Vectors::mcols(myGTFdata)) %>% dplyr::mutate(dplyr::across(names(S4Vectors::mcols(myGTFdata)),~ stringi::stri_replace_all_fixed(.x,c("=","&",","),c("%3D","%26","%2C"),vectorize_all=FALSE)))

	rtracklayer::export.gff(myGTFdata,foutput)

	# currentTime<-as.numeric(Sys.time())
	# elapsedTime<-currentTime - startTime;
	# message(paste0(gtfFile, "converted to ",foutput, " in ",formatC(elapsedTime %/% 60 %% 60, format = "d", flag = "0"),"m:",formatC(elapsedTime %% 60, width = 2, format = "d", flag = "0"),"s"))

	rm(list=setdiff(ls(), "foutput"))	
	gc(reset=T)
	
	return(foutput)
	
	
}

#' Test consistency and order of a gff file
#'
#' This function tests the consistency and order of a gff file.
#'
#'  The following list indicates the code and description of the issues detected in GFF3 files
#'  \describe{
#'  \item{NCOLUMNS_EXCEEDED}{Input file contains lines with more than 9 fields}
#'  \item{NCOLUMNS_INFERIOR}{Input file contains lines with less than 9 fields}
#'  \item{NO_IDs}{ID attribute not found in any feature}
#'  \item{DUPLICATED_IDs}{There are duplicated IDs}
#'  \item{ID_IN_MULTIPLE_CHR}{The same ID has been found in more than one chromosome}
#'  \item{NO_PARENTs}{Parent attribute not found in any feature}
#'  \item{MISSING_PARENT_IDs}{There are  missing Parent IDs}
#'  \item{PARENT_IN DIFFERENT CHR}{There are features whose Parent is located in a different chromosome}
#'  \item{PARENT_DEFINED_BEFORE_ID}{Feature ids referenced in Parent attribute before being defined as ID}
#'  \item{NOT_GROUPED_BY_CHR}{Features are not grouped by chromosome}
#'  \item{NOT_SORTED_BY_COORDINATE}{Features are not sorted by start coordinate}
#'  \item{NOT_VALID_WARNING}{File cannot be recognized as valid  GFF3. Parsing warnings.}
#'  \item{NOT_VALID_ERROR}{File cannot be recognized as valid GFF3. Parsing errors.}
#'  }
#'  The following list indicates the code and description of the issues detected in GTF files
#'  \describe{
#'  \item{NCOLUMNS_EXCEEDED}{Input file contains lines with more than 9 fields}
#'  \item{NCOLUMNS_INFERIOR}{Input file contains lines with less than 9 fields}
#'  \item{NO_GENE_ID_ATTRIBUTE}{gene_id attribute not found in any feature}
#'  \item{MISSING_GENE_IDs}{There are features without gene_id attribute}
#'  \item{NO_GENE_FEATURES}{Gene features are not included in this GTF file}
#'  \item{DUPLICATED_GENE_IDs}{There are duplicated gene_ids }
#'  \item{GENE_ID_IN_MULTIPLE_CHR}{The same gene_id has been found in more than one chromosome}
#'  \item{NO_TRANSCRIPT_ID_ATTRIBUTE}{transcript_id attribute not found in any feature There are no elements with transcript_id attribute}
#'  \item{MISSING_TRANSCRIPT_IDs}{There are  features without transcript_id attribute}
#'  \item{NO_TRANSCRIPT_FEATURES}{Transcript features are not included in this GTF file}
#'  \item{DUPLICATED_TRANSCRIPT_IDs}{There are  duplicated transcript_ids }
#'  \item{TRANSCRIPT_ID_IN_MULTIPLE_CHR}{The same transcript_id has been found in more than one chromosome}
#'  \item{DUPLICATED_GENE_AND_TRANSCRIPT_IDs}{Same id has been defined as gene_id and transcript_id}
#'  \item{NOT_GROUPED_BY_CHR}{Features are not grouped by chromosome}
#'  \item{NOT_SORTED_BY_COORDINATE}{Features are not sorted by start coordinate}
#'  \item{NOT_VALID_WARNING}{File cannot be recognized as valid GTF. Parsing warnings.}
#'  \item{NOT_VALID_ERROR}{File cannot be recognized as valid GTF. Parsing errors.}
#'  }
#' 
#' @param inFile Path to the input gff file
#' @param fileType Version of the input file (GTF/GFF3). Default AUTO: determined from the file name.
#' @return A data frame of detected issues, including a short code name, a description and estimated severity each. In no issues are detected the function will return an empty data frame.
#' 
#' @export
#' @examples
#' test_gff3<-system.file("extdata", "eden.gff3", package="Rgff")
#' check_gff(test_gff3)


check_gff<-function(inFile, fileType=c("AUTO","GFF3","GTF")){
	fileType <- match.arg(fileType)

	if (!base::file.exists(inFile)) { 
		stop(paste0("Input file ",inFile," not found."))
	}

	if(fileType == "AUTO"){
		detectedExt<-tools::file_ext(inFile)
		if(tolower(detectedExt) == "gtf"){
			fileExt="gtf"
			fileType = "GTF"
		} else if( (tolower(detectedExt) == "gff") || (tolower(detectedExt) == "gff3")) {
			fileExt=tolower(detectedExt)
			fileType = "GFF3"			
		} else {
			message(paste0("Unknown file type", detectedExt,": Defaulting to gff3"))
			fileExt="gff3"		
			fileType = "GFF3"				
		}
	}
	
	if(fileType == "GFF3"){		
		return(check_gff3(inFile))
	} else {
		return(check_gtf(inFile))
	}
}


#' Test consistency and order of a GFF3 file
#'
#' This function test consistency and order of a GFF3 file.
#' 
#' @param gffFile Path to the input GFF3 file
#' @return A data frame with the errors detected
#'
#' @importFrom rlang .data
#' @keywords internal

check_gff3<-function(gffFile){
	if (!base::file.exists(gffFile)) { 
		stop(paste0("Input GFF3 file ",gffFile," not found."))
	}
	
	options(dplyr.summarise.inform = FALSE)

	
	errorsDetected<-data.frame(ERRORCODE=character(), MESSAGE=character(), SEVERITY=character())
	errorsDetected<-tryCatch({
		nFields<-table(utils::count.fields(gffFile,sep="\t", comment.char = "#",quote=""))	

		if(max(as.numeric(names(nFields))) > 9){
			errText<-paste0("Input file contains lines with more than 9 fields")
			message(errText) 
			errorsDetected[nrow(errorsDetected)+1,]<-c("NCOLUMNS_EXCEEDED",errText,"HIGH")
		} else if(min(as.numeric(names(nFields))) < 9){
			errText<-paste0("Input file contains lines with less than 9 fields")
			message(errText) 
			errorsDetected[nrow(errorsDetected)+1,]<-c("NCOLUMNS_INFERIOR",errText,"HIGH")
		}
			
		cnames <- c("seqname","source","feature","start","end","score","strand","frame","attribute")
		myGFFdata <- utils::read.delim(gffFile, header=F,  sep="\t", comment.char = "#", quote="", blank.lines.skip=T) %>%    `colnames<-`(cnames)  

		myIDPattern<-"^(.*;)?( *ID *= *([^;]+))(.*)$"
		selected_id_lines<-grep(myIDPattern,myGFFdata$attribute,perl=T,ignore.case=T)
		if(length(selected_id_lines)==0){
			errText<-paste0("ID attribute not found in any feature")
			message(errText) 
			errorsDetected[nrow(errorsDetected)+1,]<-c("NO_IDs",errText,"HIGH")
		} else {
			resID <- data.frame(ChrID= myGFFdata[selected_id_lines,"seqname"], FeatureID=myGFFdata[selected_id_lines,"feature"], ID = gsub(myIDPattern,"\\3",myGFFdata[selected_id_lines,"attribute"],perl=T,ignore.case=T), PosID=selected_id_lines)

			### CHECK REPEATED IDs IN DIFFERENT FEATURE TYPES OR CHROMOSOMES

			idCounts <- resID %>% dplyr::group_by(.data$FeatureID,.data$ID) %>% dplyr::summarise(N=dplyr::n(), FirstID=min(.data$PosID))

			if(length(idCounts$ID) != length(unique(idCounts$ID))){
				errText<-paste0("There are ",length(idCounts$ID)-length(unique(idCounts$ID)), " duplicated IDs")
				message(errText) 
				errorsDetected[nrow(errorsDetected)+1,]<-c("DUPLICATED_IDs",errText,"HIGH")
			}

			if(any(duplicated(unique(resID[,c("ID","ChrID")])$ID))){
				errText<-paste0("The same  ",length(which(duplicated(unique(resID[,c("ID","ChrID")])$ID)))," ids have been found in more than one chromosome")
				message(errText) 
				errorsDetected[nrow(errorsDetected)+1,]<-c("ID_IN_MULTIPLE_CHR",errText,"HIGH")		
			}
		}
		

		myParentPattern<-"^(.*;)?( *Parent *= *([^;]+))(.*)$"
		selected_parent_lines<-grep(myParentPattern,myGFFdata$attribute,perl=T,ignore.case=T)

		if(length(selected_parent_lines)==0){
			errText<-paste0("Parent attribute not found in any feature")
			message(errText) 
			errorsDetected[nrow(errorsDetected)+1,]<-c("NO_PARENTs",errText,"MEDIUM")
		} else {
			resParent <- data.frame(ChrChild=myGFFdata[selected_parent_lines,"seqname"],Feature=myGFFdata[selected_parent_lines,"feature"], ID = gsub(myParentPattern,"\\3",myGFFdata[selected_parent_lines,"attribute"],perl=T,ignore.case=T), PosChild=selected_parent_lines)

			## INCLUDE MULTIPLE PARENT OCCURRENCES
			withMulti<-grep(",",resParent$ID)
			if(length(withMulti)>0){
				newItems<-as.data.frame(do.call("rbind",lapply(withMulti,function(x){t(sapply(strsplit(resParent$ID[x],",")[[1]],function(y){c(resParent$ChrChild[x],resParent$Feature[x],y,as.numeric(resParent$PosChild[x]))}))})))
				colnames(newItems) <- c("ChrChild","Feature", "ID","PosChild")
				newItems$PosChild<-as.numeric(newItems$PosChild)
				resParent<-rbind(resParent[-withMulti,],newItems)
			}

			resParent<-resParent %>% dplyr::group_by(.data$Feature,.data$ID,.data$ChrChild) %>% dplyr::summarise(Nchilden=dplyr::n(), FirstChild=min(.data$PosChild)) %>% dplyr::arrange(.data$ID)


			### MISSING PARENTS
			mergedDF<-merge(x = resParent, y = resID, by = "ID", all.x = TRUE)
			missingParents<-which(is.na(mergedDF$FeatureID))
			if(length(missingParents) >0){
				errText<-paste0("There are ",length(missingParents), " missing Parent IDs")
				message(errText) 
				errorsDetected[nrow(errorsDetected)+1,]<-c("MISSING_PARENT_IDs",errText,"HIGH")
			}

			if(any(mergedDF$ChrChild != mergedDF$ChrID, na.rm=T)){
				errText<-paste0("There are ",length(which(mergedDF$ChrChild != mergedDF$ChrID)), " features whose Parent located in a different chromosome")
				message(errText) 
				errorsDetected[nrow(errorsDetected)+1,]<-c("PARENT_IN DIFFERENT CHR",errText,"HIGH")
			}


			### PARENT BEFORE ID / FEATURE ORDER
			
			parentBeforeID<-which(mergedDF$FirstChild-mergedDF$Pos < 0)
			if(length(parentBeforeID) >0){				
				errText<-paste0(length(parentBeforeID), " feature ids referenced in Parent attribute before being defined as ID")
				message(errText) 
				errorsDetected[nrow(errorsDetected)+1,]<-c("PARENT_DEFINED_BEFORE_ID",errText,"MEDIUM")
			}
		}

		### CHROMOSOME ORDER
		chrOrder<-unique(myGFFdata$seqname)
		chrFactor<-factor(myGFFdata$seqname, levels=chrOrder)
		if(is.unsorted(chrFactor)){
			errText<-paste0("Features are not grouped by chromosome")
			message(errText) 
			errorsDetected[nrow(errorsDetected)+1,]<-c("NOT_GROUPED_BY_CHR",errText,"MEDIUM")
		}

		### COORDINATE (START) ORDER
		startSorted<-myGFFdata %>% dplyr::group_by(.data$seqname) %>% dplyr::summarise(UNSORTED=is.unsorted(.data$start, strictly=FALSE))
		if(any(startSorted$UNSORTED)){
			errText<-paste0("Features are not sorted by start coordinate in ",sum(startSorted$UNSORTED), " chromosomes")
			message(errText) 
			errorsDetected[nrow(errorsDetected)+1,]<-c("NOT_SORTED_BY_COORDINATE",errText,"MEDIUM")
		}
		if(nrow(errorsDetected) == 0){
			message(paste0(gffFile, " OK: No errors detected.")) 
		}
		return(errorsDetected)
	
	},
	warning=function(cond){
		errText<-paste0("File cannot be recognized as valid GFF3. Parsing warnings.")
		message(errText) 
		errorsDetected[nrow(errorsDetected)+1,]<-c("NOT_VALID_WARNING",errText,"HIGH")
		return(errorsDetected)
	},	
	error=function(cond){
		errText<-paste0("File cannot be recognized as valid GFF3. Parsing errors.")
		message(errText) 
		errorsDetected[nrow(errorsDetected)+1,]<-c("NOT_VALID_ERROR",errText,"HIGH")
		return(errorsDetected)
	})

	rm(list=setdiff(ls(), "errorsDetected"))	
	gc(reset=T)
	
	return(errorsDetected)
}




#' Test consistency and order of a GTF file
#'
#' This function tests consistency and order of a GTF file.
#' Optionally, can fix some of the consistency errors and sort the GTF file.
#' 
#' @param gtfFile Path to the input GTF file
#' @return A data frame with the errors detected
#'
#' @importFrom rlang .data
#' 
#' @keywords internal

check_gtf<-function(gtfFile){
	if (!base::file.exists(gtfFile)) { 
		stop(paste0("Input file ",gtfFile," not found."))
	}

	errorsDetected<-data.frame(ERRORCODE=character(), MESSAGE=character(), SEVERITY=character())
	errorsDetected<-tryCatch({
		nFields<-table(utils::count.fields(gtfFile,sep="\t", comment.char = "#",quote=""))	

		if(max(as.numeric(names(nFields))) > 9){
			errText<-paste0("Input file contains lines with more than 9 fields")
			message(errText) 
			errorsDetected[nrow(errorsDetected)+1,]<-c("NCOLUMNS_EXCEEDED",errText,"HIGH")
		} else if(min(as.numeric(names(nFields))) < 9){
			errText<-paste0("Input file contains lines with less than 9 fields")
			message(errText) 
			errorsDetected[nrow(errorsDetected)+1,]<-c("NCOLUMNS_INFERIOR",errText,"HIGH")
		}

		cnames <- c("seqname","source","feature","start","end","score","strand","frame","attribute")
		myGTFdata <- utils::read.delim(gtfFile, header=F,  sep="\t", comment.char = "#", quote="", blank.lines.skip=T) %>%    `colnames<-`(cnames)  


		### ELEMENTS WITHOUT GENE_ID
		geneIdPattern<-".* *gene_id +\"?([^\"]+)\"? *(;|$).*"
		selected_geneid_lines<-grep(geneIdPattern,myGTFdata$attribute,perl=T,ignore.case=T)
		
		if(length(selected_geneid_lines) == 0){
			errText<-paste0("gene_id attribute not found in any feature")
			message(errText) 
			errorsDetected[nrow(errorsDetected)+1,]<-c("NO_GENE_ID_ATTRIBUTE",errText,"MEDIUM")
		} else {
			if(length(selected_geneid_lines) < nrow(myGTFdata)){
				errText<-paste0("There are ",nrow(myGTFdata)-length(selected_geneid_lines), " features without gene_id attribute")
				message(errText) 
				errorsDetected[nrow(errorsDetected)+1,]<-c("MISSING_GENE_IDs",errText,"LOW")	
			}
			resGene <- data.frame(Chr=myGTFdata[selected_geneid_lines,"seqname"],Feature=myGTFdata[selected_geneid_lines,"feature"], geneID = gsub(geneIdPattern,"\\1",myGTFdata[selected_geneid_lines,"attribute"],perl=T,ignore.case=T), PosGene=selected_geneid_lines)
			addGenes<-!any(resGene$Feature=="gene")
			if(addGenes){
				### MISSING GENE FEATURES
				errText<-"Gene features are not included in this GTF file"
				message(errText) 
				errorsDetected[nrow(errorsDetected)+1,]<-c("NO_GENE_FEATURES",errText,"LOW")	
			} else {
				# DUPLICATED GENE_IDs IN gene FEATURES
				if(any(duplicated(resGene[resGene$feature=="gene",]$geneID))){
					errText<-paste0("There are ",length(which(duplicated(resGene[resGene$feature=="gene",]$geneID)))," duplicated gene_ids ")
					message(errText) 
					errorsDetected[nrow(errorsDetected)+1,]<-c("DUPLICATED_GENE_IDs",errText,"HIGH")
				}
				
				# DUPLICATED GENE_IDs IN DIFFERENT CHROMOSOMES
				if(any(duplicated((resGene  %>% dplyr::distinct(.data$geneID, .data$Chr))$geneID))){
					errText<-paste0("The same gene_id has been found in more than one chromosome")
					message(errText) 
					errorsDetected[nrow(errorsDetected)+1,]<-c("GENE_ID_IN_MULTIPLE_CHR",errText,"HIGH")
				}
			}
		}

		### ELEMENTS (EXCLUDING GENES)  WITHOUT TRASNCRIPT_ID
		transcriptIdPattern<-".* *transcript_id +\"?([^\"]+)\"? *(;|$).*"
		selected_transcriptid_lines<-grep(transcriptIdPattern,myGTFdata$attribute,perl=T,ignore.case=T)
		
		if(length(selected_transcriptid_lines) == 0){
			errText<-"transcript_id attribute not found in any feature"
			message(errText) 
			errorsDetected[nrow(errorsDetected)+1,]<-c("NO_TRANSCRIPT_ID_ATTRIBUTE",errText,"MEDIUM")
		} else {		
			if(length(selected_transcriptid_lines) < nrow(myGTFdata[myGTFdata$feature != "gene",]) ){
				errText<-paste0("There are  ",myGTFdata[myGTFdata$feature != "gene",]-length(selected_transcriptid_lines), " features without transcript_id attribute")
				message(errText) 
				errorsDetected[nrow(errorsDetected)+1,]<-c("MISSING_TRANSCRIPT_IDs",errText,"LOW")	
			} 
			resTranscript <- data.frame(Chr=myGTFdata[selected_transcriptid_lines,"seqname"], Feature=myGTFdata[selected_transcriptid_lines,"feature"], transcriptID = gsub(transcriptIdPattern,"\\1",myGTFdata[selected_transcriptid_lines,"attribute"],perl=T,ignore.case=T), PosTranscript=selected_transcriptid_lines)
			addTranscript<-!any(resTranscript$Feature=="transcript")
			if(addTranscript){
				### MISSING TRANSCRIPT FEATURES
				errText<-"Transcript features are not included in this GTF file"
				message(errText) 
				errorsDetected[nrow(errorsDetected)+1,]<-c("NO_TRANSCRIPT_FEATURES",errText,"LOW")	
			} else {
				# DUPLICATED TRANSCRIPT_IDs IN transcript FEATURES
				if(any(duplicated(resTranscript[resTranscript$feature=="transcript",]$transcriptID))){
					errText<-paste0("There are ",length(which(duplicated(resTranscript[resTranscript$feature=="transcript",]$transcriptID)))," duplicated transcript_ids ")
					message(errText) 
					errorsDetected[nrow(errorsDetected)+1,]<-c("DUPLICATED_TRANSCRIPT_IDs",errText,"HIGH")
				}
				
				# DUPLICATED TRANSCRIPT_IDs IN DIFFERENT CHROMOSOMES
				if(any(duplicated((resTranscript  %>% dplyr::distinct(.data$transcriptID, .data$Chr))$transcriptID))){
					errText<-paste0("The same transcript_id has been found in more than one chromosome")
					message(errText) 
					errorsDetected[nrow(errorsDetected)+1,]<-c("TRANSCRIPT_ID_IN_MULTIPLE_CHR",errText,"HIGH")
				}
			}
			if(length(selected_geneid_lines) != 0){
				# COMMON TRANSCRIPT_ID AND GENE_ID FEATURES
				commonGeneTranscript<-intersect(unique(resGene$geneID),unique(resTranscript$transcriptID))
				if(length(commonGeneTranscript)>0){
					errText<-paste0("The same ",length(commonGeneTranscript), " ids have been defined as gene_id and transcript_id")
					message(errText) 
					errorsDetected[nrow(errorsDetected)+1,]<-c("DUPLICATED_GENE_AND_TRANSCRIPT_IDs",errText,"MEDIUM")
				}
			}
		}


		fLevels<-unique(myGTFdata$feature)
		tIndex<-which(fLevels=="transcript")
		if(length(tIndex)>0){
			fLevels<-c(fLevels[tIndex],fLevels[-tIndex])
		}
		
		gIndex<-which(fLevels=="gene")
		if(length(gIndex)>0){
			fLevels<-c(fLevels[gIndex],fLevels[-gIndex])
		}


		### CHROMOSOME ORDER
		chrOrder<-unique(myGTFdata$seqname)
		chrFactor<-factor(myGTFdata$seqname, levels=chrOrder)
		if(is.unsorted(chrFactor)){
			errText<-paste0("Features are not grouped by chromosome")
			message(errText) 
			errorsDetected[nrow(errorsDetected)+1,]<-c("NOT_GROUPED_BY_CHR",errText,"MEDIUM")
		}

		### COORDINATE (START) ORDER
		startSorted<-myGTFdata %>% dplyr::group_by(.data$seqname) %>% dplyr::summarise(UNSORTED=is.unsorted(.data$start, strictly=FALSE))
		if(any(startSorted$UNSORTED)){
			errText<-paste0("Features are not sorted by start coordinate in ",sum(startSorted$UNSORTED), " chromosomes")
			message(errText) 
			errorsDetected[nrow(errorsDetected)+1,]<-c("NOT_SORTED_BY_COORDINATE",errText,"MEDIUM")
		}

		if(nrow(errorsDetected) == 0){
				message(paste0(gtfFile, " OK: No errors detected.")) 
		}

		rm(list=setdiff(ls(), "errorsDetected"))	
		gc(reset=T)
		return(errorsDetected)
		
		
	},
	warning=function(cond){
		errText<-paste0("File cannot be recognized as valid GTF. Parsing warnings.")
		message(errText) 
		errorsDetected[nrow(errorsDetected)+1,]<-c("NOT_VALID_WARNING",errText,"HIGH")
		return(errorsDetected)
	},	
	error=function(cond){
		errText<-paste0("File cannot be recognized as valid GTF. Parsing errors.")
		message(errText) 
		errorsDetected[nrow(errorsDetected)+1,]<-c("NOT_VALID_ERROR",errText,"HIGH")
		return(errorsDetected)
	})
	rm(list=setdiff(ls(), "errorsDetected"))	
	gc(reset=T)
	
	return(errorsDetected)
}
