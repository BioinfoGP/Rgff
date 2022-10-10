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
		stop(paste0("Input GFF3 file ",gffFile," not found."))
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
	withr::with_options(c(scipen = 999),utils::write.table(DF,file=foutput,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t", append=FALSE))

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
		stop(paste0("Input GFF3 file ",gffFile," not found."))
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
	message(paste0("Creating ",foutput," file for ",gffFile))


	cnames <- c("seqname","source","feature","start","end","score","strand","frame","attribute")
	myGFFdata <- utils::read.delim(gffFile, header=F,  sep="\t", comment.char = "#", quote="", blank.lines.skip=T) %>%    `colnames<-`(cnames)  

	idPattern<-"(^|;) *ID *= *([^;]+) *(;|$)"
	idPatternFull<-"^(.*;)* *ID *= *([^;]+) *(;|$).*"


	parentPattern<-"(^|;) *Parent *= *([^;]+) *(;|$)"
	parentPatternFull<-"^(.*;)* *Parent *= *([^;]+) *(;|$).*"


	selected_parentid_lines<-grep(parentPattern,myGFFdata[,"attribute"],perl=T,ignore.case=T)
	selected_elemid_lines<-grep(idPattern,myGFFdata[,"attribute"],perl=T,ignore.case=T)

	parentids<-rep(NA,nrow(myGFFdata))
	elemids<-rep(NA,nrow(myGFFdata))

	parentids[selected_parentid_lines]<-gsub(parentPatternFull,"\\2",myGFFdata[selected_parentid_lines,"attribute"],perl=T,ignore.case=T)
	elemids[selected_elemid_lines]<-gsub(idPatternFull,"\\2",myGFFdata[selected_elemid_lines,"attribute"],perl=T,ignore.case=T)

	elemids[elemids == ""] <- NA


	fLevels<-unique(myGFFdata$feature)
	featureCounter<-rep("",nrow(myGFFdata))
	for (cf in fLevels){
		fLines<-which(myGFFdata$feature == cf)
		featureCounter[fLines]<-c(1:length(fLines))
	}


	elemids[is.na(elemids)] <- paste0(myGFFdata[is.na(elemids),]$feature,":",myGFFdata[is.na(elemids),]$feature,featureCounter[is.na(elemids)])


	
	if(length(selected_parentid_lines) >0 && length(selected_parentid_lines) < length(parentids)){
		parentids[-selected_parentid_lines] <- NA
	}
	
	myGFFdata <- myGFFdata %>% dplyr::mutate(ID=elemids, ParentID=parentids)

	

	
	
	if(any(duplicated(unique(myGFFdata[,c("ID","seqname")])$ID))){
		toRemoveDF<-unique(myGFFdata[,c("ID","seqname")])
		toRemoveDF<-toRemoveDF[duplicated(toRemoveDF$ID),]
		
		repeatedIDs<-which(myGFFdata$ID %in% unique(toRemoveDF$ID))
	
		message(paste0("Warning, there are ",length(repeatedIDs), " rows with repeated IDs in different chromosomes. These elements and any reference to them will be removed in the conversion"))
		repeatedIDList<-unique(elemids[repeatedIDs])
		parentids[parentids %in% repeatedIDList] <- NA
		myGFFdata <- myGFFdata[-repeatedIDs,]
		elemids <- elemids[-repeatedIDs]
		parentids <- parentids[-repeatedIDs]
		
		 
	}

	
	FeatureIdAndParent<-data.frame(Feature=myGFFdata$feature,ID=elemids, ParentID=parentids)
	
	uniqueFeatureIdAndParent <- unique(FeatureIdAndParent %>% dplyr::group_by(.data$ID,.data$ParentID ))
	
	withMulti<-grep(",",uniqueFeatureIdAndParent$ParentID)
	if(length(withMulti)>0){	

		nParents<-stringi::stri_count_fixed(uniqueFeatureIdAndParent$ParentID[withMulti],",")
		newFeat<-uniqueFeatureIdAndParent[withMulti,][rep(seq_len(nrow(uniqueFeatureIdAndParent[withMulti,])), nParents+1), ]
		newParents<-unlist(strsplit(as.character(uniqueFeatureIdAndParent[withMulti,]$ParentID),","))
		newFeat$ParentID<-newParents
		uniqueFeatureIdAndParent<-rbind(uniqueFeatureIdAndParent[-withMulti,],newFeat)
	}


	mergedDF<-data.frame(Chr=myGFFdata$seqname,Start=myGFFdata$start,End=myGFFdata$end,Strand=myGFFdata$strand,Feature=myGFFdata$feature, ID=elemids)
	
	mergedDF <- mergedDF %>% dplyr::left_join(y = uniqueFeatureIdAndParent, by = c("ID","Feature"),na_matches="never")

	counter<-""
	while(any(!is.na(mergedDF[,paste0("ParentID",counter)]) & (mergedDF[,paste0("ParentID",counter)]!=""))){
		if(counter==""){
			nextCounter<-1
		} else{
			nextCounter<-counter+1
		}

		# print(paste0("Checking level ",nextCounter))

		varJoin<-paste0("ParentID",counter)
		joinExp= c("ID")
		names(joinExp) <- varJoin
		mergedDF <- mergedDF %>% dplyr::left_join(y = uniqueFeatureIdAndParent, by = joinExp, suffix=c("",as.character(nextCounter)),na_matches="never")
		if(FALSE){
			mergedDF<-merge(x = mergedDF, y = fullData, by.x = paste0("ParentID",counter), by.y = "ID", all.x = TRUE,suffixes=c("",as.character(nextCounter)))
		}
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

	pathList<-as.data.frame(t(apply(mergedDF,1,function(x){c(x["Chr"],as.numeric(x["Start"]),as.numeric(x["End"]),x["Strand"],paste(x[5:length(x)][!is.na(x[5:length(x)]) & x[5:length(x)]!=""],collapse=";"))})))


	message("Copying paths file") 
	
	withr::with_options(c(scipen = 999),utils::write.table(pathList,file=foutput,col.names=FALSE,quote=FALSE,sep="\t", row.names=FALSE,append=FALSE))

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
		searchLine=paste(paste0("\t",groupBy,";"),collapse="|")
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
	
	return (DF[order(DF$Chr,DF$GeneID,as.numeric(DF$Start)),])
}


#' Summarizes the number of elements of each type in each chromosome of a GFF file
#'
#' This function summarizes the number of features of each type in each chromosome  
#' of a GFF file and returns the statistics
#'
#' @param inFile Path to the input GFF file
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


#' Summarizes the number of features of each type in a GFF file
#'
#' This function summarizes the number of features of each type in  
#' a GFF file and returns the statistics
#'
#' @param inFile Path to the input GFF file
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
	safDF<-unique(safDF[order(safDF[,2], as.numeric(safDF[,3]), -as.numeric(safDF[,4]), safDF[,1] ),])

	withr::with_options(c(scipen = 999),utils::write.table(safDF,foutput,sep="\t",quote=FALSE,row.names = FALSE))

	rm(list=setdiff(ls(), "foutput"))	
	gc(reset=T)
	
	return(foutput)
}	


#' Plots or exports an image of the feature tree from a GFF file
#'
#' This function plots the feature tree from a GFF file or, if an output file name is provided, 
#' exports an image of in the desired format ("png", "pdf" or "svg"). 
#' Packages "DiagrammeR", "DiagrammeRsvg" and "rsvg" must be installed to use this function.
#'
#' @param inFile Path to the input GFF file 
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
		stop(paste0("Input GFF file ",inFile," not found."))
	}


	cleanInFile<-inFile
	
	detectedExt<-tools::file_ext(inFile)
	if((tolower(detectedExt) == "gz")){
		cleanInFile<-paste0(sub(pattern = "(.*)\\..{0,4}$", replacement = "\\1", inFile))
		message("Gzipped file detected.")
		detectedExt<-tools::file_ext(cleanInFile)
	}

	if(fileType == "AUTO"){
		if(tolower(detectedExt) == "gtf"){
			fileExt="gtf"
			fileType = "GTF"
		} else if( (tolower(detectedExt) == "gff") || (tolower(detectedExt) == "gff3")) {
			fileExt=tolower(detectedExt)
			fileType = "GFF3"			
		} else {
			# message(paste0("Unknown file type ", detectedExt,": Defaulting to GFF3"))
			# fileExt="gff3"
			# fileType = "GFF3"	
			message(paste0("Unknown file extension ", detectedExt,": Use a standard extension (\"gtf\",\"gff\",\"gff3\") or use the fileType argument to indicate the gff version (\"GTF\",\"GFF\")"))
			return ("")
			
		}
	} else if (fileType == "GTF") {
		fileExt="gtf"
		fileType = "GTF"
	} else {
		fileExt="gff3"
		fileType = "GFF3"				
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

			foutput=paste0(cleanInFile,".",exportFormat)
		}
		
		data.tree::ToDiagrammeRGraph(myTree) %>% DiagrammeR::export_graph(file_name=foutput, file_type=exportFormat)

		rm(list=setdiff(ls(), "foutput"))	
		gc(reset=T)
	
		return(foutput) 
	}

}



#' Analyses the feature type hierarchy of a GFF file 
#'
#' Based on the feature type hierarchy a GFF file, this function creates and returns 
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

	cleanInFile<-inFile
	
	detectedExt<-tools::file_ext(inFile)
	if((tolower(detectedExt) == "gz")){
		cleanInFile<-paste0(sub(pattern = "(.*)\\..{0,4}$", replacement = "\\1", inFile))
		message("Gzipped file detected.")
	}
	
	pairsFile=paste0(cleanInFile,".pairs")
	if (!base::file.exists(pairsFile)) {
		message ("Creating pairs file...")
		get_pairs_from_gff(inFile,fileType=fileType)
	}

	if(is.data.frame(pairsFile)){
		myPairData<-pairsFile
	} else if (is.vector(pairsFile) && is.character(pairsFile) && length(pairsFile)==1) {
		if(file.exists(pairsFile)){
			myPairData<-utils::read.table(pairsFile,header=T,sep="\t", stringsAsFactors=FALSE,colClasses=c("character","character","numeric"))
		} else {
			stop("Only Data frame of pairs or a valid path of pairs file are allowed")
		}
	}
	if(!all(c("ELEMENT","PARENT") %in% names(myPairData))){
		stop("Missing ELEMENT or PARENT columns in pair data")
	} 
	
	myPairData<-myPairData %>% dplyr::relocate(.data$PARENT, .before = .data$ELEMENT)
	
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
		stop(paste0("Input GFF3 file ",gffFile," not found."))
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
	withr::with_options(c(scipen = 999),utils::write.table(orderedGff, fileConn, append=TRUE, sep="\t", quote=F, row.names = FALSE, col.names=FALSE))

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

	different_features<-unique(myGTFdata$feature)
	if(length(different_features) > 100){
		validFeatures=c("gene", "transcript", "exon",  "CDS", "Selenocysteine", "start_codon", "stop_codon", "UTR", "5UTR", "3UTR", "inter", "inter_CNS", "intron_CNS")
		message(paste0("There are too many different feature types (",length(different_features),"). Only the following features (if present) will be shown."))
		message(paste(validFeatures, collapse=","))
		myGTFdata<-myGTFdata[myGTFdata$feature %in% validFeatures,]
	}

	featureList<- myGTFdata %>% dplyr::group_by(.data$feature) %>% dplyr::summarise(N = dplyr::n())
	parentElem<-sapply(featureList$feature,function(x){ if(x == "transcript") {"gene"} else { if ( x == "gene") {""} else {"transcript"}}})

	DF<-data.frame(ELEMENT=featureList$feature, PARENT= parentElem, N=featureList$N)


	# CHECK FOR MISSING TRANSCRIPT OR GENE PARENTS
	transcriptPattern<-" *transcript_id +\"?([^\"]+)\"? *(;|$)"
	transcriptPatternFull<-".* *transcript_id +\"?([^\"]+)\"? *(;|$).*"

	

	addTranscripts<-!any(grepl(transcriptPattern,myGTFdata[myGTFdata$feature == "transcript","attribute"],perl=T,ignore.case=T))
	if(addTranscripts){
		selected_transcriptid_lines<-grep(transcriptPattern,myGTFdata[,"attribute"],perl=T,ignore.case=T)
		transcriptIDs<-unique(gsub(transcriptPatternFull,"\\1",myGTFdata[selected_transcriptid_lines,"attribute"],perl=T,ignore.case=T))
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

	
	withr::with_options(c(scipen = 999),utils::write.table(DF,file=foutput,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t", append=FALSE))

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

	different_features<-unique(myGTFdata$feature)
	if(length(different_features) > 100){
		validFeatures=c("gene", "transcript", "exon",  "CDS", "Selenocysteine", "start_codon", "stop_codon", "UTR", "5UTR", "3UTR", "inter", "inter_CNS", "intron_CNS")
		message(paste0("There are too many different feature types (",length(different_features),"). Only the following features (if present) will be shown."))
		message(paste(validFeatures, collapse=","))
		myGTFdata<-myGTFdata[myGTFdata$feature %in% validFeatures,]
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


	# CHECK FOR MISSING GENE PARENTS

	allGenes<-sort(unique(myGTFdata[myGTFdata$feature != "gene",]$gene_id))
	if("" %in% allGenes){
		allGenes<-allGenes[-which(allGenes == "")]
	} 


	featuredGenes<-sort(unique(myGTFdata[myGTFdata$feature == "gene",]$gene_id))
	if("" %in% featuredGenes){
		featuredGenes<-featuredGenes[-which(featuredGenes == "")]
	} 

	addGenes= (length(allGenes) > length(featuredGenes))
	
	if(addGenes){
		missingGenes<-allGenes[which(!(allGenes %in% featuredGenes))]
		geneDF <- myGTFdata[myGTFdata$gene_id %in% missingGenes,] %>% dplyr::group_by(.data$gene_id,.data$seqname,.data$strand) %>% dplyr::summarise(start=min(.data$start), end=max(.data$end)) %>% dplyr::mutate(ASSIGNED=0)
	}

	# CHECK FOR MISSING TRANSCRIPT PARENTS

	allTranscripts<-sort(unique(myGTFdata[myGTFdata$feature != "gene" & myGTFdata$feature != "transcript",]$transcript_id))
	if("" %in% allTranscripts){
		allTranscripts<-allTranscripts[-which(allTranscripts == "")]
	} 

	featuredTranscripts<-sort(unique(myGTFdata[myGTFdata$feature == "transcript",]$transcript_id))
	if("" %in% featuredTranscripts){
		featuredTranscripts<-featuredTranscripts[-which(featuredTranscripts == "")]
	} 

	addTranscripts= (length(allTranscripts) > length(featuredTranscripts))

	if(addTranscripts){
		missingTranscripts<-allTranscripts[which(!(allTranscripts %in% featuredTranscripts))]
		transcriptDF<- myGTFdata[myGTFdata$feature != "gene" & (myGTFdata$transcript_id %in% missingTranscripts),] %>% dplyr::group_by(.data$transcript_id,.data$gene_id,.data$seqname,.data$strand)  %>% dplyr::summarise(start=min(.data$start), end=max(.data$end)) %>% dplyr::mutate(ASSIGNED=0)

	}

	## CREATE PATHS DATA  

	pathsData<-data.frame(seqname=character(0),start=numeric(0),end=numeric(0),strand=character(0),path=character(0), featureType=numeric(0))
	myGTFdataGene<-myGTFdata[myGTFdata$feature == "gene",]
	if(nrow(myGTFdataGene)>0){
		pathsData<-rbind(pathsData,data.frame(seqname=myGTFdataGene$seqname,start=as.numeric(myGTFdataGene$start),end=as.numeric(myGTFdataGene$end), strand=myGTFdataGene$strand,path=paste0("gene;",myGTFdataGene$gene_id),featureType=rep(0,nrow(myGTFdataGene))))
	}
	myGTFdataTranscript<-myGTFdata[myGTFdata$feature == "transcript",]	
	if(nrow(myGTFdataTranscript)>0){
		pathsData<-rbind(pathsData,data.frame(seqname=myGTFdataTranscript$seqname,start=as.numeric(myGTFdataTranscript$start),end=as.numeric(myGTFdataTranscript$end), strand=myGTFdataTranscript$strand,path=paste0("transcript;",myGTFdataTranscript$transcript_id,";gene;",myGTFdataTranscript$gene_id),featureType=rep(1,nrow(myGTFdataTranscript))))
	}

	myGTFdataOthers<-myGTFdata[myGTFdata$feature != "transcript" & myGTFdata$feature != "gene",]	
	if(nrow(myGTFdataOthers)>0){
		pathsData<-rbind(pathsData,pathsDataO<-data.frame(seqname=myGTFdataOthers$seqname,start=as.numeric(myGTFdataOthers$start),end=as.numeric(myGTFdataOthers$end), strand=myGTFdataOthers$strand,path=paste0(myGTFdataOthers$feature,";",myGTFdataOthers$feature,"_",myGTFdataOthers$fcounter,";transcript;",myGTFdataOthers$transcript_id,";gene;",myGTFdataOthers$gene_id),featureType=rep(2,nrow(myGTFdataOthers))))
	}


	if (addGenes){
		pathsData<-rbind(pathsData,data.frame(seqname=geneDF$seqname,start=as.numeric(geneDF$start),end=as.numeric(geneDF$end), strand=geneDF$strand,path=paste0("gene;",geneDF$gene_id),featureType=rep(0,nrow(geneDF))))
	}

	if (addTranscripts){
		pathsData<-rbind(pathsData,data.frame(seqname=transcriptDF$seqname,start=as.numeric(transcriptDF$start),end=as.numeric(transcriptDF$end), strand=transcriptDF$strand,path=paste0("transcript;",transcriptDF$transcript_id,";gene;",transcriptDF$gene_id),featureType=rep(1,nrow(transcriptDF))))
	}
	
	# SORT AND WRITE TO FILE
	orderedPathsData<- pathsData %>% dplyr::arrange(.data$seqname, .data$start, dplyr::desc(.data$end), .data$featureType)
	withr::with_options(c(scipen = 999),utils::write.table(orderedPathsData[,c(1:5)], foutput, sep="\t", quote=F, row.names = FALSE, col.names=FALSE))

	

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
	withr::with_options(c(scipen = 999),utils::write.table(orderedGTFdata, fileConn, append=TRUE, sep="\t", quote=F, row.names = FALSE, col.names=FALSE))

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
#' This function produces a sorted GFF file from an unsorted GFF file.
#' The default order is by Chromosome, Start, End (reverse) and feature (based on the precedency in feature tree)
#'
#' @param inFile Path to the input GFF file
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
	
	cleanInFile<-inFile
	
	detectedExt<-tools::file_ext(inFile)
	if((tolower(detectedExt) == "gz")){
		cleanInFile<-paste0(sub(pattern = "(.*)\\..{0,4}$", replacement = "\\1", inFile))
		message("Gzipped file detected.")
		detectedExt<-tools::file_ext(cleanInFile)
	}

	if(fileType == "AUTO"){
		if(tolower(detectedExt) == "gtf"){
			fileExt="gtf"
			fileType = "GTF"
		} else if( (tolower(detectedExt) == "gff") || (tolower(detectedExt) == "gff3")) {
			fileExt=tolower(detectedExt)
			fileType = "GFF3"			
		} else {
			# message(paste0("Unknown file type ", detectedExt,": Defaulting to GFF3"))
			# fileExt="gff3"
			# fileType = "GFF3"	
			message(paste0("Unknown file extension ", detectedExt,": Use a standard extension (\"gtf\",\"gff\",\"gff3\") or use the fileType argument to indicate the gff version (\"GTF\",\"GFF\")"))
			return ("")
			
		}
	} else if (fileType == "GTF") {
		fileExt="gtf"
		fileType = "GTF"
	} else {
		fileExt="gff3"
		fileType = "GFF3"				
	}

	if(missing(outFile)){
		foutput=paste0(sub(pattern = "(.*)\\..{0,4}$", replacement = "\\1", cleanInFile),".sorted.",fileExt)
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
#' @param inFile Path to the input GFF file
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


	cleanInFile<-inFile
	
	detectedExt<-tools::file_ext(inFile)
	if((tolower(detectedExt) == "gz")){
		cleanInFile<-paste0(sub(pattern = "(.*)\\..{0,4}$", replacement = "\\1", inFile))
		message("Gzipped file detected.")
		detectedExt<-tools::file_ext(cleanInFile)
	}

	if(fileType == "AUTO"){
		if(tolower(detectedExt) == "gtf"){
			fileExt="gtf"
			fileType = "GTF"
		} else if( (tolower(detectedExt) == "gff") || (tolower(detectedExt) == "gff3")) {
			fileExt=tolower(detectedExt)
			fileType = "GFF3"			
		} else {
			# message(paste0("Unknown file type ", detectedExt,": Defaulting to GFF3"))
			# fileExt="gff3"
			# fileType = "GFF3"	
			message(paste0("Unknown file extension ", detectedExt,": Use a standard extension (\"gtf\",\"gff\",\"gff3\") or use the fileType argument to indicate the gff version (\"GTF\",\"GFF\")"))
			return ("")
			
		}
	} else if (fileType == "GTF") {
		fileExt="gtf"
		fileType = "GTF"
	} else {
		fileExt="gff3"
		fileType = "GFF3"				
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
		foutput=paste0(cleanInFile,fSuffix,".saf")
	} else {
		foutput=outFile
	}


	if(base::file.exists(foutput) && forceOverwrite==FALSE){
		message("Output SAF file exists. Use forceOverwrite=TRUE to overwrite the existing file.")
		return(foutput)
	}
	
	PATHSfile=paste0(cleanInFile,".paths")
	if (!base::file.exists(PATHSfile)) {
		message ("Creating PATHS file... please be patient...")
		if(fileType == "GFF3"){
			gff3_to_paths(inFile, outFile=PATHSfile)
		} else {
			gtf_to_paths(inFile, outFile=PATHSfile)			
		}
	}

	

	for (i in c(1:length(featureList))){
		if(blockList[i] == ""){
			featureParam=c(featureList[i])
			blockParam=c()
		} else {
			blockParam=c(blockList[i])
			featureParam=c(featureList[i])
		}
		if(i==1){
			safDF<-unique(saf_from_paths(PATHSfile,groupBy=featureParam,block=blockParam))
		} else {

			safDF<-rbind(safDF,unique(saf_from_paths(PATHSfile,groupBy=featureParam,block=blockParam)))
		}
	}	
	safDF<-unique(safDF[order(safDF[,2], as.numeric(safDF[,3]), -as.numeric(safDF[,4]), safDF[,1] ),])

	withr::with_options(c(scipen = 999),utils::write.table(safDF,foutput,sep="\t",quote=FALSE,row.names = FALSE))

	rm(list=setdiff(ls(), "foutput"))	
	gc(reset=T)
	
	return(foutput)
}	


#' Creates a features pair file from a GFF file
#'
#' This function creates a features pair file from a GFF file
#'
#' @param inFile Path to the input GFF file
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


	cleanInFile<-inFile
	
	detectedExt<-tools::file_ext(inFile)
	if((tolower(detectedExt) == "gz")){
		cleanInFile<-paste0(sub(pattern = "(.*)\\..{0,4}$", replacement = "\\1", inFile))
		message("Gzipped file detected.")
		detectedExt<-tools::file_ext(cleanInFile)
	}

	if(fileType == "AUTO"){
		if(tolower(detectedExt) == "gtf"){
			fileExt="gtf"
			fileType = "GTF"
		} else if( (tolower(detectedExt) == "gff") || (tolower(detectedExt) == "gff3")) {
			fileExt=tolower(detectedExt)
			fileType = "GFF3"			
		} else {
			# message(paste0("Unknown file type ", detectedExt,": Defaulting to GFF3"))
			# fileExt="gff3"
			# fileType = "GFF3"	
			message(paste0("Unknown file extension ", detectedExt,": Use a standard extension (\"gtf\",\"gff\",\"gff3\") or use the fileType argument to indicate the gff version (\"GTF\",\"GFF\")"))
			return ("")
			
		}
	} else if (fileType == "GTF") {
		fileExt="gtf"
		fileType = "GTF"
	} else {
		fileExt="gff3"
		fileType = "GFF3"				
	}

	if(missing(outFile)){
		foutput=paste0(cleanInFile,".pairs")
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
		return(get_pairs_from_gff3(inFile, foutput, forceOverwrite))
	} else {
		return(get_pairs_from_gtf(inFile, foutput, forceOverwrite))
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
		cleanGtfFile<-gtfFile
		detectedExt<-tools::file_ext(gtfFile) 
		if((tolower(detectedExt) == "gz")){
			cleanGtfFile<-paste0(sub(pattern = "(.*)\\..{0,4}$", replacement = "\\1", gtfFile))
		}	
		
		foutput=paste0(cleanGtfFile,".gff3")
	} else {
		foutput=outFile
	}

	if(base::file.exists(foutput) && forceOverwrite==FALSE){
		message("Output GFF3 file exists. Use forceOverwrite=TRUE to overwrite the existing file")
		return(foutput)
	}


	# startTime <- as.numeric(Sys.time())

	myGTFdata<-rtracklayer::import.gff2(gtfFile)


	if("" %in% myGTFdata$gene_id){
		empty_gene_ids<-which(myGTFdata$gene_id == "")
		warning(paste0("There are ",length(empty_gene_ids)," empty gene_id values and will be ignored"))
		myGTFdata[empty_gene_ids,]$gene_id<-NA
	}

	if("" %in% myGTFdata$transcript_id){
		empty_transcript_ids<-which(myGTFdata$transcript_id == "")
		warning(paste0("There are ",length(empty_transcript_ids)," empty transcript_id values and will be ignored"))
		myGTFdata[empty_transcript_ids,]$transcript_id<-NA
	}


	############# CHECK MISSING TRANSCRIPTS
	allTranscripts<-sort(unique(myGTFdata[myGTFdata$type != "gene" & myGTFdata$type != "transcript",]$transcript_id))
	featuredTranscripts<-sort(unique(myGTFdata[myGTFdata$type == "transcript",]$transcript_id))
	addTranscripts= (length(allTranscripts) > length(featuredTranscripts))

	if(addTranscripts){
		missingTranscripts<-allTranscripts[which(!(allTranscripts %in% featuredTranscripts))]


		TranscriptsCoords <- as.data.frame(myGTFdata[myGTFdata$type != "gene" & myGTFdata$transcript_id %in% missingTranscripts]) %>% dplyr::group_by(.data$transcript_id, .data$gene_id ,.data$seqnames, .data$strand, .data$source) %>% dplyr::summarise(start = min(.data$start), end = max(.data$end))
		TranscriptsCoords$type<-"transcript"

		

		TranscriptsCoordsGR<-GenomicRanges::makeGRangesFromDataFrame(TranscriptsCoords,
						 keep.extra.columns=TRUE,
						 ignore.strand=FALSE,
						 seqinfo=NULL,
						 seqnames.field=c("seqnames"),
						 start.field="start",
						 end.field="end",
						 strand.field="strand")

		myGTFdata<-S4Vectors::append(myGTFdata,TranscriptsCoordsGR)

	}

	############# CHECK MISSING GENES
	allGenes<-sort(unique(myGTFdata[myGTFdata$type != "gene",]$gene_id))
	
	featuredGenes<-sort(unique(myGTFdata[myGTFdata$type == "gene",]$gene_id))
	addGenes= (length(allGenes) > length(featuredGenes))
	
	if(addGenes){
		missingGenes<-allGenes[which(!(allGenes %in% featuredGenes))]

		GenesCoords <- as.data.frame(myGTFdata[myGTFdata$gene_id %in% missingGenes,]) %>% dplyr::group_by(.data$gene_id,.data$seqnames, .data$strand, .data$source) %>% dplyr::summarise(start = min(.data$start), end = max(.data$end))
		GenesCoords$type<-"gene"
		


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


	# SORT DATA
	tIndex<-which(levels(myGTFdata$type)=="transcript")
	newLevels<-c("transcript",levels(myGTFdata$type)[-tIndex])
	gIndex<-which(newLevels=="gene")
	newLevels<-c("gene",newLevels[-gIndex])

	
	orderedGff<- as.data.frame(myGTFdata) %>% dplyr::arrange(.data$seqnames, .data$start, dplyr::desc(.data$end), factor(.data$type, levels = newLevels))


	myGTFdata<-GenomicRanges::makeGRangesFromDataFrame(orderedGff,
					 keep.extra.columns=TRUE,
					 ignore.strand=FALSE,
					 seqnames.field=c("seqnames"),
					 start.field="start",
					 end.field="end",
					 strand.field="strand")
			


	# ADD MISSING IDs
	fLevels<-unique(myGTFdata$type)
	featureCounter<-rep("",length(myGTFdata))
	cFeatures<-fLevels[fLevels != "transcript" & fLevels != "gene"]

	for (cf in cFeatures){
		fLines<-which(myGTFdata$type == cf)
		featureCounter[fLines]<-c(1:length(fLines))
	}
	S4Vectors::mcols(myGTFdata)$featureCounter<-featureCounter

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
		warning(paste0("There are ",nrow(n_occur[n_occur$Freq > 1,])," repeated IDs after the conversion. ID: ",n_occur[n_occur$Freq > 1,]$IDs))
	}
	

	S4Vectors::mcols(myGTFdata)$ID<-IDs
	ParentIDs<-apply(as.data.frame(S4Vectors::mcols(myGTFdata)), MARGIN = 1, function(x) { if (as.character(x["type"])=="gene"){NA} else if(as.character(x["type"]) =="transcript"){as.character(x["gene_id"])} else {as.character(x["transcript_id"])}})
	S4Vectors::mcols(myGTFdata)$Parent<-ParentIDs

	S4Vectors::mcols(myGTFdata)$gene_id<-NULL
	S4Vectors::mcols(myGTFdata)$transcript_id<-NULL
	
	# Replace special characters
	S4Vectors::mcols(myGTFdata)<-as.data.frame(S4Vectors::mcols(myGTFdata)) %>% dplyr::mutate(dplyr::across(names(S4Vectors::mcols(myGTFdata)),~ stringi::stri_replace_all_fixed(.x,c("=","&",","),c("%3D","%26","%2C"),vectorize_all=FALSE)))

	withr::with_options(c(scipen = 999),rtracklayer::export.gff(myGTFdata,foutput))


	rm(list=setdiff(ls(), "foutput"))	
	gc(reset=T)
	
	return(foutput)
	
	
}

#' Test consistency and order of a GFF file
#'
#' This function tests the consistency and order of a GFF file.
#'
#'  The following list indicates the code and description of the issues detected in GFF3 files
#'  \describe{
#'  \item{NCOLUMNS_EXCEEDED}{Input file contains lines with more than 9 fields}
#'  \item{NCOLUMNS_INFERIOR}{Input file contains lines with less than 9 fields}
#'  \item{TOO_MANY_FEATURE_TYPES}{Input file contains too many (more than 100) different feature types}
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
#'  \item{TOO_MANY_FEATURE_TYPES}{Input file contains too many (more than 100) different feature types}
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
#' @param inFile Path to the input GFF file
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
	cleanInFile<-inFile
	
	detectedExt<-tools::file_ext(inFile) 
	if((tolower(detectedExt) == "gz")){
		cleanInFile<-paste0(sub(pattern = "(.*)\\..{0,4}$", replacement = "\\1", inFile))
		message("Gzipped file detected.")
		detectedExt<-tools::file_ext(cleanInFile)
	}

	if(fileType == "AUTO"){
		if(tolower(detectedExt) == "gtf"){
			fileExt="gtf"
			fileType = "GTF"
		} else if( (tolower(detectedExt) == "gff") || (tolower(detectedExt) == "gff3")) {
			fileExt=tolower(detectedExt)
			fileType = "GFF3"			
		} else {
			# message(paste0("Unknown file type ", detectedExt,": Defaulting to GFF3"))
			# fileExt="gff3"
			# fileType = "GFF3"	
			message(paste0("Unknown file extension ", detectedExt,": Use a standard extension (\"gtf\",\"gff\",\"gff3\") or use the fileType argument to indicate the gff version (\"GTF\",\"GFF\")"))
			return ("")
			
		}
	} else if (fileType == "GTF") {
		fileExt="gtf"
		fileType = "GTF"
	} else {
		fileExt="gff3"
		fileType = "GFF3"				
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

		different_features<-unique(myGFFdata$feature)
		if(length(different_features) > 100){
			errText<-paste0("There are too many different feature types (",length(different_features),")")
			message(errText) 
			errorsDetected[nrow(errorsDetected)+1,]<-c("TOO_MANY_FEATURE_TYPES",errText,"MEDIUM")		
		}


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


		### TOO MANY FEATURE TYPES
		different_features<-unique(myGTFdata$feature)
		if(length(different_features) > 100){
			errText<-paste0("There are too many different feature types (",length(different_features),")")
			message(errText) 
			errorsDetected[nrow(errorsDetected)+1,]<-c("TOO_MANY_FEATURE_TYPES",errText,"MEDIUM")		
		}

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
