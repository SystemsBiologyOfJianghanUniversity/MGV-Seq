#----------------------------------------------#
#--------parsing parameters
#----------------------------------------------#
args=(commandArgs(TRUE))

if(length(args)!=2)
{
	stop("Rscript genotypeFromAmpliseq.R allele_file amplicon.bed")
}

f = args[1] #as the input files
f_target_Bed = args[2] #as the target region (in bed format)

if(!file.exists(f))
{
	stop(paste(f, "does not exist!",sep=""))
}

if(!file.exists(f_target_Bed))
{
	stop(paste(f_target_Bed," does not exists!",sep=""))
}
#----------------------------------------#
#---Preparation
#----------------------------------------#
library(igraph,quietly=F)
#Read in amplicon information.

ampList = try(read.table(f_target_Bed,as.is=T),silent=T)
if(class(ampList)=="try-error")
{
	stop(paste("Cannot read in ",f_target_Bed,sep=""))
}


#----------------------------------------#
#---Function
#----------------------------------------#

getStutter <- function(a)
{
	a$stutter = a$V4/a$V4[1]
	return(a)
}


getDiffNT <- function(tarAllele,refAllele)
{
	ignoreNT = c("N","-")
	tarAllele = strsplit(tarAllele,"")[[1]]
	refAllele = strsplit(refAllele,"")[[1]]
	index=which(refAllele!=tarAllele & (!refAllele %in% ignoreNT) & (!tarAllele %in% ignoreNT))
	tmp = as.data.frame(cbind(index,tarAllele[index]),stringsAsFactors=F)
	tmp$index = as.numeric(tmp$index)
	return(tmp)
	
}
getPos <-function(tmpGeno)
{
	return(sort(as.numeric(strsplit(tmpGeno,"[ATCG]")[[1]])))
}
maskConsecutive <- function(tmpGeno)
{
	tarNumBPConsecutiveSNP=2
	if(tmpGeno=="") return(c(0,tmpGeno))
	tmpPos = getPos(tmpGeno)
	if(length(tmpPos)==1) return(c(1,tmpGeno))
	tmpNumDiff = length(tmpPos)
	if(min(diff(tmpPos))<=tarNumBPConsecutiveSNP)
	{
		tmpNT = strsplit(tmpGeno,"\\d+")[[1]][-1]
		index = which(diff(tmpPos)<=tarNumBPConsecutiveSNP)
		toRemove = unique(c(index,index+1))
		tmpGeno = paste(tmpPos[-toRemove],tmpNT[-toRemove],sep="",collapse="")
		tmpNumDiff = length(tmpPos[-toRemove])
	}
	return(c(tmpNumDiff,tmpGeno))
}

getAlleleDiffNT_reorder <- function(a)
{
	#----compare with main Allele 
	tmpNTS = lapply(a$V3,getDiffNT,a$V3[1])
	tmpNTS_str = lapply(tmpNTS,function(a){paste(a$index,a$V2,sep="")})
	
	a$numDiffNT  = unlist(lapply(tmpNTS,nrow))
	a$DiffNT = gsub(" ","",unlist(lapply(tmpNTS_str,paste,collapse="")))
	
	#----reorder
	tab = sort(tapply(a$V4,a$DiffNT,sum),decreasing=T)
	tarMainGeno = names(tab)[1]
	if(tarMainGeno!="")
	{
		tarMainAllele = a$V3[a$DiffNT==tarMainGeno][1]
		tmpNTS = lapply(a$V3,getDiffNT,tarMainAllele)
		tmpNTS_str = lapply(tmpNTS,function(a){paste(a$index,a$V2,sep="")})
		a$numDiffNT  = unlist(lapply(tmpNTS,nrow))
		a$DiffNT = gsub(" ","",unlist(lapply(tmpNTS_str,paste,collapse="")))
	
	}
	
	return(a)
}



refineAlleles <- function(a)
{
	
	
	if(max(a$numDiffNT2)<2) return(NULL)
	tmpNTS = do.call(rbind,strsplit(a$V3,""))
	sigDiffNTList = unique(a$DiffNT2)
	sigDiffNTList = sigDiffNTList[sigDiffNTList!=""]
	if(length(sigDiffNTList)==1) return(NULL)
	tmpDat = NULL
	for(i in 1:(length(sigDiffNTList)-1))
	{
		tarAllele1 = sigDiffNTList[i]
		tarNTS1 = getPos(tarAllele1)
		
		for(j in (i+1):length(sigDiffNTList))
		{
			tarAllele2 = sigDiffNTList[j]
			tarNTS2 = getPos(tarAllele2)
			if(sum(tarNTS1 %in% tarNTS2)>0)
			{
				tmpDat = rbind(tmpDat,c(tarAllele1,tarAllele2))
			}
			
		}
	}
	if(is.null(tmpDat)) return(NULL)
	g = graph.data.frame(as.data.frame(tmpDat),directed=F)
	g_cluster = clusters(g)
	tmpRes = list()
	for(i in 1:g_cluster$no)
	{
		tmpList = names(g_cluster$membership)[g_cluster$membership==i]
		tmpPos = sort(unique(unlist(lapply(tmpList,getPos))))
		if(length(tmpPos)==1) next
		tmpNTS1 = tmpNTS[,tmpPos]
		tmpNTS1[tmpNTS1=="-"]="N"
		tmpNTS1_str = apply(tmpNTS1,1,paste,collapse="")
		tmp = list(tmpList,tmpNTS1_str)
		tmpRes[[length(tmpRes)+1]] = tmp
	}
	if(length(tmpRes)==0) return(NULL)
	return(tmpRes)
}
getRefinedAlleles <- function(tarMinorGeno,tmpRefinedAlleles)
{
	tmpNTS_str = NULL
	if(!is.null(tmpRefinedAlleles))
	{
		tmpNTS_str = unlist(lapply(tmpRefinedAlleles,function(b){if(tarMinorGeno %in% b[[1]]){b[[2]]}}))
	}
	return(tmpNTS_str)
}

	
#-----extract both main and minor alleles
genotypeMNP <- function(a)
{

	
	tarStrandFDR=0.05
	tarBiasFold=10
	tarBiasFoldDiff=5
	isRefined=T
	tarMainBiasFold=10
	tarNumBPConsecutiveSNP=2
	errorRateDat = data.frame(V1=c(1,2),V2=c(1.03E-02,9.94E-04))
	#----NT threshold
	maxNT = max(errorRateDat$V1)
	
	#----- Total number of reads
	numTotal = sum(a$V4)
	
	if(numTotal < 100) return(list(NULL,NULL))
	
	numTotalForward = sum(a$V5)
	numTotalReverse = numTotal - numTotalForward
	
	
	
	#----- check number of reads on forward strand
	if(!"V5" %in% names(a))
	{
		cat("Number of forward reads should be provided.\n")
		return(NULL)
	}
	#------ check if diffNT has been reported
	if(!"numDiffNT" %in% names(a) | !"DiffNT" %in% names(a))
	{
		a = getAlleleDiffNT_reorder(a)
	}
	
	#------- check if diffNT2 has been reported
	if(!"numDiffNT2" %in% names(a) | !"DiffNT2" %in% names(a))
	{
		a$numDiffNT2 = a$numDiffNT
		a$DiffNT2 = a$DiffNT
		a[,c("numDiffNT2","DiffNT2")]=do.call(rbind,lapply(a$DiffNT,maskConsecutive,tarNumBPConsecutiveSNP))
		a$numDiffNT2 = as.numeric(a$numDiffNT2)
	}
	
	#---- Reorder alleles if it's necessary
	tab = sort(tapply(a$V4,a$DiffNT2,sum),decreasing=T)
	if(names(tab)[1]!="")
	{
		a = getAlleleDiffNT_reorder(a)
		a$numDiffNT2 = a$numDiffNT
		a$DiffNT2 = a$DiffNT
		a[,c("numDiffNT2","DiffNT2")]=do.call(rbind,lapply(a$DiffNT,maskConsecutive,tarNumBPConsecutiveSNP))
		a$numDiffNT2 = as.numeric(a$numDiffNT2)
	}
	tab = sort(tapply(a$V4,a$DiffNT2,sum),decreasing=T)
	tmpMinorGeno = names(tab)[-1]
	
	#--- Refine alleles
	if(isRefined==T)
	{
		tmpRefinedAlleles = refineAlleles(a)    
	}else
	{
		tmpRefinedAlleles = NULL
	}
	
	
	#---- check mail allele
	b = a[a$numDiffNT2==0,]
	majorAllele = b[1,]
	
	# adjust the major allele's read count
	numTotalMajor = sum(b$V4)
	numTotalMajorForward = sum(b$V5)
	numTotalMajorReverse = numTotalMajor-numTotalMajorForward
	
	
	#  Strand bias of major allele
	majorStrandBias = numTotalMajorForward/numTotalMajorReverse
	majorAllele$V4=numTotalMajor
	majorAllele$V5=numTotalMajorForward
	
	majorStrandBias = numTotalMajorForward/numTotalMajorReverse
	
	absMajorStrandBias = 10^abs(log10(majorStrandBias))
	if(absMajorStrandBias > tarMainBiasFold) 
	{
		return(list(NULL,NULL))
	}
	
	
	
	#---- check each minor alleles
	tmpPHA = NULL
	if(length(tmpMinorGeno)>0)
	{
		for(k in 1:length(tmpMinorGeno))
		{
			tarMinorGeno = tmpMinorGeno[k]
			tmpNumDiffNT = length(strsplit(tarMinorGeno,"[ATCG]")[[1]])
			
			tmpNTS_str = getRefinedAlleles(tarMinorGeno,tmpRefinedAlleles)
			if(!is.null(tmpNTS_str))
			{
				tmpIndex = which(a$DiffNT2==tarMinorGeno)
				tmpIndex2 = grep("N",tmpNTS_str)
				if(length(tmpIndex2)>0) tmpIndex = tmpIndex[!tmpIndex %in% tmpIndex2]
				if(length(tmpIndex)==0) next
				numTotalMinor = sum(a$V4[tmpIndex])
				numTotalMinorForward = sum(a$V5[tmpIndex])
				
				
			}else
			{
				numTotalMinor = sum(a$V4[a$DiffNT2==tarMinorGeno])
				numTotalMinorForward = sum(a$V5[a$DiffNT2==tarMinorGeno])
				
			}
			numTotalMinorReverse = numTotalMinor-numTotalMinorForward
			
			minorStrandBias = numTotalMinorForward/numTotalMinorReverse
			absminorStrandBias = 10^abs(log10(minorStrandBias))
			
			if(absminorStrandBias<=tarBiasFold & minorStrandBias/majorStrandBias <= tarBiasFoldDiff & minorStrandBias/majorStrandBias >= 1/tarBiasFoldDiff)
			{
				
				numTotalNoTargetMinor = numTotal-numTotalMinor
				numTotalNoTargetMinorForward = numTotalForward-numTotalMinorForward
				
				numTotalNoTargetMinorReverse = numTotalNoTargetMinor-numTotalNoTargetMinorForward
				tmpMat = matrix(c(numTotalMinorForward,numTotalMinorReverse,numTotalNoTargetMinorForward,numTotalNoTargetMinorReverse),nrow=2)
				ft = fisher.test(tmpMat)
				
				
				P_error = 1
				if(tmpNumDiffNT<maxNT)
				{
					P_error = 1-pbinom(numTotalMinor-1,numTotal,errorRateDat$V2[errorRateDat$V1==tmpNumDiffNT])
				}else
				{
					P_error = 1-pbinom(numTotalMinor-1,numTotal,errorRateDat$V2[errorRateDat$V1==maxNT])
				}
				
				
				tmpRes = a[a$DiffNT2==tarMinorGeno,c("V1","V4","V5","numDiffNT","DiffNT","numDiffNT2","DiffNT2")][1,]
				tmpRes$V4 = numTotalMinor
				tmpRes$V5 = numTotalMinorForward
				tmpRes$numTotal = numTotal
				tmpRes$numTotalForward = numTotalForward
				tmpRes$strandP = ft$p.value
				tmpRes$p_error = P_error
				tmpPHA = rbind(tmpPHA,tmpRes)
				
			}
		}
	}
	
	if(!is.null(tmpPHA))
	{
		tmpPHA$strandFDR = p.adjust(tmpPHA$strandP)
		tmpPHA = tmpPHA[tmpPHA$strandFDR >= tarStrandFDR,]
		if(nrow(tmpPHA)==0)
		{
			tmpPHA = NULL
		}
		
	}
	return(list(majorAllele,tmpPHA))
}


#--------------------------------------#
#-----Genotype
#--------------------------------------#
	#----output file
	outF = paste(f,"_minor",sep="")#for minor allele
    outF2 = paste(f,"_major",sep="")#for major allele
    
	#----process input alleles
	x <- read.table(f,as.is=T)
	x = x[x$V1 %in% ampList$V4,] #only blight amplicons
	x$V3 = gsub(",I\\S+","",x$V3)
	tmpAmpTotal = tapply(x$V4,x$V1,sum)

	y = split(x,x$V1)
	z <- lapply(y,getAlleleDiffNT_reorder)
	
	z_dat = do.call(rbind,z)
	z_dat$numDiffNT2=-1
	z_dat$DiffNT2=""
	z_dat[,c("numDiffNT2","DiffNT2")]=do.call(rbind,lapply(z_dat$DiffNT,maskConsecutive))
	row.names(z_dat)=1:nrow(z_dat)

	x = z_dat
	y = split(x,x$V1)
    
    tmpRes = lapply(y,genotypeMNP)
    tmpMajorAllele = do.call(rbind,lapply(tmpRes,"[[",1))
    tmpMinorAllele = do.call(rbind,lapply(tmpRes,"[[",2))
   
 
    
    if(!is.null(tmpMajorAllele))
    {
      write.table(tmpMajorAllele[,c(1,3)],outF2,quote=F,row.names=F,col.names=c("Amplicon","MajorAllele"),sep="\t")
	  cat("Major alleles are listed in ", outF2,"\n")
      
      tmpNumValidAmp = nrow(tmpMajorAllele)
      if(!is.null(tmpMinorAllele))
      {
        tmpMinorAllele$numTotal = tmpAmpTotal[tmpMinorAllele$V1]
        tmpMinorAllele$fdr <- p.adjust(tmpMinorAllele$p_error)
		tmpMinorAllele = tmpMinorAllele[tmpMinorAllele$fdr<0.005,]
        write.table(tmpMinorAllele[,c("V1","DiffNT2")],outF,quote=F,row.names=F,col.names=c("Amplicon","MinorAllele"),sep="\t")
		cat("Minor alleles are listed in ",outF,"\n")
        
      }else
	  {
		  cat("No minor alleles detected!")
	  }
	}else
	{
		cat("No valid markers!")
	}   
 


	