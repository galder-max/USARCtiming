################################################
## ALL FUNCTIONS (PLOTS, TIMING, ETC.)
################################################
## 1) functions to summarise CNAs
################################################
## returns fraction LOH along the genome
getLOH.fraction <- function(cn) ## cn is a matrix from an ascat profile
{
    loh <- (cn[,"Tumour.BCN"]==0)
    allG <-     sum(cn[,"End.Position"]/1000-cn[,"Start.Position"]/1000,na.rm=T)
    LOH <-     sum(cn[loh,"End.Position"]/1000-cn[loh,"Start.Position"]/1000,na.rm=T)
    LOH/allG
}

## returns fraction homozygous deletion (HD) along the genome
getHD.fraction <- function(cn)
{
    hd=(as.character(cn$nMin1_A)=="0" & as.character(cn$nMaj1_A)=="0")
    allG <-     sum(cn[,"endpos"]/1000-cn[,"startpos"]/1000,na.rm=T)
    HD <-     sum(cn[hd,"endpos"]/1000-cn[hd,"startpos"]/1000,na.rm=T)
    HD/allG
}


## returns fraction subclonal CNAs along the genome
getSubclonal.fraction <- function(cn)
{
    cond <- p.adjust(cn[,"pval"],method="fdr")<0.05
    allG <-     sum(cn[,"endpos"]/1000-cn[,"startpos"]/1000,na.rm=T)
    SUBCLONAL <-     sum(cn[cond,"endpos"]/1000-cn[cond,"startpos"]/1000,na.rm=T)
    SUBCLONAL/allG
}

## returns total size of genome with HD
getHD.size <- function(cn)
{
    hd=(as.character(cn$nMin1_A)=="0" & as.character(cn$nMaj1_A)=="0")
    sum(cn[hd,"endpos"]-cn[hd,"startpos"])
}

## returns total size of genome with subclonal fraction around 50% (indicative of missed WGD)
getSize50CCF <- function(cn)
{
    ccf50 <- cn$frac1_A<.55 & cn$frac1_A>.45
    sum(cn[ccf50,"endpos"]/1000000-cn[ccf50,"startpos"]/1000000)
}
####################################


####################################
## 2) Filter SNVs by context
####################################
## get reference genome (to filter for (C>T)pG SNVs)
getRefGenome <- function(fasta=FASTA)
{
    dna <- readDNAStringSet(fasta,format="fasta")
    dna <- dna[1:25,]
    return(dna)
}

## filter (C>T)pG SNVs in mut (vcf-like matrix)
filterCtoTatCpG <- function(mut,colPosition=3,refCol=3,altCol=4)
{
    c2 <- colPosition
    keepCT <- mut[,refCol]=="C" & mut[,altCol]=="T"
    keepGA <- mut[,refCol]=="G" & mut[,altCol]=="A"
    keepCTGA <-  keepCT |keepGA
    contexts <- sapply(which(keepCTGA),function(x) as(dna[[as.character(mut[x,1])]][mut[x,c2]+(-1:1)],"character"))
    ##contextsrev <- sapply(contexts,function(string) as(reverseComplement(as(string,"DNAString")),"character"))
    keepCT <- contexts%in%paste0(c("C","G","T","A"),"CG")
    keepGA <- contexts%in%paste0("CG",c("C","G","T","A"))
    mut[keepCTGA,][keepCT | keepGA,]
}
####################################


####################################
## 3) Timing of WGD(s)
####################################
## dummy return of timing of WGD when no WGD is called
getT1 <- function(muts,cna)
{
    return(1)
}

## return of timing and mutation counts of WGD
getT2 <- function(muts,cna)
{
    cna <- cbind(cna,cna[,"Tumour.TCN"]-cna[,"Tumour.BCN"],cna[,"Tumour.BCN"],
                 cna[,c("Start.Position","End.Position")])
    colnames(cna)[(ncol(cna)-3):ncol(cna)] <- c("nMaj1_A","nMin1_A","startpos","endpos")
    muts <- cbind(muts,muts[,"no.chrs.bearing.mut"])
    colnames(muts)[ncol(muts)] <- "multiplicity"
    require(GenomicRanges)
    ## 2+0
    cna. <- cna[cna[,"nMaj1_A"]==2 & cna[,"nMin1_A"]==0,]
    gr1 <- GRanges(cna.[,"Chromosome"],IRanges(cna.[,"startpos"],cna.[,"endpos"]))
    gr2 <- GRanges(muts[,"chr"],IRanges(muts[,"end"],muts[,"end"]))
    ov <- findOverlaps(gr2,gr1)
    mm <- muts[queryHits(ov),]
    counts2 <- sum(mm$multiplicity==2,na.rm=T)
    counts1 <- sum(mm$multiplicity==1,na.rm=T)/2
    ## 2+2
    cna. <- cna[cna[,"nMaj1_A"]==2 & cna[,"nMin1_A"]==2,]
    gr1 <- GRanges(cna.[,"Chromosome"],IRanges(cna.[,"startpos"],cna.[,"endpos"]))
    gr2 <- GRanges(muts[,"chr"],IRanges(muts[,"end"],muts[,"end"]))
    ov <- findOverlaps(gr2,gr1)
    mm <- muts[queryHits(ov),]
    counts2.22 <- sum(mm$multiplicity>=2,na.rm=T)
    counts1.22 <- sum(mm$multiplicity==1,na.rm=T)/2
    list(counts2=counts2,
         halfcounts1=counts1,
         timingWGD1=counts2/(counts1+counts2),
         counts2.22=counts2.22,
         halfcounts1.22=counts1.22,
         timingWGD1.22=counts2.22/(counts1.22+counts2.22),
         timingWGD1.tot=(counts2.22+counts2)/(counts1.22+counts2.22+counts1+counts2))
}

## return of timing and mutation counts of second WGD
getT4 <- function(muts,cna)
{
    cna <- cbind(cna,cna[,"Tumour.TCN"]-cna[,"Tumour.BCN"],cna[,"Tumour.BCN"],
                 cna[,c("Start.Position","End.Position")])
    colnames(cna)[(ncol(cna)-3):ncol(cna)] <- c("nMaj1_A","nMin1_A","startpos","endpos")
    muts <- cbind(muts,muts[,"no.chrs.bearing.mut"])
    colnames(muts)[ncol(muts)] <- "multiplicity"
    require(GenomicRanges)
    ## 4+0
    cna. <- cna[cna[,"nMaj1_A"]==4 & cna[,"nMin1_A"]==0,]
    gr1 <- GRanges(cna.[,"Chromosome"],IRanges(cna.[,"startpos"],cna.[,"endpos"]))
    gr2 <- GRanges(muts[,"chr"],IRanges(muts[,"end"],muts[,"end"]))
    ov <- findOverlaps(gr2,gr1)
    mm <- muts[queryHits(ov),]
    counts3 <- sum(mm$multiplicity==2,na.rm=T)/2
    counts4 <- sum(mm$multiplicity>=4,na.rm=T)+counts3
    counts2 <- (sum(mm$multiplicity==2,na.rm=T)+counts3)/2
    counts1 <- sum(mm$multiplicity==1,na.rm=T)/4
    ## 4+4
    cna. <- cna[cna[,"nMaj1_A"]==4 & cna[,"nMin1_A"]==4,]
    gr1 <- GRanges(cna.[,"Chromosome"],IRanges(cna.[,"startpos"],cna.[,"endpos"]))
    gr2 <- GRanges(muts[,"chr"],IRanges(muts[,"end"],muts[,"end"]))
    ov <- findOverlaps(gr2,gr1)
    mm <- muts[queryHits(ov),]
    counts3.44 <- sum(mm$multiplicity==3,na.rm=T)/2
    counts4.44 <- sum(mm$multiplicity>=4,na.rm=T)+counts3.44
    counts2.44 <- (sum(mm$multiplicity==2,na.rm=T)+counts3.44)/2
    counts1.44 <- sum(mm$multiplicity==1,na.rm=T)/4
    list(counts4=counts4,
         halfcounts2=counts2,
         quartercounts1=counts1,
         timingWGD1=counts4/(counts1+counts2+counts4),
         timingWGD2=(counts2+counts4)/(counts1+counts2+counts4),
         counts4.44=counts4.44,
         halfcounts2.44=counts2.44,
         quartercounts1.44=counts1.44,
         timingWGD1.44=counts4.44/(counts1.44+counts2.44+counts4.44),
         timingWGD2.44=(counts2.44+counts4.44)/(counts1.44+counts2.44+counts4.44),
         timingWGD1.tot=(counts4.44+counts4)/(counts1.44+counts2.44+counts4.44+counts1+counts2+counts4),
         timingWGD2.tot=(counts2.44+counts4.44+counts2+counts4)/(counts1.44+counts2.44+counts4.44+counts1+counts2+counts4))
}

## get timing depending on mode of the major allele
getTiming <- function(muts,cna,mode)
{
    switch(mode,
           "1"=getT1(muts,cna),
           "2"=getT2(muts,cna),
           "3"=getT4(muts,cna),
           "4"=getT4(muts,cna),
           "5"=getT4(muts,cna))
}

## customised try catch
mytry <- function(x,ret=list(timingWGD1=NA),silent=T,...)
{
    res <- try(x,silent=silent,...)
    if(inherits(res,"try-error")) return(ret)
    res
}
###############################################################################



###############################################################################
## 4) Timing of drivers SNVs and indels
###############################################################################
## Drivers
## subset SNV mutations falling in driver genes
getDriversSNV <- function(x)
{
    require(GenomicRanges)
    gr1 <- GRanges(as.character(x[,1]),IRanges(x[,2],x[,2]))
    gr2 <- GRanges(gsub("chr","",drivers[,2]),IRanges(as.numeric(drivers[,3]),as.numeric(drivers[,4])))
    ov <- findOverlaps(gr1,gr2)
    ret <- x[queryHits(ov),]
    ret <- cbind(ret,as.character(drivers[subjectHits(ov),1]))
}

## subset indel mutations falling in driver genes
getDriversINDEL <- function(x)
{
    require(GenomicRanges)
    gr1 <- GRanges(as.character(x[,1]),IRanges(x[,2],x[,2]))
    gr2 <- GRanges(gsub("chr","",drivers[,2]),IRanges(as.numeric(drivers[,3]),as.numeric(drivers[,4])))
    ov <- findOverlaps(gr1,gr2)
    ret <- x[queryHits(ov),]
    ret <- cbind(ret,as.character(drivers[subjectHits(ov),1]))
}

## (not used) subset CNA segments involving driver genes with LOH, HD or Amplification
getDriversCNA <- function(x)
{
    require(GenomicRanges)
    gr1 <- GRanges(as.character(x[,2]),IRanges(x[,3],x[,4]))
    gr2 <- GRanges(gsub("chr","",drivers[,2]),IRanges(as.numeric(drivers[,3]),as.numeric(drivers[,4])))
    ov <- findOverlaps(gr1,gr2)
    ret <- x[queryHits(ov),]
    ret <- cbind(ret,drivers[subjectHits(ov),1])
    keepAmp <- (ret[,"Tumour.TCN"]-ret[,"Tumour.BCN"])>=5
    keepHD <- (ret[,"Tumour.TCN"]-ret[,"Tumour.BCN"])==0
    keepLOH <- ret[,"Tumour.BCN"]==0
    ret <- cbind(ret,rep(NA,nrow(ret)))
    ret[keepLOH,ncol(ret)] <- "LOH"
    ret[keepHD,ncol(ret)] <- "HD"
    ret[keepAmp,ncol(ret)] <- "Amplified"
    ret[,]
}

## (not used) annotate mutations as clonal or subclonal and for which cluster they are assigned to (from DPCLust output)
getClusters <- function(muts)
{
    clusts <- sort(tapply(1:nrow(muts),muts[,"cluster"],function(x) mean(muts[x,"ccf"],na.rm=T)),decreasing=T)
    clonal <- which(clusts>=.95)
    if(length(clonal)==0) clonal <- 1
    subclonal <- (1:length(clusts))[-c(clonal)]
    list(clonal=names(clusts)[clonal],
         subclonal=names(clusts)[subclonal])
}

## time driver mutations using inferred multiplicities and cna states
## no WGD
timeDrive.1 <- function(ind)
{
    snv <- mutsDrive[[ind]][[1]]
    if(nrow(snv)==0) return(NA)
    ##clonality <- sapply(snv[,ncol(snv)-1],function(x) if(x%in%clusters[[ind]][[1]]) "clonal" else "subclonal")
    cna <- mutsDrive[[ind]][[2]]
    cna <- cbind(cna,cna[,"Tumour.TCN"]-cna[,"Tumour.BCN"],cna[,"Tumour.BCN"],
                 cna[,c("Start.Position","End.Position")])
    colnames(cna)[(ncol(cna)-3):ncol(cna)] <- c("nMaj1_A","nMin1_A","startpos","endpos")
    timing <- sapply(1:nrow(snv),function(x)
    {
        gene <- as.character(snv[x,ncol(snv)])
        ww <- which(cna[,ncol(cna)-5]==gene)[1]
        states <- cna[ww,c("nMaj1_A","nMin1_A")]
        if(any(is.na(states))) return(NA)
        if(states[1]>1)
        {
            if(snv[x,"no.chrs.bearing.mut"]>1) return("early")
            if(states[2]==0) return("late")
        }
        return(NA)
    })
    list(SNV=cbind(snv,timing),
         CNA=cna[cna[,ncol(cna)]=="HD" & !is.na(cna[,ncol(cna)]),])
}

## timing logic for multiplicity relative to WGD
timeWGD <- function(states, ## copy number states
                    mult) ## mutation multiplicity
{
    if(states[1]>1)
    {
        if(mult>1) return("beforeWGD1")
        if(states[2]==0 | states[2]>1) return("afterWGD1")
    }
    if(mult>1) return("beforeWGD1")
    return(NA)
}

## timing logic for multiplicity relative to second WGD
timeWGD2 <- function(states,
                     mult)
{
    if(states[1]>1 & states[1]<4)
    {
        if(mult>1) return("beforeWGD2")
        if(states[2]==0 | states[2]>1) return("afterWGD2")
    }
    if(states[1]>=4)
    {
        if(states[2]==0)
        {
            if(mult==1) return("afterWGD2")
            if(mult<4) return("afterWGD1beforeWGD2")
            if(mult>=4) return("beforeWGD1")
        }
        if(states[2]==4)
        {
            if(mult==1) return("afterWGD2")
            if(mult<4) return("afterWGD1beforeWGD2")
            if(mult>=4) return("beforeWGD1")
        }
    }
    if(mult>1) return("beforeWGD2")
    return(NA)
}

## time driver mutations using inferred multiplicities and cna states
## 1xWGD
timeDrive.2 <- function(ind)
{
    snv <- mutsDrive[[ind]][[1]]
    if(nrow(snv)==0) return(NA)
    ##clonality <- sapply(snv[,ncol(snv)-1],function(x) if(x%in%clusters[[ind]][[1]]) "clonal" else "subclonal")
    cna <- mutsDrive[[ind]][[2]]
    cna <- cbind(cna,cna[,"Tumour.TCN"]-cna[,"Tumour.BCN"],cna[,"Tumour.BCN"],
                 cna[,c("Start.Position","End.Position")])
    colnames(cna)[(ncol(cna)-3):ncol(cna)] <- c("nMaj1_A","nMin1_A","startpos","endpos")
    timing <- sapply(1:nrow(snv),function(x)
    {
        gene <- as.character(snv[x,ncol(snv)])
        ww <- which(cna[,ncol(cna)-5]==gene)[1]
        states <- cna[ww,c("nMaj1_A","nMin1_A")]
        if(all(!is.na(states))) return(timeWGD(states,snv[x,"no.chrs.bearing.mut"]))
        NA
    })
    list(SNV=cbind(snv,timing),
         CNA=cna[cna[,ncol(cna)]=="HD" & !is.na(cna[,ncol(cna)]),])
}

## time driver mutations using inferred multiplicities and cna states
## >1xWGD
timeDrive.4 <- function(ind)
{
    snv <- mutsDrive[[ind]][[1]]
    if(nrow(snv)==0) return(NA)
    ##clonality <- sapply(snv[,ncol(snv)-1],function(x) if(x%in%clusters[[ind]][[1]]) "clonal" else "subclonal")
    cna <- mutsDrive[[ind]][[2]]
    cna <- cbind(cna,cna[,"Tumour.TCN"]-cna[,"Tumour.BCN"],cna[,"Tumour.BCN"],
                 cna[,c("Start.Position","End.Position")])
    colnames(cna)[(ncol(cna)-3):ncol(cna)] <- c("nMaj1_A","nMin1_A","startpos","endpos")
    timing <- sapply(1:nrow(snv),function(x)
    {
        gene <- as.character(snv[x,ncol(snv)])
        ww <- which(cna[,ncol(cna)-5]==gene)[1]
        states <- cna[ww,c("nMaj1_A","nMin1_A")]
        if(all(!is.na(states))) return(timeWGD2(states,snv[x,"no.chrs.bearing.mut"]))
        NA
    })
    list(SNV=cbind(snv,timing),
         CNA=cna[cna[,ncol(cna)]=="HD" & !is.na(cna[,ncol(cna)]),])
}

## time driver mutation depending on mode
timeDrive <- function(x)
{
    mode <- modeCN[x]
    return(if(mode%in%c("0","1"))
               timeDrive.1(x)
           else if(mode=="2")
               timeDrive.2(x)
           else
               timeDrive.4(x))
}
###############################################################################


###############################################################################
## 5) plotting functions for timing
###############################################################################
## (not used) get confidence intervals for timing of WGDx2 using binomial distrib
getCIs.WGD2 <- function(count1,count2,count4,ci=c(0.025,.975))
{
    tot <- round(count2+count4+count1)
    qbinom(ci,size=tot,prob=(count4+count2)/tot)/tot
}

## (not used) get confidence intervals for timing of WGDx1 using binomial distrib
getCIs.WGD1 <- function(count1,count2,ci=c(0.025,.975),count4=0)
{
    tot <- round(count1+count2+count4)
    qbinom(ci,size=tot,prob=count2/tot)/tot
}


## (old version/over-written) plotting function for the timing
plotTiming <- function(timings)
{
    tWGD1s <- sapply(timings,function(x)
    {
        if(length(x)==1) return(NA)
        if("timingWGD1"%in%names(x))
        {
            return(x$timingWGD1.tot)
        }
        return(NA)
    })
    tWGD2s <- sapply(timings,function(x)
    {
        if(length(x)==1) return(NA)
        if("timingWGD2"%in%names(x))
        {
            return(x$timingWGD2.tot)
        }
        return(NA)
    })
    CIWGD1 <- lapply(timings,function(x)
    {
        if(length(x)==1) return(NA)
        if("timingWGD1"%in%names(x))
        {
            if(!"timingWGD2"%in%names(x))
                return(getCIs.WGD1(x$halfcounts1+x$halfcounts1.22,x$counts2+x$counts2.22))
            else
                return(getCIs.WGD1(x$halfcounts2+x$halfcounts2.44,x$counts4+x$counts4.44,
                                   count4=x$quartercounts1+x$quartercounts1.44))
        }
        return(NA)
    })
    CIWGD2 <- lapply(timings,function(x)
    {
        if(length(x)==1) return(NA)
        if("timingWGD2"%in%names(x))
        {
            return(getCIs.WGD2(x$quartercounts1+x$quartercounts1.44,
                               x$halfcounts2+x$halfcounts2.44,
                               x$counts4+x$counts4.44))
        }
        return(NA)
    })
    ord <- order(tWGD1s,tWGD2s,decreasing=F)
    ciwgd1.low <- sapply(CIWGD1,function(x) if(length(x)==2) x[1] else NA)
    ciwgd1.high <- sapply(CIWGD1,function(x) if(length(x)==2) x[2] else NA)
    ciwgd2.low <- sapply(CIWGD2,function(x) if(length(x)==2) x[1] else NA)
    ciwgd2.high <- sapply(CIWGD2,function(x) if(length(x)==2) x[2] else NA)
    plot(1:length(tWGD1s),tWGD1s[ord],
         frame=F,
         xlab="Samples",
         ylab="Relative Mutational Timing",
         ##ylim=c(min(ciwgd1.low,na.rm=T),1),
         ylim=c(-.5,1.2),
         xlim=c(0,length(tWGD1s)+10),
         xaxt="n",
         yaxt="n",
         pch=21,
         col=rgb(.4,.4,.4,.4))
    ##axis(side=1,at=seq(0,70,10),seq(0,70,10))
    text(1:length(tWGD1s)+1,rep(-.1,length(tWGD1s)),labels=samps[ord],cex=.55,srt=45,pos=2)
    axis(side=2,at=seq(0,1,.2),seq(0,1,.2))
    points(1:length(tWGD1s),tWGD2s[ord],
           pch=19,col=rgb(.2,.2,.2,.4))
    abline(h=c(0,.25,0.5,.75,1),lty=2,lwd=.5,col=rgb(.5,.5,.5,.5))
    segments(1:length(tWGD1s),ciwgd1.low[ord],1:length(tWGD1s),
             ciwgd1.high[ord],
             lwd=.7,col=rgb(.5,.5,.5,.5))
    segments(1:length(tWGD2s),ciwgd2.low[ord],1:length(tWGD2s),
             ciwgd2.high[ord],lwd=.7,
             col=rgb(.5,.5,.5,.5))
    NN <- -5
    text(length(timings)-NN,.97,"MRCA")
    arrows(length(timings)-NN,.45,length(timings)-NN,.25,length=.07)
    arrows(length(timings)-NN,.55,length(timings)-NN,.75,length=.07)
    text(length(timings)-NN,.5,"Clonal")
    text(length(timings)-NN,.03,"Early life")
    ##text(length(timings)-NN,1.1,"Subclonal")
    list(ord=ord,
         tWGD1s=tWGD1s[ord],
         tWGD2s=tWGD2s[ord],
         ciwgd1.low=ciwgd1.low[ord],
         ciwgd1.high=ciwgd1.high[ord],
         ciwgd2.low=ciwgd2.low[ord],
         ciwgd2.high=ciwgd2.high[ord])
}

## adding drivers to the plot
plotPointDrive <- function(x,mD,retPlot,cols,PCH)
{
    mD <- cbind(mD,rep("clonal",nrow(mD)))
    colnames(mD)[ncol(mD)] <- "clonality"
    COLGENE <- 17
    MINIOFF <- .03
    for(i in 1:nrow(mD))
    {
        if(mD[i,"clonality"]=="clonal")
        {
            if(mD[i,"timing"]=="beforeWGD1" & !is.na(mD[i,"timing"]))
            {
                points(x,
                       retPlot$ciwgd1.low[x]-MINIOFF-OFFSETS[as.character(mD[i,COLGENE])]*MINIOFF,
                       col=COLS[as.character(mD[i,COLGENE])],
                       pch=PCH)
            }
            if(grepl("beforeWGD2",mD[i,"timing"])& !is.na(mD[i,"timing"]))
            {
                points(x,
                       retPlot$ciwgd2.low[x]-MINIOFF-OFFSETS[as.character(mD[i,COLGENE])]*MINIOFF,
                       col=COLS[as.character(mD[i,COLGENE])],
                       pch=PCH)
            }
            if(mD[i,"timing"]=="afterWGD1"& !is.na(mD[i,"timing"]))
            {
                points(x,
                       retPlot$ciwgd1.high[x]+MINIOFF+OFFSETS[as.character(mD[i,16])]*MINIOFF,
                       col=COLS[as.character(mD[i,16])],
                       pch=PCH)
            }
            if(mD[i,"timing"]=="afterWGD2"& !is.na(mD[i,"timing"]))
            {
                points(x,
                       retPlot$ciwgd2.high[x]+MINIOFF+OFFSETS[as.character(mD[i,COLGENE])]*MINIOFF,
                       col=COLS[as.character(mD[i,COLGENE])],
                       pch=PCH)
            }
        }
        if(mD[i,"clonality"]=="subclonal" & !grepl("before",mD[i,"timing"])& !is.na(mD[i,"timing"]))
        {
            points(x,
                   1.1+OFFSETS[as.character(mD[i,COLGENE])]*MINIOFF,
                   col=COLS[as.character(mD[i,COLGENE])],
                   pch=PCH)
        }
    }
}

## (not used) adding subclonal mutations
plotSubclonalMuts <- function(mutsDriveTimed,
                              retPlot,
                              cols=COLS,PCH=19)
{
    mD <- lapply(retPlot$ord,function(x) mutsDriveTimed[[x]])
    ##  plot(1:length(mD),col=rgb(0,0,0,0),xlab="",ylab="Subclonal",xaxt="n",yaxt="n",frame=F)
    sapply(1:length(mD),function(x)
    {
        cat(".")
        if(!is.null(nrow(mD[[x]][[1]])))
            if(nrow(mD[[x]][[1]])>0)
            {
                plotPointDrive(x,mD[[x]][[1]],retPlot,cols,PCH=PCH)
            }
    })
}

## (not used) heatmap of presence/absence of driver mutations per sample
heatmapDrive <- function(mutsDriveTimed,retPlot)
{
    mD <- lapply(retPlot$ord,function(x) mutsDriveTimed[[x]])
    mat <- matrix(F,nrow(drivers),length(mD))
    rownames(mat) <- drivers[,1]
    for(i in 1:length(mD))
    {
        if(!is.null(nrow(mD[[i]][[1]])))
            if(nrow(mD[[i]][[1]])>0)
            {
                mat[,i] <- drivers[,1]%in%as.character(mD[[i]][[1]][,17])
            }
    }
    LIMX <- length(retPlot$tWGD1s)+10
    plot(0,0,xlab="",ylab="",
         xaxt="n",yaxt="n",
         frame=F,
         col=rgb(0,0,0,0),
         xlim=c(0,LIMX),
         ylim=c(0,1))
    nc <- ncol(mat)
    pasc <- length(retPlot$tWGD1s)/nc
    nr <- nrow(mat)
    pasr <- 1/nr
    for(i in 1:nrow(mat))
    {
        for(j in 1:ncol(mat))
        {
            polygon(c((j-1)*pasc,(j-1)*pasc,j*pasc,j*pasc)+.5,
                    c((i-1)*pasr,i*pasr,i*pasr,(i-1)*pasr),
                    border=NA,
                    col=if(mat[i,j]) COLS[rownames(mat)[i]] else rgb(0,0,0,0))
        }
    }
}

## adding dummy plot+legend
mylegend <- function(PCH=19)
{
    plot(0,0,xlab="",ylab="",
         xaxt="n",yaxt="n",
         frame=F,
         col=rgb(0,0,0,0),
         xlim=c(0,1),
         ylim=c(0,1))
    graphics:::legend("left",pch=PCH,col=COLS,legend=names(COLS),box.col=rgb(0,0,0,0))
}

## (not used) boxplot+points for subclonal fractions
boxpoints <- function(l,col1=rgb(0,0,0,.5),pch1=19,cex1=1,...)
{
    boxplot(l,frame=F,col=rgb(0,0,0,0),...)
    jitter <- list()
    for(i in 1:length(l))
    {
        jitter[[i]] <- rnorm(length(l[[i]]),sd=0.1)
        points(i+jitter[[i]],l[[i]],col=col1,pch=pch1,cex=cex1)
    }
    ##segments(1+jitter[[1]],l[[1]],2+jitter[[2]],l[[2]],lwd=1,col=rgb(.3,.3,.3,.2))
}



## (not used) add a wrapping background for real time projection in molecular-time-based plots
addAge <- function(ages,ord,colL)
  {
    ages <- ages[ord]
    lines <- lapply(c(20,40,60,80),function(x)
                    x/ages)
    ##lapply(lines,function(x) points(1:length(x),x,type="l",col=colL))
    polygon(x=c(1,1:length(lines[[1]]),length(lines[[1]])),
            y=c(0,lines[[1]],0),col=rgb(0,0,0,.02),border=NA)
    polygon(x=c(1:length(lines[[1]]),length(lines[[1]]):1),
            y=c(lines[[1]],lines[[2]][length(lines[[2]]):1]),
            col=rgb(0,0,0,.1),border=NA)
    polygon(x=c(1:length(lines[[1]]),length(lines[[1]]):1),
            y=c(lines[[2]],lines[[3]][length(lines[[2]]):1]),
            col=rgb(0,0,0,.15),border=NA)
    polygon(x=c(1:length(lines[[1]]),length(lines[[1]]):1),
            y=c(lines[[3]],lines[[4]][length(lines[[2]]):1]),
            col=rgb(0,0,0,.2),border=NA)
    polygon(x=c(1,1:length(lines[[1]]),length(lines[[1]])),
            y=c(2,lines[[4]],2),
            col=rgb(0,0,0,.3),border=NA)
  }

## plot summary of molecular timing of WGD and drivers
## returns confidence intervals
plotTiming <- function(timings,groups,samps..,colL=rgb(.4,.4,.4,.1))
{
    tWGD1s <- sapply(timings,function(x)
    {
        if(length(x)==1) return(NA)
        if("timingWGD1"%in%names(x))
        {
            return(x$timingWGD1.tot)
        }
        return(NA)
    })
    tWGD2s <- sapply(timings,function(x)
    {
        if(length(x)==1) return(NA)
        if("timingWGD2"%in%names(x))
        {
            return(x$timingWGD2.tot)
        }
        return(NA)
    })
    CIWGD1 <- lapply(timings,function(x)
    {
        if(length(x)==1) return(NA)
        if("timingWGD1"%in%names(x))
        {
            if(!"timingWGD2"%in%names(x))
                return(getCIs.WGD1(x$halfcounts1+x$halfcounts1.22,x$counts2+x$counts2.22))
            else
                return(getCIs.WGD1(x$halfcounts2+x$halfcounts2.44,x$counts4+x$counts4.44,
                                   count4=x$quartercounts1+x$quartercounts1.44))
        }
        return(NA)
    })
    CIWGD2 <- lapply(timings,function(x)
    {
        if(length(x)==1) return(NA)
        if("timingWGD2"%in%names(x))
        {
            return(getCIs.WGD2(x$quartercounts1+x$quartercounts1.44,
                               x$halfcounts2+x$halfcounts2.44,
                               x$counts4+x$counts4.44))
        }
        return(NA)
    })
    theages <- clin$Age_at_diagnosis
    ord <- order(groups,tWGD1s,decreasing=F)
    ciwgd1.low <- sapply(CIWGD1,function(x) if(length(x)==2) x[1] else NA)
    ciwgd1.high <- sapply(CIWGD1,function(x) if(length(x)==2) x[2] else NA)
    ciwgd2.low <- sapply(CIWGD2,function(x) if(length(x)==2) x[1] else NA)
    ciwgd2.high <- sapply(CIWGD2,function(x) if(length(x)==2) x[2] else NA)
    COOL <- rgb(t(col2rgb(RColorBrewer:::brewer.pal(12,"Paired"))/255))
    ##COOL <- (rgb(t(col2rgb(sample(topo.colors(10))))/255))
    COLSCN <- COOL[groups[ord]]
    plot(1:length(tWGD1s),tWGD1s[ord],
         frame=F,
         xlab="Samples",
         ylab="Relative Mutational Timing",
         ##ylim=c(min(ciwgd1.low,na.rm=T),1),
         ylim=c(-.5,1.2),
         xlim=c(0,length(tWGD1s)+10),
         xaxt="n",
         yaxt="n",
         pch=21,
         col=rgb(.4,.4,.4,.4))
    ##axis(side=1,at=seq(0,70,10),seq(0,70,10))
    text(1:length(tWGD1s)+1,rep(-.1,length(tWGD1s)),
         labels=samps..[ord],
         cex=.55,
         srt=45,
         col=COLSCN,
         pos=2)
    axis(side=2,at=seq(0,1,.2),seq(0,1,.2))
    points(1:length(tWGD1s),tWGD2s[ord],
           pch=19,col=rgb(.2,.2,.2,.4))
    ##abline(h=c(0,.25,0.5,.75,1),lty=2,lwd=.5,col=rgb(.5,.5,.5,.5))
    addAge(ages=theages,
           ord=ord,
           colL=colL)
    segments(1:length(tWGD1s),ciwgd1.low[ord],1:length(tWGD1s),
             ciwgd1.high[ord],
             lwd=.7,col=rgb(.5,.5,.5,.5))
    segments(1:length(tWGD2s),ciwgd2.low[ord],1:length(tWGD2s),
             ciwgd2.high[ord],lwd=.7,
             col=rgb(.5,.5,.5,.5))
    legend("topright",paste0("CNsig",1:7),col=COOL[1:7],pch=15)
    NN <- -5
    text(length(timings)-NN,.97,"MRCA")
    arrows(length(timings)-NN,.45,length(timings)-NN,.25,length=.07)
    arrows(length(timings)-NN,.55,length(timings)-NN,.75,length=.07)
    text(length(timings)-NN,.5,"Clonal")
    text(length(timings)-NN,.03,"Early life")
    ##text(length(timings)-NN,1.1,"Subclonal")
    list(ord=ord,
         tWGD1s=tWGD1s[ord],
         tWGD2s=tWGD2s[ord],
         ciwgd1.low=ciwgd1.low[ord],
         ciwgd1.high=ciwgd1.high[ord],
         ciwgd2.low=ciwgd2.low[ord],
         ciwgd2.high=ciwgd2.high[ord])
}


## get confidence interval of timing using bootstrapping for second WGD
getCIs.WGD2 <- function(count1,count2,count4,ci=c(0.025,.975),Nrep=1000)
{
  c2 <- round(count2+count4+.51)
  c1 <- round(count1+0.51)
  tot <- c2+c1
  vv <- c(rep("c1",c1),rep("c2",c2))
  quantile(sapply(1:Nrep,function(x)
                  {
                    vv. <- sample(vv,rep=T)
                    sum(vv.=="c2")/tot
                  }),probs=ci)
}

## get confidence interval of timing using bootstrapping for first WGD
getCIs.WGD1 <- function(count1,count2,ci=c(0.025,.975),count4=0,Nrep=1000)
{
  c2 <- round(count2+.51)
  c1 <- round(count1+0.51)
  c4 <- round(count4+.51)
  tot <- c2+c1+c4
  vv <- c(rep("c1",c1+c4),rep("c2",c2))
  quantile(sapply(1:Nrep,function(x)
                  {
                    vv. <- sample(vv,rep=T)
                    sum(vv.=="c2")/tot
                  }),probs=ci)
}

## boxplot+points
boxpoints <- function(l,col1=rgb(0,0,0,.5),pch1=19,cex1=1,...)
{
    boxplot(l,frame=F,col=rgb(0,0,0,0),...)
    for(i in 1:length(l))
    {
        points(i+rnorm(length(l[[i]]),sd=0.1),l[[i]],col=col1,pch=pch1,cex=cex1)
    }
}

## plotting function for real-time timing of WGDs per signature group
boxpoints.segs <- function(l,
                           lSegslow,
                           lSegshigh,
                           lGroups,
                           lPch,
                           col1=rgb(0,0,0,.5),pch1=19,cex1=1,...)
  {
    ll <- sapply(l,length)
    nms <- names(l)
    l <- lapply(which(ll>1),function(x) l[[x]])
    lSegslow <- lapply(which(ll>1),function(x) lSegslow[[x]])
    lSegshigh <- lapply(which(ll>1),function(x) lSegshigh[[x]])
    lGroups <- lapply(which(ll>1),function(x) lGroups[[x]])
    lPch <- lapply(which(ll>1),function(x) lPch[[x]])
    names(l) <- nms[ll>3]
    meds <- sapply(l,median)
    ord <- order(meds,decreasing=T)
    nms <- names(l)
    l <- lapply(ord,function(x) l[[x]])
    lSegslow <- lapply(ord,function(x) lSegslow[[x]])
    lSegshigh <- lapply(ord,function(x) lSegshigh[[x]])
    lGroups <- lapply(ord,function(x) lGroups[[x]])
    lPch <- lapply(ord,function(x) lPch[[x]])
    names(l) <- nms[ord]
    boxplot(l,frame=F,col=rgb(0,0,0,0),...)
    jitter <- list()
    lSegslow <- lapply(1:length(lSegslow),function(x) lSegslow[[x]]*as.numeric(!is.na(l[[x]])))
    lSegshigh <- lapply(1:length(lSegshigh),function(x) lSegshigh[[x]]*as.numeric(!is.na(l[[x]])))
    COOL <- rgb(t(col2rgb(RColorBrewer:::brewer.pal(12,"Paired"))/255))
    for(i in 1:length(l))
      {
        jitter[[i]] <- sort(rnorm(length(l[[i]]),sd=0.1))
        ord <- order(l[[i]],decreasing=F)
        points(i+jitter[[i]],l[[i]][ord],pch=lPch[[i]][ord],cex=cex1,col=paste0(COOL[lGroups[[i]]][ord],"99"))
        segments(i+jitter[[i]],lSegslow[[i]][ord],i+jitter[[i]],lSegshigh[[i]][ord],col=paste0(COOL[lGroups[[i]]][ord],"21"))
      }
    ##segments(1+jitter[[1]],l[[1]],2+jitter[[2]],l[[2]],lwd=1,col=rgb(.3,.3,.3,.2))
  }
########################################################



########################################################
## 6) Inferring multiplicities of SNVs and indels
########################################################
## functions for pre-formatting and inference of multiplicities using dpclust preprocessing functions
## writes a file with positions and base change for SNV
createLociFile <- function(vcfdat, outfile, chrom_col, pos_col, ref_col, alt_col)
{
    loci = vcfdat[, c(chrom_col, pos_col, ref_col, alt_col)]
    loci_file = "loci.txt"
    write.table(loci, file=outfile, sep="\t", quote=F, row.names=F, col.names=F)
}

## writes a file with positions and base change for indels
createLociFile.indel <- function(vcfdat, outfile, chrom_col, pos_col, ref_col, alt_col)
{
    loci = cbind(vcfdat[, c(chrom_col, pos_col)],rep("A",nrow(vcfdat)),rep("T",nrow(vcfdat)))
    loci_file = "loci.txt"
    write.table(loci, file=outfile, sep="\t", quote=F, row.names=F, col.names=F)
}

## extract allele counts from VCF (Caveman)
getCounts <- function(vec)
{
    A1 <- as.numeric(gsub("(.*):(.*):(.*):(.*):(.*):(.*):(.*):(.*):(.*):(.*)",
                          "\\2",vec))
    A2 <- as.numeric(gsub("(.*):(.*):(.*):(.*):(.*):(.*):(.*):(.*):(.*):(.*)",
                          "\\6",vec))
    As <- A1+A2
    C1 <- as.numeric(gsub("(.*):(.*):(.*):(.*):(.*):(.*):(.*):(.*):(.*):(.*)",
                          "\\3",vec))
    C2 <- as.numeric(gsub("(.*):(.*):(.*):(.*):(.*):(.*):(.*):(.*):(.*):(.*)",
                          "\\7",vec))
    Cs <- C1+C2
    G1 <- as.numeric(gsub("(.*):(.*):(.*):(.*):(.*):(.*):(.*):(.*):(.*):(.*)",
                          "\\4",vec))
    G2 <- as.numeric(gsub("(.*):(.*):(.*):(.*):(.*):(.*):(.*):(.*):(.*):(.*)",
                          "\\8",vec))
    Gs <- G1+G2
    T1 <- as.numeric(gsub("(.*):(.*):(.*):(.*):(.*):(.*):(.*):(.*):(.*):(.*)",
                          "\\5",vec))
    T2 <- as.numeric(gsub("(.*):(.*):(.*):(.*):(.*):(.*):(.*):(.*):(.*):(.*)",
                          "\\9",vec))
    Ts <- T1+T2
    allcounts <- cbind(As,Cs,Gs,Ts)
    colnames(allcounts) <- c("A","C","G","T")
    allcounts
}

## write alle counts file for inferring multiplicities for SNV
createAlleleCountsFile <- function(vcfdat, outfile)
{
    CT <- getCounts(vcfdat[,11])
    mutCount <- CT[,"A"]
    mutCount[vcfdat[,5]=="G"] <- CT[vcfdat[,5]=="G","G"]
    mutCount[vcfdat[,5]=="C"] <- CT[vcfdat[,5]=="C","C"]
    mutCount[vcfdat[,5]=="T"] <- CT[vcfdat[,5]=="T","T"]
    refCount <- CT[,"A"]
    refCount[vcfdat[,4]=="G"] <- CT[vcfdat[,4]=="G","G"]
    refCount[vcfdat[,4]=="C"] <- CT[vcfdat[,4]=="C","C"]
    refCount[vcfdat[,4]=="T"] <- CT[vcfdat[,4]=="T","T"]
    counts_table <- mutwt2allelecounts(counts.alt=mutCount,
                                       counts.ref=refCount,
                                       allele.alt=as.character(vcfdat[,5]),
                                       allele.ref=as.character(vcfdat[,4]))
    output <- data.frame(as.character(vcfdat[,1]), vcfdat[,2],
                         counts_table, rowSums(counts_table))
    colnames(output) <- c("#CHR","POS","Count_A",
                          "Count_C",
                          "Count_G","Count_T","Good_depth")
    write.table(output, file=outfile, sep="\t", quote=F, row.names=F)
}

## write alle counts file for inferring multiplicities for indels
createAlleleCountsFile.indel <- function(vcfdat, outfile)
{
    refCount <- mutCount <- rep(NA,nrow(vcfdat))
    for(i in 1:nrow(vcfdat))
    {
        flags <- strsplit(as.character(vcfdat[i,9]),split=":")[[1]]
        refC <- which(flags=="WTR")
        altC <- which(flags=="MTR")
        sre <- c(".*",rep(":.*",length(flags)-1))
        sre[refC] <- ":(.*)"
        sre <- paste0(sre,collapse="")
        refCount[i] <- as.numeric(gsub(sre,"\\1",vcfdat[i,11]))
        sre <- c(".*",rep(":.*",length(flags)-1))
        sre[altC] <- ":(.*)"
        sre <- paste0(sre,collapse="")
        mutCount[i] <- as.numeric(gsub(sre,"\\1",vcfdat[i,11]))
    }
    counts_table <- mutwt2allelecounts(counts.alt=mutCount,
                                       counts.ref=refCount,
                                       allele.alt=rep("A",length(mutCount)),
                                       allele.ref=rep("T",length(mutCount)))
    output <- data.frame(as.character(vcfdat[,1]), vcfdat[,2],
                         counts_table, rowSums(counts_table))
    colnames(output) <- c("#CHR","POS","Count_A",
                          "Count_C",
                          "Count_G","Count_T","Good_depth")
    write.table(output, file=outfile, sep="\t", quote=F, row.names=F)
}

## formatting of allele counts for ref and alt
mutwt2allelecounts <- function(counts.alt, counts.ref, allele.alt, allele.ref)
{
    output = array(0, c(length(allele.ref), 4))
    nucleotides = c("A", "C", "G", "T")
    nucleo.index = match(allele.alt, nucleotides)
    for (i in 1:nrow(output)) {
        output[i,nucleo.index[i]] = counts.alt[i]
    }
    nucleo.index = match(allele.ref, nucleotides)
    for (i in 1:nrow(output)) {
        output[i,nucleo.index[i]] = counts.ref[i]
    }
    return(output)
}

## uses dpclust3p pre-processing to derive multiplicities
## https://github.com/Wedge-Oxford/dpclust_smchet_docker
writeDPfile <- function(dpFile="dpInput.txt",
                        sampleName,
                        battenberg_rho_psi_file,
                        VCFFILE,
                        BBFILE,
                        GENDER)
{
  vcfdat <- read.table(VCFFILE)
  loci_file = "loci.txt"
  createLociFile(vcfdat, loci_file, 1,2,4,5)
  ## Create allelecounts file
  allelecounts_file = "alleleCounts.txt"
  createAlleleCountsFile(vcfdat, allelecounts_file)
  suppressWarnings(runGetDirichletProcessInfo(loci_file=loci_file,
                                              allele_frequencies_file=allelecounts_file,
                                              cellularity_file=battenberg_rho_psi_file,
                                              subclone_file=BBFILE,
                                              gender=GENDER,
                                              SNP.phase.file="NA",
                                              mut.phase.file="NA",
                                              output_file=dpFile))
}

## uses dpclust3p pre-processing to derive multiplicities for indels
## https://github.com/Wedge-Oxford/dpclust_smchet_docker
writeDPfile.indel <- function(dpFile="dpInput.txt",
                              sampleName,
                              battenberg_rho_psi_file,
                              VCFFILE,
                              BBFILE,
                              GENDER)
{
  vcfdat <- read.table(VCFFILE)
  loci_file = "loci.indel.txt"
  createLociFile.indel(vcfdat, loci_file, 1,2,4,5)
  ## Create allelecounts file
  allelecounts_file = "alleleCounts.indel.txt"
  createAlleleCountsFile.indel(vcfdat, allelecounts_file)
  suppressWarnings(runGetDirichletProcessInfo(loci_file=loci_file,
                                              allele_frequencies_file=allelecounts_file,
                                              cellularity_file=battenberg_rho_psi_file,
                                              subclone_file=BBFILE,
                                              gender=GENDER,
                                              SNP.phase.file="NA",
                                              mut.phase.file="NA",
                                              output_file=dpFile))
}

## dpclust3p is used to see Battenberg (BB) output files
## the following function takes an ASCAT file and format it as a BB-like file
writeBBlike.ASCAT <- function(ascat,ASCATFILE)
  {
    cn <- c("chr", "startpos", "endpos", "BAF", "pval",
            "LogR", "ntot", "nMaj1_A", "nMin1_A", "frac1_A",
            "nMaj2_A", "nMin2_A", "frac2_A", "SDfrac_A",
            "SDfrac_A_BS", "frac1_A_0.025", "frac1_A_0.975", "nMaj1_B",
            "nMin1_B", "frac1_B", "nMaj2_B", "nMin2_B", "frac2_B", "SDfrac_B",
            "SDfrac_B_BS", "frac1_B_0.025", "frac1_B_0.975", "nMaj1_C", "nMin1_C",
            "frac1_C", "nMaj2_C", "nMin2_C", "frac2_C", "SDfrac_C",
            "SDfrac_C_BS", "frac1_C_0.025", "frac1_C_0.975", "nMaj1_D",
            "nMin1_D", "frac1_D", "nMaj2_D", "nMin2_D", "frac2_D", "SDfrac_D",
            "SDfrac_D_BS", "frac1_D_0.025", "frac1_D_0.975", "nMaj1_E",
            "nMin1_E", "frac1_E", "nMaj2_E", "nMin2_E", "frac2_E", "SDfrac_E",
            "SDfrac_E_BS", "frac1_E_0.025", "frac1_E_0.975", "nMaj1_F", "nMin1_F",
            "frac1_F", "nMaj2_F", "nMin2_F", "frac2_F", "SDfrac_F", "SDfrac_F_BS",
            "frac1_F_0.02frac1_F_0.975")
    bb <- matrix("NA",nrow(ascat),length(cn))
    colnames(bb) <- cn
    bb[,1] <- as.character(ascat[,"Chromosome"])
    bb[,2] <- as.character(ascat[,"Start.Position"])
    bb[,3] <- as.character(ascat[,"End.Position"])
    bb[,5] <- rep("1",nrow(ascat))
    bb[,8] <- as.character(ascat[,"Tumour.TCN"]- ascat[,"Tumour.BCN"])
    bb[,9] <- as.character(ascat[,"Tumour.BCN"])
    bb[,10] <- rep("1",nrow(ascat))
    bb <- as.data.frame(bb)
    bb[,1] <- as.factor(as.character(bb[,1]))
    bb[,2] <- as.numeric(as.character(bb[,2]))
    bb[,3] <- as.numeric(as.character(bb[,3]))
    bb[,5] <- as.numeric(as.character(bb[,5]))
    bb[,8] <- as.numeric(as.character(bb[,8]))
    bb[,9] <- as.numeric(as.character(bb[,9]))
    bb[,10] <- as.numeric(as.character(bb[,10]))
    write.table(bb,file=ASCATFILE,sep="\t",col.names=T,row.names=F,quote=F)
    NULL
  }

## write purity/cellularity of the sample
writeCellularity <- function(SAMPLENAME,cellularity_file)
  {
    ww <- as.character(clin$sample)==SAMPLENAME
    cp <- cbind(1-clin$ascatContam[ww],clin$ploidy[ww])
    cp <- cbind(cp,2*(1-cp[1,1])+cp[1,2]*cp[1,1])
    colnames(cp) <- c("cellularity","ploidy","psi")
    write.table(cp,file=cellularity_file,sep="\t",col.names=T,row.names=T,quote=F)
  }

## write global CNA-based variables for samples (purity, ploidy)
writeRhoPsi <- function(SAMPLENAME,rhopsifile)
{
    ww <- as.character(clin$sample)==SAMPLENAME
    cp <- cbind(1-clin$ascatContam[ww],clin$ploidy[ww])
    cp <- cbind(cp,2*(1-cp[1,1])+cp[1,2]*cp[1,1])
    colnames(cp) <- c("cellularity","ploidy","psi")
    RP <- matrix(NA,3,5)
    colnames(RP) <- c("rho",	"psi",	"ploidy",	"distance",	"is.best")
    rownames(RP) <- c("ASCAT","FRAC_GENOME","REF_SEG")
    RP[3,4] <- "Inf"
    RP[3,5] <- "FALSE"
    RP[2,5] <- "TRUE"
    RP[2,4] <- .42
    RP[1,1:3] <- cp[c(1,3,2)]
    RP[2,1:3] <- cp[c(1,3,2)]
    write.table(RP,file=rhopsifile,sep="\t",col.names=T,row.names=T,quote=F)
}

## prepare dpclust pre-processed input to infer multiplicities of SNV mutations
getDPinput <- function(SAMPLENAME,BBKEY)
  {
    VCFFILE <- paste0("/lustre/scratch117/casm/team176/mt17/Sarcomas/vcfs/",SAMPLENAME,".vcf.gz")
    BBFILE <- paste0("/lustre/scratch117/casm/team176/mt17/Sarcomas/",SAMPLENAME,"/",
                     SAMPLENAME,"_",BBKEY,"_ASCAT.txt")
    writeBBlike.ASCAT(allCNAs[[SAMPLENAME]],BBFILE)
    battenberg_rho_psi_file <- paste0("/lustre/scratch117/casm/team176/mt17/Sarcomas/",
                                      SAMPLENAME,"/",SAMPLENAME,"_",BBKEY,"_ASCAT_rho_and_psi.txt")
    cellularity_file <- paste0("/lustre/scratch117/casm/team176/mt17/Sarcomas/",
                               SAMPLENAME,"/",SAMPLENAME,"_",BBKEY,"_ASCAT_cellularity_ploidy.txt")
    writeCellularity(SAMPLENAME,cellularity_file)
    writeRhoPsi(SAMPLENAME,battenberg_rho_psi_file)
    GENDER <- switch(as.character(clin$Gender[clin$sample==SAMPLENAME]),"M"="male","F"="female")
    DPFILE <- paste0(BBKEY,"_dpInput.txt")
    writeDPfile(sampleName=SAMPLENAME,
                dpFile=DPFILE,
                battenberg_rho_psi_file=battenberg_rho_psi_file,
                VCFFILE=VCFFILE,
                BBFILE,
                GENDER=GENDER)
    DP <- read.table(DPFILE,sep="\t",header=T)
  }

## prepare dpclust pre-processed input to infer multiplicities of indel mutations
getDPinput.indel <- function(SAMPLENAME,BBKEY)
  {
    VCFFILE <- paste0("~/Sarcoma/subCNA/INDELDRIVERS/",SAMPLENAME,".txt")
    BBFILE <- paste0("/lustre/scratch117/casm/team176/mt17/Sarcomas/",SAMPLENAME,"/",
                     SAMPLENAME,"_",BBKEY,"_ASCAT.txt")
    writeBBlike.ASCAT(allCNAs[[SAMPLENAME]],BBFILE)
    battenberg_rho_psi_file <- paste0("/lustre/scratch117/casm/team176/mt17/Sarcomas/",
                                      SAMPLENAME,"/",SAMPLENAME,"_",BBKEY,"_ASCAT_rho_and_psi.txt")
    cellularity_file <- paste0("/lustre/scratch117/casm/team176/mt17/Sarcomas/",
                               SAMPLENAME,"/",SAMPLENAME,"_",BBKEY,"_ASCAT_cellularity_ploidy.txt")
    writeCellularity(SAMPLENAME,cellularity_file)
    writeRhoPsi(SAMPLENAME,battenberg_rho_psi_file)
    GENDER <- switch(as.character(clin$Gender[clin$sample==SAMPLENAME]),"M"="male","F"="female")
    DPFILE <- paste0(BBKEY,"_indel_dpInput.txt")
    writeDPfile.indel(sampleName=SAMPLENAME,
                      dpFile=DPFILE,
                      battenberg_rho_psi_file=battenberg_rho_psi_file,
                      VCFFILE=VCFFILE,
                      BBFILE,
                      GENDER=GENDER)
    DP <- read.table(DPFILE,sep="\t",header=T)
  }

## get ref/alt from vcf
reannotateRefAlt <- function(mut,SAMPLENAME)
  {
    VCFFILE <- paste0("/lustre/scratch117/casm/team176/mt17/Sarcomas/vcfs/",SAMPLENAME,".vcf.gz")
    mutvcf <- read.table(VCFFILE)
    themuts <- paste(mut[,1],mut[,3],sep=":")
    themutsvcf <- paste(mutvcf[,1],mutvcf[,2],sep=":")
    rownames(mutvcf) <- themutsvcf
    mut <- cbind(mut,mutvcf[themuts,c(4,5)])
    colnames(mut)[(ncol(mut)-1):ncol(mut)] <- c("ref", "alt")
    mut
  }

## get indel drivers from file
getDriverIndels <- function(SAMPLENAME)
  {
    cat(".")
    truedrive <- calledD[as.character(calledD[,1])==SAMPLENAME,]
    FILEOUT <- paste0("~/Sarcoma/subCNA/INDELDRIVERS/",SAMPLENAME,".txt")
    if(nrow(truedrive)>0)
      {
        dirname <- "/nfs/cancer_ref01/nst_links/live/1263/"
        filename <- paste0(dirname,SAMPLENAME,"/",SAMPLENAME,".pindel.annot.vcf.gz")
        positions <- as.character(truedrive[,3])
        cmd <- paste0("zcat ",filename," | grep ", paste("-e",positions,collapse=" ")," > ", FILEOUT)
        system(cmd)
      }
  }
##tmpNULL <- lapply(samps,getDriverIndels)
###############################################################################



###############################################################################
## PIPELINE to load data, calculate timing and plot figures
###############################################################################
## Library path
RLIB <- "~/R/x86_64-pc-linux-gnu-library/3.4/"
.libPaths( c( .libPaths(), RLIB) )
################################################
## loading libraries
## https://github.com/Wedge-Oxford/dpclust_smchet_docker
library(DPClust)
library(dpclust3p)
library(Rsamtools)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
## ##################################
## load genome from reference file
FASTA <- "~/hs37d5.fa"
dna <- getRefGenome()
names(dna) <- c(1:22,"X","Y","MT")
################################################

################################################
## driver gene + chr:start-end
drivers <-rbind(c("TP53", "chr17", "7571720", "7590868"),
                c("RB1", "chr13", "48877883", "49056026"),
                c("CDKN2A", "chr9", "21967751", "21975132"),
                c("ATRX", "chrX", "76760356", "77041719"),
                c("PTEN", "chr10", "89623195", "89728532"),
                c("MEN1", "chr11", "64570986", "64578188"),
                c("MSH2","chr2", "47630206","47710367"))
################################################
## colours and offset for plotting
COLS <- c("TP53"=rgb(.6,.5,.3),
          "RB1"=rgb(.5,.5,.9),
          "ATRX"=rgb(.9,.5,.6),
          "CDKN2A"=rgb(.5,.2,.2),
          "MEN1"=rgb(.2,.2,.4),
          "PTEN"=rgb(.5,.8,.4),
          "MSH2"=rgb(.3,.5,.2))
OFFSETS <- c("TP53"=1,
             "RB1"=2,
             "ATRX"=3,
             "CDKN2A"=4,
             "MEN1"=5,
             "PTEN"=6,
             "MSH2"=7)
################################################


################################################
## Where all copy number files are located
allsamps <- dir("~/Sarcoma/OUTPUT/")
allsamps <- allsamps[grepl("subclones.txt",allsamps)]
PREFIXES <- sapply(allsamps,function(x)
                 {
                   if(grepl("refit",x))
                     {
                       return(gsub("(.*)_(.*)_(.*)_(.*)","\\2",x))
                     }
                   return("Original")
                 })
SAMPS <- gsub("","",gsub("(.!*)_(.*)","\\1",allsamps))
################################################

###############################################################################
## where all vcf files are located
vcfdir <- "/lustre/scratch117/casm/team176/mt17/Sarcomas/vcfs/"
###############################################################################

###############################################################################
## loading clinical parameters and additional annotation for two samples (historic)
clin <- read.csv("~/Sarcoma/subCNA/clinical.txt")
clin[clin[,1]=="PD26882a","ascatContam"] <- 0.61
clin[clin[,1]=="PD26866a","ascatContam"] <- 0.03
ploidy <- read.csv("~/Sarcoma/subCNA/ploidy.patch.txt")
clin$ploidy <- sapply(1:nrow(clin),function(x)
                      {
                        ploidy[as.character(ploidy$sample)==as.character(clin$sample[x]),"Ploidy"]
                      })
clin$ascatContam <- sapply(1:nrow(clin),function(x)
                      {
                        ploidy[as.character(ploidy$sample)==as.character(clin$sample[x]),"rho"]
                      })
###############################################################################
## loading all ASCAT segments
ascat <- read.table("~/Sarcoma/INPUT/ASCAT-seg.txt",sep="\t",header=T)
###############################################################################

###############################################################################
## separating merged profiles into sample profiles
allCNAs <- lapply(as.character(unique(clin[,1])),function(x) ascat[ascat[,1]==x,])
names(allCNAs) <- as.character(unique(clin[,1]))
###############################################################################
samps <- names(allCNAs)
###############################################################################
## get loh fraction
loh <- sapply(allCNAs,getLOH.fraction)
names(loh) <- samps
###############################################################################


###############################################################################
## (not used) one way to compute the mode of the major allele
modeCN2 <- sapply(allCNAs,function(x)
                 {
                   y <- round(x[,"Tumour.TCN"]-x[,"Tumour.BCN"])
                   y[y>=5] <- 5
                   y <- tapply(1:nrow(x),y,function(z) sum(x[z,"End.Position"]-x[z,"Start.Position"]/1000000))
                   y <- y[c("1","2","3","4","5")]
                   y <- y/sum(y,na.rm=T)
                   y[is.na(y)] <- 0
                   ##ord <- order(y,decreasing=T)
                   ##y <- y[ord]
                   ##y/sum(y)
                   ##y <- y[y>.2]
                   y <- cumsum(y)
                   ww <- which(y>=0.5)[1]
                   names(y)[if(length(ww)==0) which.max(y) else ww]
                 })

## the way to compute the mode of the major allele
modeCN <- sapply(allCNAs,function(x)
                 {
                   y <- round(x[,"Tumour.TCN"]-x[,"Tumour.BCN"])
                   y[y>=5] <- 5 ## "amplified state" all segments>5 get pooled in the same category
                   y <- tapply(1:nrow(x),y,function(z) sum(x[z,"End.Position"]-x[z,"Start.Position"]/1000000))
                   ord <- order(y,decreasing=T)
                   y <- y[ord]
                   ##y/sum(y)
                   ##y <- y[y>.2]
                   names(y)[which.max(y)]
                 })
###############################################################################
## write summary table for the mode of the major allele
##write.table(data.frame(sample=names(modeCN),modeMajor=modeCN),file="~/Sarcoma/subCNA/modeCN.txt",quote=F,col.names=T,row.names=F)
###############################################################################
## pch parameter for plotting depending on mode of major allele
pchs <- sapply(modeCN,function(x) if(x%in%1) 1 else if(x%in%2:3) 2 else 3)
pchs <- as.numeric(modeCN)
###############################################################################


###############################################################################
## writing all computed multiplicities on disk (comment this out to run)
##for(thesamp in samps)
##  {
##    cat(".")
    ##if(!file.exists(paste0("~/Sarcoma/allDPASCAT/",thesamp,".txt.gz")))
##      mytry(write.table(getDPinput(thesamp,PREFIXES[as.character(SAMPS)==thesamp]),
##                        file=gzfile(paste0("~/Sarcoma/allDPASCAT/",thesamp,".txt.gz")),
##                        sep="\t",col.names=T,row.names=F,quote=F),silent=F)
##  }
###############################################################################
## writing all computed multiplicities on disk (comment this out to run) for the two added samples
##for(thesamp in c("PD26882a","PD26866a"))
##  {
##    cat(".")
##    ##if(!file.exists(paste0("~/Sarcoma/allDPASCAT/",thesamp,".txt.gz")))
##      mytry(write.table(getDPinput(thesamp,PREFIXES[as.character(SAMPS)==thesamp]),
##                        file=gzfile(paste0("~/Sarcoma/allDPASCAT/",thesamp,".txt.gz")),
##                        sep="\t",col.names=T,row.names=F,quote=F),silent=F)
##  }
###############################################################################

###############################################################################
## getting all computed multiplicities for indels
allindels <- lapply(samps,function(x)
                    {
                      cat(".")
                      mytry(getDPinput.indel(SAMPLENAME=x,BBKEY=PREFIXES[SAMPS==as.character(x)]))
                    })
###############################################################################


###############################################################################
## read all annotated mutations (including pre-computed multiplicities)
allmuts <- lapply(samps,function(x)
                  {
                    cat(".")
                    file <- paste0("~/Sarcoma/allDPASCAT/",x,".txt.gz")
                    kk <- mytry(read.table(file,sep="\t",header=T))
                  })

##allmuts <- lapply(1:length(allmuts),function(x)
##                  {
##                    cat(".")
##                    reannotateRefAlt(allmuts[[x]],samps[x])
##                  })

## number of mutations per sample
Nmuts <- sapply(allmuts,function(x)
                {
                  ret <- nrow(x)
                  if(is.null(ret)) return(NA)
                  ret
                })
## only (C>T)pG mutations
##allmutsC2T <- lapply(allmuts,function(x)
##                     {
##                       cat(".")
##                       filterCtoTatCpG(x,refCol="ref",altCol="alt")
##                     })
##save(allmutsC2T,file="~/Sarcoma/subCNA/allmutsC2T.Rda")
load(file="~/Sarcoma/subCNA/allmutsC2T.Rda")
###############################################################################



###############################################################################
## removing the metastatic sample from the figures
noMetKeep <- clin$sample!="PD26905c"
###############################################################################
## Figures on calibration using all muts or only (C>T)pG vs. age of the patients
pdf("~/Sarcoma/subCNA/UPS.calibrationTiming.noMet.pdf")
par(mfcol=c(2,2))
xx <- clin$Age_at_diagnosis[noMetKeep]
yy <- (sapply(allmutsC2T,nrow))[noMetKeep]
keep <- yy<2000
xx <- xx[keep]
yy <- (yy[keep])
plot(xx,
     yy,
     xlab="Age at diagnosis (y)",
     ylab="Number of (C>T)pG SNVs",
     pch=19,frame=F)
abline(lm(yy~xx))
yy2 <- (sapply(allmuts,nrow))[noMetKeep]
yy2 <- (yy2[keep])
keep2 <- yy2<11000
xx <- xx[keep2]
yy2 <- yy2[keep2]
plot(xx,
     yy2,
     xlab="Age at diagnosis (y)",
     ylab="TotalnNumber of SNVs",
     pch=19,frame=F)
abline(lm(yy2~xx))
xx <- log10(sapply(allmutsC2T,nrow))[noMetKeep]
yy <- log10(sapply(allmuts,nrow))[noMetKeep]
plot(xx,
     yy,
     xlab="Number of (C>T)pG SNVs",
     ylab="Total number of SNVs",
     pch=19,
     frame=F)
abline(0,1)
##require(MASS)
##abline(lqs(yy~xx))
dev.off()
###############################################################################



###############################################################################
## get SNV and CNA in drivers for all samples
mutsDrive <- lapply(1:length(allmuts),function(x)
                    {
                      cat(".")
                      list(SNV=mytry(getDriversSNV(allmuts[[x]])),
                           CNA=mytry(getDriversCNA(allCNAs[[x]])))
                    })
###############################################################################
## time SNV in drivers
mutsDriveTimed <- lapply(1:length(mutsDrive),function(x)
                         {
                           cat(".")
                           mytry(timeDrive(x))
                         })
###############################################################################
## get indel and CNA drivers for all samples
mutsDrive <- lapply(1:length(allmuts),function(x)
                    {
                      cat(".")
                      list(SNV=mytry(getDriversINDEL(allindels[[x]])),
                           CNA=mytry(getDriversCNA(allCNAs[[x]])))
                    })
## time indel in drivers
mutsDriveTimed2 <- lapply(1:length(mutsDrive),function(x)
                         {
                           cat(".")
                           mytry(timeDrive(x))
                         })
###############################################################################


###############################################################################
## get ploidies
ploidies <- sapply(samps,function(x) clin$ploidy[as.character(clin$sample)==x])
###############################################################################
## Visualise WGD in fraction LOH vs. ploidy dimensions, shaped by mode of major allele
pdf("~/Sarcoma/subCNA/ploidy.vs.lohv4.noMet.pdf",width=7,height=5)
###############################################################################
layout(mat=cbind(c(1,1,1,1,1),c(5,4,2,3,6)),widths=c(4,1),heights=c(1,1,1,1,1))
par(mar=c(4.5,4.5,1,1))
plot(loh[noMetKeep],
     ploidies[noMetKeep],
     pch=pchs,
     xlab="Fraction of genome with LOH",
     ylab="Tumour ploidy",
     cex.lab=1.5,
     cex.axis=1.3,
     frame=F)
abline(v=seq(0,1,.2),h=1:8,lwd=.3,lty=2,col=rgb(.4,.4,.4,.4))
###############################################################################
polygon(y=c(1.6,2.5,1.5,.6),x=c(0,0,1,1),col=rgb(.3,.3,.3,.13),border=NA)
polygon(y=c(3,4,3,2),x=c(0,0,1,1),col=rgb(.5,.1,.1,.13),border=NA)
polygon(y=c(5,8,7,4),x=c(0,0,1,1),col=rgb(.1,.1,.5,.13),border=NA)
###############################################################################
text(.13,2.1,"Diploid",cex=1,col="black")
text(.13,3.3,"WGD",cex=1,col="dark red")
text(.13,5.3,"WGDx2",cex=1,col="dark blue")
###############################################################################
par(mar=c(0,0,0,0))
plot(0,0,col=rgb(0,0,0,0),frame=F,xaxt="n",yaxt="n",xlab="",ylab="")
legend("left",
       legend=c("mode of major allele",
         "1","2","3",
         "4","5"),
       box.col=rgb(0,0,0,0),
       pch=c(1,1,2,3,4,5),
       col=c(rgb(0,0,0,0),
         rgb(0,0,0,.5),
         rgb(.4,.1,.1,.5),
         rgb(.4,.1,.4,.5),
         rgb(.1,.1,.4,.5),
         rgb(.1,.5,.4,.5)))
###############################################################################
dev.off()
###############################################################################


###############################################################################
## TIMING of WGD(s)
## using only (C>T)pG
timingsWGD.CT <- lapply(1:length(allmuts),function(x)
                  {
                    cat(".")
                    mytry(getTiming(allmutsC2T[[x]],allCNAs[[x]],modeCN[x]),silent=T)
                  })

## using all mutations
timingsWGD <- lapply(1:length(allmuts),function(x)
                     {
                    cat(".")
                    mytry(getTiming(allmuts[[x]],allCNAs[[x]],modeCN[x]),silent=T)
                  })
###############################################################################

###############################################################################
## called drivers from the cohort
calledD <- read.csv("~/Sarcoma/subCNA//drivers.txt")
calledD <- calledD[calledD[,7]%in%drivers[,1],]
###############################################################################

###############################################################################
## double check all drivers are found and matched the call list
ddVCF <- lapply(samps,function(SAMPLENAME)
       {
         cat(".")
         VCFFILE <- paste0("/lustre/scratch117/casm/team176/mt17/Sarcomas/vcfs/",
                           SAMPLENAME,
                           ".vcf.gz")
         truedrive <- calledD[as.character(calledD[,1])==SAMPLENAME,]
         muts <- read.table(VCFFILE,sep="\t",header=F)
         themuts <- paste(muts[,1],muts[,2],sep=":")
         thedrive <- paste(truedrive[,2],truedrive[,3],sep=":")
         list(sum(themuts%in%thedrive),length(thedrive))
       })
###############################################################################

###############################################################################
## subset mutations in drivers to called drivers
mutsDriveTimed <- lapply(1:length(mutsDriveTimed),function(x)
                         {
                           cat(".")
                           y <- mutsDriveTimed[[x]]
                           if(length(y)==2)
                             if(!is.null(nrow(y[[1]])) & nrow(y[[1]])>0 & length(nrow(y[[1]]))>0)
                               {
                                 truedrive <- calledD[as.character(calledD[,1])==samps[x],]
                                 SNV <- y[[1]]
                                 themuts <- paste(SNV[,1],SNV[,3],sep=":")
                                 thedrive <- paste(truedrive[,2],truedrive[,3],sep=":")
                                 SNV <- SNV[themuts%in%thedrive,]
                                 y[[1]] <- SNV
                                 return(y)
                               }
                           y
                         })
###############################################################################

###############################################################################
## get annotations of samples for copy number signature groups
cnsigsind <- which(grepl("CNsig",colnames(clin)))
cnsiggroup <- apply(clin[,cnsigsind],1,which.max)
###############################################################################


########################################################
## Timing figure with metastatic sample
pdf("~/Sarcoma/subCNA/timingPlot.ASCAT.CAVEMAN.v5.pdf",width=7,height=3)
layout(mat=cbind(c(2,1,1,1,1,1),c(5,5,3,3,4,4)),widths=c(4,1),heights=c(1,1,2,2,2,2))
par(mar=c(0,4.5,0,0))
retPlot <- plotTiming(timingsWGD,cnsiggroup,samps)
plotSubclonalMuts(mutsDriveTimed,retPlot)
plotSubclonalMuts(mutsDriveTimed2,retPlot,PCH=15)
par(mar=c(0,4.5,0,0))
heatmapDrive(mutsDriveTimed,retPlot)
par(mar=c(0,0,0,0))
mylegend()
mylegend(PCH=15)
dev.off()
########################################################


###############################################################################
## Timing figure no metastatic sample
pdf("~/Sarcoma/subCNA/timingPlot.ASCAT.CAVEMAN.v5.noMet.pdf",width=7,height=3)
layout(mat=cbind(c(2,1,1,1,1,1),c(5,5,3,3,4,4)),widths=c(4,1),heights=c(1,1,2,2,2,2))
par(mar=c(0,4.5,0,0))
retPlot <- plotTiming(lapply(which(noMetKeep),function(xx) timingsWGD[[xx]]),cnsiggroup[noMetKeep],samps[noMetKeep])
plotSubclonalMuts(lapply(which(noMetKeep),function(xx)mutsDriveTimed[[x]]),retPlot)
plotSubclonalMuts(lapply(which(noMetKeep),function(xx) mutsDriveTimed2[xx]]),retPlot,PCH=15)
par(mar=c(0,4.5,0,0))
heatmapDrive(lapply(which(noMetKeep),function(xx)mutsDriveTimed[[xx]]),retPlot)
par(mar=c(0,0,0,0))
mylegend()
mylegend(PCH=15)
dev.off()
########################################################


###############################################################################
## Timing figure using only (C>T)pG SNVs
pdf("~/Sarcoma/subCNA/timingPlot.ASCAT.CAVEMAN.C2T.v1.pdf",width=7,height=3)
layout(mat=cbind(c(2,1,1,1,1,1),c(5,5,3,3,4,4)),widths=c(4,1),heights=c(1,1,2,2,2,2))
par(mar=c(0,4.5,0,0))
retPlot <- plotTiming(timingsWGD.CT,cnsiggroup,samps)
plotSubclonalMuts(mutsDriveTimed,retPlot)
plotSubclonalMuts(mutsDriveTimed2,retPlot,PCH=15)
par(mar=c(0,4.5,0,0))
heatmapDrive(mutsDriveTimed,retPlot)
par(mar=c(0,0,0,0))
mylegend()
mylegend(PCH=15)
dev.off()
########################################################



########################################################
## retrieve age of patients for real-time timing calibration
theages <- clin$Age_at_diagnosis
tgdbd <- sapply(1:length(theages),function(x) {
  if("timingWGD1.tot"%in%names(timingsWGD[[x]])) return(theages[x]-timingsWGD[[x]]$timingWGD1.tot*theages[x])
  return(NA)
})
theages <- theages[noMetKeep]
tgdbd <- tgdbd[noMetKeep]
hists <- cnsiggroup[noMetKeep]
########################################################


########################################################
## which samples have more than Nsnv. SNVs
Nsnv. <- 10000
hypermut <- ifelse(sapply(allmuts,nrow)[noMetKeep]>Nsnv.,15,19)
########################################################


########################################################
## real-time timing plot of WGDs by copy number signatures
pdf("~/Sarcoma/subCNA/timingWGD.UPS.beforediagnosis.CNsig.noMet.pdf",width=4,height=4)
boxpoints.segs(split(tgdbd[retPlot$ord],hists[retPlot$ord]),
               split((1-retPlot$ciwgd1.low)*(theages[retPlot$ord]),hists[retPlot$ord]),
               split((1-retPlot$ciwgd1.high)*(theages[retPlot$ord]),hists[retPlot$ord]),
               split(hists[retPlot$ord],hists[retPlot$ord]),
               split(hypermut[retPlot$ord],hists[retPlot$ord]),
               ylim=c(0,80),ylab="years before diagnosis",
               cex.axis=.7,cex.lab=.7)
dev.off()
########################################################




########################################################
## save data and write summary data frame of signature, timing and mode values per sample
save.image("~/Sarcoma/subCNA/imag.timing.USARC.Rda")
df4chris <- data.frame(sample=samps[retPlot$ord],
                    timingFirstWGD=retPlot$tWGD1s,
                    timingSecondWGD=retPlot$tWGD2s,
                    modeCN=modeCN[retPlot$ord],
                    age=theages[retPlot$ord],
                    timeBeforeDiagnosis=tgdbd[retPlot$ord],
                    signature=paste0("CNsig",hists[retPlot$ord]))
write.table(df4chris,col.names=T,row.names=F,quote=F,sep="\t",file="~/Sarcoma/subCNA/timing.df4chris.USARC.v3.txt")
########################################################
