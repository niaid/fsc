# scripts provided by Mani 

# this fn is somewhat obsolete, use v1v2 fn below instead.
parse.M1 <- function(day)
{
  a = read.delim("M1.exported.txt", stringsAsFactors=F)
  stopifnot(grepl(".*M1 H1N1-", a$Sample))
  a = a[!grepl("CHI-007",a$Sample),]

  a$subject = do.call(rbind, strsplit(a$Sample, "-"))[,2]
  a$day = gsub("day ", "day", do.call(rbind, strsplit(a$Sample, "-"))[,3])
  print(table(a$day))
  a = a[a$day==day,]
  print(table(a$day))

  print(colnames(a))
  a$ID315 = a$Monocytes.Single.cells.Alive.CD14..CD16..class..Freq..of.Parent / a$Monocytes.Single.cells.Alive.Count
  a$ID300 = a$Monocytes.Single.cells.Alive.CD14..CD16..interm.Freq..of.Parent / a$Monocytes.Single.cells.Alive.Count
  a$ID334 = a$Monocytes.Single.cells.Alive.CD14.CD16...non.class..Freq..of.Parent / a$Monocytes.Single.cells.Alive.Count

  b = as.matrix(a[,c('ID315','ID300','ID334')])
  rownames(b) = paste0("X",a$subject)
  stopifnot(b >= 0, b <= 1)
  b = t(b)
  #write.table(file='tmp.M1.day0.txt', b, quote=F, sep="\t")
  b
}

# ID299 is the CD45+ subset of "Monocytes/Single.cells/Alive" cell population. The three monocyte subsets are expressed as fraction of the former ID299 population in this fn, instead of the latter parent population as in the above fn. Using ID299 is more consistent with the description in the virtual flow manuscript of how all CHI cell populations are expressed as fraction of "higher-level CD45+ leukocyte or lymphocyte" gated populations. 
parse.M1.v1v2 <- function(day, debugcheck=(day=='day0'), controls=F, counts=F)
{
  v1 = read.delim("M1.exported.txt", stringsAsFactors=F) #contains ID3xxx subsets (Nr. of the fraction)
  v2 = read.delim("M1.v2.exported.txt", stringsAsFactors=F) #contains ID299 (Dr. of the fraction)
  stopifnot(v1$Sample == v2$Sample)
  a = v1
  a$ID299 = v2$ID299
  stopifnot(grepl(".*M1 H1N1-", a$Sample))
  if (!controls) {  a = a[!grepl("CHI-007",a$Sample),];  } else {  a$Sample = gsub("CHI-007","CHI007",a$Sample);  }

  if (debugcheck) #check if common cols. are same, and how different ID299 vs. its parents are. 
  {
    v1 = v1[,intersect(colnames(v1),colnames(v2))]
    v2 = v2[,intersect(colnames(v1),colnames(v2))]
    print(all.equal(v1,v2))
    source("~/codes/Rfns/helptools.R")
    ridx = arrinds(v1 != v2)[,1]
    stopifnot(v1[-ridx,]==v2[-ridx,], v1[ridx,1:2]==v2[ridx,1:2])
    print(cbind(v1[ridx,],v2[ridx,-(1:2)])) #excluding two bridge sample measurements, only sample 215.day0 has some minor discrepancies between v1 and v2
    print(summary(a$ID299 / a$Monocytes.Single.cells.Alive.Count)) #bummer: these two cell counts are not that different, so results should not change much!
  }

  a$subject = do.call(rbind, strsplit(a$Sample, "-"))[,2]
  a$day = gsub("day ", "day", do.call(rbind, strsplit(a$Sample, "-"))[,3])
  print(table(a$day))
  if (!controls) {  a = a[a$day==day,];  } else {  a = a[a$subject=="CHI007",];  }
  print(table(a$day))
  print(table(a$subject))

  print(colnames(a))
  if (counts) {  a$ID299=1;  }
  a$ID315 = a$Monocytes.Single.cells.Alive.CD14..CD16..class..Freq..of.Parent / a$ID299
  a$ID300 = a$Monocytes.Single.cells.Alive.CD14..CD16..interm.Freq..of.Parent / a$ID299
  a$ID334 = a$Monocytes.Single.cells.Alive.CD14.CD16...non.class..Freq..of.Parent / a$ID299

  b = as.matrix(a[,c('ID315','ID300','ID334')])
  rownames(b) = paste0("X",a$subject,ifelse(grepl("^day",a$day),"",paste0(".",a$day)))
  stopifnot(b >= 0, all(b <= 1) || counts)
  b = t(b)
  #write.table(file='tmp.M1.day0.txt', b, quote=F, sep="\t")
  b
}

parse.T4 <- function(day, controls=F, counts=F)
{
  hdr = read.delim("T4.exported.txt", stringsAsFactors=F, header=F, nrows=1)
  a = read.delim("T4.exported.txt", stringsAsFactors=F)
  stopifnot(grepl(".*T4  H1N1-", a$Sample), length(hdr)==ncol(a))
  if (!controls) {  a = a[!grepl("CHI-007",a$Sample),];  } else {  a$Sample = gsub("CHI-007","CHI007",a$Sample);  }

  # get the Tfh1 and Tfh1/Tfh17 cell popns, along with their CD278+PD1+ subsets.
  cidx = c(IDhlg = which(hdr == "Lymphocyte/single cells/Alive cells,Count"), #higher-level gate in the hierarchy
          ID242 = grep("Tfh1,Count$",hdr),
          ID244 = grep("Tfh1\\/.*CD278\\+.*PD\\-1\\+,Count$",hdr),
          ID247 = grep("Tfh1 Tfh17,Count$",hdr),
          ID249 = grep("Tfh1 Tfh17\\/.*CD278\\+.*PD\\-1\\+,Count$",hdr))
  stopifnot(length(cidx)==5)
  print(cbind(names(cidx), t(hdr[cidx])))
  stopifnot(make.names(hdr)[cidx] == colnames(a)[cidx])

  a$subject = do.call(rbind, strsplit(a$Sample, "-"))[,2]
  a$day = gsub("day ", "day", do.call(rbind, strsplit(a$Sample, "-"))[,3])
  print(table(a$day))
  if (!controls) {  a = a[a$day==day,];  } else {  a = a[a$subject=="CHI007",];  }
  print(table(a$day))
  print(table(a$subject))

  if (counts) {  a[,cidx['IDhlg']]=1;  }
  stopifnot(names(cidx)[1]=="IDhlg")
  b = matrix(NA, nrow(a), length(cidx)-1, dimnames=list(paste0("X",a$subject,ifelse(grepl("^day",a$day),"",paste0(".",a$day))), names(cidx)[-1]))
  for (i in names(cidx)[-1])
  {
    b[,i] = a[,cidx[i]]/a[,cidx['IDhlg']]
  }
  stopifnot(b >= 0, all(b <= 1) || counts)
  b = t(b)
  #write.table(file='tmp.T4.day0.txt', b, quote=F, sep="\t")
  b
}

# sapply(c('day0','day1','day7'), function(x) {  parse.M1T4(x);  })
parse.M1T4 <- function(day="day0")
{
  M1 = parse.M1.v1v2(day)
  T4 = parse.T4(day)
  print(dim(M1))
  print(dim(T4))
  stopifnot(colnames(M1)==colnames(T4), nrow(M1)==3, nrow(T4)==4, !(rownames(M1) %in% rownames(T4)), !duplicated(colnames(M1)))
  write.table(file=paste0('tmp.M1T4.',gsub(" ","",day),'.txt'), rbind(M1,T4), quote=F, sep="\t")
}

parse.M1T4.controls <- function()
{
  dummy.day = 'day0'
  M1 = parse.M1.v1v2(dummy.day, controls=T)
  T4 = parse.T4(dummy.day, controls=T)
  print(dim(M1))
  print(dim(T4))
  stopifnot(colnames(M1)==colnames(T4), nrow(M1)==3, nrow(T4)==4, !(rownames(M1) %in% rownames(T4)), !duplicated(colnames(M1)))
  rbind(M1,T4)
}

parse.M1T4.counts <- function(day='day0')
{
  M1 = parse.M1.v1v2(day, counts=T)
  T4 = parse.T4(day, counts=T)
  print(dim(M1))
  print(dim(T4))
  stopifnot(colnames(M1)==colnames(T4), nrow(M1)==3, nrow(T4)==4, !(rownames(M1) %in% rownames(T4)), !duplicated(colnames(M1)))
  rbind(M1,T4)
}

# day7-day0 data
day7minusday0 <- function()
{
  day0 = read.delim('M1T4.day0.txt', stringsAsFactors=F)
  day7 = read.delim('M1T4.day7.txt', stringsAsFactors=F)
  stopifnot(identical(dimnames(day0), dimnames(day7)))
  write.table(file='M1T4.day7minusday0.txt', day7-day0, sep="\t", quote=F)
}

# alternative gating based on different markers (CD64 HLA-DR instead of CD14 CD16) for defining monocyte subsets such as non-classical monocytes.
alternative.gating.checks <- function()
{
  a = read.delim("newpops/M1.v2.exported.txt", stringsAsFactors=F)
  a = a[grepl("day 0",a$Sample),]
  a$subject = do.call(rbind, strsplit(a$Sample, "-"))[,2]
  a$subject = make.names(a$subject)
  stopifnot(!duplicated(a$subject))
  a$ID349 = a$ID349/a$ID299
  a$ID353 = a$ID353/a$ID299
  a$ID357 = a$ID357/a$ID299
  b = a[,c('ID349','ID353','ID357')]
  rownames(b) = a$subject
  write.table(t(b), 'newpops/M1altgate.day0.txt', quote=F, sep="\t")
  colnames(b) = paste('obs',colnames(b),'day0',sep='.')

  source("~/codes/Rfns/helptools.R")
  a = read.delim('newpops/obspreddata.txt', row.names=1)
  a = a[,grepl("(ID315|ID300|ID334).day0",colnames(a))]
  d = rbind.fill.rownames(as.data.frame(t(a)), as.data.frame(t(b)))
  d = as.data.frame(t(d))

  if (is.null(d$pred.ID315.day0))
  {
    stopifnot(is.null(d$pred.ID300.day0), !is.null(d$pred.ID334.day0))
    d$pred.ID315.day0 = d$pred.ID300.day0 = NA
  }

  library(ggplot2)
  library(gridExtra)
  pdf(file='newpops/M1altgate.checks.pdf')
  print(ggplot(d, aes(obs.ID315.day0, obs.ID353.day0)) + geom_point())
  print(ggplot(d, aes(obs.ID300.day0, obs.ID349.day0)) + geom_point())
  pvstr = sprintf('(P=%.2g)', cor.test(d$obs.ID334.day0, d$obs.ID357.day0, method='spearman')$p.value);  pvstr=""
  print(gp1 <- ggplot(d, aes(obs.ID334.day0, obs.ID357.day0)) + geom_point() + labs(title=pvstr))
  print(ggplot(d, aes(pred.ID315.day0, obs.ID353.day0)) + geom_point())
  print(ggplot(d, aes(pred.ID300.day0, obs.ID349.day0)) + geom_point())
  pvstr = sprintf('(P=%.2g)', cor.test(d$pred.ID334.day0, d$obs.ID357.day0, method='spearman')$p.value);  pvstr=""
  print(gp2 <- ggplot(d, aes(pred.ID334.day0, obs.ID357.day0)) + geom_point() + labs(title=pvstr))
  grid.arrange(gp1, gp2, nrow=2, ncol=2)
  dev.off()

}

