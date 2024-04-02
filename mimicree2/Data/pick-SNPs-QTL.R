rm(list=ls())


args <- commandArgs(TRUE)

path = getwd()
input=c()
output = path
nl=c()
rep=1
d="gamma"

i = -1
find_annotation = FALSE
while (i < I(length(args)-1)) {
  i = i + 2
  
  if (args[i] %in% c("--input","-I")) {
    input=args[i + 1]
  }else if (args[i] %in% c("--loci","-L")) {
    nl=as.numeric(args[i + 1])
  }else if (args[i] %in% c("--output", "-O")) {
    output = args[i + 1]
  }else if (args[i] %in% c("--replicate", "-R")) {
    rep = as.numeric(args[i + 1])
  }else if (args[i] %in% c("--distribution", "-D")) {
    if (args[i + 1] %in% c("gamma","uniform")) {
      d = args[i + 1]
    } else {
      message("\n###ERROR\n# Value of parameter : invalid distribution\n")
    }
  }else{
    message("\n### ERROR ###\n\n Argument not recognized:")
    print(args[i])
    message("\n### ERROR ###\n")
    q("no")
  }
}

setwd(path)
snp=readRDS(file = input)
samp=snp[sample(1:dim(snp)[1],nl,replace = F),]
mimhap=samp[,-194]
write.table(mimhap,file=paste0(output,"/selectedloci_",nl,"_",rep,".mimhap"),quote=F,col.names=F,row.names=F,sep="\t")

effectsize=as.data.frame(cbind(as.character(paste(samp$V1)),samp$V2,samp$V4,samp$effect.size_ori,0,samp$effect.size_ori,0,samp$effect.size_ori,0))
effectsize$V1=as.character(paste(effectsize$V1))
effectsize$V4=as.numeric(paste(effectsize$V4))
effectsize$V5=as.numeric(paste(effectsize$V5))
effectsize$V6=as.numeric(paste(effectsize$V6))
effectsize$V7=as.numeric(paste(effectsize$V7))
effectsize$V8=as.numeric(paste(effectsize$V8))
effectsize$V9=as.numeric(paste(effectsize$V9))
effectsize[,c(4,6,8)]=effectsize[,c(4,6,8)]/sum(effectsize[,4])
if (d=="uniform") effectsize[,c(4,6,8)]=1/nl
effectsize[,c(4,6,8)]=effectsize[,c(4,6,8)]*sample(size = nl,x = c(-1,1),p=c(1,1),replace = T)
write.table(effectsize,file=paste0(output,"/effectsize_",d,"_",nl,"_",rep,".txt"),quote=F,col.names=F,row.names=F,sep="\t")

