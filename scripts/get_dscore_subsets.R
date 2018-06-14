#########################################################
# To Run: 
# Rscript get_dscore_subsets.R file.dlogp.txt p-cutoff dsb-cutoff 

#############################################
# Read in command line arguments 
args <- commandArgs(trailingOnly=TRUE)
dfile=args[1]
pcf=args[2]
predlowcf=args[3]

datf <- read.table(dfile, header=TRUE)

pvalcf=log(pcf, base=10)
dscorecf = qt(pcf, df=4, lower.tail=FALSE)
dnullcf = qt(.7+pcf, df=4, lower.tail=TRUE)

cancHpredL <- datf[(datf$dscore > dscorecf) & (datf$predicted > predlowcf),]
cancHpredL$subset="cancHpredL"
cancHpredH <- datf[(datf$log.p < pvalcf) & (datf$predlogp < pvalcf)  & 
                   (datf$dscore < dnullcf) & (datf$dscore > -1 * dnullcf) & 
                   (datf$predicted > predlowcf),]
cancHpredH$subset="cancHpredH"
predHcancL <- datf[(datf$dscore < -1 * dscorecf) & (datf$predicted > predlowcf),]
predHcancL$subset="cancLpredH"
cancHpredL2 <-datf[(datf$log.p < pvalcf) & (datf$predicted < predlowcf),]
cancHpredL2$subset="cancHpredL2"
    
dat2=rbind(cancHpredL, cancHpredH, predHcancL, cancHpredL2)
write.table(dat2, sub(".dlogp.txt", ".sets", dfile)), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

