# Depends on our data, e.g we might use robust linear model or linear regresstion model. 
# e.g. robust linear model.  
RLMtest = function(meth_matrix, methcol, pheno, age, gender, plate, CD8, CD4, NK, B, Mono, Neu) {mod = try(rlm(meth_matrix[, methcol] ~ pheno+age+gender+plate+CD8+CD4+NK+B+Mono+Neu, maxit=200))
	cf = try(coeftest(mod, vcov=vcovHC(mod, type="HC0")))

if (class(cf)=="try-error") {
  bad <- as.numeric(rep(NA, 3))
  names(bad)<- c("Estimate", "Std. Error", "Pr(>|z|)")
  bad
}
else{
  cf[2, c("Estimate", "Std. Error", "Pr(>|z|)")]
}
}

res <- mclapply(setNames(seq_len(ncol(mv1.t)), dimnames(mv1.t)[[2]]), RLMtest, meth_matrix=mv1.t, pheno=pheno, age=pd$age, gender=pd$gender, plate=pd$Sample_Plate, CD8=pd$CD8T, CD4=pd$CD4T, NK=pd$NK, B=pd$Bcell, Mono=pd$Mono, Neu=pd$Neu)

setattr(res, 'class', 'data.frame')
setattr(res, "row.names", c(NA_integer_,4))
setattr(res, "names", make.names(names(res), unique=TRUE))
probelistnamesB <- names(res)
result <- t(data.table(res))
result<-data.table(result)
result[, probeID := probelistnamesB]
setnames(result, c("BETA","SE", "P_VAL", "probeID")) # rename columns
setcolorder(result, c("probeID","BETA","SE", "P_VAL"))
result$padj <- p.adjust(result$P_VAL, method ="BH")


