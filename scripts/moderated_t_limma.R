# Validation only: exported Python design exactly matches canonical OLS.
args <- commandArgs(trailingOnly=TRUE)
out <- args[1]
suppressPackageStartupMessages(library(limma))
ytab <- read.csv(file.path(out,'abundance_matrix.csv'),check.names=FALSE)
dtab <- read.csv(file.path(out,'design_matrix.csv'),check.names=FALSE)
y <- as.matrix(ytab[,-1]); rownames(y) <- ytab$row_position
X <- as.matrix(dtab[,-1]); rownames(X) <- dtab$sample_id
stopifnot(identical(colnames(y),rownames(X)), all(is.finite(X)))
fit <- lmFit(y,design=X,method='ls')
ordinary <- fit$coefficients[,2]/(fit$stdev.unscaled[,2]*fit$sigma)
moderated <- eBayes(fit,trend=FALSE,robust=FALSE,proportion=0.01,stdev.coef.lim=c(0.1,4))
tab <- data.frame(row_position=ytab$row_position,coefficient=fit$coefficients[,2],
 ordinary_t_limma=ordinary,moderated_t=moderated$t[,2],moderated_p=moderated$p.value[,2],
 moderated_q=p.adjust(moderated$p.value[,2],method='BH'),residual_variance=fit$sigma^2,
 posterior_variance=moderated$s2.post,df_residual=fit$df.residual,
 df_total=moderated$df.total,average_abundance=fit$Amean)
write.csv(tab,file.path(out,'limma_statistics.csv'),row.names=FALSE)
write.csv(data.frame(df_prior=moderated$df.prior,s2_prior=moderated$s2.prior),file.path(out,'limma_hyperparameters.csv'),row.names=FALSE)
saveRDS(moderated,file.path(out,'limma_fit.rds'))
writeLines(c(R.version.string,paste('limma',as.character(packageVersion('limma'))),capture.output(sessionInfo())),file.path(out,'R_session_info.txt'))
