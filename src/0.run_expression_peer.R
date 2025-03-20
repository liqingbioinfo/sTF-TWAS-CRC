# conda activate beforetwas

library("peer")

# 285 samples
exp <- read.table("/data/example_expression_matrix.tsv",header=T,row.names=1)
exp <- as.matrix(exp)
cov <- read.table("./data/covariatesPEER.txt",header=T,row.names=1)
cov <- as.matrix(cov)

# run peer
model=PEER()
PEER_setPhenoMean(model,exp)
PEER_setCovariates(model,cov)
PEER_setNk(model,45) #45 factors for 250≤ N<350
PEER_update(model)
PEER_getNk(model)
factor=PEER_getX(model)
dim(factor)
# 285 53
residuals = PEER_getResiduals(model)

rownames(residuals) <- rownames(exp)
colnames(residuals) <- colnames(exp)

residuals_update <- apply(residuals, 2, function(x) {qnorm( rank(x,ties.method="r") / (length(x)+1)  )})
write.table(residuals_update, file = "example_expression_matrix_rmCovPEER.tsv",row.names=T, quote=F,sep = '\t')

