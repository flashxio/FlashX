# m is a D x I matrix. Each column is a training instance and each instance has D features.
# labels contains I labels, one for each training instance
# k specifies #columns of the output matrix
# The output is a D x k matrix.
LOL <- function(m, labels, k) {
	counts <- fm.table(labels)
	num.labels <- length(counts$val)
	num.features <- dim(m)[1]
	nv <- k - (num.labels - 1)
	gr.sum <- fm.groupby(m, 1, fm.as.factor(labels, 2), fm.bo.add)
	gr.mean <- fm.mapply.row(gr.sum, counts$Freq, fm.bo.div, FALSE)
	diff <- fm.get.cols(gr.mean, 1) - fm.get.cols(gr.mean, 2)
	svd <- fm.svd(m, nv=0, nu=nv)
	fm.cbind(diff, svd$u)
}
