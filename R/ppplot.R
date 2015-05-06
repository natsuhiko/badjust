ppplot <-
function(x, col=1, ylim=c(0,10), add=F, main=""){
    Log=T
    alpha=0.05
    bf=F
    fdr=F
    GC=F
    blue=rgb(0:100/100,0:100/100,1)	    
    if(Log){
        FUNC = function(x){-log10(x)}
    }else{
	FUNC = function(x){x}
    }
    if(length(col)==1){col=rep(col,length(x))}
        
    x[x==0&!is.na(x)] = min(x[x>0],na.rm=T)
    col = col[!is.na(FUNC(x))&x<=1]
    x = x[!is.na(FUNC(x))&x<=1]
    col = col[order(x)]
    x = sort(x)
    n = length(x)
    
    # expected p-values
    q = 1:n/(n+1)

    # 95% concentration bands besed on Beta distribution
    a = qbeta(1-0.025, 1:n, n-(1:n)+1)
    b = qbeta(0.025, 1:n, n-(1:n)+1)
    
    # axis transformation
    if(add==F){
    if(Log){
        plot(c(0,FUNC(q[1])), c(0,FUNC(x[1])/ifelse(GC,median(log10(x))/log10(.5),1)), type="n", 
             ylim=ylim, ylab=expression(paste("Observed -", log[10], italic(P))), xlab=expression(paste("Expected -", log[10], italic(P))), yaxt="n")
    }else{
        plot(c(0,1),c(0,1),type="n", ylab="Observed p-value", 
             xlab="Expected p-value")
    }
    axis(2,las=1)
    

    # lines of concentration bands
    #points(FUNC(q), FUNC(a), lty=3, type="l")
    #points(FUNC(q), FUNC(b), lty=3, type="l")
    polygon(c(FUNC(q),rev(FUNC(q))),c(FUNC(a),rev(FUNC(b))),border=NA,col=blue[30])
    }
    # plotting pp-points
    points(FUNC(q), FUNC(x)/ifelse(GC,median(log10(x))/log10(.5),1), col=col, pch=20)

    # a diagonal line
    abline(0,1)
    lgc=round(median(log10(x))/log10(.5)*1000)/1000
    title(main)
    title(sub=substitute(lambda[GC]==lgc,list(lgc=lgc)))
    if(fdr){abline(-log10(fdr),1,col="gray")}
    if(bf){abline(-log10(bf/n),0,col="gray")}
}
