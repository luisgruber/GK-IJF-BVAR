library(colorspace)

dire <- "figures/"
if(!dir.exists(dire)) dir.create(dire)

data <- readRDS("results/figure_4.rds")

myrows <- c("MP", "HS0", "DL0", "NG0", "R2D20", "SSVS0")
summaries <- c("median", "IQR")
laymat0 <- c(rep(1,1), rep(2,4), rep(3,4))
laymat0 <- matrix(c(rep(laymat0,8),rep(c(4,rep(5,8)),2)), ncol = 9, byrow = TRUE)
laymat01 <- sort(c(rep(max(laymat0)+c(1,2), 4)))
laymat01 <- matrix(laymat01,nrow = 8, ncol=8,byrow = TRUE)
laymat01 <- rbind(laymat01, 8,8)
laymat00 <- cbind(laymat0, laymat01)
laymat <- laymat00
for(j in 2:length(myrows)){

  laymat <- rbind(laymat, laymat00+(j*max(laymat00)-max(laymat00)))
}
laymat <- rbind(max(laymat)+c(1,rep(2,8), rep(3,8)),max(laymat)+c(1,rep(2,8), rep(3,8)),laymat)
M <- 21
p <- 2
Kp <- p*M+1


pdf(paste0(dire, "figure_4.pdf"), height = 6/8*20, width = 20)
varnames <- dimnames(data$MP_LIT$median)[[2]]
alpha <- 1

layout(laymat)
for(r in myrows){

  if(r=="MP"){
    models <- c("MP_LIT", "SHM")
    model_quote <- models

  }else if(r=="HS0"){
    models <- c("HS", "HS_star")
    model_quote <- c("HS", bquote(HS^"*"))
  }else if(r=="DL0"){
    models <- c("DL", "DL_star")
    model_quote <- c(bquote(DL[a]), bquote(DL[a]^"*"))
  }else if(r=="NG0"){
    models <- c("NG", "NG_star")
    model_quote <- c(bquote(NG[a]), bquote(NG[a]^"*"))
  }else if(r=="R2D20"){
    models <- c("R2D2", "R2D2_star")
    model_quote <- c(bquote(R2D2[a]), bquote(R2D2[a]^"*"))
  }else if(r=="SSVS0"){
    models <- c("SSVS_p", "SSVS_star")
    model_quote <- c(bquote(SSVS[p]), bquote(SSVS[p]^"*"))

  }

  par(mar=c(0,7.5,2,0), xaxs = "i", yaxs = "i") #, xaxs = "i", yaxs = "i"
  plot.new()
  axis(2, at = (c(0:(M-1))+0.5)/(M), labels = varnames,
       tick = FALSE, las = 2, cex.axis=1, pos = 1)

  for(fun in summaries){
    sums <- array(NA,c(2,Kp,M))
    myi <- 0
    for(i in models){

      myi <- myi+1
      sums[myi,,] <- data[[i]][[fun]]# apply(coefs, 2:3, eval(str2lang(fun)))
    }
    for(i in seq.int(length(models))){
      if(fun=="median"){
        zlim <- c(-1,1)*max(abs(sums))
        themin <- min(sums)
        themax <- max(sums)

        colspace <- colorspace::diverge_hcl(1001, alpha = alpha, palette = "Blue-Red")
      }else if(fun=="IQR"){
        zlim <- c(0, max(sums) )
        colspace <- colorspace::sequential_hcl(1001, alpha = alpha, rev = TRUE,
                                               palette = "Reds 2")
      }
      colbreaks <- seq(zlim[1]*1.001, zlim[2]*1.001, len=1002)
      par(mar=c(0,0,2,2))
      image(t(apply(sums[i,,],1,rev)), zlim = zlim ,xaxt = "n", yaxt="n", col = colspace,
            bty="n", main = "")

      mtext(parse(text = model_quote[i]), side = 3, cex = 1.3)
      if((Kp-1)>M){
        abline(v = seq(M,Kp-1, M)/(Kp-1) - 0.5/(Kp-1), xpd = FALSE)
      }

    }

    if(fun=="median") plot.new()

    par(mar=c(1.5,0,0.5,2), xaxs = "r", yaxs = "r")
    plot.new()

    if(fun=="median"){
      if(themin>zlim[1]){
        colbarlim <- length(which(colbreaks<themin))
        adjustedcolspace <- colspace[-c(1:(colbarlim-1))]

        llabels <- round(c(colbreaks[colbarlim],0,tail(colbreaks,1)),2)
        at <- sort(c(0,1,(0-colbreaks[colbarlim])/(tail(colbreaks,1)-colbreaks[colbarlim])))
      }else if(themax<zlim[2]){

        colbarlim <- length(which(colbreaks<themax))
        adjustedcolspace <- colspace[1:(colbarlim+1)]

        llabels <- round(c(head(colbreaks,1),0,colbreaks[colbarlim]),2)
        at <- sort(c(0,1,abs(head(colbreaks,1))/(colbreaks[colbarlim]-head(colbreaks,1))))
      }else adjustedcolspace <- colspace



    }else if(fun=="IQR"){
      adjustedcolspace <- colspace
      llabels <- round(zlim,2)
      at <- c(0,1)
    }

    len <- length(adjustedcolspace)
    xxx <- seq(0,1,length=len+1)
    rect(xxx[1:len], rep(0, len), xxx[-1],
         rep(1, len), col = adjustedcolspace, border = adjustedcolspace)

    axis(1, labels = llabels, at = at,
         tick = FALSE, las=1, cex.axis=1.5, line=-.8)

  }
}
par(mar=c(0,0,0,0), xaxs = "i", yaxs = "i")
plot.new()
par(mar=c(0,0,2.5,2), xaxs = "i", yaxs = "i")
plot.new()
abline(h=1)
mtext("Posterior median", cex = 1.8)
plot.new()
abline(h=1)
mtext("Posterior interquartile range", cex = 1.8)
dev.off()

