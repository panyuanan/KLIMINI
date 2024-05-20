
spp<-b_p


N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  #获取模型的 R2
coef(m.fit)*N

bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'#出现频率低于中性群落模型预测的部分
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'#出现频率高于中性群落模型预测的部分
library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 
#grid.text(x=unit(0,'npc')-unit(-1,'lines'), y=unit(0,'npc')-unit(-15,'lines'),label='Mean Relative Abundance (log)', gp=gpar(fontface=2)) 
#grid.text(round(coef(m.fit)*N),x=unit(0,'npc')-unit(-5,'lines'), y=unit(0,'npc')-unit(-15,'lines'),gp=gpar(fontface=2)) 
#grid.text(label = "Nm=",x=unit(0,'npc')-unit(-3,'lines'), y=unit(0,'npc')-unit(-15,'lines'),gp=gpar(fontface=2))
#grid.text(round(Rsqr,2),x=unit(0,'npc')-unit(-5,'lines'), y=unit(0,'npc')-unit(-16,'lines'),gp=gpar(fontface=2))
#grid.text(label = "Rsqr=",x=unit(0,'npc')-unit(-3,'lines'), y=unit(0,'npc')-unit(-16,'lines'),gp=gpar(fontface=2))
draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)

#loop
b_c<-bac[1:2,]
f_c<-fungi[1:2,]
b_v1<-b_p1[1:3,]
f_v1<-f_p1[1:3,]
b_v3<-b_p1[1:3,]
f_v3<-f_p1[1:3,]
b_v10<-b_p10[1:3,]
f_v10<-f_p10[1:3,]
b_i1<-b_p1[4:5,]
f_i1<-f_p1[4:5,]
b_i3<-b_p3[4:6,]
f_i3<-f_p3[4:6,]
b_i10<-b_p10[4:6,]
f_i10<-f_p10[4:6,]
dataframes_list <- list(b_c,b_p,b_v1,b_v3,b_v10,b_i1,b_i3,b_i10,b_ns,b_ne,
                        f_c,f_p,f_p1,f_p3,f_p10,f_i1,f_i3,f_i10,f_ns,f_ne)

# Initialize an empty dataframe to store results
results_df <- data.frame(Dataframe = character(), Rsqr = numeric(), Coef_mN = numeric(), stringsAsFactors = FALSE)

# Loop over each dataframe
for (i in seq_along(dataframes_list)) {
  # Extract the dataframe
  spp <- dataframes_list[[i]]
  
  # Calculations
  N <- mean(apply(spp, 1, sum))
  p.m <- apply(spp, 2, mean)
  p.m <- p.m[p.m != 0]
  p <- p.m / N
  spp.bi <- 1 * (spp > 0)
  freq <- apply(spp.bi, 2, mean)
  freq <- freq[freq != 0]
  
  C <- merge(p, freq, by = 0)
  C <- C[order(C[, 2]), ]
  C <- as.data.frame(C)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))), ]
  p <- C.0[, 2]
  freq <- C.0[, 3]
  names(p) <- C.0[, 1]
  names(freq) <- C.0[, 1]
  
  d <- 1 / N
  m.fit <- nlsLM(freq ~ pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE), start = list(m = 0.1))
  
  m.ci <- confint(m.fit, 'm', level = 0.95)
  freq.pred <- pbeta(d, N * coef(m.fit) * p, N * coef(m.fit) * (1 - p), lower.tail = FALSE)
  pred.ci <- binconf(freq.pred * nrow(spp), nrow(spp), alpha = 0.05, method = "wilson", return.df = TRUE)
  Rsqr <- 1 - (sum((freq - freq.pred) ^ 2)) / (sum((freq - mean(freq)) ^ 2))
  
  # Store results in the dataframe
  results_df <- rbind(results_df, data.frame(Dataframe = paste0("df", i), Rsqr = Rsqr, Coef_mN = coef(m.fit) * N))
}

# Print the results dataframe
print(results_df)