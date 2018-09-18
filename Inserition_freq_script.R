Insertion_matrix_mariner<-read.table("Mariner_insertion_freq")
par(cex.axis=1.25, lwd=3)
hist(Insertion_matrix_mariner$FREQUENCY, ylim=c(0,60))
axis(side = 1, font = 2)
axis(side = 2, font = 2)

Insertion_matrix_maverick<-read.table("Maverick_insertion_freq")
par(cex.axis=1.25, lwd=3)
hist(Insertion_matrix_maverick$FREQUENCY, ylim=c(0,350))
axis(side = 1, font = 2)
axis(side = 2, font = 2)



#### plot rather density function of frequencies


 dMAverick<-density(Insertion_matrix_maverick[,1])
 dmariner<-density(Insertion_matrix_mariner[,1])

 plot(range( dmariner$x, dMAverick$x), range(dmariner$y,dMAverick$y), type = "n", xlab = "x",ylab = "Density",cex = 1,lwd=2,font.axis=2)
 lines(dmariner, col = "red",lwd=2)
 lines(dMAverick, col = "blue",lwd=2)
  abline(v=c(0.05),lty=2)
  