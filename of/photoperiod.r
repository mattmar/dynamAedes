
plot(NA,xlab='Day of the year',ylab='day length [h]',xlim=c(3,21),ylim=c(1,365))
col<-colorRampPalette(c("red","yellow","blue"),alpha=0.5,bias=0.8)(65)
abline(v=11,lwd=5)
for(lat in 1:65) {
	jd <- JD(seq(ISOdate(2019,1,1),ISOdate(2019,12,31),by='day'))
	points(y=1:365,x=daylength(lat,5,jd,1)[,3],col=col[lat])
}
y=1
for(lat in 1:65) {
	y=y+5
	text(y=y,x=daylength(lat,5,jd,1)[,3][1],label=lat,cex=0.8)
}
