x=c(0,2,4,5,6,7, 8,10, 14:17, 24,25, 32)
y= lead(x)-x
z= y<2
b= ifelse(lag(z)==TRUE & z==FALSE , TRUE, z)
data.frame(x,y,z,b)

scan.kat3<- function(VCF, max.threshold=200, consecutive=3)
{
  kat<- VCF %>%
    group_by(chr) %>%
    mutate(lead= lead(start.position)-start.position,
           bmth= lead < max.threshold
          # bmth = ifelse(lag(bmth)==TRUE & bmth==FALSE , TRUE, bmth),
          # bmth = (lag(bmth, consecutive-1) == TRUE & bmth == TRUE) | (lead(bmth, consecutive-1) == TRUE & bmth == TRUE)
           )
  bmth<-kat$bmth
  
  burst.gr<-0
  vec=rep(NA, length(bmth))
  last=F
  
  for(i in 1:length(bmth))
  {
    if(bmth[i]==FALSE | is.na(bmth[i]))
    {
      vec[i]<-0
      last<- FALSE
    }
    
    if(bmth[i]==TRUE & last==FALSE & !is.na(bmth[i]))
    {
      burst.gr<- burst.gr +1
      last<-TRUE
      vec[i]<- burst.gr
    }
    
    if(bmth[i]==TRUE & last==TRUE & !is.na(bmth[i]))
    {
      vec[i]<- burst.gr
    }
    
  }
  return(as.tbl(data.frame(VCF, bursts=as.factor(vec))))
  
}