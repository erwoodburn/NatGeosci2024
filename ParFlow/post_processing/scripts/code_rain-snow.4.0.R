#ERW 12-1-23 
#Read in East River 7 year forcing (wy15-21), determine rain vs snow partitions

library("MASS")
library(graphics)


#First Read indicator and pfb binary code
dir_in=sprintf("/pscratch/sd/w/woodburn/er_deltat_fall22/forcing")
dir_in
setwd(dir_in)

source("PFB-ReadFcn.R")
 
 nx=150
 ny=170
 nz=5

 T_frz=273.66
 
#Second load watershed mask, which requires a few steps because it's flipped x and y
#Load subbasin mask
subbasin=matrix(0,nrow=170,ncol=150)
subbasin_f=matrix(0,nrow=150,ncol=170)
subbasin_fa=matrix(0,nrow=150,ncol=170)
subbasin=as.matrix(read.table(file="subbasin_indicator.nh.txt"))

mask_count=0

count_f=170+1
for (j in 1:170)  {	
  count_f=count_f-1
  for (i in 1:150)  {
    subbasin_f[i,count_f]=subbasin[j,i]

      if (subbasin[j,i]>0) {
        mask_count=mask_count+1
    }

  }
}

cat(mask_count)

count_fa=170+1
for (j in 1:170)  {	
  count_fa=count_fa-1
  for (i in 1:150)  {
    subbasin_fa[i,j]=subbasin_f[i,count_fa]
  }
}
dim(subbasin_fa)





#for (wy in 2015:2022)  {
for (wy in seq(2016,2022))  {
	

dir_in=sprintf("/pscratch/sd/w/woodburn/er_deltat_fall22/forcing/pt_pfb_only_wy%i",wy)
dir_in
setwd(dir_in)


 precip_sum_all=array(0,dim=c(8760))
 snowfall_sum_all=array(0,dim=c(8760))
 rainfall_sum_all=array(0,dim=c(8760))
 
 
 #Main time loop here
 t_end=8759
 if ((wy==2016)||(wy==2020)) {
 	t_end=8783
 }

 for (t in 1:t_end) {
	 
   
   fin=sprintf("/pscratch/sd/w/woodburn/er_deltat_fall22/forcing/pt_pfb_only_wy%i/NLDAS.APCP.%06d.pfb",wy,t) 
   precip=readpfb(fin,F)
   
   fin_temp=sprintf("/pscratch/sd/w/woodburn/er_deltat_fall22/forcing/pt_pfb_only_wy%i/NLDAS.Temp.%06d.pfb",wy,t)    
   T_air=readpfb(fin_temp,F)

      
   for (j in 1:ny)  {	
     for (i in 1:nx)  {
       if (subbasin_fa[i,j]>0) {	
  
         precip_sum_all[t]=precip_sum_all[t]+precip[i,j,1]
         
         rainfall_incr=0
         snowfall_incr=0
         
         if (T_air[i,j,1]<T_frz){
           snowfall_sum_all[t]=snowfall_sum_all[t]+precip[i,j,1]
         }
         if ((T_air[i,j,1]>=T_frz)&&(T_air[i,j,1]<=(T_frz+2))){
           rainfall_incr=(-54.632+(0.2*T_air[i,j,1]))*precip[i,j,1]
           rainfall_sum_all[t]=rainfall_sum_all[t]+rainfall_incr
           snowfall_incr=precip[i,j,1]-rainfall_incr
           snowfall_sum_all[t]=snowfall_sum_all[t]+snowfall_incr
         }
         if (T_air[i,j,1]>T_frz+2){
           rainfall_sum_all[t]=rainfall_sum_all[t]+precip[i,j,1]
         }
         
       }#endif mask
     }}#end nx ny
   
   print(t)
   
 }#endt
     
dir_out_p=sprintf("/pscratch/sd/w/woodburn/er_deltat_fall22/forcing/precip_sum.4.0_wy%i.txt",wy)
dir_out_s=sprintf("/pscratch/sd/w/woodburn/er_deltat_fall22/forcing/snowfall_sum.4.0_wy%i.txt",wy)
dir_out_r=sprintf("/pscratch/sd/w/woodburn/er_deltat_fall22/forcing/rainfall_sum.4.0_wy%i.txt",wy)

#These are annual sums, in mm/s; convert to mm per year given these were at hourly intervals 
#Divide by mask count because need a spatial average (currently spatial sum)
 write.matrix(precip_sum_all[1:t_end]*3600.0/mask_count, file=dir_out_p)
 write.matrix(snowfall_sum_all[1:t_end]*3600.0/mask_count, file=dir_out_s)
 write.matrix(rainfall_sum_all[1:t_end]*3600.0/mask_count, file=dir_out_r)

a=sum(precip_sum_all[1:t_end])
b=sum(snowfall_sum_all[1:t_end])
c=sum(rainfall_sum_all[1:t_end])
abc=paste(a,b,c,sep=" ")
abc

dir_out_all=sprintf("/pscratch/sd/w/woodburn/er_deltat_fall22/forcing/annual_sums.precip_snowfall_rainfall.4.0_wy%i.txt",wy,wy)
write.matrix(abc, file = dir_out_all)

} #end wy
