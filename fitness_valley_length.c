#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>


# define endcluster 10000


#define MAXITER 30
# define max 3000

int L,kk,maxstart,lmin[endcluster+2],suman,min,minl,maxl,valley,minstart,lserint,
nakopl[12],test,tr,length[11],jj1,tip,y,suml,prohod[1000],nakop[10],
viddaleno[10][100],lll[10],min1[11],lminend,jj,fact[100],aa,
//trajectory[1000][10],projshly[1000][1000],pair[endcluster+2],three[endcluster+2],
numberlink,dead,suma,kill,link[101][101],
chyslo[10000],degree[10000],connect[10000][10000],array[10000][10000],
arrayvpered[1000][1000],arraynazad[1000][1000],pozycija,
pobud,numlet,mindis,mind,distance[10002],previous,
jpoper,vvv,tak[100];
 double dysl,dyslmin,testn,serer, lsered[endcluster+1],seredmin,seredmindown,prob[1000],
 probtime[1000][endcluster+1],probmintime[1000][endcluster+1],denom[1000],denommin[1000],
 probtimeser[1000],probmintimeser[1000],dystime[1000],dysmintime[1000],
 seredne[endcluster+2],probmin[10],lethalsituation,alive,ttime,ttt,gamamax,
 waga[101];
 float doom,RR,deadfit,mu,gama,mustart,comb,minf1,minf2,
 minf3,birthfit,mutfit,mf3,mf1,mf2,prbin,factor,timeser,munew,riz,dysp,sum,
 //prob[TIME+2],
 //concentrationmax[TIME+2],
 fitness[max+1],maxfitnessstart,fitnessstart[max+1],
 //maxf[tmax+1],
 maxfit,
 fitnessstart1[max+1],fitnessstart2[max+1],
 fitnessstart3[max+1],
 //frequency[101][TIME+1],
 frequencytest[1000];
 //ftest[101][tmax]
 
 int running,lser,dist,letal,t,maxmutation,fr,rozpbin,binttt,count,found, clusteryes,maxfound,maxlevy1,maxlevy2,
 maxlevy3,vybirmu,lich,lichf,runing,fint,position,
 numer,maxp,example[max][21],maxtime,spr,ggg, chystyj,roz1,roz2,cilyj,xx,yy,ss,r,
 rr,NN,
 numm, x,level,
 x1,x2,yy1,y2,porjadokkut,g1,g2,man,kut,levy,timme,serlevy,serlevy2,numser,numlevy,ra;

//  float w0[M+2],w[M+2],summa,w0norm[M+2],e,d,d0;

 float q,shuffle,perev[max+2],perev2[max+2],perev3[max+2];

 int zzz,test1,run,povtor,zrobleno[max+1][21],zmina,bb, cc[100][100],bin,perc,gi,
 i1,i0,i2,i3,i4,i5,i6,i7,i8,j2,j3,j4,j5,j6,j7,j8,l1,l2,l3,l4,l5,l6,l7,l8,labelmax,
maxim,u,numerbonn,ll, ffff,cluster, kr,cl,  gy,rapira,add,vuzl,
  i,nn,neww,hh,enrich;

 
 float sigma,zag1,fin,ymax;
 
 // for random number generators

 long int text;
 int ppp,povt,j,x,p,prun,g,bah,startx1,fff, gm;     

 int dod,trial1,trial2,sumalink,conf,zagal,nnn,ii,h,krok,  zag4, zag2,zag3,f,ff;
 
 
 main ()  
{

// for random numbers generations

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7

#define RNMX (1.0-EPS)
int M,jjj;
long k;
static long idum2=123456789;
static long iy=0;
static long iv[NTAB];
float temp;
static long idum;


  // start dlja vypadkovyh chysel
 idum = -1000000;
 void rrr()
{
if (idum <= 0) 
{ 
if (-(idum) < 1) idum=1; 
else idum = -(idum);
idum2=(idum);
for (jjj=NTAB+7;jjj>=0;jjj--) {
k=(idum)/IQ1;
idum=IA1*(idum-k*IQ1)-k*IR1;
if (idum < 0) idum += IM1;
if (jjj < NTAB) iv[jjj] = idum;
}
iy=iv[0];
}
k=(idum)/IQ1; 
idum=IA1*(idum-k*IQ1)-k*IR1; 
if (idum < 0) idum += IM1; 
k=idum2/IQ2;
idum2=IA2*(idum2-k*IQ2)-k*IR2; 
if (idum2 < 0) idum2 += IM2;
jjj=iy/NDIV; 
iy=iv[j]-idum2; 
iv[jjj] = idum; 
if (iy < 1) iy += IMM1;
if ((temp=AM*iy) > RNMX) {
shuffle=RNMX;}
else {
shuffle=temp;}
}


char name7[50];
 char name1[50];
 

FILE  *data7, *data1;
int tri,a,b,vybir,s;      

 srand((unsigned)time(NULL));
 a = 1;
 b = 4; 
 
 
 printf("Enter sequence length L: \n" );
  // L = getchar();
    scanf("%d", &zzz); L=zzz; 
 
 
 sprintf (name7,"ValleylengthL=%d.dat",L);                            
 data7 = fopen(name7,"w+");
 
// kilkist bakterij
for(position=L;position<=L;position=position+1)
{                                               
numberlink=pow(2,position);
void factorial(int aa)
{
fact[aa]=1;     
for(x=1;x<=aa;x++)
{
fact[aa]=fact[aa]*x; 
}
}

void binom(int aa, int bb)
{
fact[aa]=1;     
for(x=1;x<=aa;x++)
{
fact[aa]=fact[aa]*x; 
}
fact[bb]=1;     
for(x=1;x<=bb;x++)
{
fact[bb]=fact[bb]*x; 
}
fact[aa-bb]=1;     
for(x=1;x<=aa-bb;x++)
{
fact[aa-bb]=fact[aa-bb]*x; 
}
cc[aa][bb]=fact[aa]/(fact[bb]*fact[aa-bb]);
}
         
for(i=0;i<1000;++i)
{
prob[i]=0.0;
}

factorial(position);

for(i=0;i<=10;++i)
{binom(position,i);
//prob[i]=0.0;
probmin[i]=0.0;
//printf("distance %d cc %d \n",i,cc[position][i]);
}
//printf("different %d Factorial  %d cc %d \n",numberlink,fact[position],cc[position][2]);

 
//for(gama=0.00;gama<=0.951;gama=gama+0.05)
{
lethalsituation=0;
alive=0;
clusteryes=0;

for(j=0;j<=100;++j)
{
waga[j]=0;
}

 for(i=0;i<1000;++i)
{
denom[i]=0.0;denommin[i]=0.0;
}

for (cluster=1;cluster<=endcluster;++cluster) 
{
    printf("replica %d \n",cluster);
    for(i=0;i<1000;++i)
{
probtime[i][cluster]=0.0;
probmintime[i][cluster]=0.0;
}

    
//printf("cluster %d gamma %f position %d\n",cluster,gama,position);

//fprintf(data11,"%d ",1);
jpoper=1;  
  ////////////////////////////////////     procedura zadannja landscape

for(i=0;i<=100;++i)
 {                                    
tak[i]=0;
}


for(j=1;j<=position;++j)
{                 
example[1][j]=-1;
}
nn=1;

// startova tochka vsih trajektorij

for(kk=1;kk<=position;++kk)
{
   numer=nn;
// zminjujemo po kozhnij pozyciji 0 na 1
// u vsih poperednih konfiguracijah vzhe zbudovanyh                     
   for(j=nn;j>=1;j=j-1)
   { 
   numer=numer+1;   
           for(i=1;i<=position;++i)
            {
            example[numer][i]=example[j][i];
            }
 // vsi inshi chysla zalyshajutsja takymy zh, a na pozyciji kk zmina na 1           
            example[numer][kk]=1;                                     
   }  
 nn=numer;                  
}

maxmutation=0;

for(i=0;i<=numer;++i)
{
fitness[i]=0;
fitnessstart[i]=0;
test=0;
for(j=1;j<=position;++j)
{
if(example[i][j]==1)
{test=test+1;}
}
if(test==position)
{maxmutation=i;}
}
//printf("maxmutation %d\n",maxmutation);

/*
for(j=1;j<=fact[position];++j)
{ 
for(i=1;i<=position;++i)
{trajectory[j][i]=0;}                                    
trajectory[j][0]=1;
}
*/


// graph of connections
for(i=1;i<=numer;++i)
{
chyslo[i]=0;
degree[i]=0;
for(j=1;j<=numer;++j)
{array[i][j]=0;
arrayvpered[i][j]=0;
arraynazad[i][j]=0;
connect[i][j]=0;}
for(j=1;j<=position;++j)
{
    if(example[i][j]==1)
       {chyslo[i]=chyslo[i]+1;}                     
}
}

for(i=1;i<=numer;++i)
{
for(j=i+1;j<=numer;++j)
{                     
if( (chyslo[i]-chyslo[j])==1 ^(chyslo[j]-chyslo[i])==1)
{
    pozycija=0;
    for(jj=1;jj<=position;++jj)
    {
    if (example[i][jj]==1 && example[j][jj]==1)
    {pozycija=pozycija+1;}
    }
//printf("i %d j %d chyslo %d, pozycija %d \n",i,j,chyslo[i],pozycija);    
if( pozycija==(chyslo[i]-1) && (chyslo[i]>chyslo[j]) )
{
degree[i]=degree[i]+1;
degree[j]=degree[j]+1;
connect[i][j]=1;
connect[j][i]=1;
array[i][degree[i]]=j;
array[j][degree[j]]=i;
//printf("je zvjazok \n");   
}
if( pozycija==(chyslo[i]) && (chyslo[i]<chyslo[j]) )
{
degree[i]=degree[i]+1;
degree[j]=degree[j]+1;
connect[i][j]=1;
connect[j][i]=1;
array[i][degree[i]]=j;
array[j][degree[j]]=i;
//printf("je zvjazok \n");   
}
}
}
}


distance[0]=0;
distance[1]=0;

for(j=1;j<=degree[1];++j)
{
jj=array[1][j];                         
distance[jj]=1;
}

for(i=1;i<numer;++i)
{
prohod[i]=0;
}



for(previous=1;previous<=100;++previous)
{
for(i=2;i<=numer;++i)
{
mindis=2000;
mind=0;
for(j=1;j<=degree[i];++j)
{
jj=array[i][j];
{
if(distance[jj]<mindis && distance[jj]>0)
{mindis=distance[jj];
mind=jj;}
}                         
}
// znahodymo minimum z vidstanej 
if(distance[i]==0 && distance[mind]>0)
{
distance[i]=distance[mind]+1;
//printf("i %d mind %d distance %d \n",i,mind,distance[mind]);
}
}
}
     


for(i=0;i<=position;++i)
{
nakop[i]=0;
for(j=0;j<100;++j)
{
viddaleno[i][j]=0;
}}


for(i=1;i<=numer;++i)
{
y=1;yy=1;
// kilkist galuzhen dlkja ruhu vpered
prohod[i]=position-distance[i];
nakop[distance[i]]=nakop[distance[i]]+1;
viddaleno[distance[i]][nakop[distance[i]]]=i;
 for(h=1;h<=degree[i];++h)
 {
 if(distance[array[i][h]]>distance[i])
 {
 arrayvpered[i][y]=array[i][h];
 y=y+1;
 }
 if(distance[array[i][h]]<distance[i])
 {
 arraynazad[i][yy]=array[i][h];
 yy=yy+1;
 }          
}                            
}



again: for(j=1;j<=numer;++j)
{                 
perev[j]=0;
}
//budujemo try versiji  
mf1=0.0;
minf1=10.0;
for(i=1;i<=numer;++i)
{
 vy: ss=rand()%(numer)+1; 
 if(perev[ss]==1)
 { goto vy;}                    
 else
 {
 rrr();
// kozhnomu iz vsemozhlyvyh genotypiv vypadkovo zadajemo fitnes
fitnessstart[ss]=shuffle; 
perev[ss]=1;
}
}
fitnessstart[1]=0.5;
fitnessstart[maxmutation]=1.0;


lmin[cluster]=0;
lsered[cluster]=0.0;

for(dist=1;dist<=10;++dist)
{lll[dist]=0;min1[dist]=0;}

for(dist=0;dist<=10;++dist)
{length[dist]=0;
nakopl[dist]=0;
min1[dist]=0;}


suml=0;
dist=1;
valley=0;
min=10;
lser=0;


for(i=0;i<=10;++i){length[i]=0;}


maxstart=0;
maxfitnessstart=0.0;
for(i=1;i<100;++i)
{
 if(viddaleno[1][i]>0)
 {                 
//printf("vuzol %d fit %f \n",viddaleno[1][i],fitnessstart[viddaleno[1][i]]);                  
if(fitnessstart[viddaleno[1][i]]>maxfitnessstart)
{maxfitnessstart=fitnessstart[viddaleno[1][i]];
maxstart=viddaleno[1][i];
}}
}

//printf("max   vuzol %d fit %f \n",maxstart,maxfitnessstart);

for(i=0;i<position;++i){length[i]=0;}
{
dist=1;
min=10;
//po vsih vuzlah, vidstan jakyh vid 0 rivna dist
for(j=1;j<=cc[position][dist];++j)
{                                                                                      
if(viddaleno[dist][j]==maxstart)
{
jj1=viddaleno[dist][j];
minstart=jj1;
//printf("jj1 %d fit %f max %d \n",jj1,fitnessstart[jj1],maxstart);
if(distance[jj1]<position)
{
///3!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
           
       for(i1=1;i1<=prohod[jj1];++i1)
        {
        j2= arrayvpered[jj1][i1];
        if(distance[j2]==position-1)
        {
      
           if(fitnessstart[minstart]<0.5)
       {valley=valley+1;
      length[valley]=1;
      }
      
          if(fitnessstart[j2]<fitnessstart[minstart])
       {
     valley=valley+1;
     length[valley]=1;
    // printf("male j2 %d fitnes %f valley %d lenth %d \n",j2, fitnessstart[j2],valley, length[valley]);                                                                   
       }
 if(fitnessstart[j2]>fitnessstart[minstart])
       {
     minstart=j2;
     valley=0;
    // printf("novyj start %d fitnes %f\n",j2, fitnessstart[j2]);                                                                   
       }
        maxl=0;
        for(i=0;i<=position;++i)
        {if(length[i]>0){maxl=i;}}
       // if(fitnessstart[minstart]<0.5){maxl=maxl+1;}
        nakopl[maxl]=nakopl[maxl]+1;
        if(maxl<min)
        {min=maxl;}
       suml=suml+1;
    //   printf("sum %d j2 %d  max %d nakop %d \n",suml,j2,maxl,nakopl[maxl]);                                                           
        valley=0;
        for(i=0;i<10;++i){length[i]=0;}
        minstart=jj1;    
     }
     
     
     
 //4!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
            for(i2=1;i2<=prohod[j2];++i2)
            {
            j3= arrayvpered[j2][i2];
             if(distance[j3]==position-1)
        {
      
           if(fitnessstart[minstart]<0.5)
       {valley=valley+1;
      length[valley]=1;
      }
      
          
          if(fitnessstart[j2]<fitnessstart[minstart])
       {
     valley=valley+1;
     length[valley]=1;
  //   printf("male j2 %d fitnes %f lenth %d \n",j2, fitnessstart[j2],valley);                                                                   
       }
 if(fitnessstart[j2]>fitnessstart[minstart])
       {
     minstart=j2;
     valley=0;
    // printf("novyj start %d fitnes %f\n",j2, fitnessstart[j2]);                                                                   
       }
          if(fitnessstart[j3]<fitnessstart[minstart])
       {
     valley=valley+1;
     length[valley]=1;
    // printf("male j3 %d fitnes %f lenth %d \n",j3, fitnessstart[j3],valley);                                                                   
       }
 if(fitnessstart[j3]>fitnessstart[minstart])
       {
     minstart=j3;
     valley=0;
    // printf("novyj start %d fitnes %f\n",j3, fitnessstart[j3]);                                                                   
       }
        maxl=0;
        for(i=0;i<=position;++i)
        {if(length[i]>0){maxl=i;}}
        //if(fitnessstart[minstart]<0.5){maxl=maxl+1;}
        nakopl[maxl]=nakopl[maxl]+1;
        if(maxl<min)
        {min=maxl;}
       suml=suml+1;
      // printf("sum %d j3 %d  max %d\n",suml,j2,maxl);
        valley=0;
        for(i=0;i<10;++i){length[i]=0;}
        minstart=jj1;    
     }
     
     
                     
///// 5!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111                 
            
               for(i3=1;i3<=prohod[j3];++i3)
            {
            j4= arrayvpered[j3][i3];
              if(distance[j4]==position-1)
        {

               if(fitnessstart[jj1]<0.5)
       {valley=valley+1;
      length[valley]=1;
      }
      

          if(fitnessstart[j2]<fitnessstart[minstart])
       {
     valley=valley+1;
     length[valley]=1;
       }
 if(fitnessstart[j2]>fitnessstart[minstart])
       {
     minstart=j2;
     valley=0;
       }
       
          if(fitnessstart[j3]<fitnessstart[minstart])
       {
     valley=valley+1;
     length[valley]=1;
       }
 if(fitnessstart[j3]>fitnessstart[minstart])
       {
     minstart=j3;
     valley=0;
       }
       
          if(fitnessstart[j4]<fitnessstart[minstart])
       {
     valley=valley+1;
     length[valley]=1;
       }
 if(fitnessstart[j4]>fitnessstart[minstart])
       {
     minstart=j4;
     valley=0;
       }

        maxl=0;
        for(i=0;i<=position;++i)
        {            
        //  if(i==4 && length[i]>0){printf("yes %d fit %f !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",length[i],fitnessstart[minstart]);}                 
       if(length[i]>0){maxl=i;}}
       // if(fitnessstart[minstart]<0.5){maxl=maxl+1;}
         if(maxl<min)
        {min=maxl;} 
        suml=suml+1;
        nakopl[maxl]=nakopl[maxl]+1;
    // if(nakopl[maxl]==24) { printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! length %d  min %d\n",maxl,min);}
     // if(maxl>3){ printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!sum %d  maxl %d min %d\n",suml,maxl, min );}
        valley=0;
        for(i=0;i<10;++i){length[i]=0;}
        minstart=jj1;    
     }
        
 //6!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11            
                 
                       for(i4=1;i4<=prohod[j4];++i4)
                       {
                        j5= arrayvpered[j4][i4];

                           if(distance[j5]==position-1)
                           {

                      if(fitnessstart[minstart]<0.5)
       {valley=valley+1;
      length[valley]=1;
      }
      
                                                     if(fitnessstart[j2]<fitnessstart[minstart])
                                                       {
                                                        valley=valley+1;
                                                        length[valley]=1;
                                   //                      printf("male j2 %d fitnes %f lenth %d \n",j2, fitnessstart[j2],valley);    
                                                        }
                                                        if(fitnessstart[j2]>fitnessstart[minstart])
                                                        {
                                                         minstart=j2;
                                                         valley=0;
                                     //                    printf("new peak j2 %d fitnes %f  \n",j2, fitnessstart[j2]); 
                                                         }
                                                             if(fitnessstart[j3]<fitnessstart[minstart])
                                                             {
                                                              valley=valley+1;
                                                              length[valley]=1;
                                       //                         printf("male j3 %d fitnes %f lenth %d \n",j3, fitnessstart[j3],valley);   
                                                              }
                                                              if(fitnessstart[j3]>fitnessstart[minstart])
                                                              {
                                                              minstart=j3;
                                                              valley=0;
                                         //                     printf("new peak j3 %d fitnes %f  \n",j3, fitnessstart[j3]); 
                                                              }
                                                              if(fitnessstart[j4]<fitnessstart[minstart])
                                                              {
                                                              valley=valley+1;
                                                              length[valley]=1;
                                           //                     printf("male j4 %d fitnes %f lenth %d \n",j4, fitnessstart[j4],valley);   
                                                              }
                                                              if(fitnessstart[j4]>fitnessstart[minstart])
                                                              {
                                                              minstart=j4;
                                                              valley=0;
                                             //                 printf("new peak j4 %d fitnes %f  \n",j4, fitnessstart[j4]); 
                                                              }
                                                              if(fitnessstart[j5]<fitnessstart[minstart])
                                                              {
                                                              valley=valley+1;
                                                              length[valley]=1;
                                               //                printf("male j5 %d fitnes %f lenth %d \n",j5, fitnessstart[j5],valley);   
                                                              }
                                                              if(fitnessstart[j5]>fitnessstart[minstart])
                                                              {
                                                              minstart=j5;
                                                              valley=0;
                                                 //              printf("new peak j5 %d fitnes %f  \n",j5, fitnessstart[j5]);   
                                                              
                                                              }
                                                              maxl=0;
                                                              for(i=0;i<=position;++i)
                                                             { if(length[i]>0){maxl=i;}}
                                                            //  if(fitnessstart[minstart]<0.5){maxl=maxl+1;}
                                                               if(maxl<min)
                                                               {min=maxl;}
                                                              suml=suml+1;
                                                              nakopl[maxl]=nakopl[maxl]+1;                                          
                                                             
                                                             // printf("sum %d j5 %d  max %d min %d\n",suml,j5,maxl,min);
                                                              valley=0;
                                                              for(i=0;i<10;++i){length[i]=0;}
                                                              minstart=jj1;    
                                                              }
            
            
            
// 7!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

                      for(i5=1;i5<=prohod[j5];++i5)
                       {
                        j6= arrayvpered[j5][i5];

                       if(distance[j6]==position-1)
                           {

//printf("jj1 %d minstart %d maxstart %d fit %f n",jj1, minstart, maxstart, fitnessstart[minstart]);
                            if(fitnessstart[minstart]<0.5)
       {valley=valley+1;
      length[valley]=1;
      }
      
                                                    if(fitnessstart[j2]<fitnessstart[minstart])
                                                       {
                                                        valley=valley+1;
                                                        length[valley]=1;
                                                        }
                                                        if(fitnessstart[j2]>fitnessstart[minstart])
                                                        {
                                                         minstart=j2;
                                                         valley=0;
                                                         }
                                                             if(fitnessstart[j3]<fitnessstart[minstart])
                                                             {
                                                              valley=valley+1;
                                                              length[valley]=1;
                                                              }
                                                              if(fitnessstart[j3]>fitnessstart[minstart])
                                                              {
                                                              minstart=j3;
                                                              valley=0;
                                                              }
                                                              if(fitnessstart[j4]<fitnessstart[minstart])
                                                              {
                                                              valley=valley+1;
                                                              length[valley]=1;
                                                              }
                                                              if(fitnessstart[j4]>fitnessstart[minstart])
                                                              {
                                                              minstart=j4;
                                                              valley=0;
                                                              }
                                                              if(fitnessstart[j5]<fitnessstart[minstart])
                                                              {
                                                              valley=valley+1;
                                                              length[valley]=1;
                                                              }
                                                              if(fitnessstart[j5]>fitnessstart[minstart])
                                                              {
                                                              minstart=j5;
                                                              valley=0;
                                                              }
                                                               if(fitnessstart[j6]<fitnessstart[minstart])
                                                              {
                                                              valley=valley+1;
                                                              length[valley]=1;
                                                              }
                                                              if(fitnessstart[j6]>fitnessstart[minstart])
                                                              {
                                                              minstart=j6;
                                                              valley=0;
                                                              }
                                                              maxl=0;
                                                              for(i=0;i<=position;++i)
                                                              {if(length[i]>0){maxl=i;}}
                                                            //  if(fitnessstart[minstart]<0.5){maxl=maxl+1;}
                                                              nakopl[maxl]=nakopl[maxl]+1;
                                                              
                                                               if(maxl<min)
                                                               {min=maxl;}
                                                              suml=suml+1;
                   
    //if(nakopl[maxl]==720) { printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! length %d  min %d\n",maxl,min);}
                                                             // printf("sum %d j5 %d  max %d min %d\n",suml,j5,maxl,min);
                                                              valley=0;
                                                              
                                                              valley=0;
                                                              for(i=0;i<10;++i){length[i]=0;}
                                                              minstart=jj1;    
                                                              }
            
    // 8!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

                      for(i6=1;i6<=prohod[j6];++i6)
                       {
                        j7= arrayvpered[j6][i6];

                       if(distance[j7]==position-1)
                           {

     if(fitnessstart[minstart]<0.5)
       {valley=valley+1;
      length[valley]=1;
      }
      
    
                                                    if(fitnessstart[j2]<fitnessstart[minstart])
                                                       {
                                                        valley=valley+1;
                                                        length[valley]=1;
                                                        }
                                                        if(fitnessstart[j2]>fitnessstart[minstart])
                                                        {
                                                         minstart=j2;
                                                         valley=0;
                                                         }
                                                             if(fitnessstart[j3]<fitnessstart[minstart])
                                                             {
                                                              valley=valley+1;
                                                              length[valley]=1;
                                                              }
                                                              if(fitnessstart[j3]>fitnessstart[minstart])
                                                              {
                                                              minstart=j3;
                                                              valley=0;
                                                              }
                                                              if(fitnessstart[j4]<fitnessstart[minstart])
                                                              {
                                                              valley=valley+1;
                                                              length[valley]=1;
                                                              }
                                                              if(fitnessstart[j4]>fitnessstart[minstart])
                                                              {
                                                              minstart=j4;
                                                              valley=0;
                                                              }
                                                              if(fitnessstart[j5]<fitnessstart[minstart])
                                                              {
                                                              valley=valley+1;
                                                              length[valley]=1;
                                                              }
                                                              if(fitnessstart[j5]>fitnessstart[minstart])
                                                              {
                                                              minstart=j5;
                                                              valley=0;
                                                              }
                                                               if(fitnessstart[j6]<fitnessstart[minstart])
                                                              {
                                                              valley=valley+1;
                                                              length[valley]=1;
                                                              }
                                                              if(fitnessstart[j6]>fitnessstart[minstart])
                                                              {
                                                              minstart=j6;
                                                              valley=0;
                                                              }
                                                               if(fitnessstart[j7]<fitnessstart[minstart])
                                                              {
                                                              valley=valley+1;
                                                              length[valley]=1;
                                                              }
                                                              if(fitnessstart[j7]>fitnessstart[minstart])
                                                              {
                                                              minstart=j7;
                                                              valley=0;
                                                              }
                                                              maxl=0;
                                                              for(i=0;i<=position;++i)
                                                              {if(length[i]>0){maxl=i;}}
                                                            //  if(fitnessstart[minstart]<0.5){maxl=maxl+1;}
                                                              suml=suml+1;
                                                              nakopl[maxl]=nakopl[maxl]+1;
                                                               if(maxl<min)
                                                                 {min=maxl;}
                                                              
                                                              //printf("sum %d j7 %d  max %d\n",suml,j7,maxl);
                                                              valley=0;
                                                              for(i=0;i<10;++i){length[i]=0;}
                                                              minstart=jj1;    
                                                              }
            
     // 9!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

                      for(i7=1;i7<=prohod[j7];++i7)
                       {
                        j8= arrayvpered[j7][i7];

                       if(distance[j8]==position-1)
                           {
                                                        if(fitnessstart[minstart]<0.5)
       {valley=valley+1;
      length[valley]=1;
      }
      
                                                   
                                                    if(fitnessstart[j2]<fitnessstart[minstart])
                                                       {
                                                        valley=valley+1;
                                                        length[valley]=1;
                                                        }
                                                        if(fitnessstart[j2]>fitnessstart[minstart])
                                                        {
                                                         minstart=j2;
                                                         valley=0;
                                                         }
                                                             if(fitnessstart[j3]<fitnessstart[minstart])
                                                             {
                                                              valley=valley+1;
                                                              length[valley]=1;
                                                              }
                                                              if(fitnessstart[j3]>fitnessstart[minstart])
                                                              {
                                                              minstart=j3;
                                                              valley=0;
                                                              }
                                                              if(fitnessstart[j4]<fitnessstart[minstart])
                                                              {
                                                              valley=valley+1;
                                                              length[valley]=1;
                                                              }
                                                              if(fitnessstart[j4]>fitnessstart[minstart])
                                                              {
                                                              minstart=j4;
                                                              valley=0;
                                                              }
                                                              if(fitnessstart[j5]<fitnessstart[minstart])
                                                              {
                                                              valley=valley+1;
                                                              length[valley]=1;
                                                              }
                                                              if(fitnessstart[j5]>fitnessstart[minstart])
                                                              {
                                                              minstart=j5;
                                                              valley=0;
                                                              }
                                                               if(fitnessstart[j6]<fitnessstart[minstart])
                                                              {
                                                              valley=valley+1;
                                                              length[valley]=1;
                                                              }
                                                              if(fitnessstart[j6]>fitnessstart[minstart])
                                                              {
                                                              minstart=j6;
                                                              valley=0;
                                                              }
                                                               if(fitnessstart[j7]<fitnessstart[minstart])
                                                              {
                                                              valley=valley+1;
                                                              length[valley]=1;
                                                              }
                                                              if(fitnessstart[j7]>fitnessstart[minstart])
                                                              {
                                                              minstart=j7;
                                                              valley=0;
                                                              }
                                                               if(fitnessstart[j8]<fitnessstart[minstart])
                                                              {
                                                              valley=valley+1;
                                                              length[valley]=1;
                                                              }
                                                              if(fitnessstart[j8]>fitnessstart[minstart])
                                                              {
                                                              minstart=j8;
                                                              valley=0;
                                                              }
                                                              maxl=0;
                                                              
                                                              for(i=0;i<=position;++i)
                                                              {if(length[i]>0){
                                                             
                                                                            maxl=i;}}
                      //       if(fitnessstart[minstart]<0.5){maxl=maxl+1;}
                                                               nakopl[maxl]=nakopl[maxl]+1;
                                                              suml=suml+1;
                                                               if(maxl<min)
                                                               {min=maxl;}
                                                              //printf("sum %d j7 %d  max %d\n",suml,j7,maxl);
                                                              valley=0;
                                                              for(i=0;i<10;++i){length[i]=0;}
                                                              minstart=jj1;    
                                                              }
            
//8
}
//7
}
//6 
} 
//5 
}        
//4
}        
//3
}        
//2
}                                                 
//1
}
}
}
}
                             


                             
//if(version==1)
{
lmin[cluster]=min;
}

//lmindown[cluster]=mindown;
//printf("min %d lmin %d \n", min,lmin[cluster]);
//printf("mindown %d lmin %d \n", mindown,lmindown[cluster]);
suman=0;
probmin[min]=probmin[min]+1.0;
//printf("min %d prob %f \n", min,probmin[min]);

for(i=0;i<=position;++i)
{
suman=suman+nakopl[i];
if(nakopl[i]>0)
{
//printf("length %d nakop %d suman %d\n",i,nakopl[i],suman);
}
}
lsered[cluster]=0;
for(i=0;i<=position;++i)
{
lsered[cluster]=lsered[cluster]+1.0*i*nakopl[i]/suman;
}
//printf(" lsered %f lserint %d \n", lsered[cluster],lserint);
lserint=100*(lsered[cluster]+0.005);
printf(" minimum %d laveraged %f \n", min,lsered[cluster]);
prob[lserint]=prob[lserint]+1.0;


//cluster
}


serer=0;
seredmin=0;
seredmindown=0;
dysl=0;
dyslmin=0;
for(cl=1;cl<=endcluster;++cl)
{
serer=serer+1.0*lsered[cl];
seredmin=seredmin+1.0*lmin[cl];
}

for(cl=1;cl<=endcluster;++cl)
{
dysl=dysl+(serer/endcluster-1.0*lsered[cl])*(serer/endcluster-1.0*lsered[cl]);
dyslmin=dyslmin+(seredmin/endcluster-1.0*lmin[cl])*(seredmin/endcluster-1.0*lmin[cl]);
}

fprintf(data7,"%d  %f  \n",position,1.0*seredmin/endcluster);




// gama
}


serer=0;



/*
sprintf (name1,"ProblavL%d.dat",position); 
data1 = fopen(name1,"w+");
for(j=0;j<1000;j=j+1)
{
  if(prob[j]>0.0)
  {                   
fprintf(data1, "%f %f  \n", 0.01*j,prob[j]/endcluster);
}
}
fclose(data1);   



sprintf (name1,"ProblminL%d.dat",position);

data1 = fopen(name1,"w+");
for(j=0;j<position;j=j+1)
{
 if(probmin[j]>0.0)
 {                         
fprintf(data1, "%d %f  \n", j,probmin[j]/endcluster);
}
}
fclose(data1);   
*/
// position
}



 fclose(data7);


system ("pause"); 
return 0; 


}