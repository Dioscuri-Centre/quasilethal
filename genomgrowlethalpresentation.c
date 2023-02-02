#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
# define endcluster 10000
#define MAXITER 30
# define N 21
# define Nmin 20
# define Nmax 1000
# define increase 50
# define max 3000

//unsigned long long int tt,nopp;
int L,spivpad,zzz,startcluster,build,numberlink,dead,suma,kill,link[101][101],state[Nmax+1],
//story[N+1][tttmax],
jpoper,vvv,tak[100];
 double serer, lethalsituation,alive,ttime,ttt,gamamax,waga[101],
 dystt,dysttt,tser[endcluster+1],tttser[endcluster+1],tserend,tttserend,tt,delta;
 float dystime,doom,RR,nakop,deadfit,mu,gama,mustart,comb,birthfit,mutfit,mf3,mf1,mf2,prbin,factor,timeser,munew,riz,dysp,sum,
 fitness[max+1],fitnessstart[max+1], maxfit, frequencytest[1000];
 
 int letal,jj,t,maxmutation,fr,rozpbin,binttt,count,found, test,clusteryes,maxfound,maxlevy1,maxlevy2,
 maxlevy3,vybirmu,lich,lichf,runing,fint,position,
 numer,maxp,genom[Nmax+1][21],genomtest[Nmax+1],example[max][21],stepx[6],stepy[6],DA, 
 maxtime,spr, chystyj,roz1,roz2,cilyj,xx,yy,ss,r,
 rr,NN, numm, aa,x,level,x1,x2,yy1,y2,porjadokkut,g1,g2,man,kut,levy,timme,serlevy,serlevy2,numser,numlevy,ra;

 float q,ggg,shuffle,perev[max+2],perev2[max+2],perev3[max+2];
 int test1,run,povtor,zrobleno[max+1][21],zmina,bb, bin, kk,betta,testx[Nmax+2],testy[Nmax+2],p1,p2,perc,gi,labelmax,
maxim,u,numerbonn,ll, ffff,cluster, kr,cl,  gy,rapira,add,vuzl, i,nn,neww,hh,enrich;

 
 float sigma,zag1,fin,ymax;
 
 // for random number generators

 long int text;
 int ppp,povt,cc,j,x,g,fff, gm;     

 int trial1,stopp[N+2],trial2,sumalink,conf,zagal,nnn,ii,h,krok,  zag4, zag2,zag3,f,ff;
 
 
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



//# include <C:\Users\WellDone\Documents\Programs\Epidemy\random-generator.rc>
//# include <C:\Users\WellDone\Documents\Programs\Networks-Scale-Free\eigenvalues.rc>

//# include <C:\Users\Виктория\Documents\Programs\Scale-free-networks\random-generator.rc>
//# include <C:\Users\Виктория\Documents\Programs\Scale-free-networks\eigenvalues.rc>
//# include <C:\Users\Виктория\Documents\Programs\Epidemy\random-generator.rc>
 char name9[50];
 char name10[50];
int buildclaster;
buildclaster=build;
FILE  *data9, *data10;
int tri,a,b,vybir,s;      

 srand((unsigned)time(NULL));
 startcluster=1;
  
// kilkist bakterij
for(NN=Nmin;NN<=N;NN=NN+increase)
{
printf("Enter sequence length L: \n" );
  // L = getchar();
    scanf("%d", &zzz); L=zzz; 
printf("Enter mutation probability mu: \n" );                        
  scanf("%f", &ggg); mu=ggg;

printf("Enter probability to mutate into lethal state gamma: \n" );                        
  scanf("%f", &ggg); gama=ggg;

position=L;
numberlink=pow(2,position);
mustart=mu;
sprintf (name9,"T_L=%dgamma%f.dat",position,gama); 
data9 = fopen(name9,"w+");

 //for(gama=0.00;gama<=0.001;gama=gama+0.75)
{          
sprintf (name10,"HistogramL=%dgamma=%1.2f.dat",position,gama); 
data10 = fopen(name10,"w+");

tserend=0;
tttserend=0;

for (cluster=1;cluster<=endcluster;++cluster) 
{
beginning: tser[cluster]=0;
tttser[cluster]=0;


  // creating fitness landscape
for(j=1;j<=position;++j)
{                 
example[1][j]=-1;
}

nn=1;
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
again: for(j=1;j<=numer;++j)
{                 
perev[j]=0;
}

  

for(i=1;i<=numer;++i)
{
 vy: ss=rand()%(numer)+1; 
 if(perev[ss]==1)
 { goto vy;}                    
 else
 {
// fitness is prescribed randomly to each genotype
vvv=rand()%(1000)+1;
fitnessstart[ss]=0.001*vvv;
perev[ss]=1;
}
}

fitnessstart[1]=0.5;
fitnessstart[maxmutation]=1.0;

for(j=0;j<=1000;++j)
{                 
frequencytest[j]=0;
}


// probability of mutation
for(i=1;i<=NN;++i)
{
state[i]=1;
//story[i][0]=1;
// zadajemo "kolektyv" iz N genomiv
for(j=1;j<=position;++j)
{
genom[i][j]=-1;
}
fitness[i]=fitnessstart[1];
fint=100*(fitness[i]+0.005); 
//frequency[fint][tt]=frequency[fint][tt]+1.0;
frequencytest[fint]=frequencytest[fint]+1.0;
}

///////////////////////////////////// kinec zadannja landscape////////////////////////////////////////////////////////////////////    
printf("cluster %d gama %f  \n",cluster,gama);    
    
tt=0.0;
ttt=0;
ttime=0;
zmina=0;
lich=0;
lichf=0;
found=0;
suma=NN;
RR=0;


while(found<1)
{
  letal=0;           
  kill=0;
   rrr();
// probability to die
doom=shuffle;
if(doom<=(1.0*suma/(Nmax)))
{
ss=rand()%(suma)+1;
fint=100*(fitness[ss]+0.005); 
frequencytest[fint]=frequencytest[fint]-1.0;
kill=1;
}

// probability to replicate
vybirth: xx=rand()%(suma)+1;
 rrr();
 birthfit=shuffle;
 if(birthfit<=fitness[xx])
 {
if(kill==0)
{
ss=suma+1;
suma=suma+1;}
     rrr(); 
     mutfit=shuffle;
     comb=(1-mu-gama);
     if(mutfit<=mu)
     {  
   // "mutation"
            i=rand()%(position)+1;
         for(j=1;j<=position;++j)
         {              
         genomtest[j]=genom[xx][j];
         }  
         genomtest[i]=(-1)*genom[xx][i];
        zmina=1;
     }     
     if(mutfit>mu && mutfit<=(mu+comb))
       {
  //  "no mutation";
       for(hh=1;hh<=position;++hh)
         {              
         genomtest[hh]=genom[xx][hh];
         }
        }
        
     if(mutfit>(mu+comb))
     {
     // "mutation into lethal state"                     
      letal=1;
      fitness[ss]=0.0;      
fint=100*(fitness[ss]+0.005);
frequencytest[0]=frequencytest[0]+1.0;                  
     }   
     
     
      for(j=1;j<=position;++j)
           {                     
    //  printf("j %d genom %d \n",j,genomtest[j]);
      } 
    // printf("ttt %f end! \n",ttt); 
     if(letal==0)
     {   
    
     for(j=1;j<=numer;++j)
     {
           spivpad=0;
           for(jj=1;jj<=position;++jj)
           {                     
           if(genomtest[jj]==example[j][jj]){spivpad=spivpad+1;} }
       if(spivpad==position)
     {
    // printf("vidpovidaje %d \n",j); 
      
      state[ss]=j;
      fitness[ss]=fitnessstart[j];
      for (i=1;i<=position;++i)
      {genom[ss][i]=example[j][i];}
      fint=100*(fitness[ss]+0.005);
       frequencytest[fint]=frequencytest[fint]+1.0;

if(j==maxmutation)
{

alive=alive+1.0*(suma-frequencytest[0]);
tttserend=tttserend+ttt;
tttser[cluster]=ttt;
clusteryes=clusteryes+1;
found=1;
printf("adaptation time  %f \n",1.0*(ttt));
//printf("nakop %f\n",tttserend);

//fprintf(data10, "%f \n",1.0*(tt)/(suma-frequencytest[0]));
fprintf(data10, "%f \n",1.0*ttt);



break;
zmina=0;
  }
 }      
}
}
}



  else goto vybirth;  

//fprintf(data6, "%f %d \n",ttime,suma);

ttt=ttt+1.0/(suma-frequencytest[0]);
ttime=ttime+1.0;

// end
}


//fclose(data6);

ending: x=5;
        //cluster
}

fclose(data10);


dysttt=0;
for(cl=1;cl<=endcluster;++cl)
{
dysttt=dysttt+(tttserend/endcluster-tttser[cl])*(tttserend/endcluster-tttser[cl]);
}

fprintf(data9," %f %f %f \n",gama,tttserend/endcluster,sqrt(dysttt/(endcluster*endcluster)));
printf(" gama %f time %f \n",gama,tttserend/endcluster);


// gama
}

fclose(data9);
system ("pause"); 
return 0; 
 
//NN
}
 




}