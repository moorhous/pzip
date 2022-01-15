#include <stdio.h>
#include <stdlib.h>
#define N 255
#define PL(x,y) (colperm[plane[n1*rowperm[x]+y]])
int npts,*plane,*mols,nfree,sz,*rowperm,*colperm;
int initbit,nlatinchar,nbits,*perm;
char *latinchar;
unsigned n,n1,informat,recmode,allones,buf;
FILE *infile;

unsigned readk(unsigned m)  // Read the binary expression for k (from {0,1,2,...,m})
{int sm,tail;
unsigned t,newbuf;
sm=szreqd(m);
if ((tail=sm-nfree) > 0)
   {fread(&newbuf,nbits/8,1,infile);
   t=(buf<<tail)|(newbuf>>(nbits-tail));
   buf=newbuf&(allones>>tail);
   nfree=nbits-tail;
   return t;}
else
   {nfree=-tail;
   t=buf>>nfree;
   if (nfree == 0)
      {fread(&buf,nbits/8,1,infile);
      nfree=nbits;}
   else buf&=(allones>>(nbits-nfree));
   return t;}
}

int main(argc,argv)
int argc;
char *argv[];
{int e,i,ii,j,jj,k,kk;
unsigned m,r;
char instr[120];

if (argc != 2)
   {fprintf(stderr,
      "Usage: punzip filename\nwww.uwyo.edu/moorhouse/pzip.html\n");
   return 0;}

allones=~0;
nbits=1;
m=2;
while ((m>>nbits) == 1)
   {nbits++;
   m*=2;}
fprintf(stderr,"punzip identifies your %d-bit processor\n",nbits);
if ((nbits%8) != 0)
   {fprintf(stderr,"But %d is not a multiple of 8!?!\n",nbits);
   return 1;}
sprintf(instr,"%s.pz",argv[1]);
if ((infile=fopen(instr,"rb")) == NULL)
   {sprintf(instr,"%s",argv[1]);
   if ((infile=fopen(instr,"rb")) == NULL)
      {fprintf(stderr,"Cannot find file %s.pz or %s\n",instr,instr);
      return 1;}
   }
nfree=0;
recmode=readk(1);
informat=readk(2);
if (informat == 1) initbit=readk(1);
n=readk(N);
n1=n+1;
npts=n*n1+1;
plane=malloc(n1*npts*sizeof(int));
mols=malloc((n-1)*n*n*sizeof(int));
perm=malloc(npts*sizeof(int));
rowperm=malloc(npts*sizeof(int));
colperm=malloc(npts*sizeof(int));
latinchar=malloc(n*sizeof(char));

if (informat < 2)
   {if (recmode)
      {for (k=0; k<npts; k++) perm[k]=1;
      for (i=0; i<npts; i++)
         {if (i == npts-1)
            {e=0;
   	    while (e<npts && perm[e]==0) e++;}
         else
            {r=readk(npts-i-1);
            ii=0;
            for (e=0; e<npts && ii<=r; e++) ii+=perm[e];
   	    if ((--e)>=npts || ii<r)
   	       {fprintf(stderr,
	          "^GError in decompression algorithm: e=%d, ii=%d, r=%d\n",e,ii,r);
   	       return 1;}
            perm[e]=0;}
         rowperm[e]=i;}
      for (k=0; k<npts; k++) perm[k]=1;
      for (i=0; i<npts; i++)
         {if (i == npts-1)
            {e=0;
   	    while (e<npts && perm[e]==0) e++;}
         else
            {r=readk(npts-i-1);
            ii=0;
            for (e=0; e<npts && ii<=r; e++) ii+=perm[e];
   	    if ((--e)>=npts || ii<r)
   	       {fprintf(stderr,
	          "^GError in decompression algorithm: e=%d, ii=%d, r=%d\n",e,ii,r);
   	       return 1;}
            perm[e]=0;}
         colperm[e]=i;}
      }
   else for (k=0; k<npts; k++) rowperm[k]=colperm[k]=k;}

for (i=0; i<n-1; i++) for (j=0; j<n; j++) for (k=0; k<n; k++)
   {if ((recmode && informat==2) || (j>0 && (i>0 || k>0)))
      {for (e=0; e<n; e++) perm[e]=1;
      for (jj=0; jj<j; jj++) perm[mols[n*(n*i+jj)+k]]=0;
      for (kk=0; kk<k; kk++) perm[mols[n*(n*i+j)+kk]]=0;
      for (jj=0; jj<j; jj++) for (kk=0; kk<n; kk++) if (kk != k)
         {ii=0;
         while (ii<i && mols[n*(n*ii+jj)+kk]!=mols[n*(n*ii+j)+k]) ii++;
         if (ii < i) perm[mols[n*(n*i+jj)+kk]]=0;}
      m=0;
      for (e=0; e<n; e++) m+=perm[e];
      if (m < 1)
         {fprintf(stderr,"Error in decompression algorithm: m=%d\n",m);
         return 1;}
      if (m == 1)
         {e=0;
         while (perm[e] == 0) e++;}
      else
         {r=readk(m-1);
         ii=0;
         for (e=0; e<n && ii<=r; e++) ii+=perm[e];
         if ((--e)>=n || ii<r)
            {fprintf(stderr,
	       "Error in decompression algorithm: e=%d, ii=%d, r=%d\n",e,ii,r);
            return 1;}
         }
      }
   else 
      {if (j == 0) e=k;
      else /* i==0 && k==0 */ e=j;}
   mols[n*(n*i+j)+k]=e;}
if (informat == 2)
   {if (recmode) for (k=0; k<n; k++) latinchar[k]=getc(infile);
   else
      {for (i=0; i<10 && i<n; i++) latinchar[i]='0'+i;
      if (n > 10)
         {for (i=0; i<26 && 10+i<n; i++) latinchar[10+i]='a'+i;
         if (n > 36)
            {for (i=0; i<26 && 36+i<n; i++) latinchar[36+i]='A'+i;
	    if (n > 62)
               {for (i=0; i<15 && 62+i<n; i++) latinchar[62+i]='!'+i;
	       if (n > 77)
                  {for (i=0; i<6 && 77+i<n; i++) latinchar[77+i]='['+i;
	          if (n > 83)
                     {for (i=0; i<4 && 83+i<n; i++) latinchar[83+i]='{'+i;
	             if (n > 87)
	                {fprintf(stderr,
	                   "Only 87 characters currently implemented for Latin squares\n");
	                return 0;}
		     }
	          }
	       }
	    }
	 }
      }
   }
fclose(infile);

if (informat < 2)
   {j=0;
   for (i=0; i<n1; i++)
      {plane[n1*i]=0;
      for (k=1; k<n1; k++) plane[n1*i+k]=++j;}
   for (i=0; i<n; i++)
      {plane[n1*(n1+i)]=1;
      for (k=0; k<n; k++) plane[n1*(n1+i)+k+1]=n1+i+n*k;}
   for (i=0; i<n-1; i++) for (ii=0; ii<n; ii++)
      {plane[n1*(n*(i+1)+n1+ii)]=i+2;
      for (j=0; j<n; j++) plane[n1*(n*(i+1)+n1+mols[n*(n*i+j)+ii])+j+1]=n1+n*j+ii;}
   }
switch (informat)
   {case 0:
      for (i=0; i<npts; i++)
         {for (j=0; j<npts; j++) perm[j]=0;
	 for (k=0; k<n1; k++) perm[PL(i,k)]=1;
         for (j=0; j<npts; j++) printf("%d",perm[j]);
	 putchar('\n');}
      break;
   case 1:
      for (i=0; i<npts; i++)
         {for (j=0; j<npts; j++) perm[j]=0;
	 for (k=0; k<n1; k++) perm[PL(i,k)]=1;
	 k=0;
         for (j=0; j<npts; j++) if (perm[j])
	    {if (++k > 0) putchar(' ');
	    printf("%d",j+initbit);}
	 putchar('\n');}
      break;
   case 2:
      for (i=0; i<n-1; i++)
         {for (j=0; j<n; j++)
            {for (k=0; k<n; k++) putchar(latinchar[mols[n*(n*i+j)+k]]);
            putchar('\n');}
         if (i < n-2) putchar('\n');}
      break;
   default:
      fprintf(stderr,"Unrecognized format\n");
      break;
   }
}

int szreqd(unsigned m)  // Number of bits required to store an integer in {0,1,2,...,m}
{int i;
i=nbits-1;
while (i>=0 && ((m>>i)&1)==0) i--;
return i+1;}

int showlist(int pos)
{int i,j;
printf("\nPosition %d\n",pos);
for (i=0; i<npts; i++)
   {for (j=0; j<npts; j++) perm[j]=0;
   for (j=0; j<n1; j++) perm[PL(i,j)]=1;
   for (j=0; j<npts; j++) if (perm[j]) printf(" %d",j);
   printf("\n");}
return 0;}
