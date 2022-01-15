#include <stdio.h>
#include <stdlib.h>
#define N 255
#define PL(x,y) (colperm[plane[n1*rowperm[x]+y]])
#define SHOWMOLS 0
// Change to 1 if you want to see MOLS arising from input incidence matrix or list of line sets
int npts,*plane,*mols,nfree,sz,*rowperm,*colperm;
int initbit,nlatinchar,nbits,*perm,*prm,*sq;
char *latinchar;
unsigned n,n1,informat,recmode,allones,buf;
FILE *infile,*outfile;

int main(argc,argv)
int argc;
char *argv[];
{int e,i,ii,j,jj,k,kk;
unsigned m,r;
char outstr[120];

if (argc != 5)
   {fprintf(stderr,
      "Usage: pzip n format recovery_mode filename\n   where format is in {i,l,m}");
   fprintf(stderr," (incidence matrix, list, mols);\n  recovery_mode is in {e,i}");
   fprintf(stderr," (exact, isomorph only).\n\nwww.uwyo.edu/moorhouse/pzip.html\n");
   return 0;}

allones=~0;
nbits=1;
m=2;
while ((m>>nbits) == 1)
   {nbits++;
   m*=2;}
fprintf(stderr,"pzip identifies your %d-bit processor\n",nbits);
if ((nbits%8) != 0)
   {fprintf(stderr,"But %d is not a multiple of 8!?!\n",nbits);
   return 1;}

sscanf(argv[1],"%d",&n);
n1=n+1;
npts=n*n1+1;
plane=malloc(n1*npts*sizeof(int));
mols=malloc((n-1)*n*n*sizeof(int));
perm=malloc(npts*sizeof(int));
prm=malloc(npts*sizeof(int));
rowperm=malloc(npts*sizeof(int));
colperm=malloc(npts*sizeof(int));
latinchar=malloc(n*sizeof(char));
sq=malloc(n*n*sizeof(int));
for (i=0; i<npts; i++) rowperm[i]=colperm[i]=i;  // Omit???
infile=fopen(argv[4],"r");
if (n > N)
   {fprintf(stderr,"Error: order %d exceeds %d\n",n,N);
   return 1;}
switch (argv[3][0])
   {case 'e':
      recmode=1;
      break;
   case 'i':
      recmode=0;
      break;
   default:
      fprintf(stderr,"Unrecognized recovery_mode %s\n",argv[3]);
      return 1;
   }
switch (argv[2][0])
   {case 'i':
      informat=0;
      for (i=0; i<npts; i++)
         {j=0;
	 for (k=0; k<npts; k++)
	    {e=getc(infile)-'0';
	    while (e<0 || e>1)
	       {if (feof(infile) || ferror(infile))
	          {fprintf(stderr,"End of input file %s\n",argv[4]);
	          return 1;}
	       e=getc(infile)-'0';}
	    if (e)
	       {if (j > n)
	          {fprintf(stderr,
		     "Error: exceeding %d lines through point %d\n",n1,i);
	          return 1;}
	       plane[n1*i+(j++)]=k;}
	    }
	 if (j <= n)
	    {fprintf(stderr,"Error: fewer than %d lines through point %d\n",n1,i);
	    return 1;}
	 }
      if ((e=list2latin()) < 0)
         {fprintf(stderr,"Error converting list to mols\n");
	 return 1;}
      break;
   case 'l':
      informat=1;
      k=1000;
      kk=-1;
      for (i=0; i<npts; i++) for (j=0; j<n1; j++)
         {fscanf(infile,"%d",&jj);
	 if (feof(infile) || ferror(infile))
	    {fprintf(stderr,"End of input file %s\n",argv[4]);
	    return 1;}
         plane[n1*i+j]=jj;
	 if (jj < k) k=jj;
	 if (jj > kk) kk=jj;}
      if (k==1 && kk==npts)
         {initbit=1;
	 for (i=0; i<npts; i++) for (j=0; j<n1; j++) plane[n1*i+j]--;}
      else
         {if (k==0 && kk==npts-1) initbit=0;
	 else
	    {fprintf(stderr,
	       "Error: unexpected range [%d,%d] for labels of lines\n",
	       k,kk);
            return 1;}
	 }
      if ((e=list2latin()) < 0)
         {fprintf(stderr,"Error converting list to mols\n");
	 return 1;}
      break;
   case 'm':
      informat=2;
      nlatinchar=0;
      i=j=k=0;
      while (i < n-1)
         {if ((e=readLatinChar()) < 0)
	    {if (e == -1) fprintf(stderr,"End of input file %s\n",argv[4]);
	    else fprintf(stderr,"Exceeding %d distinct chars in input mols\n",n);
	    return 1;}
         mols[n*(n*i+j)+k]=e;
         if (++k >= n)
            {k=0;
   	    if (++j >= n)
   	       {j=0;
   	       i++;}
            }
         }
      break;
   default:
      fprintf(stderr,"Error: Unrecognized input format %s\n",argv[2]);
      return 1;
   }

if (SHOWMOLS) showmols();

if (checkmols())
   {fprintf(stderr,"Error: input failed the test as %d MOLS(%d)\n",n-1,n);
   return 1;}
fprintf(stderr,"Input verified as %d MOLS(%d)\n",n-1,n);

if (recmode==0 && informat==2)
   {for (i=0; i<n-1; i++)  // Standardise first row of mols[i][][]
      {for (k=0; k<n; k++) perm[mols[n*n*i+k]]=k;
      for (j=0; j<n; j++) for (k=0; k<n; k++) mols[n*(n*i+j)+k]=perm[mols[n*(n*i+j)+k]];}
   // Standardise first column of mols[0][][]
   for (j=0; j<n; j++) rowperm[mols[n*j]]=j;
   for (i=0; i<n-1; i++) for (k=0; k<n; k++)
      {for (j=0; j<n; j++) perm[j]=rowperm[mols[n*(n*i+j)+k]];
      for (j=0; j<n; j++) mols[n*(n*i+j)+k]=perm[j];}
   }

nfree=nbits;
buf=0;
sz=0;
sprintf(outstr,"%s.pz",argv[4]);
outfile=fopen(outstr,"wb");

stork(recmode,1);
stork(informat,2);
if (informat == 1) stork(initbit,1);
stork(n,N);
if (recmode && (informat<2))
   {for (k=0; k<npts; k++) perm[k]=1;
   for (i=0; i<npts; i++)
      {e=rowperm[i];
      r=0;
      for (ii=0; ii<e; ii++) r+=perm[ii];
      stork(r,npts-i-1);
      perm[e]=0;}
   for (k=0; k<npts; k++) perm[k]=1;
   for (i=0; i<npts; i++)
      {e=colperm[i];
      r=0;
      for (ii=0; ii<e; ii++) r+=perm[ii];
      stork(r,npts-i-1);
      perm[e]=0;}
   }

for (i=0; i<n-1; i++) for (j=0; j<n; j++) for (k=0; k<n; k++)
   if ((recmode && informat==2) || (j>0 && (i>0 || k>0)))
   {for (e=0; e<n; e++) perm[e]=1;
   for (jj=0; jj<j; jj++) perm[mols[n*(n*i+jj)+k]]=0;
   for (kk=0; kk<k; kk++) perm[mols[n*(n*i+j)+kk]]=0;
   for (jj=0; jj<j; jj++) for (kk=0; kk<n; kk++) if (kk != k)
      {ii=0;
      while (ii<i && mols[n*(n*ii+jj)+kk]!=mols[n*(n*ii+j)+k]) ii++;
      if (ii < i) perm[mols[n*(n*i+jj)+kk]]=0;}
   m=0;
   for (e=0; e<n; e++) m+=perm[e];
   e=mols[n*(n*i+j)+k];
   if (perm[e] == 0)
      {fprintf(stderr,"Error in compression algorithm\n");
      return 1;}
   r=0;
   for (ii=0; ii<e; ii++) r+=perm[ii];
   stork(r,m-1);}
 
if (nfree < nbits)
   {fwrite(&buf,nbits/8,1,outfile);
   sz+=nbits/8;}
if (recmode && informat==2)
   {for (k=0; k<n; k++) putc(latinchar[k],outfile);
   sz+=n;}
fclose(outfile);
fprintf(stderr,"Done; %d bytes written to %s\n",sz,outstr);}

int stork(unsigned k,unsigned m)  // Append the binary expression for k (from {0,1,2,...,m}) to the output
{int sm,tail;
if (m == 0) return 0;
sm=szreqd(m);
if ((tail=sm-nfree) >= 0)
   {buf|=(k>>tail);
   fwrite(&buf,nbits/8,1,outfile);
   sz+=nbits/8;
   nfree=nbits-tail;
   if (tail > 0) buf=k<<nfree;
   else buf=0;}
else
   {nfree=-tail;
   buf|=(k<<nfree);}
return 0;}

int szreqd(unsigned m)  // Number of bits required to store an integer in {0,1,2,...,m}
{int i;
i=nbits-1;
while (i>=0 && ((m>>i)&1)==0) i--;
return i+1;}

int swap_row(int p1,int p2)
{int t;
t=rowperm[p1];
rowperm[p1]=rowperm[p2];
rowperm[p2]=t;
return 1;}

int swap_col(int l1,int l2,int pos)
{int i,u1,u2;
if (l1<0 || l1>=npts || l2<0 || l2>=npts || l1==l2)
   {fprintf(stderr,"Error: cannot swap columns %d and %d\n",l1,l2);
   return 1;}

u1=u2=1;
i=0;
while (i<npts && u1 && u2)
   {if (colperm[i] == l1)
      {colperm[i]=l2;
      u2=0;}
   else if (colperm[i] == l2)
      {colperm[i]=l1;
      u1=0;}
   i++;}
while (i<npts && u2)
   {if (colperm[i] == l1)
      {colperm[i]=l2;
      return 0;}
   i++;}
while (i<npts && u1)
   {if (colperm[i] == l2)
      {colperm[i]=l1;
      return 0;}
   i++;}
if (i>=npts || u1 || u2) fprintf(stderr,"Error in position %d\n",pos);
return (i>=npts || u1 || u2);}

int inc(int p,int l)
{int i;
if (p<0 || p>=npts || l<0 || l>=npts)
   {fprintf(stderr,"Error: inc(%d,%d) undefined\n",p,l);
   return -1;}
for (i=0; i<n1; i++) if (PL(p,i) == l) return 1;
return 0;}

int readLatinChar()
{int i;
char c;
c=getc(infile);
if (feof(infile) || ferror(infile)) return -1;
while (c==' ' || c=='\n' || c=='\t' || c=='\r')
   {c=getc(infile);
   if (feof(infile) || ferror(infile)) return -1;}
i=0;
while (i<nlatinchar && latinchar[i]!=c) if (++i >= n) return -2;
if (i == nlatinchar) latinchar[nlatinchar++]=c;
return i;}

int list2latin()
{int i,j,k,ii,jj,kk;
// Permute rows and columns into standard form.
//
// Row 0
for (i=0; i<npts; i++) rowperm[i]=colperm[i]=i;
for (k=0; k<n1; k++) if ((j=PL(0,k)) > k) swap_col(k,j,1);

// Col 0
i=1;
while (inc(i,0)) i++;
ii=npts-1;
while (inc(ii,0) == 0) ii--;
while (i < ii)
   {swap_row(i,ii);
   while (inc(i,0)) i++;
   while (inc(ii,0) == 0) ii--;}
if (i!=n1 || ii!=n) printf("Oops: i=%d, ii=%d\n",i,ii);

// Rows 1,2,...,n
for (k=0; k<npts; k++) perm[k]=0;
for (i=1; i<n1; i++) for (k=0; k<n1; k++) if ((j=PL(i,k)) > 0) perm[j]=i;
i=n1;
while (i<npts && (ii=perm[i])==(jj=((i-1)/n))) i++;
while (i < npts)
   {j=i+1;
   while (j<npts && perm[j]!=jj) j++;
   if (j >= npts) return -1;
   perm[i]=jj;
   perm[j]=ii;
   swap_col(i,j,2);
   // i++;
   while (i<npts && (ii=perm[i])==(jj=((i-1)/n))) i++;}

// Cols 1,2,...,n
for (k=0; k<npts; k++) perm[k]=0;
for (j=n1; j<npts; j++)
   {k=0;
   while ((i=PL(j,k)) > n) k++;
   perm[j]=i;}
i=n1;
while (i<npts && (ii=perm[i])==(jj=((i-1)/n))) i++;
while (i < npts)
   {j=i+1;
   while (j<npts && perm[j]!=jj) j++;
   if (j >= npts) return -1;
   perm[i]=jj;
   perm[j]=ii;
   swap_row(i,j);
   // i++;
   while (i<npts && (ii=perm[i])==(jj=((i-1)/n))) i++;}

// Rows n+1,n+2,...,2n
for (i=n1; i<=2*n; i++) for (k=0; k<n1; k++) if ((j=PL(i,k)) > n) perm[j]=i;
i=n1;
while (i<npts && (ii=perm[i])==(jj=(n1+((n+i-1)%n)))) i++;
while (i < npts)
   {j=i+1;
   while (j<npts && perm[j]!=jj) j++;
   if (j >= npts) return -1;
   perm[i]=jj;
   perm[j]=ii;
   swap_col(i,j,3);
   i++;
   while (i<npts && (ii=perm[i])==(jj=(n1+((n+i-1)%n)))) i++;}

// Cols n+1,n+2,...,2n
for (j=n1; j<npts; j++)
   {k=0;
   while (k<n1 && ((i=PL(j,k))<n1 || i>2*n)) k++;
   perm[j]=i;}
i=n1;
while (i<npts && (ii=perm[i])==(jj=(n1+((n+i-1)%n)))) i++;
while (i < npts)
   {j=i+1;
   while (j<npts && perm[j]!=jj) j++;
   if (j >= npts) return -1;
   perm[i]=jj;
   perm[j]=ii;
   swap_row(i,j);
   i++;
   while (i<npts && (ii=perm[i])==(jj=(n1+((n+i-1)%n)))) i++;}

// Stage 7 normalisation
for (i=0; i<n; i++)
   {k=0;
   while (((PL(n1+n+i,k)+n-1)%n) > 0) k++;
   perm[(PL(n1+n+i,k)-n-1)/n]=i;}
for (j=1; j<n; j++) if ((jj=perm[j]) > j)
   {k=1;
   while (k<n && perm[k]!=j) k++;
   if (k == n) return -1;
   perm[j]=j;
   perm[k]=jj;
   for (kk=0; kk<n; kk++) swap_col(n1+n*j+kk,n1+n*k+kk,4);}
for (i=1; i<n1; i++)
   {k=0;
   while ((j=PL(i,k)) == 0) k++;
   perm[i]=(j-1)/n;}
i=1;
while (i <= n)
   {if ((ii=perm[i]) == i) i++;
   else
      {perm[i]=perm[ii];
      perm[ii]=ii;
      swap_row(i,ii);}
   }

// Create the n-1 MOLS(n)
for (i=0; i<n*(n-1); i++)
   {ii=i/n;
   for (k=0; k<n1; k++) if ((j=PL(n1+n+i,k)-n-1) >= 0)
      {jj=j/n;
      mols[n*(n*ii+jj)+(j%n)]=i%n;}
   }
return 0;}

int showlist(int pos)
{int i,j;
printf("\nPosition %d\n",pos);
for (i=0; i<npts; i++)
   {for (j=0; j<npts; j++) perm[j]=0;
   for (j=0; j<n1; j++) perm[PL(i,j)]=1;
   for (j=0; j<npts; j++) if (perm[j]) printf(" %d",j);
   printf("\n");}
return 0;}

int checkmols()
{int i,j,k,ii,jj,kk,e;
for (i=0; i<n-1; i++)  // Check that mols[i][][] is a Latin square
   {for (j=0; j<n; j++)
      {for (k=0; k<n; k++) perm[k]=0;
      for (k=0; k<n; k++)
         {if ((e=mols[n*(n*i+j)+k])<0 || e>=n)
   	    {fprintf(stderr,
	       "Error: square %d has (%d,%d)-entry equal to %d outside of the range 0..%d\n",
   	       i,j,k,e,n-1);
   	    return 1;}
   	 if (perm[e])
   	    {fprintf(stderr,"Error: Entry %d occurs twice in row %d of square %d\n",e,j,i);
   	    return 1;}
         perm[e]=1;}
      }
   for (k=0; k<n; k++)
      {for (j=0; j<n; j++) perm[j]=0;
      for (j=0; j<n; j++)
         {e=mols[n*(n*i+j)+k];
   	 if (perm[e])
   	    {fprintf(stderr,"Error: Entry %d occurs twice in col %d of square %d\n",e,k,i);
   	    return 1;}
         perm[e]=1;}
      }
   }
for (i=0; i<n-2; i++) for (ii=i+1; ii<n-1; ii++)  // Check that squares i and ii are orthogonal
   {for (j=0; j<n; j++) for (k=0; k<n; k++) sq[n*j+k]=0;
   for (j=0; j<n; j++) for (k=0; k<n; k++)
      {jj=mols[n*(n*i+j)+k];
      kk=mols[n*(n*ii+j)+k];
      if (sq[n*jj+kk])
         {fprintf(stderr,
   	    "Error: Squares %d and %d are not orthogonal; entry pair (%d,%d) occurs more than once.\n",
   	    i,ii,jj,kk);
	 fprintf(stderr,"The second such occurrence is in position (%d,%d)\n",j,k);
         return 1;}
      sq[n*jj+kk]=1;}
   }
return 0;}

int showmols()
{int i,j,k;
for (i=0; i<10 && i<n; i++) latinchar[i]='0'+i;
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
for (i=0; i<n-1; i++)
   {if (i > 0) putchar('\n');
   for (j=0; j<n; j++)
      {for (k=0; k<n; k++) putchar(latinchar[mols[n*(n*i+j)+k]]);
      putchar('\n');}
   }
fflush(stdout);
return 0;}
