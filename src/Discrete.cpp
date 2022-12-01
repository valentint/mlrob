#include <R.h>
#include <Rmath.h>
#define  MAX 1000
using namespace std;

double mid_point(double n1, double n2);

double entCI(double input[], int cMatrix[], double partition,int nrows, int begin, int end);

double logd(double x);

double ent(double s[], int cs[], int sCount);

double fullent(double s[], int cs[], int begin, int end);

double p(int cls, int cs[], int sCount);

double p(int cls, int cs[], int begin, int end);

void doMaxMinC(int cs[],int sCount, int& maxC, int& minC);

void doMaxMinC(int cs[],int& maxC, int& minC,int begin, int end);

double delta(double input[], int cInput[], double partition, int nrows, int begin, int end);

void do_the_rest(double input[],int cInput[],double t[],double t_sel[],int nrows,int begin,int end,int& index);

int min(double num[], int start, int limit);

double min2(double num[], int start, int limit);


extern "C" {

void discrete(double *mIn, int *cIn, double *t,int *ccsize, double *points, int *nparti)
{

	int index=0,d, csize = *ccsize;
   double t_sel[MAX];

   for(int j=0;j<csize;j++) t_sel[j]=0.0;

   do_the_rest(mIn,cIn,t,t_sel,csize,0,csize,index);

   for(d=0;d<csize;d++) points[d] = t_sel[d];

   *nparti = index + 1;

}



void Points(double *matrix, int *csize, double *tMatrix)
{

	int j;
   for(j=0;j<*csize-1;j++)
   {
      if(matrix[j+1]==matrix[j])
        {
         do
           {
            tMatrix[j]=0.0;
            j++;
           } while(matrix[j+1]==matrix[j]);
         }
         tMatrix[j] = mid_point(matrix[j], matrix[j+1]);
   }
   tMatrix[(*csize)-1]=0.0;

}


}


double mid_point(double n1, double n2) 
{
  return ((n1+n2)/2.0);
}


double entCI(double input[], int cMatrix[], double partition, int nrows,int begin, int end)
{
  double* s1=new double[nrows];
  double* s2=new double[nrows];
  double entropy;
  int* cs1= new int[nrows];
  int* cs2=new int[nrows];
  int s1Count=0, s2Count=0, sCount=0;
    while(input[begin]<partition)
    {
        cs1[s1Count]=cMatrix[begin];
        s1[s1Count++]=input[begin++];
    }
    while(begin<end)
    {
        cs2[s2Count]=cMatrix[begin];
        s2[s2Count++]=input[begin++];
    }
    sCount=s1Count+s2Count;
    entropy=(s1Count/double(sCount))*ent(s1,cs1,s1Count)
               +(s2Count/double(sCount))*ent(s2,cs2,s2Count);
    return entropy;
delete [] s1;
delete [] s2;
delete [] cs1;
delete [] cs2;
}

double ent(double s[], int cs[], int sCount)
{
    double sum=0.0;
    int maxC, minC;
    doMaxMinC(cs, sCount,maxC,minC);
    for(int i=minC;i<=maxC;i++)
    {
       if(p(i,cs,sCount)!=0.0)
       {
         sum+=(p(i,cs,sCount)*logd(p(i,cs,sCount)));
       }
    }
    if(sum==0.0)
      return sum;
    sum=-sum;
    return sum;
}

double fullent(double s[], int cs[], int begin, int end)
{
    double sum=0.0;
    int maxC, minC;
    doMaxMinC(cs,maxC,minC,begin,end);
    for(int i=minC;i<=maxC;i++)
    {
       if(p(i,cs,begin,end)!=0.0)
       {
         sum+=(p(i,cs,begin,end)*logd(p(i,cs,begin,end)));
       }
    }
    if(sum==0.0)
      return sum;
    sum=-sum;
    return sum;
}

void doMaxMinC(int cs[],int sCount, int& maxC, int& minC)
{
    maxC=minC=cs[0];
    for(int i=1;i<sCount;i++)
       if(cs[i]>maxC)
          maxC=cs[i];
       else if(cs[i]<minC)
          minC=cs[i];
}

void doMaxMinC(int cs[],int& maxC, int& minC,int begin, int end)
{
    maxC=minC=cs[begin];
    for(int i=begin+1;i<end;i++)
       if(cs[i]>maxC)
          maxC=cs[i];
       else if(cs[i]<minC)
          minC=cs[i];
}

double p(int cls, int cs[], int sCount)
{
    int clsCount=0;
    for(int i=0;i<sCount;i++)
       if(cs[i]==cls)
          clsCount++;
    if(clsCount==0)
        return 0.0;
    return (clsCount/double(sCount));
}

double p(int cls, int cs[], int begin, int end)
{
	int clsCount=0;
    for(int i=begin;i<end;i++)
       if(cs[i]==cls)
          clsCount++;
    if(clsCount==0)
        return 0.0;
    return (clsCount/double(end-begin));
}

double logd(double x)
{
    return (log(x)/log(2.0));
}

double delta(double input[], int cInput[], double partition, int nrows, int begin, int end)
{
    int nbegin=begin;
    double* s1=new double[nrows];
    double* s2=new double[nrows];
    double delta;
    int* cs1=new int[nrows];
    int* cs2=new int[nrows];
    int s1Count=0, s2Count=0;
    while(input[begin]<partition)
    {
        cs1[s1Count]=cInput[begin];
        s1[s1Count++]=input[begin++];
    }
    while(begin<end)
    {
        cs2[s2Count]=cInput[begin];
        s2[s2Count++]=input[begin++];
    }
    int c=0, c1=0, c2=0, maxC[3], minC[3];
    doMaxMinC(cInput,maxC[0],minC[0],nbegin,end);
    doMaxMinC(cs1,s1Count,maxC[1],minC[1]);
    doMaxMinC(cs2,s2Count,maxC[2],minC[2]);

    for(int i=0;i<3;i++)
    {
        for(int j=minC[i];j<=maxC[i];j++)
        {
            if(i==0 && p(j,cInput,nbegin,end)!=0.0)
               c++;
            else if(i==1 && p(j,cs1,s1Count)!=0.0)
               c1++;
            else if(i==2 && p(j,cs2,s2Count)!=0.0)
               c2++;
        }
    }
    delta = logd(pow(3.0,double(c))-2)-(c*fullent(input,cInput,nbegin,end)-
                     c1*ent(s1,cs1,s1Count)-c2*ent(s2,cs2,s2Count));
    return delta;
  delete [] s1;
  delete [] s2;
  delete [] cs1;
  delete [] cs2;
}

void do_the_rest(double input[],int cInput[],double t[],double t_sel[],int nrows,int begin,int end,int& index)
{
   double* entropy=new double [nrows];
   for(int i=begin;i<end;i++)
   {
        if(t[i]==0.0)
	{
          entropy[i]=0.0;
	}
        else
        {
          entropy[i]=entCI(input,cInput,t[i],nrows,begin,end);
        }
   }
   int idxMinEnt = min(entropy,begin,end);
   if(idxMinEnt>=0)
   {
     double gain = fullent(input,cInput,begin,end)-entropy[idxMinEnt];
     double right_side = (logd((end-begin)-1)/(end-begin))
                   +(delta(input,cInput,t[idxMinEnt],nrows,begin,end)/(end-begin));

	 if(gain>right_side && idxMinEnt>=0)
	 {
	   t_sel[index]=t[idxMinEnt];
	   index++;
	   int nEnd = idxMinEnt+1;
	     do_the_rest(input,cInput,t,t_sel,nrows,begin,nEnd,index);
	     do_the_rest(input,cInput,t,t_sel,nrows,nEnd,end,index);
	 }
   }
   else
   {
       return;
   }
}

int min(double num[],int start, int limit)
{
   double min;
   int i=start;
   while(num[i]==0.0)
   {
      i++;
   }
   min=num[i];
   int result=i;
   for(int j=i;j<limit;j++)
   {
      if(num[j]<min && num[j]!=0.0)
      {
         min=num[j];
         result=j;
      }
   }
   if(result==limit)
	   return -1;
   return (result);
}

double min2(double num[],int start, int limit)
{
   double min;
   min=num[start];
   for(int j=start+1;j<limit;j++)
   {
      if(num[j]<min)
      {
         min=num[j];
      }
   }
   return (min);
}

