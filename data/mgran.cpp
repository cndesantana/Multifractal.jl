/* Programa para cálculo do espectro
multifractal de dados porosimetria utilizando
o método de Chhabra e Jansen	*/

#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX(mx,x) (x>mx?x:mx)
#define MIN(mx,x) (x<mx?x:mx)

struct Hstc {
    double sl,sd,r,in,ea,eb;
};

void help();
struct Hstc Fitting(double *, double *, int );
double *dvector(long);
long *lvector(long);
double **dmatrix(int,int);
void free_dmatrix(double **m,int nrh,int nch);
void free_dvector(double *);
void free_lvector(long *);
void nrerror(char *);
void Chext(char *, char *);
double calc_SumM(double *x,double *y,double Ei,double Ef,int N);
double E;

int main(int argc, char **argv)
{
    FILE *Fin,*Fout,*FoutFa;
    char NFout[200],NFoutFa[200];
    long N,i;
    double lx,ly,Qi,Qf,dq,Dq,q,EDq,mq;
    double MaxY=-999,MinY=999,SomaY=0;
    double MaxX=-999,MinX=999,SomaX=0;
    struct Hstc FDq, FAq,FFq;
    double *x;
    double *y;
    double *e,**Md,*Ma,*Mf;
    int Np,k,I;
    double RmDq,RmFa;					// R2 minimo a ser considerado
    char Type;					// Define o tipo de saida S para spectro e R p/resumo

    double Dqmx=-999.9,EDqmx,RDqmx;
    double Dqmn= 999.9,EDqmn,RDqmn;
    double qMin,qMax;
    double Do,RDo,EDo;
    double D1,RD1,ED1;
    double D2,RD2,ED2;
    double ao,EAo,RAo;
    double Amn=999,EAmn,RAmn;		// Alfa minimo, erro e r2
    double Amx=-999,EAmx,RAmx;	// Alfa maximo, erro e r2
    double Aqmx=-999, Aqmn=999;	// q para o alfa maximo e mínimo.
    double Fmx=-999, Fmn=999;   // f(alfa) min e maximo
    int Io=0;					// Initial partition


    D2=D1=RD1=RD2=ED1=ED2=-1;	// -1 indicates that for the especific q (2 or 1) the R was not calculated

    /******************* leitura dos arquivos *****************/
    if(argc==2)
    {
        printf("Data\t'+q\t'-q\tDmin\tEDmin\tRDmin\tDmax\tEDmax\tRDmax\tDo\tEDo\tRDo\tD1\tED1\tRD1\tD2\tED2\tRD2\t'qAMin\t'qAMax\tAo\tEAo\tRAo\tAmax\tEAmax\tRAmax\tAmin\tEAmin\tRAmin\tFamin\tFamax\n");
        exit(0);
    }

    if(argc!=12)
    {
        help();
        exit(0);
    }

    strcpy(NFout,argv[1]);
    Chext(NFout,argv[2]);

    strcpy(NFoutFa,argv[1]);
    Chext(NFoutFa,argv[3]);

    Fin=fopen(argv[1],"r");
    Fout=fopen(NFout,"w+");
    FoutFa=fopen(NFoutFa,"w+");
    Qi=atof(argv[4]);
    Qf=atof(argv[5]);
    dq=atof(argv[6]);
    Np=atoi(argv[7]);
    RmDq=atof(argv[8]);
    RmFa=atof(argv[9]);
    Type=argv[10][0];
    Io=atoi(argv[11]);		// Initial partition

    /* Fix the size of the file, the maximum and minimum */
    N=0;
    do
    {
        N++;
        fscanf (Fin, "%lf %lf\n", &lx, &ly);
    }
    while (feof(Fin)==0);

    x=dvector(N);
    y=dvector(N);
    Md=dmatrix((int)((Qf-Qi)/dq)+1,Np+1);
    Ma=dvector(Np+1);
    Mf=dvector(Np+1);

    e=dvector(Np+1);
    /* Load data  */
    SomaY=0;
    SomaX=0;
    fseek(Fin,(long)0,SEEK_SET);
    for(i=0;i<N;i++){
        fscanf (Fin, "%lf %lf\n", &x[i],&y[i]);
        MaxY=MAX(MaxY,y[i]);			// Calculates the maximum and minimum
        MinY=MIN(MinY,y[i]);
        MaxX=MAX(MaxX,x[i]);
        MinX=MIN(MinX,x[i]);
        SomaY+=y[i];
    }
    for(i=0;i<N;i++)
        x[i]=(x[i]-MinX)/(MaxX-MinX);			// Escala fica entre zero e um

    // Begins the "thing"
    I=Io;			// Initial partition, for I=1 the mi(Epson) finalize with Epson=1/2
    for(q=Qi;q<=Qf;q+=dq)
    {
        //printf("%g \n",q);
        for(k=0;k<Np;k++)Ma[k]=Mf[k]=Md[(int)((q-Qi)/dq)][k]=0;	// inicialize fitting vectors
        for(k=I;k<Np;k++)						// Loop for partition numbers
        {
            //printf("%d \n",k);
            double Nor=0,m;
            int Pr;
            Pr=pow(2,k);
            E=(double)1.0/Pr;						// Size of each partition
            e[k-I]=log10(1.0/Pr);

            for(i=1;i<=Pr;i++)						// To estimate f(alfa)
            {
                m=calc_SumM(x,y,(i-1)*E,i*E,N)/SomaY;
                if(m)Nor+=pow(m,q);
            }

            for(i=1;i<=Pr;i++)			// loop for scan over the partition
            {
                m=calc_SumM(x,y,(i-1)*E,i*E,N)/SomaY;
                if(m==0) continue;		// Evita divergencias de medidas nulas
                if(q>(1-dq/2)&&q<(1+dq/2))
                    Md[(int)((q-Qi)/dq)][k-I]+= (m*log10(m)/Nor);	// if q==1
                else	 Md[(int)((q-Qi)/dq)][k-I]+=pow(m,q);

                mq=pow(m,q);					// To estimate f(alfa)
                mq/=Nor;
                Ma[k-I]+=mq*log10(m);
                Mf[k-I]+=mq*log10(mq);
            }
            if(q==0)
                int y=0;


            if(! (q>(1-dq/2) && q<(1+dq/2)) )
                Md[(int)((q-Qi)/dq)][k-I]= log10(Md[(int)((q-Qi)/dq)][k-I]); // if q!=1
        }
        if(q==0)
            int y=0;
        FAq=Fitting(e,Ma,Np);
        FFq=Fitting(e,Mf,Np);
        FDq=Fitting(e,Md[(int)((q-Qi)/dq)],Np);
        if((q>(1-dq/2) && q<(1+dq/2)))Dq=FDq.sl;
        else {Dq=FDq.sl/(q-1);FDq.sd/=fabs(q-1);}
        if(FAq.r>=RmFa && FFq.r>=RmFa)
        {
            if(Type=='S')
            {
                fprintf(FoutFa,"%lf %lf %le %lf %lf %le\n",FAq.sl,FAq.sd,FAq.r,FFq.sl,FFq.sd,FFq.r);
            }
            else
            {
                if(FAq.sl>Amx) {Amx=FAq.sl;EAmx=FAq.sd;RAmx=FAq.r;Aqmx=q;}
                if(FAq.sl<Amn) {Amn=FAq.sl;EAmn=FAq.sd;RAmn=FAq.r;Aqmn=q;}
                if(FFq.sl<Fmn) {Fmn=FFq.sl;}
                if(FFq.sl>Fmx) {Fmx=FFq.sl;}
                if((q>(0-dq/2) && q<(0+dq/2))){ao=FAq.sl;EAo=FAq.sd;RAo=FAq.r;}
            }
        }
        if(FDq.r>=RmDq)
        {
            if(Type=='S')
            {
                fprintf(Fout,"%f %lf %lf %lf %le\n",q,Dq,Dq*(q-1),FDq.sd,FDq.r);
            }
            else
            {
                EDq=((q>(1-dq/2) && q<(1+dq/2))?FDq.ea:fabs(FDq.ea/(q-1)));

                if(Dq>Dqmx){Dqmx=Dq;qMax=q;EDqmx=EDq;RDqmx=FDq.r;}
                if(Dq<Dqmn){Dqmn=Dq;qMin=q;EDqmn=EDq;RDqmn=FDq.r;}
                if((q>(0-dq/2) && q<(0+dq/2))){Do=Dq;RDo=FDq.r;EDo=EDq;}
                if((q>(1-dq/2) && q<(1+dq/2))){D1=Dq;RD1=FDq.r;ED1=EDq;}
                if((q>(2-dq/2) && q<(2+dq/2))){
                    D2=Dq;RD2=FDq.r;ED2=EDq;}
            }
        }
    }
    if(Type=='R'){
        fprintf(Fout,"Scale\t");								// Gera Medidas Tau(q) para conferencia
        for(q=Qi;q<=Qf;q+=dq)fprintf(Fout,"\'%f\t",q);
        fprintf(Fout,"\n");
        for(k=0;k<Np;k++)
        {
            fprintf(Fout,"%lf\t",e[k]);
            for(q=Qi;q<=Qf;q+=dq)
                fprintf(Fout,"%lf\t",Md[(int)((q-Qi)/dq)][k]);
            fprintf(Fout,"\n");
        }
        fprintf(Fout,"\n");
        printf("%s\t%3.1f\t%3.1f\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",argv[1],qMin,qMax,Dqmn,EDqmn,RDqmn,Dqmx,EDqmx,RDqmx);
        printf("%lf\t%lf\t%lf\t",Do,EDo,RDo);
        printf("%lf\t%lf\t%lf\t",D1,ED1,RD1);
        printf("%lf\t%lf\t%lf\t",D2,ED2,RD2);
        printf("%3.1f\t%3.1f\t",Aqmn,Aqmx);
        printf("%lf\t%lf\t%lf\t",ao,EAo,RAo);
        printf("%lf\t%lf\t%lf\t",Amx,EAmx,RAmx);
        printf("%lf\t%lf\t%lf\t",Amn,EAmn,RAmn);
        printf("%lf\t%lf\n",Fmn,Fmx);

    }

    free_dvector(x);
    free_dvector(y);
    free_dvector(e);
    free_dvector(Mf);
    free_dvector(Ma);
    free_dmatrix(Md,(int)((Qf-Qi)/dq)+1,Np+1);
    fclose(Fin);
    fclose(Fout);
    fclose(FoutFa);
}


// Calculates mesure into partition
double calc_SumM(double *x,double *y,double Ei,double Ef,int N)
{
    int i;
    double ret=0;
    for(i=0;i<N;i++)
        if(x[i]>Ei && x[i]<=Ef)
            ret+=y[i];
    return ret;
}


double *dvector(long nh)
{
    double *v;
    unsigned long tam;
    tam=(nh+1)*sizeof(double);

    v=(double *)malloc((unsigned long) tam);
    if (!v) nrerror("allocation failure in dvector()");
    return v;
}
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
double **dmatrix(int nrh,int nch)
{
    int i;
    double **m;

    /* allocate pointers to rows */
    m=(double **) malloc((unsigned) (nrh+1)*sizeof(double*));
    if (!m) nrerror("allocation failure 1 in dmatrix()");

    /* allocate rows and set pointers to them */
    for(i=0;i<=nrh;i++) {
        m[i]=(double *) malloc((unsigned) (nch+1)*sizeof(double));
        if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
    }
    /* return pointer to array of pointers to rows */
    return m;
}


/* free a double matrix allocated by dmatrix() */
void free_dmatrix(double **m,int nrh,int nch)
{
    int i;

    for(i=nrh;i>=0;i--) free((char*) (m[i]));
    free((char*) (m));
}

long *lvector(long nh)
{
    long *v;
    unsigned long tam;
    tam=(nh+1)*sizeof(long);

    v=(long *)malloc((unsigned long) tam);
    if (!v) nrerror("allocation failure in dvector()");
    return v;
}


void nrerror(char error_text[])
{
    fprintf(stderr,"Numerical Recipes run-time error...\n");
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"...now exiting to system...\n");
    exit(0);
}


void help()
{
    fprintf(stderr,"usage: MultiPerfisChhabra {input file} {output extes. Dq} {output extes. Fa} {initial q} {final q} {dq} {Np}[R2 min dq][R2 min fa][R/S] \n");
    fprintf(stderr,"       where R/S means: S for f(alfa) spectra and R for resum\n");

}


void free_dvector(double *v)
{
    free(v);
}
void free_lvector(long *v)
{
    free(v);
}

struct Hstc Fitting(double *vx, double *vy, int N)
{
    struct Hstc Ht;
    double s,sx=0.0,sy=0.0,sx2=0.0,sxy=0.0,sy2=0.0,a,b,r,rx,ry,w,sa,sb;
    int i;

    /*	CALCULATE THE SUM */
    for(i=0;i<N;i++)
    {
        sx+=vx[i];sy+=vy[i];
        sxy+=vx[i]*vy[i];
        sx2+=vx[i]*vx[i]; sy2+=vy[i]*vy[i];
    }
    /*  	CALCULATE THE SLOPE */

    /*  DETERMINE THE COEFICIENTS	*/
    s=sx2-sx*sx/N;
    a=(sxy-sx*sy/N)/s;
    b=(sy-a*sx)/N;
    w=sy2+a*a*sx2+N*b*b;
    w=w-2.0*a*sxy-2.0*b*sy+2.0*a*b*sx;
    if (w<0.0) w=0.0;
    else w=sqrt(w/(N-2));
    rx=sx2-sx*sx/N;
    ry=sy2-sy*sy/N;

    // Slope error
    sa=(sy2+N*b*b+a*a*sx2-2*(b*sy-a*b*sx+a*sxy))/(N-2);
    sb=sqrt( (sx2*sa)/(N*sx2-sx*sx) );
    sa=sqrt( (N*sa)/(N*sx2-sx*sx) );

    if(fabs(ry)<1.0e-10)
    {
        if(fabs(a)<1.0e-10) r=1.0;
        else r=30000.0;
    }
    else r=a*a*rx/ry;

    Ht.sl=a;
    Ht.sd=w;
    Ht.r=r;
    Ht.in=b;
    Ht.ea=sa;
    Ht.eb=sb;

    return Ht;
}

void Chext(char *name, char *ext)
{
    int idx;
    for(idx=strlen(name);idx>=0 && name[idx]!='.';idx--);
    if(idx>=0)
        name[idx+1]=0;
    strcat(name,ext);
}

