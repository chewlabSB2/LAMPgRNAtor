#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <unistd.h>
#include <time.h>
#include <sys/stat.h>
#include <ctype.h>
#include <common.h>

/// generate a read; int length: the length of reads
void generate(char *seq,char out[],int pos,int length)
{
    int i;
    for(i=0;i<length;i++)
    {
        out[i]=seq[pos+i];
    }
    out[i]='\0';
}

///check the GC-content; int length: the length of read
float gc(char seq[],int length)
{
    int i,number;
    float gc;

    number=0;
    for(i=0;i<length;i++)
    {
        if(seq[i]=='C')
        {
            number++;
            continue;
        }
    
        if(seq[i]=='G')
        {
            number++;
        }
    }

    gc=1.0*number/length*100;
    return gc;
}

///translate A...G to int
int translate(char a)
{
    if(a=='A')
        return 0;
    if(a=='T')
        return 1;
    if(a=='C')
        return 2;
    return 3;
}

//caculate tm
float tm(char seq[],float deltah[],float deltas[],int length)
{
    int i,pos;
    float total_deltah,total_deltas,result;

    total_deltah=0;
    total_deltas=0;
    for(i=0;i<length-1;i++)
    {
        pos=translate(seq[i]);
        pos=pos*4+translate(seq[i+1]);
        total_deltah+=deltah[pos];
        total_deltas+=deltas[pos];
    }

    total_deltah=(-1.0)*total_deltah;
    total_deltas=(-1.0)*total_deltas;
    if((seq[0]=='A')||(seq[0]=='T'))
    {
        total_deltah+=2.3;
        total_deltas+=4.1;
    }
    else
    {
        total_deltah+=0.1;
        total_deltas-=2.8;
    }
        if((seq[length-1]=='A')||(seq[length-1]=='T'))
        {
                total_deltah+=2.3;
                total_deltas+=4.1;
        }
        else
        {
                total_deltah+=0.1;
                total_deltas-=2.8;
        }
    result=1000.0*total_deltah/(total_deltas-0.51986*(length-1)-36.70381)-273.15;
    return result;
}

///caculate stability, int strand: 0 is 5' and 1 is 3'
void stability(char seq[],float stab[],int length,float Svalue[])
{
    int i,pos;
    
    pos=0;
    for(i=0;i<6;i++)
        pos=pos*4+translate(seq[i]);
    Svalue[0]=stab[pos];
//the other part
        pos=0;
        for(i=0;i<6;i++)
        pos=pos*4+translate(seq[i+length-6]);
    Svalue[1]=stab[pos];
}

//whether species chars in reads
int words(char *seq,int position,int length)
{
    int i;
    
    for(i=0;i<length;i++)
    {
        if(seq[position+i]=='N')
        {
            return 0;
        }
    }
    return 1;
}

int check_long_ploy(char primer[],int length)
{
    int i,same;
    char ref;

    same=1;
    ref=primer[0];
    for(i=1;i<length;i++)
    {
        if(primer[i]==ref)
            same++;
        else
        {
            if(same>=6)
                return 0;
            same=1;
            ref=primer[i];
        }
    }
    if(same>=6)
        return 0;
    return 1;
}

int dimer(char primer[],int length)
{
//same
    if(primer[length-1]==primer[length-2]&&primer[length-1]==primer[length-3]&&primer[length-1]==primer[length-4])
        return 0;
    if(primer[length-1]=='A')
    {
        if(primer[length-6]!='T')
            return 1;
    }
    else if(primer[length-1]=='T')
    {
        if(primer[length-6]!='A')
            return 1;
    }
    else if(primer[length-1]=='C')
    {
        if(primer[length-6]!='G')
            return 1;
    }
    else
    {
        if(primer[length-6]!='C')
            return 1;
    }

    if(primer[length-2]=='A')
        {        
                if(primer[length-5]!='T')
                        return 1;        
        }
        else if(primer[length-2]=='T')
        {        
                if(primer[length-5]!='A')
                        return 1;
        }                
        else if(primer[length-2]=='C')
        {
                if(primer[length-5]!='G')
                        return 1;  
        }
        else
        {
                if(primer[length-5]!='C')   
                        return 1;  
        }

    if(primer[length-3]=='A')
        {        
                if(primer[length-4]!='T')
                        return 1;        
        }
        else if(primer[length-3]=='T')
        {        
                if(primer[length-4]!='A')
                        return 1;
        }                
        else if(primer[length-3]=='C')
        {
                if(primer[length-4]!='G')
                        return 1;  
        }
        else
        {
                if(primer[length-4]!='C')   
                        return 1;  
        }
    return 0;
}
///function: int length: the length of genome
void candidate_primer(char *seq,int flag[],FILE *Inner,FILE *Outer,FILE *Loop,int Num[],float stab[],float deltah[],float deltas[],int length,double stackEntropies[],double stackEnthalpies[],double stackint2Entropies[],double stackint2Enthalpies[],double dangleEntropies3[],double dangleEnthalpies3[],double dangleEntropies5[],double dangleEnthalpies5[],double hairpinLoopEntropies[],double interiorLoopEntropies[],double bulgeLoopEntropies[],double hairpinLoopEnthalpies[],double interiorLoopEnthalpies[],double bulgeLoopEnthalpies[],double tstackEntropies[],double tstackEnthalpies[],double tstack2Entropies[],double tstack2Enthalpies[],char *triloopEntropies1,char *triloopEnthalpies1,char *tetraloopEntropies1,char *tetraloopEnthalpies1,double *triloopEntropies2,double *triloopEnthalpies2,double *tetraloopEntropies2,double *tetraloopEnthalpies2,int numTriloops,int numTetraloops,double atpS[],double atpH[])
{
    int i,circle,check,inner_plus,inner_minus,outer_plus,outer_minus,loop_plus,loop_minus,plus,minus;
    char primer[30],rev[30],*file;
    double secondary;
    float GC_content,Tm,Svalue[2];

    for(circle=0;circle<=(length-15);circle++)
    {
        for(i=15;i<=25;i++)  //read length is from 18 to 25
        {
            if(circle+i>length)
                continue;
            check=words(seq,circle,i);
            if(check==0)
                                break;
            memset(primer,'\0',30);
            generate(seq,primer,circle,i);

            //  check=check_long_ploy(primer,i);
            //  if(check==0)
            //      break;
            GC_content=gc(primer,i);
            if(GC_content<30||GC_content>70)
                continue;
            Tm=tm(primer,deltah,deltas,i);
            if(Tm>68||Tm<55)
                continue;

            reverse(primer,rev,i);
            plus=dimer(primer,i);
            minus=dimer(rev,i);
            if(flag[8]&&plus)
            {
                secondary=thal(primer,primer,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,1);
                if(secondary>Tm-10)
                    plus=0;
            }
            if(flag[8]&&plus)
            {
                secondary=thal(primer,primer,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,2);
                if(secondary>Tm-10)
                    plus=0;
            }
            if(flag[8]&&plus)
            {
                secondary=thal(primer,primer,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,4);
                if(secondary>Tm-10)
                    plus=0;
            }

            if(flag[8]&&minus)
            {
                secondary=thal(rev,rev,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,1);
                if(secondary>Tm-10)
                    minus=0;
            }
            
            if(flag[8]&&minus)
            {
                secondary=thal(rev,rev,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,2);
                if(secondary>Tm-10)
                    minus=0;
            }
            if(flag[8]&&minus)
            {
                secondary=thal(rev,rev,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,4);
                if(secondary>Tm-10)
                    minus=0;
            }
            if(plus+minus==0)
                continue;
            
            inner_plus=0;
            inner_minus=0;
            outer_plus=0;
            outer_minus=0;
            if(flag[7])
            {
                loop_plus=0;
                loop_minus=0;
            }
            if(plus==1)
            {
                stability(primer,stab,i,Svalue);
                
                //inner
                if(Svalue[0]>=4&&Svalue[1]>=3&&i<=22&&Tm>=64&&GC_content>=40)
                    inner_plus++; //GC-rich
                if(Svalue[0]>=4&&Svalue[1]>=3&&i>=20&&Tm>=60&&Tm<=63&&GC_content<=65)
                    inner_plus=inner_plus+2; //AT-rich
                if(Svalue[0]>=4&&Svalue[1]>=3&&i>=20&&i<=22&&Tm>=64&&Tm<=66&&GC_content>=40&&GC_content<=65)
                    inner_plus=inner_plus+4;

                //outer
                if(Svalue[0]>=3&&Svalue[1]>=4&&i<=20&&Tm>=59&&Tm<=63&&GC_content>=40)
                    outer_plus++;
                if(Svalue[0]>=3&&Svalue[1]>=4&&i>=18&&Tm<=58&&GC_content<=65)
                    outer_plus+=2;
                if(Svalue[0]>=3&&Svalue[1]>=4&&i>=18&&i<=20&&Tm>=59&&Tm<=61&&GC_content>=40&&GC_content<=65)
                    outer_plus+=4;
                
                //loop
                if(flag[7]&&Svalue[0]>=3&&Svalue[1]>=4&&i<=22&&Tm>=64&&GC_content>=40)
                    loop_plus++; //GC-rich
                if(flag[7]&&Svalue[0]>=3&&Svalue[1]>=4&&i>=20&&Tm>=60&&Tm<=63&&GC_content<=65)
                    loop_plus=loop_plus+2; //AT-rich
                if(flag[7]&&Svalue[0]>=3&&Svalue[1]>=4&&i>=20&&i<=22&&Tm>=64&&Tm<=66&&GC_content>=40&&GC_content<=65)
                    loop_plus=loop_plus+4;
            }
            
            if(minus==1)
            {
                stability(rev,stab,i,Svalue);
                //inner
                if(Svalue[0]>=4&&Svalue[1]>=3&&i<=22&&Tm>=64&&GC_content>=40)
                    inner_minus++; //GC-rich
                if(Svalue[0]>=4&&Svalue[1]>=3&&i>=20&&Tm>=60&&Tm<=63&&GC_content<=65)
                    inner_minus=inner_minus+2; //AT-rich
                if(Svalue[0]>=4&&Svalue[1]>=3&&i>=20&&i<=22&&Tm>=64&&Tm<=66&&GC_content>=40&&GC_content<=65)
                    inner_minus=inner_minus+4;

                //outer
                if(Svalue[0]>=3&&Svalue[1]>=4&&i<=20&&Tm>=59&&Tm<=63&&GC_content>=40)
                    outer_minus++;
                if(Svalue[0]>=3&&Svalue[1]>=4&&i>=18&&Tm<=58&&GC_content<=65)
                    outer_minus+=2;
                if(Svalue[0]>=3&&Svalue[1]>=4&&i>=18&&i<=20&&Tm>=59&&Tm<=61&&GC_content>=40&&GC_content<=65)
                    outer_minus+=4;
                
                //loop
                if(flag[7]&&Svalue[0]>=3&&Svalue[1]>=4&&i<=22&&Tm>=64&&GC_content>=40)
                    loop_minus++; //GC-rich
                if(flag[7]&&Svalue[0]>=3&&Svalue[1]>=4&&i>=20&&Tm>=60&&Tm<=63&&GC_content<=65)
                    loop_minus=loop_minus+2; //AT-rich
                if(flag[7]&&Svalue[0]>=3&&Svalue[1]>=4&&i>=20&&i<=22&&Tm>=64&&Tm<=66&&GC_content>=40&&GC_content<=65)
                    loop_minus=loop_minus+4;
            }

            if(inner_plus||inner_minus)
            {
                fprintf(Inner,"pos:%d\tlength:%d\t+:%d\t-:%d\t%0.2f\n",circle,i,inner_plus,inner_minus,Tm);
                Num[0]++;
            }
            
            if(outer_plus||outer_minus)
            {
                fprintf(Outer,"pos:%d\tlength:%d\t+:%d\t-:%d\t%0.2f\n",circle,i,outer_plus,outer_minus,Tm);
                Num[1]++;
            }
            
            if(flag[7]&&(loop_plus||loop_minus))
            {
                fprintf(Loop,"pos:%d\tlength:%d\t+:%d\t-:%d\t%0.2f\n",circle,i,loop_plus,loop_minus,Tm);
                Num[2]++;
            }
        }
    }
    return;
}

void SINGLE(const char *reference, const char *prefix, char *par_path, char *store_path_py) {
    double stackEntropies[625],stackEnthalpies[625],stackint2Entropies[625],stackint2Enthalpies[625],dangleEntropies3[125],dangleEnthalpies3[125],dangleEntropies5[125],dangleEnthalpies5[125];
    double hairpinLoopEntropies[30],interiorLoopEntropies[30],bulgeLoopEntropies[30],hairpinLoopEnthalpies[30],interiorLoopEnthalpies[30],bulgeLoopEnthalpies[30],tstackEntropies[625],tstackEnthalpies[625],tstack2Entropies[625],tstack2Enthalpies[625];
    char *triloopEntropies1,*triloopEnthalpies1,*tetraloopEntropies1,*tetraloopEnthalpies1;
    double *triloopEntropies2,*triloopEnthalpies2,*tetraloopEntropies2,*tetraloopEnthalpies2,atpS[25],atpH[25];
    char *store_path,*stab_path,*tm_path,*curren_path,*input, *seq, *temp;
    FILE *fp,*Inner,*Outer,*Loop;
    int i,flag[10],length,numTriloops,numTetraloops,Num[3];
    float deltah[16],deltas[16],stab[4096],temp1,temp2;
    struct stat statbuf;
    time_t start,end;

    start=time(NULL);
    for(i=0;i<10;i++){
        flag[i]=1;
    }

    // -dir (Directory)
    length=strlen(store_path_py);
    store_path=(char *)malloc(length+1);
    memset(store_path,'\0',length+1);
    strcpy(store_path,store_path_py);

    // -in (reference)
    length=strlen(reference);
    seq=(char *)malloc(length*sizeof(char));
    strcpy(seq, reference);

    length=strlen(store_path)+strlen(prefix)+10;
    temp=(char *)malloc(length);
    memset(temp,'\0',length);
    strcpy(temp,store_path);
    strcat(temp,"Inner/");
    mkdir(temp,0755);
    strcat(temp,prefix);
    Inner=fopen(temp,"w");
    if(Inner==NULL)
    {
        printf("Error! Can't create the %s file!\n",temp);
        exit(1);
    }

    memset(temp,'\0',length);
    strcpy(temp,store_path);
    strcat(temp,"Outer/");
    mkdir(temp,0755);
    strcat(temp,prefix);
    Outer=fopen(temp,"w");
    if(Outer==NULL)
    {
        printf("Error! Can't create the %s file!\n",temp);
        exit(1);
    }

    if(flag[7]==1)
    {
        memset(temp,'\0',length);       
        strcpy(temp,store_path);      
        strcat(temp,"Loop/");
        mkdir(temp,0755);
        strcat(temp,prefix);
        Loop=fopen(temp,"w");
        if(Loop==NULL)
        {
            printf("Error! Can't create the %s file!\n",temp);
            exit(1);       
        }
    }
    free(temp);

    //stability parameter file
    length=strlen(par_path);
    stab_path=(char *)malloc(length+30);
    memset(stab_path,'\0',length+30);
    strcpy(stab_path,par_path);
    strcat(stab_path,"stab_parameter.txt");
    
    //tm parameter file
    tm_path=(char *)malloc(length+30);
    memset(tm_path,'\0',length+30);
    strcpy(tm_path,par_path);
    strcat(tm_path,"tm_nn_parameter.txt");

    if(flag[8])
    {
        getStack(stackEntropies,stackEnthalpies,par_path);
        getStackint2(stackint2Entropies,stackint2Enthalpies,par_path);
        getDangle(dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,par_path);
        getLoop(hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,par_path);
        getTstack(tstackEntropies,tstackEnthalpies,par_path);
        getTstack2(tstack2Entropies,tstack2Enthalpies,par_path);

        numTriloops=get_num_line(par_path,0);
        triloopEntropies1=(char *)malloc(numTriloops*5);
        triloopEnthalpies1=(char *)malloc(numTriloops*5);
        triloopEntropies2=(double *)malloc(numTriloops*sizeof(double));
        triloopEnthalpies2=(double *)malloc(numTriloops*sizeof(double));
        getTriloop(triloopEntropies1,triloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,par_path);
        
        numTetraloops=get_num_line(par_path,1);
        tetraloopEntropies1=(char *)malloc(numTetraloops*6);
        tetraloopEnthalpies1=(char *)malloc(numTetraloops*6);
        tetraloopEntropies2=(double *)malloc(numTetraloops*sizeof(double));
        tetraloopEnthalpies2=(double *)malloc(numTetraloops*sizeof(double));
        getTetraloop(tetraloopEntropies1,tetraloopEnthalpies1,tetraloopEntropies2,tetraloopEnthalpies2,par_path);
        tableStartATS(6.9,atpS);
        tableStartATH(2200.0,atpH);
    }

    length=strlen(seq);
    
    //input Tm parameter
    fp=fopen(tm_path,"r");  //read the paramter of deltah and deltas
    if(fp==NULL)
    {
        printf("Error: can't open the %s file!\n",tm_path);
        exit(1);
    }
    while(fscanf(fp,"%d\t%f\t%f",&i,&temp1,&temp2)!=EOF)
    {
        deltah[i]=temp1;
        deltas[i]=temp2;
    }
    fclose(fp);

    //input stability parameter
    fp=fopen(stab_path,"r");  //read the parameters of stability
    if(fp==NULL)
    {
        printf("Error: can't open the %s file!\n",stab_path);
        exit(1);
    }
    while(fscanf(fp,"%d\t%f",&i,&temp1)!=EOF)
    {
        stab[i]=temp1;
    }
    fclose(fp);

    end=time(NULL);
    printf("Single LAMPgRNAtor - It takes %d seconds to prepare.\n",(int)difftime(end,start));
    start=time(NULL);
    Num[0]=0;
    Num[1]=0;
    Num[2]=0;
    candidate_primer(seq,flag,Inner,Outer,Loop,Num,stab,deltah,deltas,length,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH);
    printf("Single LAMPgRNAtor - There ara %d candidate primers used as F3/F2/B2/B3.\n",Num[1]);
    fclose(Outer);
    printf("Single LAMPgRNAtor - There are %d candidate primers used as F1c/B1c.\n",Num[0]);
    fclose(Inner);
    if(flag[7]==1)
    {
        printf("Single LAMPgRNAtor - There are %d candidate primers used as LF/LB.\n",Num[2]);
        fclose(Loop);
    }
    
    //check
    if(Num[1]<4)
        printf("Single LAMPgRNAtor - Warning: there don't have enough primers(>=4) used as F3/F2/B2/B3.\n");
    if(Num[0]<2)
        printf("Single LAMPgRNAtor - Warning: there don't have enough primers(>=2) used as F1c/B1c.\n");
    if(flag[7]==1 && Num[2]<1)
        printf("Single LAMPgRNAtor - Warning: there don't have enough primers(>=1) used as LF/LB. But you can design LAMP primers without loop primer.\n");

    end=time(NULL);
    printf("Single LAMPgRNAtor - It takes %d seconds to identify candidate single primer regions.\n",(int)difftime(end,start));

    free(stab_path);
    free(tm_path);
    if(flag[8])
    {
        free(triloopEntropies1);
        free(triloopEnthalpies1);
        free(tetraloopEntropies1);
        free(tetraloopEnthalpies1);
        free(triloopEntropies2);
        free(triloopEnthalpies2);
        free(tetraloopEntropies2);
        free(tetraloopEnthalpies2);
    }
}