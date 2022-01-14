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

struct Node
{
    int pos;
    int gi;
    int plus;  //as a flag, 1 is OK, 0 is no
    int minus; //as a flag
    struct Node *next;
};

struct INFO
{
    char name[301];
    int turn;
    struct INFO *next;
};

struct Primer
{
    int pos;
    int len;
    int plus;
    int minus;
    int total;  //how many GIs can be used
    float Tm;
    struct Primer *self;
    struct Primer *loop;
    struct Primer *notloop;
    struct Primer *next;
    struct Node *common;
    struct Node *special;
};

//get the file size
int file_size2(char *filename)
{
    int size;
        struct stat statbuf;

        stat(filename,&statbuf);
        size=statbuf.st_size;
        return size;
}

void generate_primer(char *seq,char primer[],int start,int length,int flag) //flag=0:plus
{
        int i;
        
        if(flag==0)
        {
                for(i=0;i<length;i++)
                        primer[i]=seq[start+i];
                primer[i]='\0';
        }
        else
        {
                for(i=0;i<length;i++)
                {
                        if((seq[start+length-1-i]=='A')||(seq[start+length-1-i]=='a'))
                                primer[i]='T';
                        else if((seq[start+length-1-i]=='T')||(seq[start+length-1-i]=='t'))
                                primer[i]='A';
                        else if((seq[start+length-1-i]=='C')||(seq[start+length-1-i]=='c'))
                                primer[i]='G';
                        else 
                                primer[i]='C';
                }
                primer[length]='\0';
        }
}

int check_same(struct Primer *one[],struct Primer *two[],struct Primer *p_F3,struct Primer *p_F2,struct Primer *p_F1c,struct Primer *p_B1c,struct Primer *p_B2,struct Primer *p_B3)
{
    int i;
    for(i=0;i<10;i++)
    {
        if(one[i]==NULL)
            return 0;
        if(one[i]==p_F3)
        {
            if(two[i]==p_F2)
                return 1;
            if(two[i]==p_F1c)
                return 1;
            if(two[i]==p_B1c)
                return 1;
            if(two[i]==p_B2)
                return 1;
            if(two[i]==p_B3)
                return 1;
        }
        else if(one[i]==p_F2)
        {
            if(two[i]==p_F1c)
                                return 1;
                        if(two[i]==p_B1c)
                                return 1;
                        if(two[i]==p_B2)
                                return 1;
                        if(two[i]==p_B3)
                                return 1;
        }
        else if(one[i]==p_F1c)
        {
            if(two[i]==p_B1c)
                                return 1;
                        if(two[i]==p_B2)
                                return 1;
                        if(two[i]==p_B3)
                                return 1;
        }
        else if(one[i]==p_B1c)
        {
            if(two[i]==p_B2) 
                                return 1;
                        if(two[i]==p_B3)
                                return 1;
        }
        else if(one[i]==p_B2)
        {
            if(two[i]==p_B3)
                                return 1;
        }
    }
    return 0;
}

int add_same(struct Primer *one[],struct Primer *two[],struct Primer *p_F3,struct Primer *p_F2,struct Primer *p_F1c,struct Primer *p_B1c,struct Primer *p_B2,struct Primer *p_B3,int error[],int replace)
{
    if(error[0]==0)
    {
        one[replace]=p_F3;
        if(error[1]==1)
            two[replace]=p_F2;
        else if(error[1]==2)
            two[replace]=p_F1c;
        else if(error[1]==3)
            two[replace]=p_B1c;
        else if(error[1]==4)
            two[replace]=p_B2;
        else
            two[replace]=p_B3;
    }
    else if(error[0]==1)
    {
        one[replace]=p_F2;
        if(error[1]==2)
                        two[replace]=p_F1c;
                else if(error[1]==3)
                        two[replace]=p_B1c;
                else if(error[1]==4)
                        two[replace]=p_B2;
                else
                        two[replace]=p_B3;
    }
    else if(error[0]==2)
        {
                one[replace]=p_F1c;
        if(error[1]==3)
                        two[replace]=p_B1c;
                else if(error[1]==4)
                        two[replace]=p_B2;
                else
                        two[replace]=p_B3;
        }
    else if(error[0]==3)
        {
        one[replace]=p_B1c;
        if(error[1]==4)
                        two[replace]=p_B2;
                else
                        two[replace]=p_B3;
        }
    else
        {
                one[replace]=p_B2;
        two[replace]=p_B3;
        }
    replace++;
    if(replace==10)
        replace=0;
    return replace;
}

void DESIGNER(const char* prefix, const char* reference, int minLen, const char* output, char *par_path, char *store_path, int expect) {
    double stackEntropies[625],stackEnthalpies[625],stackint2Entropies[625],stackint2Enthalpies[625],dangleEntropies3[125],dangleEnthalpies3[125],dangleEntropies5[125],dangleEnthalpies5[125];
    double hairpinLoopEntropies[30],interiorLoopEntropies[30],bulgeLoopEntropies[30],hairpinLoopEnthalpies[30],interiorLoopEnthalpies[30],bulgeLoopEnthalpies[30],tstackEntropies[625],tstackEnthalpies[625],tstack2Entropies[625],tstack2Enthalpies[625];
    double *triloopEntropies2,*triloopEnthalpies2,*tetraloopEntropies2,*tetraloopEnthalpies2,atpS[25],atpH[25];
    int i,j,length,*result,new,*best_par,flag[13],common_num[1],turn,success,numTriloops,numTetraloops,max_loop,min_loop,error[2],replace,GC;
    char *path_fa,*inner,*outer,*loop;
    char *temp,*seq,F3[26],F2[26],F1c[26],B1c[26],B2[26],B3[26],LF[26],LB[26],*triloopEntropies1,*triloopEnthalpies1,*tetraloopEntropies1,*tetraloopEnthalpies1;
    FILE *fp;
    struct Primer *headL,*headS,*headLoop,*p_F3,*p_F2,*p_F1c,*p_B1c,*p_B2,*p_B3,*result_Loop[2],*one[10],*two[10];
    struct Node *p_node,*p_temp;
    struct INFO *p_list;
    time_t start,end,begin;

    struct Primer *read_par();
    struct INFO *read_list();
    void next_one();
    void add();
    int check_add();
    int check_gc();
    int check_structure();
    int design_loop();

    start=time(NULL);
    begin=start;
        
    for(i=0;i<13;i++){
        flag[i]=1;
    }

    flag[5]=0; //Common
    flag[6]=0; //Specific
    flag[12]=0;
    //expect=50;

    // -in
    length=strlen(reference);
    seq=(char *)malloc(length*sizeof(char));
    strcpy(seq, reference);
    
    if(flag[7])
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

    j=strlen(store_path)+strlen(prefix)+12;
    outer=(char *)malloc(j);
    memset(outer,'\0',j);
    strcpy(outer,store_path);

    inner=(char *)malloc(j);
    memset(inner,'\0',j);
    strcpy(inner,store_path);

    if(flag[10]==1)
    {
        loop=(char *)malloc(j);
        memset(loop,'\0',j);
        strcpy(loop,store_path);
    }
    
    strcat(outer,"Outer/");
    strcat(outer,prefix);
    strcat(inner,"Inner/");
    strcat(inner,prefix);
    if(flag[10]==1)
    {
        strcat(loop,"Loop/");
        strcat(loop,prefix);
    }

    common_num[0]=1;

    //read parameters
    headS=read_par(outer,0,0);
    headL=read_par(inner,0,0);
    if(flag[10])
    {
        headLoop=read_par(loop,0,0);//don't use special info
        p_F3=headLoop;
        while(p_F3->next!=NULL)
            p_F3=p_F3->next;
        max_loop=p_F3->pos;
        min_loop=headLoop->pos;
    }

    //the next one 
    next_one(headL,headL,0);
    next_one(headS,headS,0);
    next_one(headS,headL,1);
    next_one(headL,headS,1);
    if(flag[10]==1)
    {
        next_one(headS,headLoop,2);
        next_one(headL,headLoop,2);
    }

    if(flag[5]) result=(int *)malloc(common_num[0]*sizeof(int));
    best_par=(int *)malloc(expect*sizeof(int)); //include LF/LB
    for(i=0;i<expect;i++)
        best_par[i]=-1; 
    fp=fopen(output,"w");
    if(fp==NULL)
    {
        printf("Error: can't create the %s file!\n",output);
        exit(1);
    }
    end=time(NULL);
    printf("Design LAMPgRNAtor - It takes %0.0f seconds to prepare data.\n",difftime(end,start));

    //design LAMP primers
    start=time(NULL);
    if(flag[12]==0)
    {
        for(i=0;i<10;i++)
        {
            one[i]=NULL;
            two[i]=NULL;
        }
        replace=0;
    }
    turn=0;
    for(j=common_num[0];j>=1;j--)
    {
        for(p_F3=headS;p_F3;p_F3=p_F3->next)   //F3
        {
            if(flag[10]&&(p_F3->pos+200)<min_loop)
                continue;
            if(flag[10]&&(p_F3->pos-200)>max_loop)
                break;
            if(p_F3->total<j)
                continue;
            if(p_F3->plus==0)
                continue;
            success=check_add(p_F3->pos,best_par,expect); //in best_par are sorted primers by common
            if(success==0)
                continue;
            new=0;  //whether find a new primer, if find, adjust F3 position
            for(p_F2=p_F3->self;p_F2;p_F2=p_F2->next)  //F2
                {
                flag[0]=p_F2->plus&p_F3->plus;
                if(flag[0]==0)
                    continue;
                if(p_F2->total<j)
                    continue;
                if(p_F2->pos-(p_F3->pos+p_F3->len)>20)
                    break;
                for(p_F1c=p_F2->notloop;p_F1c;p_F1c=p_F1c->next)   //F1c
                {
                    flag[1]=flag[0]&p_F1c->minus;
                    if(flag[1]==0)
                        continue;
                    if(p_F1c->total<j)
                        continue;
                    //Tm
                    if(p_F1c->Tm-p_F3->Tm<3)
                        continue;
                    if(p_F1c->Tm-p_F2->Tm<3)
                        continue;

                    if(p_F1c->pos-p_F2->pos-1<40)
                        continue;
                    if(p_F1c->pos-p_F2->pos-1>60)
                        break;
                    for(p_B1c=p_F1c->self;p_B1c;p_B1c=p_B1c->next)   //B1c
                    {
                        flag[2]=flag[1]&p_B1c->plus;
                        if(flag[2]==0)
                            continue;
                        if(p_B1c->total<j)
                            continue;
                    //Tm
                        if(p_B1c->Tm-p_F3->Tm<3)
                            continue;
                        if(p_B1c->Tm-p_F2->Tm<3)
                            continue;

                        if(p_B1c->pos-p_F1c->pos<minLen)
                            continue;
                        if(p_B1c->pos-p_F1c->pos>85)
                            break;
                        for(p_B2=p_B1c->notloop;p_B2;p_B2=p_B2->next)   //B2
                        {
                            flag[3]=flag[2]&p_B2->minus;
                            if(flag[3]==0)
                                continue;
                            if(p_B2->total<j)
                                continue;
                            //Tm
                            if(p_F1c->Tm-p_B2->Tm<3)
                                continue;
                            if(p_B1c->Tm-p_B2->Tm<3)
                                continue;

                            if((p_B2->pos+p_B2->len-1)-(p_B1c->pos+p_B1c->len-1)-1<40)
                                continue;
                            if((p_B2->pos+p_B2->len-1)-(p_B1c->pos+p_B1c->len-1)-1>60)
                                break;
                            if(p_B2->pos+p_B2->len-1-p_F2->pos-1<120)
                                continue;
                            if(p_B2->pos+p_B2->len-1-p_F2->pos-1>180)
                                break;
                            
                            //check whether has enough positions for loop
                            if(flag[10]&&((p_F1c->pos-p_F2->pos-p_F2->len)<18)&&(p_B2->pos-p_B1c->pos-p_B1c->len<18))
                                continue;
                            for(p_B3=p_B2->self;p_B3;p_B3=p_B3->next)  //B3
                            {
                                if(p_B3->pos<p_B2->pos+p_B2->len)
                                    continue;
                                if(p_B3->pos-(p_B2->pos+p_B2->len)>20)
                                    break;
                                if(p_B3->total<j)
                                    continue;
                                
                                flag[4]=flag[3]&p_B3->minus;
                                if(flag[4]==0)
                                    continue;
                                GC=check_gc(seq,p_F3->pos,(p_B3->pos+p_B3->len));
                                if((GC&flag[4])!=GC)
                                    continue;
                                
                                //Tm
                                if(p_F1c->Tm-p_B3->Tm<3)
                                    continue;
                                if(p_B1c->Tm-p_B3->Tm<3)
                                    continue;

                                if(flag[12])
                                    new=1;
                                else
                                {
                                    success=check_same(one,two,p_F3,p_F2,p_F1c,p_B1c,p_B2,p_B3);
                                    if(success)
                                    {
                                        if(new)
                                            break;
                                        else
                                            continue;
                                    }
                                }
                                
                                //check second structure
                                generate_primer(seq,F3,p_F3->pos,p_F3->len,0);
                                generate_primer(seq,F2,p_F2->pos,p_F2->len,0);
                                generate_primer(seq,F1c,p_F1c->pos,p_F1c->len,1);
                                generate_primer(seq,B1c,p_B1c->pos,p_B1c->len,0);
                                generate_primer(seq,B2,p_B2->pos,p_B2->len,1);
                                generate_primer(seq,B3,p_B3->pos,p_B3->len,1);
                                if(flag[7])
                                {
                                    success=check_structure(F3,F2,F1c,B1c,B2,B3,GC,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,error);
                                    if(success==0)
                                    {
                                        replace=add_same(one,two,p_F3,p_F2,p_F1c,p_B1c,p_B2,p_B3,error,replace);
                                        if(new)
                                            break;
                                        else
                                            continue;
                                    }
                                }
                                
                                //design loop
                                if(flag[10])
                                {
                                    success=design_loop(p_F3,p_F2,p_F2->loop,p_F1c,p_B1c,p_B1c->loop,p_B2,p_B3,result_Loop,result,flag,common_num[0],seq,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,GC);
                                    if(success==0)
                                    {
                                        if(new)
                                            break;
                                        else
                                            continue;
                                    }
                                }
                                best_par[turn]=p_F3->pos;

                            //output
                                turn++;
                                fprintf(fp,"The %d LAMP primers:\n",turn);
                                fprintf(fp,"  F3: pos:%d,length:%d bp, primer(5'-3'):%s\n",p_F3->pos,p_F3->len,F3);
                                fprintf(fp,"  F2: pos:%d,length:%d bp, primer(5'-3'):%s\n",p_F2->pos,p_F2->len,F2);
                                fprintf(fp,"  F1c: pos:%d,length:%d bp, primer(5'-3'):%s\n",p_F1c->pos,p_F1c->len,F1c);
                                fprintf(fp,"  B1c: pos:%d,length:%d bp, primer(5'-3'):%s\n",p_B1c->pos,p_B1c->len,B1c);
                                fprintf(fp,"  B2: pos:%d,length:%d bp, primer(5'-3'):%s\n",p_B2->pos,p_B2->len,B2);
                                fprintf(fp,"  B3: pos:%d,length:%d bp, primer(5'-3'):%s\n",p_B3->pos,p_B3->len,B3);
                                if(flag[10])
                                {
                                    if(result_Loop[0]==NULL)
                                        fprintf(fp,"  LF: NULL\n");
                                    else
                                    {
                                        generate_primer(seq,LF,result_Loop[0]->pos,result_Loop[0]->len,1);
                                        fprintf(fp,"  LF: pos:%d,length:%d bp, primer(5'-3'):%s\n",result_Loop[0]->pos,result_Loop[0]->len,LF);
                                    }
                                    if(result_Loop[1]==NULL)
                                        fprintf(fp,"  LB: NULL\n");
                                    else
                                    {
                                        generate_primer(seq,LB,result_Loop[1]->pos,result_Loop[1]->len,0);
                                        fprintf(fp,"  LB: pos:%d,length:%d bp, primer(5'-3'):%s\n",result_Loop[1]->pos,result_Loop[1]->len,LB);
                                    }
                                }
                                
                                end=time(NULL);
                                printf("Design LAMPgRNAtor - It takes %0.0f seconds to design the %d-th LAMP primer set successfully.\n",difftime(end,start),turn);
                                start=time(NULL);
                                new=1;
                                break;
                            } //B2
                            if(turn>=expect||new)
                                break;
                                                }  //B1c
                        if(turn>=expect||new)
                            break;
                                        } //F1c
                    if(turn>=expect||new)
                        break;
                                } //F2
                if(turn>=expect||new)
                    break;
                        }  //B3
            if(turn>=expect)
                break;
                }  //F3
        if(turn>=expect)
            break;
        }
    fclose(fp);
    free(best_par);
    free(seq);
    free(inner);
    free(outer);

    //free struct list
    while(headL)
    {
        p_node=headL->common;
        while(p_node)
        {
            p_temp=p_node->next;
            free(p_node);
            p_node=p_temp;
        }
        p_node=headL->special;
        while(p_node)  
            {
                p_temp=p_node->next;  
                free(p_node);  
                p_node=p_temp;  
            }
        
        p_F3=headL->next;
        free(headL);
        headL=p_F3;
    }

    while(headS)
    {
        p_node=headS->common;  
        while(p_node)  
        {
            p_temp=p_node->next;  
            free(p_node);  
            p_node=p_temp;  
        }
        p_node=headS->special;
        while(p_node)
        {               
            p_temp=p_node->next;
            free(p_node);
            p_node=p_temp;
        }

        p_F3=headS->next;
        free(headS);
        headS=p_F3;
    }

    if(flag[7])
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

    if(flag[10])
    {
        free(loop);
        while(headLoop)
        {
            p_node=headLoop->common;
            while(p_node)
            {
                    p_temp=p_node->next;
                    free(p_node);
                    p_node=p_temp;
            }
            p_node=headLoop->special;
            while(p_node)
            {
                    p_temp=p_node->next;
                    free(p_node);
                    p_node=p_temp;
            }

            p_F3=headLoop->next;
            free(headLoop);
            headLoop=p_F3;
        }
    }
    end=time(NULL);
    printf("Design LAMPgRNAtor - It takes %0.0f seconds to free memory.\n",difftime(end,start));
    printf("Design LAMPgRNAtor - It takes total %0.0f seconds to finish this design.\n",difftime(end,begin));
}

void merge(int *store,int *temp,int *data,int start,int end,int middle)
{
        int i,j,k;
        i=start;
        j=middle+1;
        k=start;
        while(i<=middle&&j<=end)
        {
                if(data[store[i]*6]<data[store[j]*6])
                {
                        temp[k]=store[i];
                        k++;
                        i++;
                        continue;
                }
                if(data[store[j]*6]<data[store[i]*6])
                {
                        temp[k]=store[j];
                        k++;
                        j++;
                        continue;
                }
        //len
                if(data[store[i]*6+1]<data[store[j]*6+1])
                {
                        temp[k]=store[i];
                        k++;
                        i++;
                        continue;
                }
                if(data[store[j]*6+1]<data[store[i]*6+1])
                {
                        temp[k]=store[j];
                        k++;
                        j++;
                        continue;
                }
        //gi
                if(data[store[i]*6+2]<data[store[j]*6+2])
                {
                        temp[k]=store[i];
                        k++;
                        i++;
                        continue;
                }
                if(data[store[j]*6+2]<data[store[i]*6+2])
                {
                        temp[k]=store[j];
                        k++;
                        j++;
                        continue;
                }
        //position
                if(data[store[i]*6+3]<data[store[j]*6+3])
                {
                        temp[k]=store[i];
                        k++;
                        i++;
                }
                else
                {
                        temp[k]=store[j];
                        k++;
                        j++;
                        continue;
                }
        }
        while(i<=middle)
        {
                temp[k]=store[i];
                k++;
                i++;
        }
        while(j<=end)
        {
                temp[k]=store[j];
                k++;
                j++;
        }
        for(i=start;i<=end;i++)
                store[i]=temp[i];
}

void sort_merge(int *store,int *temp,int *data,int start,int end)
{
        int middle;

        if(start<end)
        {
                middle=(start+end)/2;
                sort_merge(store,temp,data,start,middle);
                sort_merge(store,temp,data,(middle+1),end);
                merge(store,temp,data,start,end,middle);
        }
}

////function read primer informatin and align information 
struct Primer *read_par(char *path,int common_flag,int special_flag)
{
    char *in,line[100];
    int pos,len,gi,position,plus,minus,size,i,flag,total,*store,*temp,*data;
    float Tm;
    struct Primer *new_primer,*p_primer,*head;
    struct Node *new_node,*p_node;
    FILE *fp;

///read the  primer file
    if(access(path,0)==-1)
    {
        printf("Error! Don't have the %s file!\n",path);
        exit(1);
    }
        fp=fopen(path,"r");
        if(fp==NULL)
        {
                printf("Error: can't open the %s file!\n",path);
                exit(1);
        }
    
    size=sizeof(struct Primer);
    i=0;
        while(fscanf(fp,"pos:%d\tlength:%d\t+:%d\t-:%d\t%f\n",&pos,&len,&plus,&minus,&Tm)!=EOF)
        {
        new_primer=(struct Primer *)malloc(size);
        new_primer->pos=pos;
        new_primer->len=len;
        new_primer->total=1;
        new_primer->plus=plus;
        new_primer->minus=minus;
        new_primer->Tm=Tm;
        new_primer->next=NULL;
        new_primer->self=NULL;
        new_primer->loop=NULL;
        new_primer->notloop=NULL;
        new_primer->common=NULL;
        new_primer->special=NULL;

        if(i==0)
        {
            head=new_primer;
            p_primer=new_primer;
            i++;
        }
        else
        {
            p_primer->next=new_primer;
            p_primer=new_primer;
        }
        }
    fclose(fp);
        if(i==0)
        {
                printf("Sorry! Don't have any candidate single primers in %s!\n",path);
                exit(1);
        }

//parameter of common
    if(common_flag==1)
    {
        i=strlen(path);
        in=(char *)malloc(i+20);
            memset(in,'\0',i+20);
            strcpy(in,path);
            strcat(in,"-common.txt"); //suffix of parameter
        if(access(in,0)==-1)
        {
            printf("Error! Don't have the %s file!\n",in);
            exit(1);
        }

            fp=fopen(in,"r");
            if(fp==NULL)
            {
                    printf("Error: can't open the %s file!\n",in);
                    exit(1);
            }
        
        total=0;
        while(fgets(line,100,fp)!=NULL)
            total++;
        rewind(fp);
        store=(int *)malloc(total*sizeof(int));
        temp=(int *)malloc(total*sizeof(int));
        data=(int *)malloc(6*total*sizeof(int));
        total=0;
        while(fscanf(fp,"%d\t%d\t%d\t%d\t%d\t%d\n",&pos,&len,&gi,&position,&plus,&minus)!=EOF)
        {
            data[6*total]=pos;
            data[6*total+1]=len;
            data[6*total+2]=gi;
            data[6*total+3]=position;
            data[6*total+4]=plus;
            data[6*total+5]=minus;
            store[total]=total;
            total++;
        }
        fclose(fp);
        sort_merge(store,temp,data,0,(total-1));

        p_primer=head;
        size=sizeof(struct Node);
        flag=0;
        i=0;
        while(p_primer&&i<total)
        {
        //pos
            if(data[store[i]*6]<p_primer->pos)
            {
                i++;
                continue;
            }
            if(data[store[i]*6]>p_primer->pos)
            {
                p_primer=p_primer->next;
                flag=0;
                continue;
            }
        //len
            if(data[store[i]*6+1]<p_primer->len)
                        {
                                i++;
                                continue;
                        }
                        if(data[store[i]*6+1]>p_primer->len)
                        {
                                p_primer=p_primer->next;
                flag=0;
                                continue;
                        }
            new_node=(struct Node *)malloc(size);
            new_node->gi=data[store[i]*6+2];
            new_node->pos=data[store[i]*6+3];
            new_node->plus=data[store[i]*6+4];
            new_node->minus=data[store[i]*6+5];
                        new_node->next=NULL;
            if(flag==0)
            {
                flag++;
                p_primer->common=new_node;
                p_node=new_node;
            }
            else
            {
                            p_node->next=new_node;
                p_node=new_node;
            }
            i++;
            }
        free(in);
        free(data);
        free(store);
        free(temp);
    }
//paramter for special
    if(special_flag==1)
    {
        i=strlen(path);
        in=(char *)malloc(i+20);
        memset(in,'\0',i+20);
            strcpy(in,path);
            strcat(in,"-specific.txt"); //suffix of parameter
        if(access(in,0)==-1)
        {
            printf("Error! Don't have the %s file!\n",in);
            exit(1);
        }

            fp=fopen(in,"r");
            if(fp==NULL)
            {
                    printf("Error: can't open the %s file!\n",in);
                    exit(1);
            }
        total=0;
                while(fgets(line,100,fp)!=NULL)
                        total++;
                rewind(fp);
                store=(int *)malloc(total*sizeof(int));
                temp=(int *)malloc(total*sizeof(int));
                data=(int *)malloc(6*total*sizeof(int));
                total=0;
                while(fscanf(fp,"%d\t%d\t%d\t%d\t%d\t%d\n",&pos,&len,&gi,&position,&plus,&minus)!=EOF)
                {
                        data[6*total]=pos;
                        data[6*total+1]=len;
                        data[6*total+2]=gi;
                        data[6*total+3]=position;
                        data[6*total+4]=plus;
                        data[6*total+5]=minus;
                        store[total]=total;
                        total++;
                }
                fclose(fp);
                sort_merge(store,temp,data,0,(total-1));

                p_primer=head;
                size=sizeof(struct Node);
                flag=0;
        i=0;
                while(p_primer&&i<total)
                {
                //pos
                        if(data[store[i]*6]<p_primer->pos)
                        {
                                i++;
                                continue;
                        }
                        if(data[store[i]*6]>p_primer->pos)
                        {
                                p_primer=p_primer->next;
                flag=0;
                                continue;
                        }
                //len
                        if(data[store[i]*6+1]<p_primer->len)
                        {
                                i++;
                                continue;
                        }
                        if(data[store[i]*6+1]>p_primer->len)
                        {
                                p_primer=p_primer->next;
                flag=0;
                                continue;
                        }
                        new_node=(struct Node *)malloc(size);
                        new_node->gi=data[store[i]*6+2];
                        new_node->pos=data[store[i]*6+3];
                        new_node->plus=data[store[i]*6+4];
                        new_node->minus=data[store[i]*6+5];
                        new_node->next=NULL;
            if(flag==0)
            {
                flag++;
                p_primer->special=new_node;
                p_node=new_node;
            }
            else
            {
                            p_node->next=new_node;
                            p_node=new_node;
            }
                        i++;
                }
        free(in);
                free(data);
                free(store);
                free(temp);
    }
    return head;
}

struct INFO *read_list(char *path,int common_num[])
{
    char *in,name[301];
    int turn,i,size;
    struct INFO *new_primer,*p_primer,*head;
    FILE *fp;

    i=strlen(path);
    in=(char *)malloc(i+20);
    memset(in,'\0',i+20);
    strcpy(in,path);
    strcat(in,"-common_list.txt");
        if(access(in,0)==-1)
        {
                printf("Error! Don't have the %s file!\n",in);
                exit(1);
        }
        fp=fopen(in,"r");
        if(fp==NULL)
        {
                printf("Error: can't open the %s file!\n",in);
                exit(1);
        }
        
        size=sizeof(struct INFO);
        i=0;
        memset(name,'\0',301);
        while(fscanf(fp,"%s\t%d\n",name,&turn)!=EOF)
        {
                new_primer=(struct INFO *)malloc(size);
                new_primer->turn=turn;
                strcpy(new_primer->name,name);
                new_primer->next=NULL;

                if(i==0)
                {
                        head=new_primer;
                        p_primer=new_primer;
                        i++;
                }
                else
                {
                        p_primer->next=new_primer;
                        p_primer=new_primer;
                }
                memset(name,'\0',301);
        }
        fclose(fp);
    common_num[0]=turn;
    free(in);
    return head;
}

//function: the next one
void next_one(struct Primer *first, struct Primer *second,int flag) //0:self,1:notloop;2:loop
{
    struct Primer *one,*two,*start;
    int pos=-1;

    one=first;
    start=second;
    two=start;

    while(one)
    {
        if(pos!=one->pos)
        {
            while(start)
            {
                if(start->pos+18<one->pos)
                    start=start->next;
                else
                    break;
            }
            pos=one->pos;
        }
        //move second
        two=start;
        while(two)
        {
            if(two->pos<one->pos+one->len)
                two=two->next;
            else
            {
                if(flag==0)
                    one->self=two;
                else if(flag==1)
                    one->notloop=two;
                else
                    one->loop=two;
                break;
            }
        }
        one=one->next;
    }           
}

int check_structure(char F3[],char F2[],char F1c[],char B1c[],char B2[],char B3[],int GC,double stackEntropies[],double stackEnthalpies[],double stackint2Entropies[],double stackint2Enthalpies[],double dangleEntropies3[],double dangleEnthalpies3[],double dangleEntropies5[],double dangleEnthalpies5[],double hairpinLoopEntropies[],double interiorLoopEntropies[],double bulgeLoopEntropies[],double hairpinLoopEnthalpies[],double interiorLoopEnthalpies[],double bulgeLoopEnthalpies[],double tstackEntropies[],double tstackEnthalpies[],double tstack2Entropies[],double tstack2Enthalpies[],char *triloopEntropies1,char *triloopEnthalpies1,char *tetraloopEntropies1,char *tetraloopEnthalpies1,double *triloopEntropies2,double *triloopEnthalpies2,double *tetraloopEntropies2,double *tetraloopEnthalpies2,int numTriloops,int numTetraloops,double atpS[],double atpH[],int error[])
{
    int i,j,threshold;
    double TH;
    char *list[6],rev1[26],rev2[26];

//prepare
    if(GC==1||GC==4)
        threshold=49;
    else
        threshold=45;
    list[0]=F3;
    list[1]=F2;
    list[2]=F1c;
    list[3]=B1c;
    list[4]=B2;
    list[5]=B3;
    for(i=0;i<6;i++)
    {
        for(j=i+1;j<6;j++)
        {
            TH=thal(list[i],list[j],stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,1);
            if(TH>threshold)
            {
                error[0]=i;
                error[1]=j;
                return 0;
            }

            TH=thal(list[i],list[j],stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,2);
            if(TH>threshold)
            {
                error[0]=i;
                error[1]=j;
                                return 0;
            }
            TH=thal(list[i],list[j],stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,3);
            if(TH>threshold)
            {
                error[0]=i;
                error[1]=j;
                                return 0;
            }

            reverse(list[j],rev1,strlen(list[j]));
            reverse(list[i],rev2,strlen(list[i]));
            TH=thal(rev1,rev2,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,2);
            if(TH>threshold)
            {
                error[0]=i;
                error[1]=j;
                                return 0;
            }
            TH=thal(rev1,rev2,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,3);
            if(TH>threshold)
            {
                error[0]=i;
                error[1]=j;
                                return 0;
            }
        }
    }
    return 1;
}           

int check_structure_loop(char F3[],char F2[],char F1c[],char B1c[],char B2[],char B3[],char LF[],char LB[],int GC,double stackEntropies[],double stackEnthalpies[],double stackint2Entropies[],double stackint2Enthalpies[],double dangleEntropies3[],double dangleEnthalpies3[],double dangleEntropies5[],double dangleEnthalpies5[],double hairpinLoopEntropies[],double interiorLoopEntropies[],double bulgeLoopEntropies[],double hairpinLoopEnthalpies[],double interiorLoopEnthalpies[],double bulgeLoopEnthalpies[],double tstackEntropies[],double tstackEnthalpies[],double tstack2Entropies[],double tstack2Enthalpies[],char *triloopEntropies1,char *triloopEnthalpies1,char *tetraloopEntropies1,char *tetraloopEnthalpies1,double *triloopEntropies2,double *triloopEnthalpies2,double *tetraloopEntropies2,double *tetraloopEnthalpies2,int numTriloops,int numTetraloops,double atpS[],double atpH[])
{
        int i,threshold;
        double TH;
        char *list[8],rev1[26],rev2[26];

//prepare
    if(GC==1||GC==4)
        threshold=49;
    else
        threshold=45;
        list[0]=F3;
        list[1]=F2;
        list[2]=LF;
        list[3]=F1c;
        list[4]=B1c;
    list[5]=LB;
    list[6]=B2;
    list[7]=B3;

    if(list[2]!=NULL)
    {
        reverse(list[2],rev1,strlen(list[2]));
            for(i=0;i<=1;i++)
        {
            TH=thal(list[i],list[2],stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,1);
            if(TH>threshold)
                return 0;
            TH=thal(list[i],list[2],stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,2);
            if(TH>threshold)
                                return 0;
                        TH=thal(list[i],list[2],stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,3);
            if(TH>threshold)
                                return 0;

            reverse(list[i],rev2,strlen(list[i]));
            TH=thal(rev1,rev2,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,2);
            if(TH>threshold)
                                return 0;
            TH=thal(rev1,rev2,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,3);
            if(TH>44+threshold)
                                return 0;
        }

        for(i=3;i<8;i++)
        {
            if(i==5)
                continue;
            TH=thal(list[2],list[i],stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,1);
            if(TH>threshold)
                                return 0;
            TH=thal(list[2],list[i],stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,2);
            if(TH>threshold)
                                return 0;
                        TH=thal(list[2],list[i],stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,3);
            if(TH>threshold)
                                return 0;

            reverse(list[i],rev2,strlen(list[i]));
            TH=thal(rev2,rev1,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,2);
            if(TH>threshold) 
                                return 0;
            TH=thal(rev2,rev1,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,3);
            if(TH>threshold)
                                return 0;
        }
        }
    if(list[5]!=NULL)
    {
        reverse(list[5],rev1,strlen(list[5]));
        for(i=0;i<5;i++)
        {
            if(i==2)
                continue;
            TH=thal(list[i],list[5],stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,1);
            if(TH>threshold)
                                return 0;
            TH=thal(list[i],list[5],stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,2);
            if(TH>threshold)
                                return 0;
                        TH=thal(list[i],list[5],stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,3);
            if(TH>threshold)
                                return 0;

            reverse(list[i],rev2,strlen(list[i]));
            TH=thal(rev1,rev2,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,2);
            if(TH>threshold)
                                return 0;
            TH=thal(rev1,rev2,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,3);
            if(TH>threshold)
                return 0;
        }
        for(i=6;i<8;i++)
        {
            TH=thal(list[5],list[i],stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,1);
            if(TH>threshold)
                                return 0;
            TH=thal(list[5],list[i],stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,2);
            if(TH>threshold)
                                return 0;
            TH=thal(list[5],list[i],stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,3);
            if(TH>threshold)
                                return 0;

            reverse(list[i],rev2,strlen(list[i]));
            TH=thal(rev2,rev1,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,2);
            if(TH>threshold)
                                return 0;
            TH=thal(rev2,rev1,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,3);
            if(TH>threshold)
                                return 0;
        }
    }
    if(list[2]!=NULL&&list[5]!=NULL)
    {
        TH=thal(list[2],list[5],stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,1);
        if(TH>threshold)
            return 0;
        TH=thal(list[2],list[5],stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,2);
        if(TH>threshold)
                        return 0;
                TH=thal(list[2],list[5],stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,3);
        if(TH>threshold)
                        return 0;

        reverse(list[2],rev1,strlen(list[2]));
        reverse(list[5],rev2,strlen(list[5]));
        TH=thal(rev2,rev1,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,2);
        if(TH>threshold)
                        return 0;
        TH=thal(rev2,rev1,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,3);
        if(TH>threshold)
                        return 0;
    }
    return 1;
}

int design_loop(struct Primer *p_F3,struct Primer *p_F2,struct Primer *p_LF,struct Primer *p_F1c,struct Primer *p_B1c,struct Primer *p_LB,struct Primer *p_B2,struct Primer *p_B3,struct Primer *result_Loop[], int *result,int flag[],int common_num,char *seq,double stackEntropies[],double stackEnthalpies[],double stackint2Entropies[],double stackint2Enthalpies[],double dangleEntropies3[],double dangleEnthalpies3[],double dangleEntropies5[],double dangleEnthalpies5[],double hairpinLoopEntropies[],double interiorLoopEntropies[],double bulgeLoopEntropies[],double hairpinLoopEnthalpies[],double interiorLoopEnthalpies[],double bulgeLoopEnthalpies[],double tstackEntropies[],double tstackEnthalpies[],double tstack2Entropies[],double tstack2Enthalpies[],char *triloopEntropies1,char *triloopEnthalpies1,char *tetraloopEntropies1,char *tetraloopEnthalpies1,double *triloopEntropies2,double *triloopEnthalpies2,double *tetraloopEntropies2,double *tetraloopEnthalpies2,int numTriloops,int numTetraloops,double atpS[],double atpH[],int GC)
{
    int success;
    struct Primer *LF,*LB;
    struct Node *c_LF,*c_LB,*s_LF,*s_LB;
    char primer_F3[26],primer_F2[26],primer_F1c[26],primer_B1c[26],primer_B2[26],primer_B3[26],primer_LF[26],primer_LB[26];

    if(flag[7])
    {
        generate_primer(seq,primer_F3,p_F3->pos,p_F3->len,0);
        generate_primer(seq,primer_F2,p_F2->pos,p_F2->len,0);
        generate_primer(seq,primer_F1c,p_F1c->pos,p_F1c->len,1);
        generate_primer(seq,primer_B1c,p_B1c->pos,p_B1c->len,0);
        generate_primer(seq,primer_B2,p_B2->pos,p_B2->len,1);
        generate_primer(seq,primer_B3,p_B3->pos,p_B3->len,1);
    }
//LF and LB 
    success=0;
    LF=p_LF;
    while(LF)
    {
        if(LF->pos+LF->len>p_F1c->pos)
            break;
        if(LF->minus==0)
        {
            LF=LF->next;
            continue;
        }
        if((LF->minus&GC)!=GC)
        {
            LF=LF->next;
            continue;
        }

        LB=p_LB;
        if(flag[7])
            generate_primer(seq,primer_LF,LF->pos,LF->len,1);
        while(LB)
        {
            if(LB->pos+LB->len>p_B2->pos)
                break;
            if(LB->plus==0)
            {
                LB=LB->next;
                continue;
            }
            if((LB->plus&GC)!=GC)
            {
                LB=LB->next;
                continue;
            }
        //check_structure
            if(flag[7])
            {
                generate_primer(seq,primer_LB,LB->pos,LB->len,0);
                success=check_structure_loop(primer_F3,primer_F2,primer_F1c,primer_B1c,primer_B2,primer_B3,primer_LF,primer_LB,GC,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH);
                if(success==0)
                            {
                                    LB=LB->next;
                                    continue;
                            }
            }
            result_Loop[0]=LF;
            result_Loop[1]=LB;
            success=1;
            break;
        }
        if(success==1)
            break;
        else
            LF=LF->next;
    }
    if(success==1)
        return success;
//only LF
    LF=p_LF;
    result_Loop[1]=NULL;
        while(LF)
        {
                if(LF->pos+LF->len>p_F1c->pos)
                        break;
        if(LF->minus==0)
        {
            LF=LF->next;
            continue;
        }
        if((LF->minus&GC)!=GC)
        {
            LF=LF->next;
            continue;
        }
    //check_structure
        if(flag[7])
        {
            generate_primer(seq,primer_LF,LF->pos,LF->len,1);
            success=check_structure_loop(primer_F3,primer_F2,primer_F1c,primer_B1c,primer_B2,primer_B3,primer_LF,NULL,GC,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH);
            if(success==0)
            {
                LF=LF->next;
                continue;
            }
        }
        result_Loop[0]=LF;
        success=1;
        break;
        }
        if(success==1)
                return success;
    
    //only LB
    LB=p_LB;
    result_Loop[0]=NULL;
        while(LB)
        {
                if(LB->pos+LB->len>p_B2->pos)
                        break;
        if(LB->plus==0)
        {
            LB=LB->next;
            continue;
        }
        if((LB->plus&GC)!=GC)
        {
            LB=LB->next;
            continue;
        }
    //check_structure
        if(flag[7])
        {
            generate_primer(seq,primer_LB,LB->pos,LB->len,0);
            success=check_structure_loop(primer_F3,primer_F2,primer_F1c,primer_B1c,primer_B2,primer_B3,NULL,primer_LB,GC,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH);
            if(success==0)
            {
                LB=LB->next;
                continue;
            }
        }
        result_Loop[1]=LB;
        success=1;
        break;
        }
    return success;
}           
        
int check_add(int F3_pos,int *best_par,int expect)
{
    int i,dis;

    for(i=0;i<expect;i++)
    {
        if(best_par[i]==-1) //the empty record
            return 1;
        dis=best_par[i]-F3_pos;
        if(abs(dis)<150)
            return 0;
    }
    return 1;
}

int check_gc(char *seq,int start,int end)//p_F3->pos,(p_B3->pos+p_B3->len))
{
    int i,total=0;
    float gc;

    for(i=start;i<end;i++)
    {
        if(seq[i]=='C'||seq[i]=='G'||seq[i]=='c'||seq[i]=='g')
            total++;
    }
    gc=total*100.0/(end-start);
    if(gc>=60)
        return 1;
    if(gc<=45)
        return 2;
    return 4;
}
    
