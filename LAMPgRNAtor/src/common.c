#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <unistd.h>
#include <time.h>
#include <sys/stat.h>
#include <ctype.h>

void reverse(char seq[],char rev[],int length)
{
        int i;
        
        for(i=0;i<length;i++)
        {
                if(seq[length-1-i]=='A')
                {
                        rev[i]='T';
                        continue;
                }
                if(seq[length-1-i]=='T')
                {
                        rev[i]='A';
                        continue;
                }
                if(seq[length-1-i]=='C')
                {
                        rev[i]='G';
                        continue;
                }
                rev[i]='C';
        }
        rev[i]='\0';
}

char str2int(char c)
{
        switch (c)
        {
                case 'A': case 'a': case '0':
                        return 0;
                case 'C': case 'c': case '1':
                        return 1;
                case 'G': case 'g': case '2':
                        return 2;
                case 'T': case 't': case '3':
                        return 3;
        }
        return 4;
}

void readLoop(FILE *file,double *v1,double *v2,double *v3)
{
        char *line,*p,*q;
        
        line=(char *)malloc(200);
        memset(line,'\0',200);
        fgets(line,200,file);

        p = line;
        while (isspace(*p)) 
                p++;
        while (isdigit(*p)) 
                p++;
        while (isspace(*p)) 
                p++;

        q = p;
        while (!isspace(*q)) 
                q++;
        *q = '\0';
        q++;
        if (!strcmp(p, "inf"))
                *v1 =1.0*INFINITY;
        else 
                sscanf(p, "%lf", v1);
        while (isspace(*q))
                q++;

        p = q;
        while (!isspace(*p))
                p++;
        *p = '\0';
        p++;
        if (!strcmp(q, "inf"))
                *v2 =1.0*INFINITY;
        else 
                sscanf(q, "%lf", v2);
        while (isspace(*p))
                p++;

        q = p;
        while (!isspace(*q) && (*q != '\0'))
                q++;
        *q = '\0';
        if (!strcmp(p, "inf"))
                *v3 =1.0*INFINITY;
        else 
                sscanf(p, "%lf", v3);
}

void getStack(double stackEntropies[],double stackEnthalpies[], char *path)
{
        int i, j, ii, jj;
        FILE *sFile, *hFile;
        char *line;

        i=strlen(path)+20;
        line=(char *)malloc(i);
        memset(line,'\0',i);
        strcpy(line,path);
        strcat(line,"stack.ds");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        sFile=fopen(line,"r");
        if(sFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }

        memset(line,'\0',i);
        strcpy(line,path);
        strcat(line,"stack.dh");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        hFile=fopen(line,"r");
        if(hFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }
        free(line);

        line=(char *)malloc(20);
        memset(line,'\0',20);
        for (i = 0; i < 5; ++i)
        {
                for (ii = 0; ii < 5; ++ii)
                {
                        for (j = 0; j < 5; ++j)
                        {
                                for (jj = 0; jj < 5; ++jj)
                                {
                                        if (i == 4 || j == 4 || ii == 4 || jj == 4) //N 
                                        {
                                                stackEntropies[i*125+ii*25+j*5+jj] = -1.0;
                                                stackEnthalpies[i*125+ii*25+j*5+jj]=1.0*INFINITY;
                                        }
                                        else 
                                        {
                                                if(fgets(line,20,sFile)==NULL)
                                                {
                                                        printf("Error! When read parameters in getStack function!\n");
                                                        exit(1);
                                                }
                                                if(strncmp(line, "inf", 3)==0)
                                                        stackEntropies[i*125+ii*25+j*5+jj]=1.0*INFINITY;
                                                else
                                                        stackEntropies[i*125+ii*25+j*5+jj] = atof(line);

                                                if(fgets(line,20,hFile)==NULL)
                                                {
                                                        printf("Error! When read parameters in getStack function!\n");
                                                        exit(1);
                                                }
                                                if(strncmp(line, "inf", 3)==0)
                                                        stackEnthalpies[i*125+ii*25+j*5+jj]=1.0*INFINITY;
                                                else
                                                        stackEnthalpies[i*125+ii*25+j*5+jj] = atof(line);

                                                if (fabs(stackEntropies[i*125+ii*25+j*5+jj])>999999999 ||fabs(stackEnthalpies[i*125+ii*25+j*5+jj])>999999999) 
                                                {
                                                        stackEntropies[i*125+ii*25+j*5+jj] = -1.0;
                                                        stackEnthalpies[i*125+ii*25+j*5+jj] =1.0*INFINITY;
                                                }
                                        }
                                }
                        }
                }
        }
        fclose(sFile);
        fclose(hFile);
        free(line);
}

void getStackint2(double stackint2Entropies[],double stackint2Enthalpies[], char *path)
{
        int i, j, ii, jj;
        FILE *sFile, *hFile;
        char *line;

        i=strlen(path)+20;
        line=(char *)malloc(i);
        memset(line,'\0',i);
        strcpy(line,path);
        strcat(line,"stackmm.ds");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        sFile=fopen(line,"r");
        if(sFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }

        memset(line,'\0',i);
        strcpy(line,path);
        strcat(line,"stackmm.dh");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        hFile=fopen(line,"r");
        if(hFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }
        free(line);

        line=(char *)malloc(20);
        memset(line,'\0',20);
        for (i = 0; i < 5; ++i)
        {
                for (ii = 0; ii < 5; ++ii)
                {
                        for (j = 0; j < 5; ++j)
                        {
                                for (jj = 0; jj < 5; ++jj)
                                {
                                        if (i == 4 || j == 4 || ii == 4 || jj == 4)
                                        {
                                                stackint2Entropies[i*125+ii*25+j*5+jj] = -1.0;
                                                stackint2Enthalpies[i*125+ii*25+j*5+jj] =1.0*INFINITY;
                                        } 
                                        else 
                                        {
                                                if(fgets(line,20,sFile)==NULL)
                                                {
                                                        printf("Error! When read parameters in getStackint2 function!\n");
                                                        exit(1);
                                                }
                                                if(strncmp(line, "inf", 3)==0)
                                                        stackint2Entropies[i*125+ii*25+j*5+jj]=1.0*INFINITY;
                                                else
                                                        stackint2Entropies[i*125+ii*25+j*5+jj] = atof(line);

                                                if(fgets(line,20,hFile)==NULL)
                                                {
                                                        printf("Error! When read parameters in getStackint2 function!\n");
                                                        exit(1);
                                                }
                                                if(strncmp(line, "inf", 3)==0)
                                                        stackint2Enthalpies[i*125+ii*25+j*5+jj]=1.0*INFINITY;
                                                else
                                                        stackint2Enthalpies[i*125+ii*25+j*5+jj] = atof(line);

                                                if(fabs(stackint2Entropies[i*125+ii*25+j*5+jj])>999999999||fabs(stackint2Enthalpies[i*125+ii*25+j*5+jj])>999999999)
                                                {
                                                        stackint2Entropies[i*125+ii*25+j*5+jj] = -1.0;
                                                        stackint2Enthalpies[i*125+ii*25+j*5+jj] =1.0*INFINITY;
                                                }
                                        }
                                }
                        }
                }
        }
        fclose(sFile);
        fclose(hFile);
        free(line);
}

void getDangle(double dangleEntropies3[],double dangleEnthalpies3[],double dangleEntropies5[],double dangleEnthalpies5[],char *path)
{
        int i, j, k;
        FILE *sFile, *hFile;
        char *line;
        
        i=strlen(path)+20;
        line=(char *)malloc(i);
        memset(line,'\0',i);
        strcpy(line,path);
        strcat(line,"dangle.ds");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        sFile=fopen(line,"r");
        if(sFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }

        memset(line,'\0',i);
        strcpy(line,path);
        strcat(line,"dangle.dh");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        hFile=fopen(line,"r");
        if(hFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }
        free(line);

        line=(char *)malloc(20);
        memset(line,'\0',20);
        for (i = 0; i < 5; ++i)
                for (j = 0; j < 5; ++j)
                        for (k = 0; k < 5; ++k) 
                        {
                                if (i == 4 || j == 4) 
                                {
                                        dangleEntropies3[i*25+k*5+j] = -1.0;
                                        dangleEnthalpies3[i*25+k*5+j] =1.0*INFINITY;
                                }
                                else if (k == 4)
                                {
                                        dangleEntropies3[i*25+k*5+j] = -1.0;
                                        dangleEnthalpies3[i*25+k*5+j] =1.0*INFINITY;
                                } 
                                else
                                {
                                        if(fgets(line,20,sFile)==NULL)
                                        {
                                                printf("Error! When read parameters in getDangle function!\n");
                                                exit(1);
                                        }
                                        if(strncmp(line, "inf", 3)==0)
                                                dangleEntropies3[i*25+k*5+j]=1.0*INFINITY;
                                        else
                                                dangleEntropies3[i*25+k*5+j]=atof(line);

                                        if(fgets(line,20,hFile)==NULL)
                                        {
                                                printf("Error! When read parameters in getDangle function!\n");        
                                                exit(1);        
                                        }
                                        if(strncmp(line, "inf", 3)==0)        
                                                dangleEnthalpies3[i*25+k*5+j]=1.0*INFINITY;           
                                        else        
                                                dangleEnthalpies3[i*25+k*5+j]=atof(line);

                                        if(fabs(dangleEntropies3[i*25+k*5+j])>999999999||fabs(dangleEnthalpies3[i*25+k*5+j])>999999999) 
                                        {
                                                dangleEntropies3[i*25+k*5+j] = -1.0;
                                                dangleEnthalpies3[i*25+k*5+j] =1.0*INFINITY;
                                        }
                                }
                        }

        for (i = 0; i < 5; ++i)
                for (j = 0; j < 5; ++j)
                        for (k = 0; k < 5; ++k) 
                        {
                                if (i == 4 || j == 4)
                                {
                                        dangleEntropies5[i*25+j*5+k] = -1.0;
                                        dangleEnthalpies5[i*25+j*5+k] =1.0*INFINITY;
                                } 
                                else if (k == 4) 
                                {
                                        dangleEntropies5[i*25+j*5+k] = -1.0;
                                        dangleEnthalpies5[i*25+j*5+k] =1.0*INFINITY;
                                }
                                else
                                {
                                        if(fgets(line,20,sFile)==NULL)
                                        {
                                                printf("Error! When read parameters in getDangle function!\n");
                                                exit(1);
                                        }
                                        if(strncmp(line, "inf", 3)==0)
                                                dangleEntropies5[i*25+j*5+k]=1.0*INFINITY;
                                        else
                                                dangleEntropies5[i*25+j*5+k]=atof(line);

                                        if(fgets(line,20,hFile)==NULL)
                                        {
                                                printf("Error! When read parameters in getDangle function!\n");        
                                                exit(1);        
                                        }
                                        if(strncmp(line, "inf", 3)==0)        
                                                dangleEnthalpies5[i*25+j*5+k]=1.0*INFINITY;           
                                        else        
                                                dangleEnthalpies5[i*25+j*5+k]=atof(line);

                                        if(fabs(dangleEntropies5[i*25+j*5+k])>999999999||fabs(dangleEnthalpies5[i*25+j*5+k])>999999999)
                                        {
                                                dangleEntropies5[i*25+j*5+k] = -1.0;
                                                dangleEnthalpies5[i*25+j*5+k] =1.0*INFINITY;
                                        }
                                }
                        }
        fclose(sFile);
        fclose(hFile);
        free(line);
}

void getLoop(double hairpinLoopEntropies[30],double interiorLoopEntropies[30],double bulgeLoopEntropies[30],double hairpinLoopEnthalpies[30],double interiorLoopEnthalpies[30],double bulgeLoopEnthalpies[30],char *path)
{
        int k;
        FILE *sFile, *hFile;
        char *line;

        k=strlen(path)+20;
        line=(char *)malloc(k);
        memset(line,'\0',k);
        strcpy(line,path);
        strcat(line,"loops.ds");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        sFile=fopen(line,"r");
        if(sFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }

        memset(line,'\0',k);
        strcpy(line,path);
        strcat(line,"loops.dh");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        hFile=fopen(line,"r");
        if(hFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }
        free(line);

        for (k = 0; k < 30; ++k)
        {
                readLoop(sFile, &interiorLoopEntropies[k], &bulgeLoopEntropies[k], &hairpinLoopEntropies[k]);
                readLoop(hFile, &interiorLoopEnthalpies[k], &bulgeLoopEnthalpies[k], &hairpinLoopEnthalpies[k]);
        }
        fclose(sFile);
        fclose(hFile);
}

void getTstack(double tstackEntropies[],double tstackEnthalpies[],char *path)
{
        int i1, j1, i2, j2;
        FILE *sFile, *hFile;
        char *line;

        i1=strlen(path)+20;
        line=(char *)malloc(i1);
        memset(line,'\0',i1);
        strcpy(line,path);
        strcat(line,"tstack_tm_inf.ds");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        sFile=fopen(line,"r");
        if(sFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }

        memset(line,'\0',i1);
        strcpy(line,path);      
        strcat(line,"tstack.dh");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }             
        hFile=fopen(line,"r");
        if(sFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);   
        }
        free(line);

        line=(char *)malloc(20);
        memset(line,'\0',20);
        for (i1 = 0; i1 < 5; ++i1)
                for (i2 = 0; i2 < 5; ++i2)
                        for (j1 = 0; j1 < 5; ++j1)
                                for (j2 = 0; j2 < 5; ++j2)
                                        if (i1 == 4 || j1 == 4)
                                        {
                                                tstackEnthalpies[i1*125+i2*25+j1*5+j2]=1.0*INFINITY;
                                                tstackEntropies[i1*125+i2*25+j1*5+j2] = -1.0;
                                        }
                                        else if (i2 == 4 || j2 == 4)
                                        {
                                                tstackEntropies[i1*125+i2*25+j1*5+j2] = 0.00000000001;
                                                tstackEnthalpies[i1*125+i2*25+j1*5+j2] = 0.0;
                                        }
                                        else
                                        {
                                                if(fgets(line,20,sFile)==NULL)
                                                {
                                                        printf("Error! When read parameters in getTstack function!\n");
                                                        exit(1);
                                                }
                                                if(strncmp(line, "inf", 3)==0)
                                                        tstackEntropies[i1*125+i2*25+j1*5+j2]=1.0*INFINITY;
                                                else
                                                        tstackEntropies[i1*125+i2*25+j1*5+j2]=atof(line);

                                                if(fgets(line,20,hFile)==NULL)
                                                {
                                                        printf("Error! When read parameters in getTstack function!\n");
                                                        exit(1);
                                                }
                                                if(strncmp(line, "inf", 3)==0)
                                                        tstackEnthalpies[i1*125+i2*25+j1*5+j2]=1.0*INFINITY;
                                                else
                                                        tstackEnthalpies[i1*125+i2*25+j1*5+j2]=atof(line);

                                                if (fabs(tstackEntropies[i1*125+i2*25+j1*5+j2])>999999999||fabs(tstackEnthalpies[i1*125+i2*25+j1*5+j2])>999999999)
                                                {
                                                        tstackEntropies[i1*125+i2*25+j1*5+j2] = -1.0;
                                                        tstackEnthalpies[i1*125+i2*25+j1*5+j2] =1.0*INFINITY;
                                                }
                                        }
        fclose(sFile);
        fclose(hFile);
        free(line);
}

void getTstack2(double tstack2Entropies[],double tstack2Enthalpies[],char *path)
{
        int i1, j1, i2, j2;
        FILE *sFile, *hFile;
        char *line;

        i1=strlen(path)+20;
        line=(char *)malloc(i1);
        memset(line,'\0',i1);
        strcpy(line,path);
        strcat(line,"tstack2.ds");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        sFile=fopen(line,"r");
        if(sFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }

        memset(line,'\0',i1);
        strcpy(line,path);      
        strcat(line,"tstack2.dh");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }             
        hFile=fopen(line,"r");
        if(sFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);   
        }
        free(line);

        line=(char *)malloc(20);
        memset(line,'\0',20);
        for (i1 = 0; i1 < 5; ++i1)
                for (i2 = 0; i2 < 5; ++i2)
                        for (j1 = 0; j1 < 5; ++j1)
                                for (j2 = 0; j2 < 5; ++j2)
                                        if (i1 == 4 || j1 == 4)
                                        {
                                                tstack2Enthalpies[i1*125+i2*25+j1*5+j2] =1.0*INFINITY;
                                                tstack2Entropies[i1*125+i2*25+j1*5+j2] = -1.0;
                                        }
                                        else if (i2 == 4 || j2 == 4)
                                        {
                                                tstack2Entropies[i1*125+i2*25+j1*5+j2] = 0.00000000001;
                                                tstack2Enthalpies[i1*125+i2*25+j1*5+j2] = 0.0;
                                        }
                                        else
                                        {
                                                if(fgets(line,20,sFile)==NULL)
                                                {
                                                        printf("Error! When read parameters in getTstack2 function!\n");
                                                        exit(1);
                                                }
                                                if(strncmp(line, "inf", 3)==0)
                                                        tstack2Entropies[i1*125+i2*25+j1*5+j2]=1.0*INFINITY;
                                                else
                                                        tstack2Entropies[i1*125+i2*25+j1*5+j2]=atof(line);

                                                if(fgets(line,20,hFile)==NULL)
                                                {
                                                        printf("Error! When read parameters in getTstack2 function!\n");
                                                        exit(1);
                                                }
                                                if(strncmp(line, "inf", 3)==0)
                                                        tstack2Enthalpies[i1*125+i2*25+j1*5+j2]=1.0*INFINITY;
                                                else
                                                        tstack2Enthalpies[i1*125+i2*25+j1*5+j2]=atof(line);


                                                if (fabs(tstack2Entropies[i1*125+i2*25+j1*5+j2])>999999999||fabs(tstack2Enthalpies[i1*125+i2*25+j1*5+j2])>999999999)
                                                {
                                                        tstack2Entropies[i1*125+i2*25+j1*5+j2] = -1.0;
                                                        tstack2Enthalpies[i1*125+i2*25+j1*5+j2] =1.0*INFINITY;
                                                }
                                        }
        fclose(sFile);
        fclose(hFile);
        free(line);
}

int get_num_line(char *path,int flag)
{
    FILE *fp;
    int i,size;
    char *line;

    i=strlen(path)+20;
        line=(char *)malloc(i);
        memset(line,'\0',i);
        strcpy(line,path);
    if(flag==0)
            strcat(line,"triloop.ds");
    else
        strcat(line,"tetraloop.ds");

        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        fp=fopen(line,"r");
        if(fp==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }

    size=0;
    while(fgets(line,i,fp)!=NULL)
        size++;
    return size;
}

void getTriloop(char *triloopEntropies1,char *triloopEnthalpies1,double *triloopEntropies2,double *triloopEnthalpies2,char *path)
{
        FILE *sFile, *hFile;
        int i,turn;
        char *line,seq[10],value[10];
        
        i=strlen(path)+20;
        line=(char *)malloc(i);
        memset(line,'\0',i);
        strcpy(line,path);
        strcat(line,"triloop.ds");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        sFile=fopen(line,"r");
        if(sFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }
    
    turn=0;
        while(fscanf(sFile,"%s\t%s\n",seq,value)!=EOF)
        {
        for (i=0;i<5;i++)
            triloopEntropies1[5*turn+i]=str2int(seq[i]);
        if(value[0]=='i')
            triloopEntropies2[turn]=1.0*INFINITY;
        else
            triloopEntropies2[turn]=atof(value);
        turn++;
        }
        fclose(sFile);

    i=strlen(path)+20;
        memset(line,'\0',i);
        strcpy(line,path);
        strcat(line,"triloop.dh");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        hFile=fopen(line,"r");
        if(hFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }
        free(line);

    turn=0;
        while(fscanf(hFile,"%s\t%s\n",seq,value)!=EOF)
        {
        for(i=0;i<5;i++)
            triloopEnthalpies1[turn*5+i]=str2int(seq[i]);
        if(value[0]=='i')
            triloopEnthalpies2[turn]=1.0*INFINITY;
        else
            triloopEnthalpies2[turn]=atof(value);
        turn++;
        }
        fclose(hFile);
}

void getTetraloop(char *tetraloopEntropies1,char *tetraloopEnthalpies1,double *tetraloopEntropies2,double *tetraloopEnthalpies2,char *path)
{
        FILE *sFile, *hFile;
        int i, turn;
        char *line,seq[10],value[10];

        i=strlen(path)+20;
        line=(char *)malloc(i);
        memset(line,'\0',i);
        strcpy(line,path);
        strcat(line,"tetraloop.ds");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        sFile=fopen(line,"r");
        if(sFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }

    turn=0;
        while(fscanf(sFile,"%s\t%s\n",seq,value)!=EOF)
        {
        for(i=0;i<6;i++)
            tetraloopEntropies1[turn*6+i]=str2int(seq[i]);
        if(value[0]=='i')
            tetraloopEntropies2[turn]=1.0*INFINITY;
        else
            tetraloopEntropies2[turn]=atof(value);
        turn++;
        }
        fclose(sFile);

        memset(line,'\0',i);
        strcpy(line,path);
        strcat(line,"tetraloop.dh");
        if(access(line,0)==-1)
        {
                printf("Error! Don't have %s file!\n",line);
                exit(1);
        }
        hFile=fopen(line,"r");
        if(hFile==NULL)
        {
                printf("Error! Can't open the %s file!\n",line);
                exit(1);
        }
        free(line);
        
    turn=0;
        while(fscanf(hFile,"%s\t%s\n",seq,value)!=EOF)
        {
        for(i=0;i<6;i++)
            tetraloopEnthalpies1[6*turn+i]=str2int(seq[i]);
        if(value[0]=='i')
            tetraloopEnthalpies2[turn]=1.0*INFINITY;
        else
            tetraloopEnthalpies2[turn]=atof(value);
        turn++;
        }
        fclose(hFile);
}

void tableStartATS(double atp_value, double atpS[])
{
        int i, j;

        for (i = 0; i < 5; ++i)
                for (j = 0; j < 5; ++j)
                        atpS[i*5+j] = 0.00000000001;
        atpS[3] = atpS[15] = atp_value;
}

void tableStartATH(double atp_value,double atpH[])
{
        int i, j;

        for (i = 0; i < 5; ++i)
                for (j = 0; j < 5; ++j)
                        atpH[i*5+j] = 0.0;
        atpH[3] = atpH[15] = atp_value;
}

void initMatrix2(int Initint[],double enthalpyDPT[],double entropyDPT[],char numSeq1[])
{
    int i,j;
    for(i=1;i<=Initint[0];++i)
        for(j=i;j<=Initint[1];++j)
            if(j-i<4 || (numSeq1[i]+numSeq1[j]!=3))
            {
                enthalpyDPT[(i-1)*Initint[2]+j-1]=1.0*INFINITY;
                entropyDPT[(i-1)*Initint[2]+j-1]=-1.0;
            }
            else
            {
                enthalpyDPT[(i-1)*Initint[2]+j-1]=0.0;
                entropyDPT[(i-1)*Initint[2]+j-1]=-3224.0;
            }
}

double Ss(int i,int j,int k,double stackEntropies[],int Initint[],char numSeq1[],char numSeq2[])
{
    if(k==2)
    {
        if(i>=j)
            return -1.0;
        if(i==Initint[0]||j==Initint[1]+1)
            return -1.0;

        if(i>Initint[0])
            i-=Initint[0];
        if(j>Initint[1])
            j-=Initint[1];
        return stackEntropies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j-1]];
    }
    else
        return stackEntropies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]];
}

double Hs(int i,int j,int k,double stackEnthalpies[],int Initint[],char numSeq1[],char numSeq2[])
{
    if(k==2)
    {
        if(i>= j)
            return 1.0*INFINITY;
        if(i==Initint[0]||j==Initint[1]+1)
            return 1.0*INFINITY;

        if(i>Initint[0])
            i-=Initint[0];
        if(j>Initint[1])
            j-=Initint[1];
        if(fabs(stackEnthalpies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j-1]])<999999999)
            return stackEnthalpies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j-1]];
        else
            return 1.0*INFINITY;
    }
    else
        return stackEnthalpies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]];
}

void maxTM2(int i,int j,double stackEntropies[],double stackEnthalpies[],double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],char numSeq1[],char numSeq2[])
{
    double T0,T1,S0,S1,H0,H1;

    S0=entropyDPT[(i-1)*Initint[2]+j-1];
    H0=enthalpyDPT[(i-1)*Initint[2]+j-1];
    T0=(H0+Initdouble[0])/(S0+Initdouble[1]+Initdouble[2]);
    if(fabs(enthalpyDPT[(i-1)*Initint[2]+j-1])<999999999)
    {
        S1=(entropyDPT[i*Initint[2]+j-2]+Ss(i,j,2,stackEntropies,Initint,numSeq1,numSeq2));
        H1=(enthalpyDPT[i*Initint[2]+j-2]+Hs(i,j,2,stackEnthalpies,Initint,numSeq1,numSeq2));
    }
    else
    {
        S1=-1.0;
        H1=1.0*INFINITY;
    }
    T1=(H1+Initdouble[0])/(S1+Initdouble[1]+Initdouble[2]);
    if(S1<-2500.0)
    {
        S1=-3224.0;
        H1=0.0;
    }
    if(S0<-2500.0)
    {
        S0=-3224.0;
        H0=0.0;
    }

    if(T1>T0)
    {
        entropyDPT[(i-1)*Initint[2]+j-1]=S1;
        enthalpyDPT[(i-1)*Initint[2]+j-1]= H1;
    }
    else
    {
        entropyDPT[(i-1)*Initint[2]+j-1]=S0;
        enthalpyDPT[(i-1)*Initint[2]+j-1]=H0;
    }
}

void calc_bulge_internal2(int i,int j,int ii,int jj,double *EntropyEnthalpy,int traceback,double stackEntropies[],double stackEnthalpies[],double stackint2Entropies[],double stackint2Enthalpies[],double interiorLoopEntropies[],double bulgeLoopEntropies[],double interiorLoopEnthalpies[],double bulgeLoopEnthalpies[],double tstackEntropies[],double tstackEnthalpies[],double atpS[],double atpH[],double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],char numSeq1[],char numSeq2[])
{
    int loopSize1,loopSize2,loopSize;
    double T1,T2,S,H;

    S=-3224.0;
    H=0.0;
    loopSize1=ii-i-1;
    loopSize2=j-jj-1;
    if(loopSize1+loopSize2>30)
    {
        EntropyEnthalpy[0]=-1.0;
        EntropyEnthalpy[1]=1.0*INFINITY;
        return;
    }

    loopSize=loopSize1+loopSize2-1;
    if((loopSize1==0&&loopSize2>0)||(loopSize2==0&&loopSize1>0))
    {
        if(loopSize2==1||loopSize1==1)
        { 
            if((loopSize2==1&&loopSize1==0)||(loopSize2==0&&loopSize1==1))
            {
                H=bulgeLoopEnthalpies[loopSize]+stackEnthalpies[numSeq1[i]*125+numSeq1[ii]*25+numSeq2[j]*5+numSeq2[jj]];
                S=bulgeLoopEntropies[loopSize]+stackEntropies[numSeq1[i]*125+numSeq1[ii]*25+numSeq2[j]*5+numSeq2[jj]];
            }
            if(traceback!=1)
            {
                H+=enthalpyDPT[(ii-1)*Initint[2]+jj-1];
                S+=entropyDPT[(ii-1)*Initint[2]+jj-1];
            }

            if(fabs(H)>999999999)
            {
                H=1.0*INFINITY;
                S=-1.0;
            }
            T1=(H+Initdouble[0])/((S+Initdouble[1])+Initdouble[2]);
            T2=(enthalpyDPT[(i-1)*Initint[2]+j-1]+Initdouble[0])/((entropyDPT[(i-1)*Initint[2]+j-1])+Initdouble[1]+Initdouble[2]);
            if((T1>T2)||((traceback&&T1>=T2)||traceback==1))
            {
                EntropyEnthalpy[0]=S;
                EntropyEnthalpy[1]=H;
            }
        }
        else
        {
            H=bulgeLoopEnthalpies[loopSize]+atpH[numSeq1[i]*5+numSeq2[j]]+atpH[numSeq1[ii]*5+numSeq2[jj]];
            if(traceback!=1)
                H+=enthalpyDPT[(ii-1)*Initint[2]+jj-1];

            S=bulgeLoopEntropies[loopSize]+atpS[numSeq1[i]*5+numSeq2[j]]+atpS[numSeq1[ii]*5+numSeq2[jj]];
            if(traceback!=1)
                S+=entropyDPT[(ii-1)*Initint[2]+jj-1];
            if(fabs(H)>999999999)
            {
                H=1.0*INFINITY;
                S=-1.0;
            }
            T1=(H+Initdouble[0])/((S+Initdouble[1])+Initdouble[2]);
            T2=(enthalpyDPT[(i-1)*Initint[2]+j-1]+Initdouble[0])/(entropyDPT[(i-1)*Initint[2]+j-1]+Initdouble[1]+Initdouble[2]);
            if((T1>T2)||((traceback&&T1>=T2)||(traceback==1)))
            {
                EntropyEnthalpy[0]=S;
                EntropyEnthalpy[1]=H;
            }
        }
    }
    else if(loopSize1==1&&loopSize2==1)
    {
        S=stackint2Entropies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j-1]]+stackint2Entropies[numSeq2[jj]*125+numSeq2[jj+1]*25+numSeq1[ii]*5+numSeq1[ii-1]];
        if(traceback!=1)
            S+=entropyDPT[(ii-1)*Initint[2]+jj-1];

        H=stackint2Enthalpies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j-1]]+stackint2Enthalpies[numSeq2[jj]*125+numSeq2[jj+1]*25+numSeq1[ii]*5+numSeq1[ii-1]];
        if(traceback!=1)
            H+=enthalpyDPT[(ii-1)*Initint[2]+jj-1];
        if(fabs(H)>999999999)
        {
            H=1.0*INFINITY;
            S=-1.0;
        }
        T1=(H+Initdouble[0])/((S+Initdouble[1])+Initdouble[2]);
        T2=(enthalpyDPT[(i-1)*Initint[2]+j-1]+Initdouble[0])/(entropyDPT[(i-1)*Initint[2]+j-1]+Initdouble[1]+Initdouble[2]);
        if((T1-T2>=0.000001)||traceback)
        {
            if((T1>T2)||((traceback&&T1>= T2)||traceback==1))
            {
                EntropyEnthalpy[0]=S;
                EntropyEnthalpy[1]=H;
            }
        }
        return;
    }
    else
    {
        H=interiorLoopEnthalpies[loopSize]+tstackEnthalpies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j-1]]+tstackEnthalpies[numSeq2[jj]*125+numSeq2[jj+1]*25+numSeq1[ii]*5+numSeq1[ii-1]];
        if(traceback!=1)
            H+=enthalpyDPT[(ii-1)*Initint[2]+jj-1];

        S=interiorLoopEntropies[loopSize]+tstackEntropies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j-1]]+tstackEntropies[numSeq2[jj]*125+numSeq2[jj+1]*25+numSeq1[ii]*5+numSeq1[ii-1]]+(-300/310.15*abs(loopSize1-loopSize2));
        if(traceback!=1)
            S+=entropyDPT[(ii-1)*Initint[2]+jj-1];
        if(fabs(H)>999999999)
        {
            H=1.0*INFINITY;
            S=-1.0;
        }

        T1=(H+Initdouble[0])/((S+Initdouble[1])+Initdouble[2]);
        T2=(enthalpyDPT[(i-1)*Initint[2]+j-1]+Initdouble[0])/((entropyDPT[(i-1)*Initint[2]+j-1])+Initdouble[1]+Initdouble[2]);
        if((T1>T2)||((traceback&&T1>=T2)||(traceback==1)))
        {
            EntropyEnthalpy[0]=S;
            EntropyEnthalpy[1]=H;
        }
    }
    return;
}

void CBI(int i,int j,double* EntropyEnthalpy,int traceback,double stackEntropies[],double stackEnthalpies[],double stackint2Entropies[],double stackint2Enthalpies[],double interiorLoopEntropies[],double bulgeLoopEntropies[],double interiorLoopEnthalpies[],double bulgeLoopEnthalpies[],double tstackEntropies[],double tstackEnthalpies[],double atpS[],double atpH[],double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],char numSeq1[],char numSeq2[])
{
    int d,ii,jj;

    for(d=j-i-3;d>=4&&d>=j-i-32;--d)
        for(ii=i+1;ii<j-d&&ii<=Initint[0];++ii)
        {
            jj=d+ii;
            if(traceback==0)
            {
                EntropyEnthalpy[0]=-1.0;
                EntropyEnthalpy[1]=1.0*INFINITY;
            }
            if(fabs(enthalpyDPT[(ii-1)*Initint[2]+jj-1])<999999999)
            {
                calc_bulge_internal2(i,j,ii,jj,EntropyEnthalpy,traceback,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,interiorLoopEntropies,bulgeLoopEntropies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);
                if(fabs(EntropyEnthalpy[1])<999999999)
                {
                    if(EntropyEnthalpy[0] <-2500.0)
                    {
                        EntropyEnthalpy[0]=-3224.0;
                        EntropyEnthalpy[1]=0.0;
                    }
                    if(traceback==0)
                    {
                        enthalpyDPT[(i-1)*Initint[2]+j-1]=EntropyEnthalpy[1];
                        entropyDPT[(i-1)*Initint[2]+j-1]=EntropyEnthalpy[0];
                    }
                }
            }
        }
    return;
}

int find_pos(char *ref,int ref_start,char *source,int length,int num)
{
    int flag,i,j;

    for(i=0;i<num;i++)
    {
        flag=0;
        for(j=0;j<length;j++)
        {
            if(ref[ref_start+j]!=source[i*length+j])
            {
                flag++;
                break;
            }
        }
        if(flag==0)
            return i;
    }
    return -1;
}

void calc_hairpin(int i,int j,double *EntropyEnthalpy,int traceback,double hairpinLoopEntropies[],double hairpinLoopEnthalpies[],double tstack2Entropies[],double tstack2Enthalpies[],char *triloopEntropies1,char *triloopEnthalpies1,char *tetraloopEntropies1,char *tetraloopEnthalpies1,double *triloopEntropies2,double *triloopEnthalpies2,double *tetraloopEntropies2,double *tetraloopEnthalpies2,int numTriloops,int numTetraloops,double atpS[],double atpH[],double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],char numSeq1[])
{
    int pos,loopSize=j-i-1;
    double T1,T2;
    
    if(loopSize < 3)
    {
        EntropyEnthalpy[0]=-1.0;
        EntropyEnthalpy[1]=1.0*INFINITY;
        return;
    }
    if(i<=Initint[0]&&Initint[1]<j)
    {
        EntropyEnthalpy[0]=-1.0;
        EntropyEnthalpy[1]=1.0*INFINITY;
        return;
    }
    else if(i>Initint[1])
    {
        i-= Initint[0];
        j-= Initint[1];
    }
    if(loopSize<=30)
    {
        EntropyEnthalpy[1]=hairpinLoopEnthalpies[loopSize-1];
        EntropyEnthalpy[0]=hairpinLoopEntropies[loopSize-1];
    }
    else
    {
        EntropyEnthalpy[1]=hairpinLoopEnthalpies[29];
        EntropyEnthalpy[0]=hairpinLoopEntropies[29];
    }

    if(loopSize>3) // for loops 4 bp and more in length, terminal mm are accounted
    {
        EntropyEnthalpy[1]+=tstack2Enthalpies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq1[j]*5+numSeq1[j-1]];
        EntropyEnthalpy[0]+=tstack2Entropies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq1[j]*5+numSeq1[j-1]];
    }
    else if(loopSize == 3) // for loops 3 bp in length at-penalty is considered
    {
        EntropyEnthalpy[1]+=atpH[numSeq1[i]*5+numSeq1[j]];
        EntropyEnthalpy[0]+=atpS[numSeq1[i]*5+numSeq1[j]];
    }

    if(loopSize==3) // closing AT-penalty (+), triloop bonus, hairpin of 3 (+) 
    {
        pos=find_pos(numSeq1,i,triloopEnthalpies1,5,numTriloops);
        if(pos!=-1)
            EntropyEnthalpy[1]+=triloopEnthalpies2[pos];

        pos=find_pos(numSeq1,i,triloopEntropies1,5,numTriloops);
        if(pos!=-1)
            EntropyEnthalpy[0]+=triloopEntropies2[pos];
    }
    else if (loopSize == 4) // terminal mismatch, tetraloop bonus, hairpin of 4
    {
        pos=find_pos(numSeq1,i,tetraloopEnthalpies1,6,numTetraloops);
        if(pos!=-1)
            EntropyEnthalpy[1]+=tetraloopEnthalpies2[pos];

        pos=find_pos(numSeq1,i,tetraloopEntropies1,6,numTetraloops);
        if(pos!=-1)
            EntropyEnthalpy[0]+=tetraloopEntropies2[pos];
    }
    if(fabs(EntropyEnthalpy[1])>999999999)
    {
        EntropyEnthalpy[1] =1.0*INFINITY;
        EntropyEnthalpy[0] = -1.0;
    }
    T1 = (EntropyEnthalpy[1] +Initdouble[0]) / ((EntropyEnthalpy[0] +Initdouble[1]+ Initdouble[2]));
    T2 = (enthalpyDPT[(i-1)*Initint[2]+j-1] +Initdouble[0]) / ((entropyDPT[(i-1)*Initint[2]+j-1]) +Initdouble[1]+ Initdouble[2]);
    if(T1 < T2 && traceback == 0) // ERROR
    {
        EntropyEnthalpy[0] =entropyDPT[(i-1)*Initint[2]+j-1];
        EntropyEnthalpy[1] =enthalpyDPT[(i-1)*Initint[2]+j-1];
    }
    return;
}

void fillMatrix2(double stackEntropies[],double stackEnthalpies[],double stackint2Entropies[],double stackint2Enthalpies[],double hairpinLoopEntropies[],double interiorLoopEntropies[],double bulgeLoopEntropies[],double hairpinLoopEnthalpies[],double interiorLoopEnthalpies[],double bulgeLoopEnthalpies[],double tstackEntropies[],double tstackEnthalpies[],double tstack2Entropies[],double tstack2Enthalpies[],char *triloopEntropies1,char *triloopEnthalpies1,char *tetraloopEntropies1,char *tetraloopEnthalpies1,double *triloopEntropies2,double *triloopEnthalpies2,double *tetraloopEntropies2,double *tetraloopEnthalpies2,int numTriloops,int numTetraloops,double atpS[],double atpH[],double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],char numSeq1[],char numSeq2[])
{
    int i, j;
    double SH[2];

    for (j = 2; j <= Initint[1]; ++j)
        for (i = j - 3 - 1; i >= 1; --i)
        {
            if (fabs(enthalpyDPT[(i-1)*Initint[2]+j-1])<999999999)
            {
                SH[0] = -1.0;
                SH[1] =1.0*INFINITY;
                maxTM2(i,j,stackEntropies,stackEnthalpies,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);
                CBI(i,j,SH,0,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,interiorLoopEntropies,bulgeLoopEntropies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);

                SH[0] = -1.0;
                SH[1] =1.0*INFINITY;
                calc_hairpin(i, j, SH, 0,hairpinLoopEntropies,hairpinLoopEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1);
                if(fabs(SH[1])<999999999)
                {
                    if(SH[0] <-2500.0) /* to not give dH any value if dS is unreasonable */
                    {
                        SH[0] =-3224.0;
                        SH[1] = 0.0;
                    }
                    entropyDPT[(i-1)*Initint[2]+j-1]= SH[0];
                    enthalpyDPT[(i-1)*Initint[2]+j-1]= SH[1];
                }
            }
        }
}

int max5(double a,double b,double c,double d,double e)
{
    if(a>b&&a>c&&a>d&&a>e)
        return 1;
    else if(b>c&&b>d&&b>e)
        return 2;
    else if(c>d&&c>e)
        return 3;
    else if(d>e)
        return 4;
    else
        return 5;
}

double Sd5(int i,int j,double dangleEntropies5[],char numSeq1[])
{
    return dangleEntropies5[numSeq1[i]*25+numSeq1[j]*5+numSeq1[j-1]];
}

double Hd5(int i,int j,double dangleEnthalpies5[],char numSeq1[])
{
    return dangleEnthalpies5[numSeq1[i]*25+numSeq1[j]*5+numSeq1[j-1]];
}

double Sd3(int i,int j,double dangleEntropies3[],char numSeq1[])
{
    return dangleEntropies3[numSeq1[i]*25+numSeq1[i+1]*5+numSeq1[j]];
}

double Hd3(int i,int j,double dangleEnthalpies3[],char numSeq1[])
{
    return dangleEnthalpies3[numSeq1[i]*25+numSeq1[i+1]*5+numSeq1[j]];
}

double Ststack(int i,int j,double tstack2Entropies[],char numSeq1[])
{
    return tstack2Entropies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq1[j]*5+numSeq1[j-1]];
}

double Htstack(int i,int j,double tstack2Enthalpies[],char numSeq1[])
{
    return tstack2Enthalpies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq1[j]*5+numSeq1[j-1]];
}

double END5_1(int i,int hs,double atpS[],double atpH[],double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],double send5[],double hend5[],char numSeq1[])
{
    int k;
    double max_tm,T1,T2,H,S,H_max,S_max;

    max_tm=-1.0*INFINITY;
    H_max=1.0*INFINITY;
    S_max=-1.0;
    for(k=0;k<=i-5;++k)
    {
        T1=(hend5[k]+Initdouble[0])/(send5[k]+Initdouble[1]+Initdouble[2]);
        T2=Initdouble[0]/(Initdouble[1]+Initdouble[2]);
        if(T1>=T2)
        {
            H=hend5[k]+atpH[numSeq1[k+1]*5+numSeq1[i]]+enthalpyDPT[k*Initint[2]+i-1];
            S=send5[k]+atpS[numSeq1[k+1]*5+numSeq1[i]]+entropyDPT[k*Initint[2]+i-1];
            if(fabs(H)>999999999||H>0||S>0)  // H and S must be greater than 0 to avoid BS
            {
                H=1.0*INFINITY;
                S=-1.0;
            }
            T1=(H+Initdouble[0])/(S+Initdouble[1]+Initdouble[2]);
        }
        else
        {
            H=atpH[numSeq1[k+1]*5+numSeq1[i]]+enthalpyDPT[k*Initint[2]+i-1];
            S=atpS[numSeq1[k+1]*5+numSeq1[i]]+entropyDPT[k*Initint[2]+i-1];
            if(fabs(H)>999999999||H>0||S>0)
            {
                H=1.0*INFINITY;
                S=-1.0;
            }
            T1=(H+Initdouble[0])/(S+Initdouble[1]+Initdouble[2]);
        }

        if(max_tm<T1)
        {
            if(S>-2500.0)
            {
                H_max=H;
                S_max=S;
                max_tm=T1;
            }
        }
    }
    if(hs==1)
        return H_max;
    return S_max;
}

double END5_2(int i,int hs,double dangleEntropies5[],double dangleEnthalpies5[],double atpS[],double atpH[],double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],double send5[],double hend5[],char numSeq1[])
{
    int k;
    double max_tm,T1,T2,H,S,H_max,S_max;

    H_max=1.0*INFINITY;
    max_tm=-1.0*INFINITY;
    S_max=-1.0;
    for(k=0;k<=i-6;++k)
    {
        T1=(hend5[k]+Initdouble[0])/(send5[k]+Initdouble[1]+Initdouble[2]);
        T2=Initdouble[0]/(Initdouble[1]+Initdouble[2]);
        if(T1>=T2)
        {
            H=hend5[k]+atpH[numSeq1[k+2]*5+numSeq1[i]]+Hd5(i,k+2,dangleEnthalpies5,numSeq1)+enthalpyDPT[(k+1)*Initint[2]+i-1];
            S=send5[k]+atpS[numSeq1[k+2]*5+numSeq1[i]]+Sd5(i,k+2,dangleEntropies5,numSeq1)+entropyDPT[(k+1)*Initint[2]+i-1];
            if(fabs(H)>999999999||H>0||S>0)
            {
                H=1.0*INFINITY;
                S=-1.0;
            }
            T1=(H+Initdouble[0])/(S+Initdouble[1]+Initdouble[2]);
        }
        else
        {
            H=atpH[numSeq1[k+2]*5+numSeq1[i]]+Hd5(i,k+2,dangleEnthalpies5,numSeq1)+enthalpyDPT[(k+1)*Initint[2]+i-1];
            S=atpS[numSeq1[k+2]*5+numSeq1[i]]+Sd5(i,k+2,dangleEntropies5,numSeq1)+entropyDPT[(k+1)*Initint[2]+i-1];
            if(fabs(H)>999999999||H>0||S>0)
            {
                H=1.0*INFINITY;
                S=-1.0;
            }
            T1=(H+Initdouble[0])/(S+Initdouble[1]+Initdouble[2]);
        }

        if(max_tm<T1)
        {
            if(S>-2500.0)
            {
                H_max=H;
                S_max=S;
                max_tm=T1;
            }
        }
    }
    if(hs==1)
        return H_max;
    return S_max;
}

double END5_3(int i,int hs,double dangleEntropies3[],double dangleEnthalpies3[],double atpS[],double atpH[],double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],double send5[],double hend5[],char numSeq1[])
{
    int k;
    double max_tm,T1,T2,H,S,H_max,S_max;

    H_max=1.0*INFINITY;
    max_tm=-1.0*INFINITY;
    S_max=-1.0;
    for(k=0;k<=i-6;++k)
    {
        T1=(hend5[k]+Initdouble[0])/(send5[k]+Initdouble[1]+Initdouble[2]);
        T2=Initdouble[0]/(Initdouble[1]+Initdouble[2]);
        if(T1>=T2)
        {
            H=hend5[k]+atpH[numSeq1[k+1]*5+numSeq1[i-1]]+Hd3(i-1,k+1,dangleEnthalpies3,numSeq1)+enthalpyDPT[k*Initint[2]+i-2];
            S=send5[k]+atpS[numSeq1[k+1]*5+numSeq1[i-1]]+Sd3(i-1,k+1,dangleEntropies3,numSeq1)+entropyDPT[k*Initint[2]+i-2];
            if(fabs(H)>999999999||H>0||S>0)
            {
                H=1.0*INFINITY;
                S=-1.0;
            }
            T1=(H+Initdouble[0])/(S+Initdouble[1]+Initdouble[2]);
        }
        else
        {
            H=atpH[numSeq1[k+1]*5+numSeq1[i-1]]+Hd3(i-1,k+1,dangleEnthalpies3,numSeq1)+enthalpyDPT[k*Initint[2]+i-2];
            S=atpS[numSeq1[k+1]*5+numSeq1[i-1]]+Sd3(i-1,k+1,dangleEntropies3,numSeq1)+entropyDPT[k*Initint[2]+i-2];
            if(fabs(H)>999999999||H>0||S>0)
            {
                H=1.0*INFINITY;
                S=-1.0;
            }
            T1=(H+Initdouble[0])/(S+Initdouble[1]+Initdouble[2]);
        }

        if(max_tm<T1)
        {
            if(S>-2500.0)
            {
                H_max=H;
                S_max=S;
                max_tm=T1;
            }
        }
    }
    if(hs==1)
        return H_max;
    return S_max;
}

double END5_4(int i,int hs,double tstack2Entropies[],double tstack2Enthalpies[],double atpS[],double atpH[],double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],double send5[],double hend5[],char numSeq1[])
{
    int k;
    double max_tm,T1,T2,H,S,H_max,S_max;

    H_max=1.0*INFINITY;
    max_tm=-1.0*INFINITY;
    S_max=-1.0;
    for(k=0;k<=i-7;++k)
    {
        T1=(hend5[k]+Initdouble[0])/(send5[k]+Initdouble[1]+Initdouble[2]);
        T2=Initdouble[0]/(Initdouble[1]+Initdouble[2]);
        if(T1>=T2)
        {
            H=hend5[k]+atpH[numSeq1[k+2]*5+numSeq1[i-1]]+Htstack(i-1,k+2,tstack2Enthalpies,numSeq1)+enthalpyDPT[(k+1)*Initint[2]+i-2];
            S=send5[k]+atpS[numSeq1[k+2]*5+numSeq1[i-1]]+Ststack(i-1,k+2,tstack2Entropies,numSeq1)+entropyDPT[(k+1)*Initint[2]+i-2];
            if(fabs(H)>999999999||H>0||S>0)
            {
                H=1.0*INFINITY;
                S=-1.0;
            }
            T1=(H+Initdouble[0])/(S+Initdouble[1]+Initdouble[2]);
        }
        else
        {
            H=atpH[numSeq1[k+2]*5+numSeq1[i-1]]+Htstack(i-1,k+2,tstack2Enthalpies,numSeq1)+enthalpyDPT[(k+1)*Initint[2]+i-2];
            S=atpS[numSeq1[k+2]*5+numSeq1[i-1]]+Ststack(i-1,k+2,tstack2Entropies,numSeq1)+entropyDPT[(k+1)*Initint[2]+i-2];
            if(fabs(H)>999999999||H>0||S>0)
            {
                H=1.0*INFINITY;
                S=-1.0;
            }
            T1=(H+Initdouble[0])/(S+Initdouble[1]+Initdouble[2]);
        }

        if(max_tm<T1)
        {
            if(S>-2500.0)
            {
                H_max=H;
                S_max=S;
                max_tm=T1;
            }
        }
    }
    if(hs==1)
        return H_max;
    return S_max;
}

void calc_terminal_bp(double temp,double dangleEntropies3[],double dangleEnthalpies3[],double dangleEntropies5[],double dangleEnthalpies5[],double tstack2Entropies[],double tstack2Enthalpies[],double atpS[],double atpH[],double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],double send5[],double hend5[],char numSeq1[])
{
    int i,max;
    double T1,T2,T3,T4,T5,G,end5_11,end5_12,end5_21,end5_22,end5_31,end5_32,end5_41,end5_42;
    
    send5[0]=send5[1]= -1.0;
    hend5[0]=hend5[1]=1.0*INFINITY;

    for(i=2;i<=Initint[0];i++)
    {
        send5[i]=-3224.0;
        hend5[i]=0;
    }

// adding terminal penalties to 3' end and to 5' end 
    for(i=2;i<=Initint[0];++i)
    {
        max=0;
        T1=(hend5[i-1]+Initdouble[0])/(send5[i-1]+Initdouble[1]+Initdouble[2]);
        end5_11=END5_1(i,1,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1);
        end5_12=END5_1(i,2,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1);
        T2=(end5_11+Initdouble[0])/(end5_12+Initdouble[1]+Initdouble[2]);
        end5_21=END5_2(i,1,dangleEntropies5,dangleEnthalpies5,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1);
        end5_22=END5_2(i,2,dangleEntropies5,dangleEnthalpies5,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1);
        T3=(end5_21+Initdouble[0])/(end5_22+Initdouble[1]+Initdouble[2]);
        end5_31=END5_3(i,1,dangleEntropies3,dangleEnthalpies3,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1);
        end5_32=END5_3(i,2,dangleEntropies3,dangleEnthalpies3,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1);
        T4=(end5_31+Initdouble[0])/(end5_32+Initdouble[1]+Initdouble[2]);
        end5_41=END5_4(i,1,tstack2Entropies,tstack2Enthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1);
        end5_42=END5_4(i,2,tstack2Entropies,tstack2Enthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1);
        T5=(end5_41+Initdouble[0])/(end5_42+Initdouble[1]+Initdouble[2]);

        max=max5(T1,T2,T3,T4,T5);
        switch(max)
        {
            case 1:
                send5[i]=send5[i-1];
                hend5[i]=hend5[i-1];
                break;
            case 2:
                G=end5_11-temp*end5_12;
                if(G<0.0)
                {
                    send5[i]=end5_12;
                    hend5[i]=end5_11;
                }
                else
                {
                    send5[i]=send5[i-1];
                    hend5[i]=hend5[i-1];
                }
                break;
            case 3:
                G=end5_21-temp*end5_22;
                if(G<0.0)
                {
                    send5[i]=end5_22;
                    hend5[i]=end5_21;
                }
                else
                {
                    send5[i]=send5[i-1];
                    hend5[i]=hend5[i-1];
                }
                break;
            case 4:
                G=end5_31-temp*end5_32;
                if(G<0.0)
                {
                    send5[i]=end5_32;
                    hend5[i]=end5_31;
                }
                else
                {
                    send5[i]=send5[i-1];
                    hend5[i]=hend5[i-1];
                }
                break;
            case 5:
                G=end5_41-temp*end5_42;
                if(G<0.0)
                {
                    send5[i]=end5_42;
                    hend5[i]=end5_41;
                }
                else
                {
                    send5[i]=send5[i-1];
                    hend5[i]=hend5[i-1];
                }
                break;
            default:
                break;
        }
    }
}

int newpush(int store[],int i,int j,int mtrx,int total,int next)
{
        int k;
        for(k=total-1;k>=next;k--)
        {
                store[(k+1)*3]=store[k*3];
                store[(k+1)*3+1]=store[k*3+1];
                store[(k+1)*3+2]=store[k*3+2];
        }
        store[next*3]=i;                  
        store[next*3+1]=j;
        store[next*3+2]=mtrx;

        return total+1;           
}

int equal(double a,double b)
{
    if(fabs(a)>999999999||fabs(b)>999999999)
        return 0;
    return fabs(a-b)<1e-5;
}

void tracebacku(int bp[],double stackEntropies[],double stackEnthalpies[],double stackint2Entropies[],double stackint2Enthalpies[],double dangleEntropies3[],double dangleEnthalpies3[],double dangleEntropies5[],double dangleEnthalpies5[],double hairpinLoopEntropies[],double interiorLoopEntropies[],double bulgeLoopEntropies[],double hairpinLoopEnthalpies[],double interiorLoopEnthalpies[],double bulgeLoopEnthalpies[],double tstackEntropies[],double tstackEnthalpies[],double tstack2Entropies[],double tstack2Enthalpies[],char *triloopEntropies1,char *triloopEnthalpies1,char *tetraloopEntropies1,char *tetraloopEnthalpies1,double *triloopEntropies2,double *triloopEnthalpies2,double *tetraloopEntropies2,double *tetraloopEnthalpies2,int numTriloops,int numTetraloops,double atpS[],double atpH[],double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],double send5[],double hend5[],char numSeq1[],char numSeq2[])
{
    int i,j,store[50],total,now,ii,jj,k,d,done;
    double SH1[2],SH2[2],EntropyEnthalpy[2];

        total=newpush(store,Initint[0],0,1,0,0);
        now=0;
        while(now<total)
        {
                i=store[3*now]; // top->i;
                j=store[3*now+1]; // top->j;
                if(store[now*3+2]==1)
                {
                        while(equal(send5[i],send5[i-1])&&equal(hend5[i],hend5[i-1])) // if previous structure is the same as this one
                                --i;
                        if(i==0)
                                continue;
                        if(equal(send5[i],END5_1(i,2,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1))&&equal(hend5[i],END5_1(i,1,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1)))
                        {
                                for(k=0;k<=i-5;++k)
                                        if(equal(send5[i],atpS[numSeq1[k+1]*5+numSeq1[i]]+entropyDPT[k*Initint[2]+i-1])&&equal(hend5[i],atpH[numSeq1[k+1]*5+numSeq1[i]]+enthalpyDPT[k*Initint[2]+i-1]))
                                        {
                                                total=newpush(store,k+1,i,0,total,now+1);                    
                                                break;
                                        }
                                        else if(equal(send5[i],send5[k]+atpS[numSeq1[k+1]*5+numSeq1[i]]+entropyDPT[k*Initint[2]+i-1])&&equal(hend5[i],hend5[k]+atpH[numSeq1[k+1]*5+numSeq1[i]]+enthalpyDPT[k*Initint[2]+i-1]))
                                        {
                                                total=newpush(store,k+1,i,0,total,now+1);
                                                total=newpush(store,k,0,1,total,now+1);
                                                break;
                                        }
                        }
                        else if(equal(send5[i],END5_2(i,2,dangleEntropies5,dangleEnthalpies5,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1))&&equal(hend5[i],END5_2(i,1,dangleEntropies5,dangleEnthalpies5,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1)))
                        {
                                for (k=0;k<=i-6;++k)
                                        if(equal(send5[i],atpS[numSeq1[k+2]*5+numSeq1[i]]+Sd5(i,k+2,dangleEntropies5,numSeq1)+entropyDPT[(k+1)*Initint[2]+i-1])&&equal(hend5[i],atpH[numSeq1[k+2]*5+numSeq1[i]]+Hd5(i,k+2,dangleEnthalpies5,numSeq1)+enthalpyDPT[(k+1)*Initint[2]+i-1]))
                                        {
                                                total=newpush(store,k+2,i,0,total,now+1);
                                                break;
                                        }
                                        else if(equal(send5[i],send5[k]+atpS[numSeq1[k+2]*5+numSeq1[i]]+Sd5(i,k+2,dangleEntropies5,numSeq1)+entropyDPT[(k+1)*Initint[2]+i-1])&&equal(hend5[i],hend5[k]+atpH[numSeq1[k+2]*5+numSeq1[i]]+Hd5(i,k+2,dangleEnthalpies5,numSeq1)+enthalpyDPT[(k+1)*Initint[2]+i-1]))
                                        {
                                                total=newpush(store,k+2,i,0,total,now+1);
                                                total=newpush(store,k,0,1,total,now+1);
                                                break;
                                        }
                        }
                        else if(equal(send5[i],END5_3(i,2,dangleEntropies3,dangleEnthalpies3,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1))&&equal(hend5[i],END5_3(i,1,dangleEntropies3,dangleEnthalpies3,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1)))
                        {
                                for (k=0;k<=i-6;++k)
                                        if(equal(send5[i],atpS[numSeq1[k+1]*5+numSeq1[i-1]]+Sd3(i-1,k+1,dangleEntropies3,numSeq1)+entropyDPT[k*Initint[2]+i-2])&&equal(hend5[i],atpH[numSeq1[k+1]*5+numSeq1[i-1]]+Hd3(i-1,k+1,dangleEnthalpies3,numSeq1)+enthalpyDPT[k*Initint[2]+i-2]))
                                        {
                                                total=newpush(store,k+1,i-1,0,total,now+1);
                                                break;
                                        }
                                        else if(equal(send5[i],send5[k]+atpS[numSeq1[k+1]*5+numSeq1[i-1]]+Sd3(i-1,k+1,dangleEntropies3,numSeq1)+entropyDPT[k*Initint[2]+i-2])&&equal(hend5[i],hend5[k]+atpH[numSeq1[k+1]*5+numSeq1[i-1]]+Hd3(i-1,k+1,dangleEnthalpies3,numSeq1)+enthalpyDPT[k*Initint[2]+i-2]))
                                        {
                                                total=newpush(store,k+1,i-1,0,total,now+1);
                                                total=newpush(store,k,0,1,total,now+1);
                                                break;
                                        }
                        }
                        else if(equal(send5[i],END5_4(i,2,tstack2Entropies,tstack2Enthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1))&&equal(hend5[i],END5_4(i,1,tstack2Entropies,tstack2Enthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1)))
                        {
                                for (k=0;k<=i-7;++k)
                                        if(equal(send5[i],atpS[numSeq1[k+2]*5+numSeq1[i-1]]+Ststack(i-1,k+2,tstack2Entropies,numSeq1)+entropyDPT[(k+1)*Initint[2]+i-2])&&equal(hend5[i],atpH[numSeq1[k+2]*5+numSeq1[i-1]]+Htstack(i-1,k+2,tstack2Enthalpies,numSeq1)+enthalpyDPT[(k+1)*Initint[2]+i-2]))
                                        {
                                                total=newpush(store,k+2,i-1,0,total,now+1);
                                                break;
                                        }
                                        else if(equal(send5[i],send5[k]+atpS[numSeq1[k+2]*5+numSeq1[i-1]]+Ststack(i-1,k+2,tstack2Entropies,numSeq1)+entropyDPT[(k+1)*Initint[2]+i-2])&&equal(hend5[i],hend5[k]+atpH[numSeq1[k+2]*5+numSeq1[i-1]]+Htstack(i-1,k+2,tstack2Enthalpies,numSeq1)+enthalpyDPT[(k+1)*Initint[2]+i-2]))
                                        {
                                                total=newpush(store,k+2,i-1,0,total,now+1);
                                                total=newpush(store,k,0,1,total,now+1);
                                                break;
                                        }
                        }
                }
                else if(store[3*now+2]==0)
                {
                        bp[i-1]=j;
                        bp[j-1]=i;
                        SH1[0]=-1.0;
                        SH1[1]=1.0*INFINITY;
                        calc_hairpin(i,j,SH1,1,hairpinLoopEntropies,hairpinLoopEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1);

                        SH2[0]=-1.0;
                        SH2[1]=1.0*INFINITY;
                        CBI(i,j,SH2,2,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,interiorLoopEntropies,bulgeLoopEntropies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);

                        if (equal(entropyDPT[(i-1)*Initint[2]+j-1],Ss(i,j,2,stackEntropies,Initint,numSeq1,numSeq2)+entropyDPT[i*Initint[2]+j-2])&&equal(enthalpyDPT[(i-1)*Initint[2]+j-1],Hs(i,j,2,stackEnthalpies,Initint,numSeq1,numSeq2)+enthalpyDPT[i*Initint[2]+j-2]))
                                total=newpush(store,i+1,j-1,0,total,now+1);
                        else if(equal(entropyDPT[(i-1)*Initint[2]+j-1],SH1[0])&&equal(enthalpyDPT[(i-1)*Initint[2]+j-1],SH1[1]));
                        else if(equal(entropyDPT[(i-1)*Initint[2]+j-1],SH2[0])&&equal(enthalpyDPT[(i-1)*Initint[2]+j-1],SH2[1]))
                        {
                                for (done=0,d=j-i-3;d>=4&&d>=j-i-32&&!done;--d)
                                        for (ii=i+1;ii<j-d;++ii)
                                        {
                                                jj=d+ii;
                                                EntropyEnthalpy[0]=-1.0;
                                                EntropyEnthalpy[1]=1.0*INFINITY;
                                                calc_bulge_internal2(i,j,ii,jj,EntropyEnthalpy,1,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,interiorLoopEntropies,bulgeLoopEntropies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);

                                                if (equal(entropyDPT[(i-1)*Initint[2]+j-1],EntropyEnthalpy[0]+entropyDPT[(ii-1)*Initint[2]+jj-1])&&equal(enthalpyDPT[(i-1)*Initint[2]+j-1],EntropyEnthalpy[1]+enthalpyDPT[(ii-1)*Initint[2]+jj-1]))
                                                {
                                                        total=newpush(store,ii,jj,0,total,now+1);
                                                        ++done;
                                                        break;
                                                }
                                        }
                        }
                }
                now++;
        }
}

double drawHairpin(int bp[],double mh,double ms,int Initint[])
{
        int i,N;

        N=0;
        if(fabs(ms)>999999999||fabs(mh)>999999999)
        {
        return 0.0;
        }
        else
        {
        for(i=1;i<Initint[0];++i)
        {
            if(bp[i-1]>0)
                N++;
                }
                return mh/(ms+(((N/2)-1)*-0.51986))-273.15;
        }
}

void initMatrix(int Initint[],double enthalpyDPT[],double entropyDPT[],char numSeq1[],char numSeq2[])
{
    int i,j;

    for(i=1;i<=Initint[0];++i)
    {
        for(j=1;j<=Initint[1];++j)
        {
            if(numSeq1[i]+numSeq2[j]!=3)
            {
                enthalpyDPT[(i-1)*Initint[2]+j-1]=1.0*INFINITY;
                entropyDPT[(i-1)*Initint[2]+j-1]=-1.0;
            }
            else
            {
                enthalpyDPT[(i-1)*Initint[2]+j-1]=0.0;
                entropyDPT[(i-1)*Initint[2]+j-1]=-3224.0;
            }
        }
    }
}

void LSH(int i,int j,double *EntropyEnthalpy,double dangleEntropies3[],double dangleEnthalpies3[],double dangleEntropies5[],double dangleEnthalpies5[],double tstack2Entropies[],double tstack2Enthalpies[],double atpS[],double atpH[],double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],char numSeq1[],char numSeq2[])
{
    double S1,H1,T1,S2,H2,T2;

    if(numSeq1[i]+numSeq2[j]!=3)
    {
        entropyDPT[(i-1)*Initint[2]+j-1]=-1.0;
        enthalpyDPT[(i-1)*Initint[2]+j-1]=1.0*INFINITY;
        return;
    }

    S1=atpS[numSeq1[i]*5+numSeq2[j]]+tstack2Entropies[numSeq2[j]*125+numSeq2[j-1]*25+numSeq1[i]*5+numSeq1[i-1]];
    H1=atpH[numSeq1[i]*5+numSeq2[j]]+tstack2Enthalpies[numSeq2[j]*125+numSeq2[j-1]*25+numSeq1[i]*5+numSeq1[i-1]];
    if(fabs(H1)>999999999)
    {
        H1=1.0*INFINITY;
        S1=-1.0;
    }
// If there is two dangling ends at the same end of duplex
    if(fabs(dangleEnthalpies3[numSeq2[j]*25+numSeq2[j-1]*5+numSeq1[i]])<999999999&&fabs(dangleEnthalpies5[numSeq2[j]*25+numSeq1[i]*5+numSeq1[i-1]])<999999999)
    {
        S2=atpS[numSeq1[i]*5+numSeq2[j]]+dangleEntropies3[numSeq2[j]*25+numSeq2[j-1]*5+numSeq1[i]]+dangleEntropies5[numSeq2[j]*25+numSeq1[i]*5+numSeq1[i-1]];
        H2=atpH[numSeq1[i]*5+numSeq2[j]]+dangleEnthalpies3[numSeq2[j]*25+numSeq2[j-1]*5+numSeq1[i]]+dangleEnthalpies5[numSeq2[j]*25+numSeq1[i]*5+numSeq1[i-1]];
        if(fabs(H2)>999999999)
        {
            H2=1.0*INFINITY;
            S2=-1.0;
        }
        T2=(H2+Initdouble[0])/(S2+Initdouble[1]+Initdouble[2]);
        if(fabs(H1)<999999999)
        {
            T1=(H1+Initdouble[0])/(S1+Initdouble[1]+Initdouble[2]);
            if(T1<T2)
            {
                S1=S2;
                H1=H2;
                T1=T2;
            }
        }
        else
        {
            S1=S2;
            H1=H2;
            T1=T2;
        }
    }
    else if(fabs(dangleEnthalpies3[numSeq2[j]*25+numSeq2[j-1]*5+numSeq1[i]])<999999999)
    {
        S2=atpS[numSeq1[i]*5+numSeq2[j]]+dangleEntropies3[numSeq2[j]*25+numSeq2[j-1]*5+numSeq1[i]];
        H2=atpH[numSeq1[i]*5+numSeq2[j]]+dangleEnthalpies3[numSeq2[j]*25+numSeq2[j-1]*5+numSeq1[i]];
        if(fabs(H2)>999999999)
        {
            H2=1.0*INFINITY;
            S2=-1.0;
        }
        T2=(H2+Initdouble[0])/(S2+Initdouble[1]+Initdouble[2]);
        if(fabs(H1)<999999999)
        {
            T1=(H1+Initdouble[0])/(S1+Initdouble[1]+Initdouble[2]);
            if(T1<T2)
            {
                S1=S2;
                H1=H2;
                T1=T2;
            }
        }
        else
        {
            S1=S2;
            H1=H2;
            T1=T2;
        }
    }
    else if(fabs(dangleEnthalpies5[numSeq2[j]*25+numSeq1[i]*5+numSeq1[i-1]])<999999999)
    {
        S2=atpS[numSeq1[i]*5+numSeq2[j]]+dangleEntropies5[numSeq2[j]*25+numSeq1[i]*5+numSeq1[i-1]];
        H2=atpH[numSeq1[i]*5+numSeq2[j]]+dangleEnthalpies5[numSeq2[j]*25+numSeq1[i]*5+numSeq1[i-1]];
        if(fabs(H2)>999999999)
        {
            H2=1.0*INFINITY;
            S2=-1.0;
        }
        T2=(H2+Initdouble[0])/(S2+Initdouble[1]+Initdouble[2]);
        if(fabs(H1)<999999999)
        {
            T1=(H1+Initdouble[0])/(S1+Initdouble[1]+Initdouble[2]);
            if(T1<T2)
            {
                S1=S2;
                H1=H2;
                T1=T2;
            }
        }
        else
        {
            S1=S2;
            H1=H2;
            T1=T2;
        }
    }

    S2=atpS[numSeq1[i]*5+numSeq2[j]];
    H2=atpH[numSeq1[i]*5+numSeq2[j]];
    T2=(H2+Initdouble[0])/(S2+Initdouble[1]+Initdouble[2]);
    if(fabs(H1)<999999999)
    {
        if(T1<T2)
        {
            EntropyEnthalpy[0]=S2;
            EntropyEnthalpy[1]=H2;
        }
        else
        {
            EntropyEnthalpy[0]=S1;
            EntropyEnthalpy[1]=H1;
        }
    }
    else
    {
        EntropyEnthalpy[0]=S2;
        EntropyEnthalpy[1]=H2;
    }
    return;
}

void maxTM(int i,int j,double stackEntropies[],double stackEnthalpies[],double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],char numSeq1[],char numSeq2[])
{
    double T0,T1,S0,S1,H0,H1;

    S0=entropyDPT[(i-1)*Initint[2]+j-1];
    H0=enthalpyDPT[(i-1)*Initint[2]+j-1];
    T0=(H0+Initdouble[0])/(S0+Initdouble[1]+Initdouble[2]); // at current position 
    if(fabs(enthalpyDPT[(i-2)*Initint[2]+j-2])<999999999&&fabs(Hs(i-1,j-1,1,stackEnthalpies,Initint,numSeq1,numSeq2))<999999999)
    {
        S1=(entropyDPT[(i-2)*Initint[2]+j-2]+Ss(i-1,j-1,1,stackEntropies,Initint,numSeq1,numSeq2));
        H1=(enthalpyDPT[(i-2)*Initint[2]+j-2]+Hs(i-1,j-1,1,stackEnthalpies,Initint,numSeq1,numSeq2));
    }
    else
    {
        S1=-1.0;
        H1=1.0*INFINITY;
    }
    T1=(H1+Initdouble[0])/(S1+Initdouble[1]+Initdouble[2]);

    if(S1<-2500.0)
    {
// to not give dH any value if dS is unreasonable
        S1=-3224.0;
        H1=0.0;
    }
    if(S0<-2500.0)
    {
// to not give dH any value if dS is unreasonable
        S0=-3224.0;
        H0=0.0;
    }
    if((T1>T0)||(S0>0&&H0>0)) // T1 on suurem 
    {
        entropyDPT[(i-1)*Initint[2]+j-1]=S1;
        enthalpyDPT[(i-1)*Initint[2]+j-1]=H1;
    }
    else if(T0>=T1)
    {
        entropyDPT[(i-1)*Initint[2]+j-1]=S0;
        enthalpyDPT[(i-1)*Initint[2]+j-1]=H0;
    }
}

void calc_bulge_internal(int i,int j,int ii,int jj,double* EntropyEnthalpy,int traceback,double stackEntropies[],double stackEnthalpies[],double stackint2Entropies[],double stackint2Enthalpies[],double interiorLoopEntropies[],double bulgeLoopEntropies[],double interiorLoopEnthalpies[],double bulgeLoopEnthalpies[],double tstackEntropies[],double tstackEnthalpies[],double atpS[],double atpH[],double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],char numSeq1[],char numSeq2[])
{
    int loopSize1,loopSize2,loopSize;
    double T1,T2,S,H;

    S=-3224.0;
    H=0;
    loopSize1=ii-i-1;
    loopSize2=jj-j-1;
    loopSize=loopSize1+loopSize2-1;
    if((loopSize1==0&&loopSize2>0)||(loopSize2==0&&loopSize1>0))// only bulges have to be considered
    {
        if(loopSize2==1||loopSize1==1) // bulge loop of size one is treated differently the intervening nn-pair must be added
        {
            if((loopSize2==1&&loopSize1==0)||(loopSize2==0&&loopSize1==1))
            {
                H=bulgeLoopEnthalpies[loopSize]+stackEnthalpies[numSeq1[i]*125+numSeq1[ii]*25+numSeq2[j]*5+numSeq2[jj]];
                S=bulgeLoopEntropies[loopSize]+stackEntropies[numSeq1[i]*125+numSeq1[ii]*25+numSeq2[j]*5+numSeq2[jj]];
            }
            H+=enthalpyDPT[(i-1)*Initint[2]+j-1];
            S+=entropyDPT[(i-1)*Initint[2]+j-1];
            if(fabs(H)>999999999)
            {
                H=1.0*INFINITY;
                S=-1.0;
            }

            T1=(H+Initdouble[0])/((S+Initdouble[1])+Initdouble[2]);
            T2=(enthalpyDPT[(ii-1)*Initint[2]+jj-1]+Initdouble[0])/((entropyDPT[(ii-1)*Initint[2]+jj-1])+Initdouble[1]+Initdouble[2]);
            if((T1>T2)||((traceback&&T1>=T2)||(traceback==1)))
            {
                EntropyEnthalpy[0]=S;
                EntropyEnthalpy[1]=H;
            }
        }
        else // we have _not_ implemented Jacobson-Stockaymayer equation; the maximum bulgeloop size is 30
        {
            H=bulgeLoopEnthalpies[loopSize]+atpH[numSeq1[i]*5+numSeq2[j]]+atpH[numSeq1[ii]*5+numSeq2[jj]];
            H+=enthalpyDPT[(i-1)*Initint[2]+j-1];

            S=bulgeLoopEntropies[loopSize]+atpS[numSeq1[i]*5+numSeq2[j]]+atpS[numSeq1[ii]*5+numSeq2[jj]];
            S+=entropyDPT[(i-1)*Initint[2]+j-1];
            if(fabs(H)>999999999)
            {
                H=1.0*INFINITY;
                S=-1.0;
            }
            T1=(H+Initdouble[0])/((S+Initdouble[1])+Initdouble[2]);
            T2=(enthalpyDPT[(ii-1)*Initint[2]+jj-1]+Initdouble[0])/(entropyDPT[(ii-1)*Initint[2]+jj-1]+Initdouble[1]+Initdouble[2]);
            if((T1>T2)||((traceback&&T1>=T2)||(traceback==1)))
            {
                EntropyEnthalpy[0]=S;
                EntropyEnthalpy[1]=H;
            }
        }
    }
    else if(loopSize1==1&&loopSize2==1)
    {
        S=stackint2Entropies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]]+stackint2Entropies[numSeq2[jj]*125+numSeq2[jj-1]*25+numSeq1[ii]*5+numSeq1[ii-1]];
        S+=entropyDPT[(i-1)*Initint[2]+j-1];

        H=stackint2Enthalpies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]]+stackint2Enthalpies[numSeq2[jj]*125+numSeq2[jj-1]*25+numSeq1[ii]*5+numSeq1[ii-1]];
        H+=enthalpyDPT[(i-1)*Initint[2]+j-1];
        if(fabs(H)>999999999)
        {
            H=1.0*INFINITY;
            S=-1.0;
        }
        T1=(H+Initdouble[0])/((S+Initdouble[1])+Initdouble[2]);
        T2=(enthalpyDPT[(ii-1)*Initint[2]+jj-1]+Initdouble[0])/(entropyDPT[(ii-1)*Initint[2]+jj-1]+Initdouble[1]+Initdouble[2]);
        if((T1-T2>=0.000001)||traceback==1)
        {
            if((T1>T2)||(traceback&&T1>=T2))
            {
                EntropyEnthalpy[0]=S;
                EntropyEnthalpy[1]=H;
            }
        }
        return;
    }
    else // only internal loops
    {
        H=interiorLoopEnthalpies[loopSize]+tstackEnthalpies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]]+tstackEnthalpies[numSeq2[jj]*125+numSeq2[jj-1]*25+numSeq1[ii]*5+numSeq1[ii-1]];
        H+=enthalpyDPT[(i-1)*Initint[2]+j-1];

        S=interiorLoopEntropies[loopSize]+tstackEntropies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]]+tstackEntropies[numSeq2[jj]*125+numSeq2[jj-1]*25+numSeq1[ii]*5+numSeq1[ii-1]]+(-300/310.15*abs(loopSize1-loopSize2));
        S+=entropyDPT[(i-1)*Initint[2]+j-1];
        if(fabs(H)>999999999)
        {
            H=1.0*INFINITY;
            S=-1.0;
        }
        T1=(H+Initdouble[0])/((S+Initdouble[1])+Initdouble[2]);
        T2=(enthalpyDPT[(ii-1)*Initint[2]+jj-1]+Initdouble[0])/((entropyDPT[(ii-1)*Initint[2]+jj-1])+Initdouble[1]+Initdouble[2]);
        if((T1>T2)||((traceback&&T1>=T2)||(traceback==1)))
        {
            EntropyEnthalpy[0]=S;
            EntropyEnthalpy[1]=H;
        }
    }
    return;
}

void fillMatrix(double stackEntropies[],double stackEnthalpies[],double stackint2Entropies[],double stackint2Enthalpies[],double dangleEntropies3[],double dangleEnthalpies3[],double dangleEntropies5[],double dangleEnthalpies5[],double interiorLoopEntropies[],double bulgeLoopEntropies[],double interiorLoopEnthalpies[],double bulgeLoopEnthalpies[],double tstackEntropies[],double tstackEnthalpies[],double tstack2Entropies[],double tstack2Enthalpies[],double atpS[],double atpH[],double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],char numSeq1[],char numSeq2[])
{
    int d,i,j,ii,jj;
    double SH[2];

    for(i=1;i<=Initint[0];++i)
    {
        for(j=1;j<=Initint[1];++j)
        {
            if(fabs(enthalpyDPT[(i-1)*Initint[2]+j-1])<999999999)
            {
                SH[0]=-1.0;
                SH[1]=1.0*INFINITY;
                LSH(i,j,SH,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,tstack2Entropies,tstack2Enthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);

                if(fabs(SH[1])<999999999)
                {
                    entropyDPT[(i-1)*Initint[2]+j-1]=SH[0];
                    enthalpyDPT[(i-1)*Initint[2]+j-1]=SH[1];
                }
                if(i>1&&j>1)
                {
                    maxTM(i,j,stackEntropies,stackEnthalpies,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);
                    for(d=3;d<=32;d++)
                    {
                        ii=i-1;
                        jj=-ii-d+(j+i);
                        if(jj<1)
                        {
                            ii-=abs(jj-1);
                            jj=1;
                        }
                        for(;ii>0&&jj<j;--ii,++jj)
                        {
                            if(fabs(enthalpyDPT[(ii-1)*Initint[2]+jj-1])<999999999)
                            {
                                SH[0]=-1.0;
                                SH[1]=1.0*INFINITY;
                                calc_bulge_internal(ii,jj,i,j,SH,0,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,interiorLoopEntropies,bulgeLoopEntropies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);

                                if(SH[0]<-2500.0)
                                {
                                    SH[0] =-3224.0;
                                    SH[1] = 0.0;
                                }
                                if(fabs(SH[1])<999999999)
                                {
                                    enthalpyDPT[(i-1)*Initint[2]+j-1]=SH[1];
                                    entropyDPT[(i-1)*Initint[2]+j-1]=SH[0];
                                }
                            }
                        }
                    }
                } // if 
            }
        } // for 
    } //for
}

void RSH(int i,int j,double EntropyEnthalpy[],double dangleEntropies3[],double dangleEnthalpies3[],double dangleEntropies5[],double dangleEnthalpies5[],double tstack2Entropies[],double tstack2Enthalpies[],double atpS[],double atpH[],double Initdouble[],char numSeq1[],char numSeq2[])
{
    double S1,S2,H1,H2,T1,T2;

    if(numSeq1[i]+numSeq2[j]!=3)
    {
        EntropyEnthalpy[0]=-1.0;
        EntropyEnthalpy[1]=1.0*INFINITY;
        return;
    }
    S1=atpS[numSeq1[i]*5+numSeq2[j]]+tstack2Entropies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]];
    H1=atpH[numSeq1[i]*5+numSeq2[j]]+tstack2Enthalpies[numSeq1[i]*125+numSeq1[i+1]*25+numSeq2[j]*5+numSeq2[j+1]];
    if(fabs(H1)>999999999)
    {
        H1=1.0*INFINITY;
        S1=-1.0;
    }
    if(fabs(dangleEnthalpies3[numSeq1[i]*25+numSeq1[i+1]*5+numSeq2[j]])<999999999&&fabs(dangleEnthalpies5[numSeq1[i]*25+numSeq2[j]*5+numSeq2[j+1]])<999999999)
    {
        S2=atpS[numSeq1[i]*5+numSeq2[j]]+dangleEntropies3[numSeq1[i]*25+numSeq1[i+1]*5+numSeq2[j]]+dangleEntropies5[numSeq1[i]*25+numSeq2[j]*5+numSeq2[j+1]];
        H2=atpH[numSeq1[i]*5+numSeq2[j]]+dangleEnthalpies3[numSeq1[i]*25+numSeq1[i+1]*5+numSeq2[j]]+dangleEnthalpies5[numSeq1[i]*25+numSeq2[j]*5+numSeq2[j+1]];
        if(fabs(H2)>999999999)
        {
            H2=1.0*INFINITY;
            S2=-1.0;
        }
        T2=(H2+Initdouble[0])/(S2+Initdouble[1]+Initdouble[2]);
        if(fabs(H1)<999999999)
        {
            T1=(H1+Initdouble[0])/(S1+Initdouble[1]+Initdouble[2]);
            if(T1<T2)
            {
                S1=S2;
                H1=H2;
                T1=T2;
            }
        }
        else
        {
            S1=S2;
            H1=H2;
            T1=T2;
        }
    }

    if(fabs(dangleEnthalpies3[numSeq1[i]*25+numSeq1[i+1]*5+numSeq2[j]])<999999999)
    {
        S2=atpS[numSeq1[i]*5+numSeq2[j]]+dangleEntropies3[numSeq1[i]*25+numSeq1[i+1]*5+numSeq2[j]];
        H2=atpH[numSeq1[i]*5+numSeq2[j]]+dangleEnthalpies3[numSeq1[i]*25+numSeq1[i+1]*5+numSeq2[j]];
        if(fabs(H2)>999999999)
        {
            H2=1.0*INFINITY;
            S2=-1.0;
        }
        T2=(H2+Initdouble[0])/(S2+Initdouble[1]+Initdouble[2]);
        if(fabs(H1)<999999999)
        {
            T1=(H1+Initdouble[0])/(S1+Initdouble[1]+Initdouble[2]);
            if(T1<T2)
            {
                S1=S2;
                H1=H2;
                T1=T2;
            }
        }
        else
        {
            S1=S2;
            H1=H2;
            T1=T2;
        }
    }

    if(fabs(dangleEnthalpies5[numSeq1[i]*25+numSeq2[j]*5+numSeq2[j+1]])<999999999)
    {
        S2=atpS[numSeq1[i]*5+numSeq2[j]]+dangleEntropies5[numSeq1[i]*25+numSeq2[j]*5+numSeq2[j+1]];
        H2=atpH[numSeq1[i]*5+numSeq2[j]]+dangleEnthalpies5[numSeq1[i]*25+numSeq2[j]*5+numSeq2[j+1]];
        if(fabs(H2)>999999999)
        {
            H2=1.0*INFINITY;
            S2=-1.0;
        }
        T2=(H2+Initdouble[0])/(S2+Initdouble[1]+Initdouble[2]);
        if(fabs(H1)<999999999)
        {
            T1=(H1+Initdouble[0])/(S1+Initdouble[1]+Initdouble[2]);
            if(T1<T2)
            {
                S1=S2;
                H1=H2;
                T1=T2;
            }
        }
        else
        {
            S1=S2;
            H1=H2;
            T1=T2;
        }
    }
    S2=atpS[numSeq1[i]*5+numSeq2[j]];
    H2=atpH[numSeq1[i]*5+numSeq2[j]];
    T2=(H2+Initdouble[0])/(S2+Initdouble[1]+Initdouble[2]);
    if(fabs(H1)<999999999)
    {
        if(T1<T2)
        {
            EntropyEnthalpy[0]=S2;
            EntropyEnthalpy[1]=H2;
        }
        else
        {
            EntropyEnthalpy[0]=S1;
            EntropyEnthalpy[1]=H1;
        }
    }
    else
    {
        EntropyEnthalpy[0]=S2;
        EntropyEnthalpy[1]=H2;
    }
    return;
}

void traceback(int i,int j,int* ps1,int* ps2,double stackEntropies[],double stackEnthalpies[],double stackint2Entropies[],double stackint2Enthalpies[],double dangleEntropies3[],double dangleEnthalpies3[],double dangleEntropies5[],double dangleEnthalpies5[],double interiorLoopEntropies[],double bulgeLoopEntropies[],double interiorLoopEnthalpies[],double bulgeLoopEnthalpies[],double tstackEntropies[],double tstackEnthalpies[],double tstack2Entropies[],double tstack2Enthalpies[],double atpS[],double atpH[],double Initdouble[],int Initint[],double enthalpyDPT[],double entropyDPT[],char numSeq1[],char numSeq2[])
{
    int d,ii,jj,done;
    double SH[2];

    ps1[i-1]=j;
    ps2[j-1]=i;
    while(1)
    {
        SH[0]=-1.0;
        SH[1]=1.0*INFINITY;
        LSH(i,j,SH,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,tstack2Entropies,tstack2Enthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);
        if(equal(entropyDPT[(i-1)*Initint[2]+j-1],SH[0])&&equal(enthalpyDPT[(i-1)*Initint[2]+j-1],SH[1]))
            break;

        done = 0;
        if(i>1&&j>1&&equal(entropyDPT[(i-1)*Initint[2]+j-1],Ss(i-1,j-1,1,stackEntropies,Initint,numSeq1,numSeq2)+entropyDPT[(i-2)*Initint[2]+j-2]))
        {
            i=i-1;
            j=j-1;
            ps1[i-1]=j;
            ps2[j-1]=i;
            done=1;
        }
        for(d=3;!done&&d<=32;++d)
        {
            ii=i-1;
            jj=-ii-d+(j+i);
            if(jj<1)
            {
                ii-=abs(jj-1);
                jj=1;
            }
            for(;!done&&ii>0&&jj<j;--ii,++jj)
            {
                SH[0]=-1.0;
                SH[1]=1.0*INFINITY;
                calc_bulge_internal(ii,jj,i,j,SH,1,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,interiorLoopEntropies,bulgeLoopEntropies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);
                if(equal(entropyDPT[(i-1)*Initint[2]+j-1],SH[0])&&equal(enthalpyDPT[(i-1)*Initint[2]+j-1],SH[1]))
                {
                    i=ii;
                    j=jj;
                    ps1[i-1]=j;
                    ps2[j-1]=i;
                    done=1;
                    break;
                }
            }
        }
    }
}

double drawDimer(int *ps1,int *ps2,double H,double S,double Initdouble[],int Initint[])
{
        int i,N;

        if(fabs(Initdouble[3])>999999999)
                return 0.0;
        else
        {
                N=0;
                for(i=0;i<Initint[0];i++)
                {
                        if(ps1[i]>0)
                                ++N;
                }
                for(i=0;i<Initint[1];i++)
                {
                        if(ps2[i]>0)
                                ++N;
                }
                N=(N/2)-1;
                return (H/(S+(N*-0.51986)+Initdouble[2]))-273.15;
        }
}

int symmetry_thermo(char seq[])
{
        int i = 0;
        int seq_len=strlen(seq);
        int mp = seq_len/2;
        if(seq_len%2==1)
                return 0;          

        while(i<mp) 
        {
                if((seq[i]=='A'&&seq[seq_len-1-i]!='T')||(seq[i]=='T'&&seq[seq_len-1-i]!='A')||(seq[seq_len-1-i]=='A'&&seq[i]!='T')||(seq[seq_len-1-i]=='T'&&seq[i]!='A'))
                        return 0;   
                if((seq[i]=='C'&&seq[seq_len-1-i]!='G')||(seq[i]=='G'&&seq[seq_len-1-i]!='C')||(seq[seq_len-1-i]=='C'&&seq[i]!='G')||(seq[seq_len-1-i]=='G'&&seq[i]!='C'))
                        return 0;
        i++;
        }
        return 1;
}

double thal(char oligo_f[],char oligo_r[],double stackEntropies[],double stackEnthalpies[],double stackint2Entropies[],double stackint2Enthalpies[],double dangleEntropies3[],double dangleEnthalpies3[],double dangleEntropies5[],double dangleEnthalpies5[],double hairpinLoopEntropies[],double interiorLoopEntropies[],double bulgeLoopEntropies[],double hairpinLoopEnthalpies[],double interiorLoopEnthalpies[],double bulgeLoopEnthalpies[],double tstackEntropies[],double tstackEnthalpies[],double tstack2Entropies[],double tstack2Enthalpies[],char *triloopEntropies1,char *triloopEnthalpies1,char *tetraloopEntropies1,char *tetraloopEnthalpies1,double *triloopEntropies2,double *triloopEnthalpies2,double *tetraloopEntropies2,double *tetraloopEnthalpies2,int numTriloops,int numTetraloops,double atpS[],double atpH[],int type)
{
    double SH[2],Initdouble[4];//0 is dplx_init_H, 1 is dplx_init_S, 2 is RC, 3 is SHleft
    int Initint[5]; //0 is len1, 1 is len2, 2 is len3, 3 is bestI, 4 is bestJ
    int i, j;
    double T1,enthalpyDPT[625],entropyDPT[625],send5[26],hend5[26],result_TH;
    int ps1[25],ps2[25];
    char numSeq1[27],numSeq2[27];
    double mh, ms;

/*** INIT values for unimolecular and bimolecular structures ***/
    if (type==4) /* unimolecular folding */
    {
        Initdouble[0]= 0.0;
        Initdouble[1] = -0.00000000001;
        Initdouble[2]=0;
    }
    else /* hybridization of two oligos */
    {
        Initdouble[0]= 200;
        Initdouble[1]= -5.7;
        if(symmetry_thermo(oligo_f) && symmetry_thermo(oligo_r))
            Initdouble[2]=1.9872* log(38/1000000000.0);
        else
            Initdouble[2]=1.9872* log(38/4000000000.0);
    }
    
    /* convert nucleotides to numbers */
    if(type==1 || type==2)
    {
        Initint[0]=strlen(oligo_f);
        Initint[1]=strlen(oligo_r);
        for(i=1;i<=Initint[0];++i)
            numSeq1[i]=str2int(oligo_f[i-1]);
        for(i=1;i<=Initint[1];++i)
            numSeq2[i]=str2int(oligo_r[Initint[1]-i]);
    }
    else if(type==3)
    {
        Initint[0]=strlen(oligo_r);
        Initint[1]=strlen(oligo_f);
        for(i=1;i<=Initint[0];++i)
            numSeq1[i]=str2int(oligo_r[i-1]);
        for(i=1;i<=Initint[1];++i)
            numSeq2[i]=str2int(oligo_f[Initint[1]-i]);
    }
    else
    {
        Initint[0]=strlen(oligo_f);
                Initint[1]=strlen(oligo_r);
        Initint[2]=Initint[1]-1;
                for(i=1;i<=Initint[0];++i)      
                        numSeq1[i]=str2int(oligo_f[i-1]);   
                for(i=1;i<=Initint[1];++i)      
                        numSeq2[i]=str2int(oligo_r[i-1]);
    }
    numSeq1[0]=numSeq1[Initint[0]+1]=numSeq2[0]=numSeq2[Initint[1]+1]=4; /* mark as N-s */

    result_TH=0;
    if (type==4) /* calculate structure of monomer */
    {
        initMatrix2(Initint,enthalpyDPT,entropyDPT,numSeq1);
        fillMatrix2(stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);
        calc_terminal_bp(310.15,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,tstack2Entropies,tstack2Enthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1);
        mh=hend5[Initint[0]];
        ms=send5[Initint[0]];
        for (i=0;i<Initint[0];i++)
            ps1[i]=0;
        if(fabs(mh)<999999999)
        {
            tracebacku(ps1,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,hairpinLoopEntropies,interiorLoopEntropies,bulgeLoopEntropies,hairpinLoopEnthalpies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,triloopEntropies1,triloopEnthalpies1,tetraloopEntropies1,tetraloopEnthalpies1,triloopEntropies2,triloopEnthalpies2,tetraloopEntropies2,tetraloopEnthalpies2,numTriloops,numTetraloops,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,send5,hend5,numSeq1,numSeq2);
            result_TH=drawHairpin(ps1,mh,ms,Initint);
            result_TH=(int)(result_TH*100+0.5)/100.0;
        }
        return result_TH;
    }
    else if(type!=4) /* Hybridization of two moleculs */
    {
        Initint[2]=Initint[1];
        initMatrix(Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);
        fillMatrix(stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,interiorLoopEntropies,bulgeLoopEntropies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);

        Initdouble[3]=-1.0*INFINITY;
    /* calculate terminal basepairs */
        Initint[3]=Initint[4]=0;
        if(type==1)
            for (i=1;i<=Initint[0];i++)
            {
                for (j=1;j<=Initint[1];j++)
                {
                    RSH(i,j,SH,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,tstack2Entropies,tstack2Enthalpies,atpS,atpH,Initdouble,numSeq1,numSeq2);
                    SH[0]=SH[0]+0.000001; /* this adding is done for compiler, optimization -O2 vs -O0 */
                    SH[1]=SH[1]+0.000001;
                    T1=((enthalpyDPT[(i-1)*Initint[2]+j-1]+ SH[1] +Initdouble[0]) / ((entropyDPT[(i-1)*Initint[2]+j-1]) + SH[0] +Initdouble[1] + Initdouble[2])) -273.15;
                    if(T1>Initdouble[3]&&((entropyDPT[(i-1)*Initint[2]+j-1]+SH[0])<0&&(SH[1]+enthalpyDPT[(i-1)*Initint[2]+j-1])<0))
                    {
                        Initdouble[3]=T1;
                        Initint[3]=i;
                        Initint[4]=j;
                    }
                }
            }
        if(type==2||type==3)
        {
         //THAL_END1
            Initint[4]=0;
            Initint[3]=Initint[0];
            i=Initint[0];
            Initdouble[3]=-1.0*INFINITY;
            for (j=1;j<=Initint[1];++j)
            {
                RSH(i,j,SH,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,tstack2Entropies,tstack2Enthalpies,atpS,atpH,Initdouble,numSeq1,numSeq2);
                SH[0]=SH[0]+0.000001; // this adding is done for compiler, optimization -O2 vs -O0,that compiler could understand that SH is changed in this cycle 
                SH[1]=SH[1]+0.000001;
                T1=((enthalpyDPT[(i-1)*Initint[2]+j-1]+SH[1]+Initdouble[0])/((entropyDPT[(i-1)*Initint[2]+j-1])+SH[0]+Initdouble[1]+Initdouble[2]))-273.15;
                if (T1>Initdouble[3]&&((SH[0]+entropyDPT[(i-1)*Initint[2]+j-1])<0&&(SH[1]+enthalpyDPT[(i-1)*Initint[2]+j-1])<0))
                {
                    Initdouble[3]=T1;
                    Initint[4]=j;
                }
            }
        }
        if(fabs(Initdouble[3])>999999999)
            Initint[3]=Initint[4]=1;
        RSH(Initint[3], Initint[4], SH,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,tstack2Entropies,tstack2Enthalpies,atpS,atpH,Initdouble,numSeq1,numSeq2);
     // tracebacking 
        for (i=0;i<Initint[0];++i)
            ps1[i]=0;
        for (j=0;j<Initint[1];++j)
            ps2[j] = 0;
        if(fabs(enthalpyDPT[(Initint[3]-1)*Initint[2]+Initint[4]-1])<999999999)
        {
            traceback(Initint[3],Initint[4],ps1,ps2,stackEntropies,stackEnthalpies,stackint2Entropies,stackint2Enthalpies,dangleEntropies3,dangleEnthalpies3,dangleEntropies5,dangleEnthalpies5,interiorLoopEntropies,bulgeLoopEntropies,interiorLoopEnthalpies,bulgeLoopEnthalpies,tstackEntropies,tstackEnthalpies,tstack2Entropies,tstack2Enthalpies,atpS,atpH,Initdouble,Initint,enthalpyDPT,entropyDPT,numSeq1,numSeq2);
            result_TH=drawDimer(ps1,ps2,(enthalpyDPT[(Initint[3]-1)*Initint[2]+Initint[4]-1]+SH[1]+Initdouble[0]),(entropyDPT[(Initint[3]-1)*Initint[2]+Initint[4]-1]+SH[0]+Initdouble[1]),Initdouble,Initint);
            result_TH=(int)(result_TH*100+0.5)/100.0;
        }
        return result_TH;
    }
    return result_TH;
}