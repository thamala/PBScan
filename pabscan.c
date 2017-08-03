/*
 
 PaBScan, program for performing selection outlier scans with population branch statistic (PBS)
 
 Code written by Tuomas Hamala
 tuomas.hamala@gmail.com
 
 Manual available at: https://github.com/thamala/PaBScan
 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include "pabscan.h"
#define empty "EMPTY"
#define merror "\nERROR: System out of memory\n"
#define version "2017.07.27"

int main(int argc, char *argv[]){
    
    int second=0, minute=0, hour=0;
    time_t timer=0;
    
    if(argc == 1){
        printf("\nPaBScan version %s\n\n", version);
        showHelp();
        return 0;
    }
    
    timer = time(NULL);
    
    openFiles(argc, argv);
    
    second = time(NULL) - timer;
    
    minute = second / 60;
    
    hour = second / 3600;
    
    if(hour > 0)
        printf("Elapsed time: %i h, %i min & %i sec\n\n", hour, minute-hour*60, second-minute*60);
    else if(minute > 0)
        printf("Elapsed time: %i min & %i sec\n\n", minute, second-minute*60);
    else if(second >= 5)
        printf("Elapsed time: %i sec\n\n", second);
    else
        printf("\n");
    
    return 0;
}

void openFiles(int argc, char *argv[]){
	
    int i=0, win=1, step=1, perm=0, div=0, indiv_n=0, snp_n=0, pop1_n=0, pop2_n=0, pop3_n=0, sims=0, min=1, *poplist;
	int v_ok=0, vp_ok=0, i_ok=0, l_ok=0, p1_ok=0, p2_ok=0, p3_ok=0, ms_ok=0, msp_ok=0, o_ok=0, min_ok=0, maf_ok=0;
    double maf=0.01, **ms;
    char **pop1, **pop2, **pop3, name[21], namec[26];
	Geno_s *sites;
    Pbs_s **estimates;
	FILE *geno_file, *pop1_file, *pop2_file, *pop3_file, *ms_file, *out_file;
    
    printf("\nPaBScan version %s\n\n", version);
	
	puts("Started with parameters:");
	
	for(i=1;i<argc;i++){
            
        if(strcmp(argv[i], "-vcf") == 0 || strcmp(argv[i], "-vcfp") == 0 || strcmp(argv[i], "-in") == 0 || strcmp(argv[i], "-likes") == 0){
            if((geno_file = fopen(argv[++i], "r")) == NULL){
                printf("\nERROR: Cannot open file %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            if(strcmp(argv[i-1], "-vcf") == 0){
                printf("\t-vcf %s\n", argv[i]);
                v_ok = 1;
            }
            else if(strcmp(argv[i-1], "-vcfp") == 0){
                printf("\t-vcfp %s\n", argv[i]);
                vp_ok = 1;
            }
            else if(strcmp(argv[i-1], "-in") == 0){
                printf("\t-in %s\n", argv[i]);
                i_ok = 1;
            }
            else{
                printf("\t-likes %s\n", argv[i]);
                l_ok = 1;
            }
        }
        
        else if(strcmp(argv[i], "-pop1") == 0){
            if((pop1_file = fopen(argv[++i], "r")) == NULL){
                printf("\nERROR: Cannot open file %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            printf("\t-pop1 %s\n", argv[i]);
            p1_ok = 1;
        }
        
        else if(strcmp(argv[i], "-pop2") == 0){
            if((pop2_file = fopen(argv[++i], "r")) == NULL){
                printf("\nERROR: Cannot open file %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            printf("\t-pop2 %s\n", argv[i]);
            p2_ok = 1;
        }
        
        else if(strcmp(argv[i], "-pop3") == 0){
            if((pop3_file = fopen(argv[++i], "r")) == NULL){
                printf("\nERROR: Cannot open file %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            printf("\t-pop3 %s\n", argv[i]);
            p3_ok = 1;
        }
        
        else if(strcmp(argv[i], "-ms") == 0 || strcmp(argv[i], "-msp") == 0){
            if((ms_file = fopen(argv[++i], "r")) == NULL){
                printf("\nERROR: Cannot open file %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            if(strcmp(argv[i-1], "-msp") == 0){
                printf("\t-msp %s\n", argv[i]);
                msp_ok = 1;
            }
            else{
                printf("\t-ms %s\n", argv[i]);
                ms_ok = 1;
            }
        }
        
        else if(strcmp(argv[i], "-out") == 0){
            strncpy(name, argv[++i], 20);
            strcpy(namec, name);
            
            strcat(namec, ".pbs");
            
            if((out_file = fopen(namec, "w")) == NULL){
                printf("\nERROR: Cannot create file %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            printf("\t-out %s\n", argv[i]);
            o_ok = 1;
        }
        
        else if(strcmp(argv[i], "-win") == 0){
            win = atoi(argv[++i]);
            if(win == 0)
                win = 1;
            printf("\t-win %i\n", win);
        }
        
        else if(strcmp(argv[i], "-step") == 0){
            step = atoi(argv[++i]);
            if(step == 0)
                step = 1;
            printf("\t-step %i\n", step);
        }
        
        else if(strcmp(argv[i], "-mc") == 0){
            perm = atoi(argv[++i]);
            printf("\t-mc %i\n", perm);
        }
        
        else if(strcmp(argv[i], "-div") == 0){
            div = atoi(argv[++i]);
            if(div > 2){
                printf("ERROR: Unknown parameter for '-div'\nUse -help' for more information\n\n");
                exit(EXIT_FAILURE);
            }
            printf("\t-div %i\n", div);
        }
        
        else if(strcmp(argv[i], "-min") == 0){
            min = atoi(argv[++i]);
            if(min == 0)
                min = 1;
            printf("\t-min %i\n", min);
            min_ok = 1;
        }
        
        else if(strcmp(argv[i], "-maf") == 0){
            maf = atof(argv[++i]);
            if(maf == 0)
                printf("\t-maf %.0f\n", maf);
            else if(maf >= 0.01)
                printf("\t-maf %.2f\n", maf);
            else if(maf >= 0.001)
                printf("\t-maf %.3f\n", maf);
            else
                printf("\t-maf %f\n", maf);
            maf_ok = 1;
        }
        
        else if(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "-help") == 0 || strcmp(argv[i], "--help") == 0){
            puts("\t-help\n");
            showHelp();
            if(argc == 2)
                exit(EXIT_FAILURE);
            i++;
        }
        
        else{
            printf("\nERROR: Unknown argument '%s'\nUse '-help' for more information\n\n", argv[i]);
            exit(EXIT_FAILURE);
        }
	}
    
    if(maf_ok == 0 && min_ok == 0)
        puts("\t-min 1\n\t-maf 0.01");
    else if(maf_ok == 0)
        puts("\t-maf 0.01");
    else if(min_ok == 0)
        puts("\t-min 1");
    
    printf("\n");
    
    if(v_ok == 1 && vp_ok == 1)
        v_ok = 0;
    
    if(ms_ok == 1 && msp_ok == 1)
        ms_ok = 0;
	
	if(ms_ok == 1 && perm > 0){
		puts("\nWarning: When using the '-ms' argument, '-mc' becomes redundant\n");
		perm = 0;
	}
    
    if((v_ok == 1 && (i_ok == 1 || l_ok == 1)) || (i_ok == 1 && (v_ok == 1 || vp_ok == 1 || l_ok == 1)) || (l_ok == 1 && (v_ok == 1 || vp_ok == 1 || i_ok == 1))){
        puts("\nERROR: More than one input format specified\nUse '-help' for more information\n");
        exit(EXIT_FAILURE);
    }
	
    if((v_ok == 1 || vp_ok == 1 || i_ok == 1 || l_ok == 1) && p1_ok == 1 && p2_ok == 1 && p3_ok == 1 && o_ok == 1){
        
        puts("Reading input files...\n");
        
        pop1 = readPop(pop1_file, &pop1_n);
        
        pop2 = readPop(pop2_file, &pop2_n);
        
        pop3 = readPop(pop3_file, &pop3_n);
        
        indiv_n = pop1_n + pop2_n + pop3_n;
        
        if((poplist = malloc((indiv_n+2)*sizeof(int))) == NULL){
            puts(merror);
            exit(EXIT_FAILURE);
        }
       
        if(v_ok == 1)
            sites = readVcf(geno_file, pop1, pop2, pop3, name, indiv_n, &snp_n, poplist, 0);
        else if(vp_ok == 1)
            sites = readVcf(geno_file, pop1, pop2, pop3, name, indiv_n, &snp_n, poplist, 1);
        else if(i_ok == 1)
            sites = readNative(geno_file, pop1, pop2, pop3, indiv_n, &snp_n, poplist);
        else
            sites = readBeagle(geno_file, pop1, pop2, pop3, indiv_n, &snp_n, poplist);
        
        printf("Read %i individuals\n", indiv_n);
        
		if(ms_ok == 1)
			ms = processMs(ms_file, name, pop1_n, pop2_n, pop3_n, &sims, div, 0, maf);
        else if(msp_ok == 1)
            ms = processMs(ms_file, name, pop1_n, pop2_n, pop3_n, &sims, div, 1, maf);
        
        if(l_ok == 1)
            estimates = estimatePBS(sites, poplist, ms, sims, snp_n, indiv_n, win, step, perm, div, min, 1, maf);
        else
            estimates = estimatePBS(sites, poplist, ms, sims, snp_n, indiv_n, win, step, perm, div, min, 0, maf);
        
        printOutput(out_file, estimates, win, perm, sims);

    }
    
	else{
        
		if(v_ok == 0 && vp_ok == 0 && i_ok == 0 && l_ok == 0)
			puts("ERROR: Genotype file missing\nUse '-help' for more information\n");
        
		if(p1_ok == 0 || p2_ok == 0 || p3_ok == 0)
			puts("ERROR: Three population files must be provided\nUse '-help' for more information\n");
        
		if(o_ok == 0)
			puts("ERROR: Output file missing\nUse '-help' for more information\n");
	}
}

char **readPop(FILE *pop_file, int *n){
    
    int i=0, char_i=0, line_i=0, maxchar=0;
    char c, **list, *line;
    
    while((c=fgetc(pop_file)) != EOF){
        char_i++;
        if(c == '\n'){
            line_i++;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }
    
    rewind(pop_file);
    
    if((list = malloc((line_i+2)*sizeof(char*))) == NULL){
        puts(merror);
        exit(EXIT_FAILURE);
    }
    
    for(i=0;i<line_i+2;i++){
        if((list[i] = malloc((maxchar+2)*sizeof(char))) == NULL){
            puts(merror);
            exit(EXIT_FAILURE);
        }
    }
    
    if((line = malloc((maxchar+5)*sizeof(char))) == NULL){
        puts(merror);
        exit(EXIT_FAILURE);
    }
    
    i = 0;
    
    while(fgets(line, maxchar+5, pop_file) != NULL){
        lineTerminator(line);
        strcpy(list[i], line);
        i++;
    }
    
    *n = i;
    
    strcpy(list[i], empty);
    
    free(line);
    
    fclose(pop_file);
    
    return(list);
}

Geno_s *readVcf(FILE *geno_file, char **pop1, char **pop2, char **pop3, char name[], int indiv_n, int *snp_i, int *poplist, int print){
	
	int i=0, j=0, line_i=0, char_i=0, maxchar=0, split_i=0, indiv_i=0, pop_i=0, id_ok=0, *tempind, **temparray;
	char c, *line, *temp=NULL, *p=NULL, *substring, **templist, ref[10], alt[10], namec[30];
    Geno_s *sites;
    FILE *pabs_file;

	while((c=fgetc(geno_file)) != EOF){
        char_i++;
        if(c == '\n'){
            line_i++;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }
	
	rewind(geno_file);
	
	if((sites = malloc((line_i+2)*sizeof(Geno_s))) == NULL){
		puts(merror);
		exit(EXIT_FAILURE);
	}
    
    if((templist = malloc((indiv_n+2)*sizeof(char*))) == NULL){
        puts(merror);
        exit(EXIT_FAILURE);
    }
    
    for(i=0;i<indiv_n+2;i++){
        if((templist[i] = malloc(6*sizeof(char))) == NULL){
            puts(merror);
            exit(EXIT_FAILURE);
        }
    }
	
	if((line = malloc((maxchar+5)*sizeof(char))) == NULL){
        puts(merror);
        exit(EXIT_FAILURE);
    }
    
    if((tempind = malloc((indiv_n+2)*sizeof(int))) == NULL){
        puts(merror);
        exit(EXIT_FAILURE);
    }
    
    if(print == 1){
        
        strcpy(namec, name);
        
        strcat(namec, ".pabscan");
        
        if((pabs_file = fopen(namec, "w")) == NULL){
            printf("\nWarning: Cannot create file %s\n\n", namec);
            print = 0;
        }
        else
            fprintf(pabs_file, "CHR\tBP\t");
    }
    
    temp = NULL;
	
	while(fgets(line, maxchar+5, geno_file) != NULL){
        
		lineTerminator(line);
		
		temp = strtok(line, "\t");
        
        if(strcmp(temp,"#CHROM") == 0 && temp != NULL){
            
            while(temp != NULL){
                
                if(split_i >= 9){
                    
                    for(i=0;strcmp(pop1[i],empty)!=0;i++){
                        if(strcmp(pop1[i], temp) == 0){
                            
                            tempind[indiv_i] = split_i;
                            poplist[indiv_i] = 1;
                            indiv_i++;
                            
                            if(print == 1)
                                fprintf(pabs_file, "%s\t", temp);
                        }
                    }
                    
                    for(i=0;strcmp(pop2[i],empty)!=0;i++){
                        if(strcmp(pop2[i], temp) == 0){
    
                            tempind[indiv_i] = split_i;
                            poplist[indiv_i] = 2;
                            indiv_i++;
                            
                            if(print == 1)
                                fprintf(pabs_file, "%s\t", temp);
                        }
                    }
                    
                    for(i=0;strcmp(pop3[i],empty)!=0;i++){
                        if(strcmp(pop3[i], temp) == 0){
                            
                            tempind[indiv_i] = split_i;
                            poplist[indiv_i] = 3;
                            indiv_i++;
                            
                            if(print == 1)
                                fprintf(pabs_file, "%s\t", temp);
                        }
                    }
                }
                
                temp = strtok(NULL, "\t");
                split_i++;
            }
            
            poplist[indiv_i] = 0;
            
            if(print == 1)
                fprintf(pabs_file, "\n");
            
            id_ok = 1;

        }
        
        else if(temp[0] != '#' && id_ok == 0){
            puts("\nERROR: VCF file is not in the specified format\n");
            exit(EXIT_FAILURE);
        }
		
        else{
            
            if((sites[*snp_i].genolist = malloc((indiv_n+2)*sizeof(int*))) == NULL){
                puts(merror);
                exit(EXIT_FAILURE);
            }
		
            strncpy(sites[*snp_i].chromo, temp, 20);
            
            if(print == 1)
                fprintf(pabs_file, "%s\t", sites[*snp_i].chromo);
            
			split_i=0;
			indiv_i=0;
			
			while(temp != NULL){
                
                if(split_i == 1){
                    
					sites[*snp_i].position = atoi(temp);
                    
                    if(print == 1)
                        fprintf(pabs_file, "%i\t", sites[*snp_i].position);
                    
                    if(*snp_i != 0){
                        if(strcmp(sites[*snp_i].chromo, sites[*snp_i-1].chromo) == 0 && sites[*snp_i].position < sites[*snp_i-1].position){
                            puts("ERROR: VCF file is not sorted\n");
                            exit(EXIT_FAILURE);
                        }
                    }
                }
                    
                for(i=0;i<indiv_n;i++){
                    
                    if(split_i == tempind[i]){
        
                        if((sites[*snp_i].genolist[indiv_i] = malloc(3*sizeof(int))) == NULL){
                            puts(merror);
                            exit(EXIT_FAILURE);
                        }
                        
                        strncpy(templist[indiv_i], temp, 4);
                        
                        indiv_i++;
                    }
                }

				temp = strtok(NULL, "\t");
                
				split_i++;
			}
			
			for(i=0;i<indiv_i;i++){
                
                sites[*snp_i].genolist[i][0] = -1;
                sites[*snp_i].genolist[i][1] = -1;
                
				p = strchr(templist[i], '/');
				
				if(p == NULL)
					p = strchr(templist[i], '|');
    
				if(p != NULL){
					substring = strtok(templist[i], ":");
		
					temp = strtok(substring, "/");
		
					if(temp == NULL)
						temp = strtok(substring, "|");
                    
                    if(isdigit(temp[0]))
                        sites[*snp_i].genolist[i][0] = atoi(temp);
                        
                    temp = strtok(NULL, "/");
        
                    if(temp == NULL)
                        temp=strtok(NULL, "|");
        
                    if(isdigit(temp[0]))
                        sites[*snp_i].genolist[i][1] = atoi(temp);
                    
                    if(print == 1){
                        
                        if(sites[*snp_i].genolist[i][0] == -1)
                            fprintf(pabs_file, "NA\t");
                        else
                            fprintf(pabs_file, "%i%i\t", sites[*snp_i].genolist[i][0], sites[*snp_i].genolist[i][1]);
                    }
				}
			}
            
            if(print == 1)
                fprintf(pabs_file, "\n");
            
			*snp_i = *snp_i+1;
		}
        
	}
    
    printf("Read %i sites\n", *snp_i);
    
    free(line);
    
    free(tempind);
    
    for(i=0;i<indiv_n+2;i++)
        free(templist[i]);
    free(templist);
    
    for(i=0;strcmp(pop1[i],empty)!=0;i++)
        free(pop1[i]);
    free(pop1[i]);
    free(pop1);
    
    for(i=0;strcmp(pop2[i],empty)!=0;i++)
        free(pop2[i]);
    free(pop2[i]);
    free(pop2);
    
    for(i=0;strcmp(pop3[i],empty)!=0;i++)
        free(pop3[i]);
    free(pop3[i]);
    free(pop3);
    
    if(print == 1)
        fclose(pabs_file);
    
    fclose(geno_file);
    
	return sites;
}

Geno_s *readNative(FILE *geno_file, char **pop1, char **pop2, char **pop3, int indiv_n, int *snp_i, int *poplist){
    
    int i=0, char_i=0, line_i=0, maxchar=0, indiv_i=0, indiv_mem=0;
    char c, *line, *temp;
    Geno_s *sites;
    
    while((c=fgetc(geno_file)) != EOF){
        char_i++;
        if(c == '\n'){
            line_i++;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }
    
    rewind(geno_file);
    
    if((sites = malloc((line_i+2)*sizeof(Geno_s))) == NULL){
        puts(merror);
        exit(EXIT_FAILURE);
    }
    
    if((line = malloc((maxchar+5)*sizeof(char))) == NULL){
        puts(merror);
        exit(EXIT_FAILURE);
    }
    
    fgets(line, maxchar+5, geno_file);
    
    lineTerminator(line);
    
    temp = strtok(line, "\t");
    
    if(strcmp(temp, "CHR") != 0){
        puts("\nERROR: Input file is not in the specified format\n");
        exit(EXIT_FAILURE);
    }
    
    else{
        
        temp = strtok(NULL, "\t");
        
        if(strcmp(temp, "BP") != 0){
            puts("\nERROR: Input file is not in the specified format\n");
            exit(EXIT_FAILURE);
        }
        
        else{
            
            temp = strtok(NULL, "\t");
            
            while(temp != NULL){
                
                indiv_mem = indiv_i;
                
                for(i=0;strcmp(pop1[i],empty)!=0;i++){
                    if(strcmp(pop1[i], temp) == 0){
                        
                        poplist[indiv_i] = 1;
                        indiv_i++;
                    }
                }
                
                for(i=0;strcmp(pop2[i],empty)!=0;i++){
                    if(strcmp(pop2[i], temp) == 0){
 
                        poplist[indiv_i] = 2;
                        indiv_i++;
                    }
                }
                
                for(i=0;strcmp(pop3[i],empty)!=0;i++){
                    if(strcmp(pop3[i], temp) == 0){
                        
                        poplist[indiv_i] = 3;
                        indiv_i++;
                    }
                }
                
                if(indiv_mem == indiv_i){
                    puts("\nERROR: Input file contains unspecified individuals\n");
                    exit(EXIT_FAILURE);
                }
                
                temp = strtok(NULL, "\t");
            }
            
            poplist[indiv_i] = 0;
        }
    }

    while(fgets(line, maxchar+5, geno_file) != NULL){
        
        indiv_i = 0;
    
        lineTerminator(line);
        
        if((sites[*snp_i].genolist = malloc((indiv_n+2)*sizeof(int*))) == NULL){
            puts(merror);
            exit(EXIT_FAILURE);
        }
        
        strncpy(sites[*snp_i].chromo, strtok(line, "\t"), 20);
        
        sites[*snp_i].position = atoi(strtok(NULL, "\t"));
        
        if(*snp_i != 0){
            if(strcmp(sites[*snp_i].chromo, sites[*snp_i-1].chromo) == 0 && sites[*snp_i].position < sites[*snp_i-1].position){
                puts("ERROR: Input file is not sorted\n");
                exit(EXIT_FAILURE);
            }
        }
        
        temp = strtok(NULL, "\t");
        
        while(temp != NULL){
            
            if((sites[*snp_i].genolist[indiv_i] = malloc(3*sizeof(int))) == NULL){
                puts(merror);
                exit(EXIT_FAILURE);
            }
            
            if(strcmp(temp, "NA") == 0){
                sites[*snp_i].genolist[indiv_i][0] = -1;
                sites[*snp_i].genolist[indiv_i][1] = -1;
            }
            
            else{
                sites[*snp_i].genolist[indiv_i][0] = temp[0] - '0';
                sites[*snp_i].genolist[indiv_i][1] = temp[1] - '0';
            }
            
            indiv_i++;
            
            temp = strtok(NULL, "\t");
        }
        
        *snp_i = *snp_i + 1;
    }
    
    printf("Read %i sites\n", *snp_i);
    
    free(line);
    
    for(i=0;strcmp(pop1[i],empty)!=0;i++)
        free(pop1[i]);
    free(pop1[i]);
    free(pop1);
    
    for(i=0;strcmp(pop2[i],empty)!=0;i++)
        free(pop2[i]);
    free(pop2[i]);
    free(pop2);
    
    for(i=0;strcmp(pop3[i],empty)!=0;i++)
        free(pop3[i]);
    free(pop3[i]);
    free(pop3);
    
    fclose(geno_file);
    
    return sites;
}
       
Geno_s *readBeagle(FILE *geno_file, char **pop1, char **pop2, char **pop3, int indiv_n, int *snp_i, int *poplist){
    
    int i=0, j=0, char_i=0, line_i=0, maxchar=0, indiv_i=0, split_i=0, geno_i=0, *tempind;
    char c, *line, *temp, *tempcp;
    Geno_s *sites;
    
    while((c=fgetc(geno_file)) != EOF){
        char_i++;
        if(c == '\n'){
            line_i++;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }
    
    rewind(geno_file);
    
    if((sites = malloc((line_i+2)*sizeof(Geno_s))) == NULL){
        puts(merror);
        exit(EXIT_FAILURE);
    }
    
    if((line = malloc((maxchar+5)*sizeof(char))) == NULL){
        puts(merror);
        exit(EXIT_FAILURE);
    }
    
    if((tempind = malloc((indiv_n+2)*3*sizeof(int))) == NULL){
        puts(merror);
        exit(EXIT_FAILURE);
    }

    fgets(line, maxchar+5, geno_file);
    
    lineTerminator(line);
    
    temp = strtok(line, "\t");
    
    if(strcmp(temp, "marker") != 0){
        puts("\nERROR: Likelihood file is not in the specified format\n");
        exit(EXIT_FAILURE);
    }
            
    while(temp != NULL){
        
        if(split_i > 2){
        
            for(i=0;strcmp(pop1[i],empty)!=0;i++){
                if(strcmp(pop1[i], temp) == 0){
                
                    strcpy(pop1[i], "F");
                    
                    for(j=0;j<3;j++){
                        tempind[geno_i] = split_i+j;
                        geno_i++;
                    }
                    
                    poplist[indiv_i] = 1;
                    indiv_i++;
                }
            }
            
            for(i=0;strcmp(pop2[i],empty)!=0;i++){
                if(strcmp(pop2[i], temp) == 0){
                    
                    strcpy(pop2[i], "F");
                    
                    for(j=0;j<3;j++){
                        tempind[geno_i] = split_i+j;
                        geno_i++;
                    }
                    
                    poplist[indiv_i] = 2;
                    indiv_i++;
                }
            }
            
            for(i=0;strcmp(pop3[i],empty)!=0;i++){
                if(strcmp(pop3[i], temp) == 0){
                    
                    strcpy(pop3[i], "F");
                    
                    for(j=0;j<3;j++){
                        tempind[geno_i] = split_i+j;
                        geno_i++;
                    }
                    
                    poplist[indiv_i] = 3;
                    indiv_i++;
                }
            }
        }

        temp = strtok(NULL, "\t");
        split_i++;
    }
        
    poplist[indiv_i] = 0;
    
    while(fgets(line, maxchar+5, geno_file) != NULL){
        
        lineTerminator(line);
        
        if((sites[*snp_i].likelist = malloc((indiv_n+2)*sizeof(double*))) == NULL){
            puts(merror);
            exit(EXIT_FAILURE);
        }

        tempcp = strtok(line, "\t");
        
        temp = strtok(NULL, "\t");
        
        split_i = 1;
        indiv_i = 0;
        geno_i = 0;
        
        while(temp != NULL){
            
            for(i=0;i<indiv_n*3;i++){
                
                if(split_i == tempind[i]){

                    if(geno_i == 0){
                
                        if((sites[*snp_i].likelist[indiv_i] = malloc(4*sizeof(double))) == NULL){
                            puts(merror);
                            exit(EXIT_FAILURE);
                        }
                    }
                    
                    sites[*snp_i].likelist[indiv_i][geno_i] = atof(temp);
                    
                    geno_i++;
                    
                    if(geno_i == 3){
                        geno_i = 0;
                        indiv_i++;
                    }
                }
            }
            
            temp = strtok(NULL, "\t");
            split_i++;
        }

        strncpy(sites[*snp_i].chromo, strtok(tempcp, "_"), 20);
        sites[*snp_i].position = atoi(strtok(NULL, "_"));
        
        if(*snp_i != 0){
            if(strcmp(sites[*snp_i].chromo, sites[*snp_i-1].chromo) == 0 && sites[*snp_i].position < sites[*snp_i-1].position){
                puts("ERROR: Input file is not sorted\n");
                exit(EXIT_FAILURE);
            }
        }
        
        *snp_i = *snp_i + 1;
    }
    
    printf("Read %i sites\n", *snp_i);
    
    free(line);
    free(tempind);
    
    for(i=0;strcmp(pop1[i],empty)!=0;i++)
        free(pop1[i]);
    free(pop1[i]);
    free(pop1);
    
    for(i=0;strcmp(pop2[i],empty)!=0;i++)
        free(pop2[i]);
    free(pop2[i]);
    free(pop2);
    
    for(i=0;strcmp(pop3[i],empty)!=0;i++)
        free(pop3[i]);
    free(pop3[i]);
    free(pop3);
    
    fclose(geno_file);
    
    return sites;
}

double **processMs(FILE *ms_file, char name[], int pop1_n, int pop2_n, int pop3_n, int *sims, int div, int print, double maf){

    int i=0, j=0, site_i=0, char_i=0, seq_i=0, line_i=0, lastline=0, maxchar=0, n=0, seqsite=0, site_ok=0, div_i[3]={0};
    double t_ms[3]={0}, pbs_ms[3]={0}, hw[3]={0}, hb[3]={0}, dxy[3]={0}, n1=0, n2=0, p1=0, p2=0;
    double **alleles, **freqlist, **ms;
	char c, *line, *temp, namec[35];
    FILE *null_file;

	while((c=fgetc(ms_file)) != EOF){
        char_i++;
		if(c == '/')
			site_i++;
        else if(c == '\n'){
			lastline++;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }
	
	rewind(ms_file);
    
    pop1_n = pop1_n * 2;
    pop2_n = pop2_n * 2;
    pop3_n = pop3_n * 2;
	
	n = pop1_n + pop2_n + pop3_n;
    
    if((ms = malloc(4*sizeof(double*))) == NULL){
        puts(merror);
        exit(EXIT_FAILURE);
    }
    
    for(i=0;i<4;i++){
        if((ms[i] = malloc((site_i/2+2)*sizeof(double))) == NULL){
            puts(merror);
            exit(EXIT_FAILURE);
        }
    }
	
	if((line = malloc((maxchar+5)*sizeof(char))) == NULL){
		puts(merror);
		exit(EXIT_FAILURE);
	}
	
    if((freqlist = malloc(4*sizeof(double*))) == NULL){
        puts(merror);
        exit(EXIT_FAILURE);
    }

	
	if((alleles = malloc((n+2)*sizeof(double*))) == NULL){
		puts(merror);
		exit(EXIT_FAILURE);
	}
    
    if(print == 1){
        
        strcpy(namec, name);

        strcat(namec, ".nulldist");
        
        if((null_file = fopen(namec, "w")) == NULL){
            printf("\nWarning: Cannot create file %s\n\n", name);
            print = 0;
        }
        
        else
            fprintf(null_file,"PBS1\tPBS2\tPBS3\n");
    }
	
	while(fgets(line, maxchar+5, ms_file) != NULL){
	
		line_i++;
	
		if(isalpha(line[0])){
			temp = strtok(line, " ");
			
			if(strcmp(temp, "segsites:") == 0){
				seqsite = atoi(strtok(NULL, " "));
				
				for(j=0;j<4;j++){
					if((freqlist[j] = malloc((seqsite+2)*sizeof(double))) == NULL){
						puts(merror);
						exit(EXIT_FAILURE);
					}
				}
				
				for(j=0;j<(n+2);j++){
					if((alleles[j] = malloc((seqsite+2)*sizeof(double))) == NULL){
						puts(merror);
						exit(EXIT_FAILURE);
					}
				}
			}
		}
		
		if((line[0] == '0' || line[0] == '1') && line_i > 4){
			for(j=0;j<strlen(line)-1;j++){
				alleles[seq_i][j] = line[j] - '0';
			}
			seq_i++;
		}
        
		
		if(((line[0] == '/' && line[1] == '/') || line_i == lastline) && line_i > 4){
			
			for(i=0;i<n;i++){
                
				for(j=0;j<seqsite;j++){
                    
                    if(i == 0){
                        freqlist[0][j] = 0;
                        freqlist[1][j] = 0;
                        freqlist[2][j] = 0;
                    }
					
                    if(i < pop1_n)
						freqlist[0][j] = freqlist[0][j] + alleles[i][j];
						
                    else if(i < pop1_n + pop2_n)
						freqlist[1][j] = freqlist[1][j] + alleles[i][j];
						
                    else
						freqlist[2][j] = freqlist[2][j] + alleles[i][j];
				}
			}
			
			for(i=0;i<seqsite;i++){
                
                for(j=0;j<3;j++) {
                    
                    switch(j){
                            
                        case 0:
                            n1 = (double)pop1_n;
                            n2 = (double)pop2_n;
                            p1 = freqlist[0][i] / (double)pop1_n;
                            p2 = freqlist[1][i] / (double)pop2_n;
                            break;
                            
                        case 1:
                            n1 = (double)pop1_n;
                            n2 = (double)pop3_n;
                            p1 = freqlist[0][i] / (double)pop1_n;
                            p2 = freqlist[2][i] / (double)pop3_n;
                            break;
                        
                        case 2:
                            n1 = (double)pop2_n;
                            n2 = (double)pop3_n;
                            p1 = freqlist[1][i] / (double)pop2_n;
                            p2 = freqlist[2][i] / (double)pop3_n;
                            break;
                    }
                    
                    //Samples not meeting the MAF cutoff are removed
                    if(((n1 * p1) + (n2 * p2)) / (n1 + n2) >= maf){
                        
                        //Calculates divergence with Weir & Cokerham's Fst
                        if(div == 1){
                            hw[j] = hw[j]+(2*(n1*n2/(n1+n2))*1/(n1+n2-2)*(n1*p1*(1-p1)+n2*p2*(1-p2)));
                            hb[j] = hb[j]+(n1*n2/(n1+n2)*pow((p1-p2),2)+(2*(n1*n2/(n1+n2))-1)*1/(n1+n2-2)*(n1*p1*(1-p1)+n2*p2*(1-p2)));
                        }
                        
                        //Calculates divergence with Dxy
                        else if(div == 2){
                            dxy[j] = dxy[j]+(p1*(1-p2)+p2*(1-p1));
                            div_i[j]++;
                        
                        }
                        
                        //As default, calculates divergence with Hudson's Fst
                        else{
                            hw[j] = hw[j]+(pow((p1-p2),2)-p1*(1-p1)/(n1-1)-p2*(1-p2)/(n2-1));
                            hb[j] = hb[j]+(p1*(1-p2)+p2*(1-p1));
                        }
                    }
                }
            }
            
            for(i=0;i<3;i++){
                if(div == 1)
                    t_ms[i] = -log(1 - (1 - hw[i] / hb[i]));
                else if(div == 2)
                    t_ms[i] = -log(1 - (dxy[i] / div_i[i]));
                else
                    t_ms[i] = -log(1 - (hw[i] / hb[i]));
                
                hw[i] = 0;
                hb[i] = 0;
                dxy[i] = 0;
                div_i[i] = 0;
            }
            
            ms[0][*sims] = (t_ms[0] + t_ms[1] - t_ms[2]) / 2;
            ms[1][*sims] = (t_ms[0] + t_ms[2] - t_ms[1]) / 2;
            ms[2][*sims] = (t_ms[2] + t_ms[1] - t_ms[0]) / 2;
            
            if(print == 1){
                for(i=0;i<3;i++)
                    fprintf(null_file, "%f\t", ms[i][*sims]);
                fprintf(null_file, "\n");
            }
			
			*sims = *sims + 1;
			seq_i = 0;
		}
	}
    
    printf("Read %i simulated samples\n", *sims);
		
	for(i=0;i<n+2;i++)
		free(alleles[i]);
	free(alleles);
    
    for(i=0;i<4;i++)
        free(freqlist[i]);
    free(freqlist);
    
	free(line);
    
    fclose(ms_file);
    
    if(print == 1)
        fclose(null_file);
	
	return(ms);
}

Pbs_s **estimatePBS(Geno_s *sites, int poplist[], double **ms, int sims, int snp_n, int indiv_n, int win, int step, int perm, int div, int min, int like, double maf){
    
    int i=0, j=0, k=0, *permsite, indiv_i=0, est_i=0, site_i=0, win_i=0, print_i=0, p_i=0;
    double mid=0, t_div[3]={0}, pbs[3]={0}, fst_ms[3]={0}, **pbs_ms;
    double t_div_perm[3]={0}, pbs_perm[3]={0}, pbs_win_perm[3]={0}, **perm_high;
    Pbs_s **estimates;
    Div_s vars[3], vars_perm[3], div_win[3], div_win_perm[3];

    //Sets seed for permutations
    if(perm > 0 && div == 0){
        printf("\nEstimating PBS with Hudson's Fst and %i permutations...\n\n", perm);
        srand(time(NULL));
    }
    
    else if(perm > 0 && div == 1){
        printf("\nEstimating PBS with Weir and Cockerham's Fst and %i permutations...\n\n", perm);
        srand(time(NULL));
    }
    
    else if(perm > 0 && div == 2){
        printf("\nEstimating PBS with Dxy and %i permutations...\n\n", perm);
        srand(time(NULL));
    }
    
    else if(sims > 0 && div == 0)
        puts("\nEstimating PBS with Hudson's Fst and simulated neutral data...\n");
    
    else if(sims > 0 && div == 1)
        puts("\nEstimating PBS with Weir & Cockerham's Fst and simulated neutral data...\n");
    
    else if(sims > 0 && div == 2)
        puts("\nEstimating PBS with Dxy and simulated neutral data...\n");
    
    else if(perm == 0 && sims == 0 && div == 1)
        puts("\nEstimating PBS with Weir & Cokerham's Fst...\n");
    
    else if(perm == 0 && sims == 0 && div == 2)
        puts("\nEstimating PBS with Dxy...\n");
    
    else
        puts("\nEstimating PBS with Hudson's Fst...\n");
    
    if((estimates = malloc(4*sizeof(Pbs_s*))) == NULL){
        puts(merror);
        exit(EXIT_FAILURE);
    }
    
    for(i=0;i<4;i++){
        
        if((estimates[i] = malloc((snp_n/win+2)*sizeof(Pbs_s))) == NULL){
            puts(merror);
            exit(EXIT_FAILURE);
        }
    }
    
    if(perm > 0){
        
        if((perm_high = malloc(4*sizeof(double*))) == NULL){
            puts(merror);
            exit(EXIT_FAILURE);
        }
        
        for(i=0;i<4;i++){
            if((perm_high[i] = malloc((snp_n/win+1)*sizeof(double))) == NULL){
                puts(merror);
                exit(EXIT_FAILURE);
            }
        }
    }
    
    if((permsite = malloc((win+1)*sizeof(int))) == NULL){
        puts(merror);
        exit(EXIT_FAILURE);
    }
    
    //Calculates PBS in windows (or single sites if win = 1) using Fst or Dxy estimates.
    //Bracketed numbers correspond to population comparisons: [0] = 1 vs. 2, [1] = 1 vs. 3, [2] = 2 vs. 3
    while(site_i < snp_n){
        
        i = site_i;
            
        win_i = 0;
        
        for(j=0;j<3;j++){
            div_win[j].hw = 0;
            div_win[j].hb = 0;
            div_win[j].dxy = 0;
        }
        
        do{
            
            vars[0] = estimateDiv(sites[i], poplist, maf, min, indiv_n, div, like, 1, 2, 0);
            vars[1] = estimateDiv(sites[i], poplist, maf, min, indiv_n, div, like, 1, 3, 0);
            vars[2] = estimateDiv(sites[i], poplist, maf, min, indiv_n, div, like, 2, 3, 0);
            
            if(isnan(vars[0].hw) == 0 && isnan(vars[1].hw) == 0 && isnan(vars[2].hw) == 0){
            
                for(j=0;j<3;j++){
                    div_win[j].hw = div_win[j].hw + vars[j].hw;
                    div_win[j].hb = div_win[j].hb + vars[j].hb;
                    div_win[j].dxy = div_win[j].dxy + vars[j].dxy;
                }
                
                if(win_i == 0){
                    strcpy(estimates[0][est_i].chromo, sites[i].chromo);
                    estimates[0][est_i].start = sites[i].position;
                }
                
                permsite[win_i] = i;
                
                win_i++;
            }
            
            else
                site_i++;
            
            i++;
            
            if(i == snp_n || strcmp(sites[i].chromo, sites[i-1].chromo) != 0){
                permsite[win_i] = -1;
                break;
            }
           
        }while(win_i < win);
        
        estimates[0][est_i].stop = sites[i-1].position;
        
        for(i=0;i<3;i++){
            if(div == 1)
                t_div[i] = -log(1 - (1 - div_win[i].hw / div_win[i].hb));
            else if(div == 2)
                t_div[i] = -log(1 - (div_win[i].dxy  / win_i));
            else
                t_div[i] = -log(1 - (div_win[i].hw / div_win[i].hb));
        }
        
        estimates[0][est_i].pbs = (t_div[0] + t_div[1] - t_div[2])/2;
        estimates[1][est_i].pbs = (t_div[0] + t_div[2] - t_div[1])/2;
        estimates[2][est_i].pbs = (t_div[2] + t_div[1] - t_div[0])/2;
        
        if(perm > 0){

            for(i=0;i<3;i++)
                perm_high[i][est_i] = -50;
            
            //Repeats the estimates 'perm' number of times, each time randomizing alleles among individuals
            for(i=0;i<perm;i++){
                
                for(j=0;j<3;j++){
                    div_win_perm[j].hw = 0;
                    div_win_perm[j].hb = 0;
                    div_win_perm[j].dxy = 0;
                }
                
                win_i = 0;

                do{
                    
                    if(permsite[win_i] == -1)
                        break;
                    
                    vars_perm[0] = estimateDiv(sites[permsite[win_i]], poplist, maf, min, indiv_n, div, like, 1, 2, 1);
                    vars_perm[1] = estimateDiv(sites[permsite[win_i]], poplist, maf, min, indiv_n, div, like, 1, 3, 1);
                    vars_perm[2] = estimateDiv(sites[permsite[win_i]], poplist, maf, min, indiv_n, div, like, 2, 3, 1);
                    
                    for(j=0;j<3;j++){
                        div_win_perm[j].hw = div_win_perm[j].hw + vars_perm[j].hw;
                        div_win_perm[j].hb = div_win_perm[j].hb + vars_perm[j].hb;
                        div_win_perm[j].dxy = div_win_perm[j].dxy + vars_perm[j].dxy;
                    }
                    
                    win_i++;

                }while(win_i < win);

                for(j=0;j<3;j++){
                    if(div == 1)
                         t_div_perm[j] = -log(1 - (1 - div_win_perm[j].hw / div_win_perm[j].hb));
                    else if(div == 2)
                        t_div_perm[j] = -log(1 - (div_win_perm[j].dxy / win_i));
                    else
                        t_div_perm[j] = -log(1 - (div_win_perm[j].hw / div_win_perm[j].hb));
                }
                
                pbs_perm[0] = (t_div_perm[0] + t_div_perm[1] - t_div_perm[2])/2;
                pbs_perm[1] = (t_div_perm[0] + t_div_perm[2] - t_div_perm[1])/2;
                pbs_perm[2] = (t_div_perm[2] + t_div_perm[1] - t_div_perm[0])/2;
                
                for(j=0;j<3;j++){
                    if(isinf(pbs_perm[j]))
                        pbs_perm[j] = -50;
                    
                    if(pbs_perm[j] > perm_high[j][est_i])
                        perm_high[j][est_i] = pbs_perm[j];
                }
            }
        }
        
        est_i++;
        
        if(step > 1){
            
            if(site_i+step-1 < snp_n){
                if(strcmp(sites[site_i+step-1].chromo, sites[site_i].chromo) != 0){
                    while(strcmp(sites[site_i-1].chromo, sites[site_i].chromo) == 0)
                        site_i++;
                }
                else
                    site_i = site_i + step - 1;
            }
            else
                site_i = site_i + step - 1;
        }
        else
            site_i++;
    }
    
    strcpy(estimates[0][est_i].chromo, empty);
    
    //Defines P-values based on the quantiles of the permutated distibution
    if(perm > 0){
        
        for(i=0;i<3;i++){
            for(j=0;j<est_i;j++){
                for(k=0;k<est_i;k++){
                    if(perm_high[i][k] > estimates[i][j].pbs)
                        p_i++;
                }
                
                estimates[i][j].p = (double)p_i / (double)est_i;
                
                p_i = 0;
            }
        }
        
        for(i=0;i<4;i++)
            free(perm_high[i]);
        free(perm_high);
    }
    
    //Defines P-values based on the quantiles of the simulated distibution
    else if(sims > 0){
        
        for(i=0;i<3;i++){
            for(j=0;j<est_i;j++){
                for(k=0;k<sims;k++){
                    if(ms[i][k] > estimates[i][j].pbs)
                        p_i++;
                }
                
                estimates[i][j].p = (double)p_i / (double)sims;
                
                p_i = 0;
            }
        }
    }
    
    if(like == 1){
        for(i=0;i<snp_n;i++){
            for(j=0;j<indiv_n;j++)
                free(sites[i].likelist[j]);
            free(sites[i].likelist);
        }
    }
    
    else{
        for(i=0;i<snp_n;i++){
            for(j=0;j<indiv_n;j++)
                free(sites[i].genolist[j]);
            free(sites[i].genolist);
        }
    }
    
    free(sites);
    
    if(sims > 0){
        for(i=0;i<4;i++)
            free(ms[i]);
        free(ms);
    }
    
    free(poplist);
    free(permsite);
    
    return(estimates);
}

Div_s estimateDiv(Geno_s sites, int poplist[], double maf, int min, int indiv_n, int div, int like, int pop1, int pop2, int perm){
    
    double n1=0, n2=0, p1=0, p2=0;
    Freq_s freqs;
    Div_s onesite;
    
    if(like == 1)
        freqs = freqLikes(sites, poplist, maf, min, indiv_n, pop1, pop2, perm);
    else
        freqs = freqCalls(sites, poplist, maf, min, indiv_n, pop1, pop2, perm);
    
    if(isnan(freqs.p1) == 0 && isnan(freqs.p2) == 0){
        
        n1 = freqs.n1;
        n2 = freqs.n2;
        p1 = freqs.p1;
        p2 = freqs.p2;
        
        onesite.hw = 0;
        onesite.hb = 0;
        onesite.dxy = 0;
    
        //Calculates divegence with Weir & Cokerham's Fst
        if(div == 1){
            onesite.hw = 2*(n1*n2/(n1+n2))*1/(n1+n2-2)*(n1*p1*(1-p1)+n2*p2*(1-p2));
            onesite.hb = n1*n2/(n1+n2)*pow((p1-p2),2)+(2*(n1*n2/(n1+n2))-1)*1/(n1+n2-2)*(n1*p1*(1-p1)+n2*p2*(1-p2));
        }

        //Calculates divegence with Dxy
        else if(div == 2)
            onesite.dxy = p1*(1-p2)+p2*(1-p1);
        
        //As default, calculates divergence with Hudson's Fst
        else{
            onesite.hw = pow((p1-p2),2)-p1*(1-p1)/(n1-1)-p2*(1-p2)/(n2-1);
            onesite.hb = p1*(1-p2)+p2*(1-p1);
        }

        return(onesite);
    }

    onesite.hw = 0.0/0.0;
    onesite.hb = 0.0/0.0;
    onesite.dxy = 0.0/0.0;
    
    return(onesite);
}

Freq_s freqCalls(Geno_s sites, int poplist[], double maf, int min, int indiv_n, int pop1, int pop2, int perm){
    
    int i=0, j=0, **genolist_new, *poplist_new, minor=-1, indiv_i=0, ind_r=0, all_r=0, temp=0;
    double alt=0, ref=0, pop1_alt=0, pop2_alt=0, pop1_ref=0, pop2_ref=0, n1=0, n2=0, n=0, minor1=0, minor2=0, p1=0, p2=0;
    Freq_s freqs;
    
    for(i=0;i<indiv_n;i++){
        
        //Minor allele is defined as the less frequent allele across all three populations
        
        if(sites.genolist[i][0] == 0)
            ref++;
        
        else if(sites.genolist[i][0] == 1)
            alt++;
        
        if(sites.genolist[i][1] == 0)
            ref++;
        
        else if(sites.genolist[i][1] == 1)
            alt++;
        
        if(poplist[i] == pop1 && sites.genolist[i][0] != -1){
            
            n1++;
            
            if(sites.genolist[i][0] == 0)
                pop1_ref++;
            else if(sites.genolist[i][0] == 1)
                pop1_alt++;
            
            if(sites.genolist[i][1] == 0)
                pop1_ref++;
            else if(sites.genolist[i][1] == 1)
                pop1_alt++;
            
        }
        
        else if(poplist[i] == pop2 && sites.genolist[i][0] != -1){
            
            n2++;
            
            if(sites.genolist[i][0] == 0)
                pop2_ref++;
            else if(sites.genolist[i][0] == 1)
                pop2_alt++;
            
            if(sites.genolist[i][1] == 0)
                pop2_ref++;
            else if(sites.genolist[i][1] == 1)
                pop2_alt++;
        }
    }
    
    n = n1 + n2;
    
    if(ref >= alt)
        minor = 1;
    else
        minor = 0;
    
    //Sites not meeting the MAF cutoff are excluded
    if(n1 >= min && n2 >= min && ((minor == 1 && (pop1_alt+pop2_alt)/(n*2) >= maf) || (minor == 0 && (pop1_ref+pop2_ref)/(n*2) >= maf))){
        
        if((genolist_new = malloc((n+3)*sizeof(int*))) == NULL){
            puts(merror);
            exit(EXIT_FAILURE);
        }
        
        for(i=0;i<n+3;i++){
            if((genolist_new[i] = malloc(3*sizeof(int))) == NULL){
                puts(merror);
                exit(EXIT_FAILURE);
            }
        }
        
        if((poplist_new = malloc((n+3)*sizeof(int))) == NULL){
            puts(merror);
            exit(EXIT_FAILURE);
        }
        
        for(i=0;i<indiv_n;i++){
            
            if((poplist[i] == pop1 || poplist[i] == pop2)  && sites.genolist[i][0] != -1){
                
                genolist_new[indiv_i][0] = sites.genolist[i][0];
                genolist_new[indiv_i][1] = sites.genolist[i][1];
                
                if(poplist[i] == pop1)
                    poplist_new[indiv_i] = pop1;
                else
                    poplist_new[indiv_i] = pop2;
                
                indiv_i++;
            }
        }
        
        //Shuffles the alleles for permutation testing
        if(perm == 1){
            
            for(i=0;i<indiv_i;i++){
                for(j=0;j<2;j++){
                    
                    ind_r = rand() % indiv_i;
                    all_r = rand() % 2;
                    temp = genolist_new[i][j];
                    genolist_new[i][j] = genolist_new[ind_r][all_r];
                    genolist_new[ind_r][all_r] = temp;
                    
                }
            }
        }
        
        for(i=0;i<indiv_i;i++){
            
            if(poplist_new[i] == pop1){
                
                if(genolist_new[i][0] == minor)
                    minor1++;
                
                if(genolist_new[i][1] == minor)
                    minor1++;
            }
            
            else{
                
                if(genolist_new[i][0] == minor)
                    minor2++;
                
                if(genolist_new[i][1] == minor)
                    minor2++;
            }
        }
        
        freqs.n1 = n1 * 2;
        freqs.n2 = n2 * 2;
        
        freqs.p1 = minor1 / freqs.n1;
        freqs.p2 = minor2 / freqs.n2;
        
        if(isnan(freqs.p1))
            freqs.p1 = 0;
        
        if(isnan(freqs.p2))
            freqs.p2 = 0;
        
        for(i=0;i<n+3;i++)
            free(genolist_new[i]);
        
        free(genolist_new);
        free(poplist_new);
        
        return freqs;
    }

    freqs.p1 = 0.0/0.0;
    freqs.p2 = 0.0/0.0;
    
    return freqs;
}

Freq_s freqLikes(Geno_s sites, int poplist[], double maf, int min, int indiv_n, int pop1, int pop2, int perm){
    
    int i=0, j=0, *poplist_new, minor=-1, indiv_i=0, ind_r=0, temp=0;
    double *likelist_new, alt=0, ref=0, pop1_alt=0, pop2_alt=0, pop1_ref=0, pop2_ref=0, n1=0, n2=0, n=0, minor1=0, minor2=0, p1=0, p2=0;
    Freq_s freqs;
    
    for(i=0;i<indiv_n;i++){
        
        //Minor allele is defined as the less frequent allele across all three populations
        
        ref = ref + (2 * sites.likelist[i][0] + sites.likelist[i][1]);
        
        alt = alt + (sites.likelist[i][1] + 2 * sites.likelist[i][2]);
        
        if(poplist[i] == pop1 && (sites.likelist[i][0] != sites.likelist[i][1] || sites.likelist[i][1] !=  sites.likelist[i][2])){
            
            n1++;
            
            pop1_ref = pop1_ref + (2 * sites.likelist[i][0] + sites.likelist[i][1]);
            
            pop1_alt = pop1_alt + (sites.likelist[i][1] + 2 * sites.likelist[i][2]);
        }
        
        else if(poplist[i] == pop2 && (sites.likelist[i][0] != sites.likelist[i][1] || sites.likelist[i][1] !=  sites.likelist[i][2])){
            
            n2++;
            
            pop2_ref = pop2_ref + (2 * sites.likelist[i][0] + sites.likelist[i][1]);
            
            pop2_alt = pop2_alt + (sites.likelist[i][1] + 2 * sites.likelist[i][2]);
        }
    }
    
    n = n1 + n2;
    
    if(ref >= alt)
        minor = 1;
    else
        minor = 0;
    
    //Sites not meeting the MAF cutoff are excluded
    if(n1 >= min && n2 >= min && ((minor == 1 && (pop1_alt+pop2_alt)/(n*2) >= maf) || (minor == 0 && (pop1_ref+pop2_ref)/(n*2) >= maf))){
        
        if((likelist_new = malloc((n+3)*sizeof(double))) == NULL){
            puts(merror);
            exit(EXIT_FAILURE);
        }
        
        if((poplist_new = malloc((n+3)*sizeof(int))) == NULL){
            puts(merror);
            exit(EXIT_FAILURE);
        }
        
        for(i=0;i<indiv_n;i++){
            
            if((poplist[i] == pop1 || poplist[i] == pop2) && (sites.likelist[i][0] != sites.likelist[i][1] || sites.likelist[i][1] !=  sites.likelist[i][2])){
                
                likelist_new[indiv_i] =  2 * sites.likelist[i][0] + sites.likelist[i][1];

                if(poplist[i] == pop1)
                    poplist_new[indiv_i] = pop1;
                else
                    poplist_new[indiv_i] = pop2;
                
                indiv_i++;
            }
        }
        
        //Shuffles the alleles for permutation testing
        if(perm == 1){
            
            for(i=0;i<indiv_i;i++){
                    
                ind_r = rand() % indiv_i;
                temp = likelist_new[i];
                likelist_new[i] = likelist_new[ind_r];
                likelist_new[ind_r] = temp;
            }
        }
        
        for(i=0;i<indiv_i;i++){
            
            if(poplist_new[i] == pop1){
                
                if(minor == 0)
                    minor1 = minor1 + likelist_new[i];
                else
                    minor1 = minor1 + (2 - likelist_new[i]);
            }
            
            else{
                
                if(minor == 0)
                    minor2 = minor2 + likelist_new[i];
                else
                    minor2 = minor2 + (2 - likelist_new[i]);
            }
        }
        
        freqs.n1 = n1 * 2;
        freqs.n2 = n2 * 2;
        
        freqs.p1 = minor1 / freqs.n1;
        freqs.p2 = minor2 / freqs.n2;
        
        if(isnan(freqs.p1))
            freqs.p1 = 0;
        
        if(isnan(freqs.p2))
            freqs.p2 = 0;
        
        free(likelist_new);
        free(poplist_new);
        
        return freqs;
    }
    
    freqs.p1 = 0.0/0.0;
    freqs.p2 = 0.0/0.0;
    
    return freqs;
}

void printOutput(FILE *out_file, Pbs_s **estimates, int win, int perm, int sims){
    
    int i=0, j=0, length=0, print_i=0;
    double mid=0;
    
    if(win > 1 && perm > 0 || sims > 0)
        fprintf(out_file, "Chromo\tBeginning\tMiddle\tEnd\tLength\tPBS1\tP1\tPBS2\tP2\tPBS3\tP3\n");
    
    else if(win > 1 && perm == 0 && sims == 0)
        fprintf(out_file, "Chromo\tBeginning\tMiddle\tEnd\tLength\tPBS1\tPBS2\tPBS3\n");
    
    else if(win == 1 && perm > 0 || sims > 0)
        fprintf(out_file, "Chromo\tPosition\tPBS1\tP1\tPBS2\tP2\tPBS3\tP3\n");
    
    else
        fprintf(out_file, "Chromo\tPosition\tPBS1\tPBS2\tPBS3\n");
    
    
    for(i=0;strcmp(estimates[0][i].chromo,empty)!=0;i++){
        
        if(isnan(estimates[0][i].pbs) == 0 && isnan(estimates[1][i].pbs) == 0 && isnan(estimates[2][i].pbs) == 0){
            if(isinf(estimates[0][i].pbs) == 0 && isinf(estimates[1][i].pbs) == 0 && isinf(estimates[2][i].pbs) == 0){
    
                //Prints position info to file
                if(win > 1){
                 
                    mid = ((double)estimates[0][i].start + (double)estimates[0][i].stop) / 2;
                     
                    length = estimates[0][i].stop - estimates[0][i].start;
                     
                    fprintf(out_file, "%s\t%i\t%.lf\t%i\t%i\t", estimates[0][i].chromo, estimates[0][i].start, mid, estimates[0][i].stop, length);
                 }
                 
                else
                    fprintf(out_file, "%s\t%i\t", estimates[0][i].chromo, estimates[0][i].start);
                
                //Prints PBS estimates to file
                for(j=0;j<3;j++){
                    fprintf(out_file, "%f\t", estimates[j][i].pbs);
                 
                    //If permutations or simulations were used, prints P-values to file
                    if(perm > 0 || sims > 0)
                        fprintf(out_file, "%f\t", estimates[j][i].p);
                }
                 
                fprintf(out_file, "\n");
                
                print_i++;
            }
        }
    }
    
    for(j=0;j<4;j++)
        free(estimates[j]);
    
    free(estimates);
    
    fclose(out_file);
     
    if(print_i == 0)
        puts("Run finished, but no records were printed to file. Are the input parameters correct?\n");
    else
        printf("Done!\nPrinted %i records to file\n", print_i);
}

void lineTerminator(char *line){
    
    int i=0;
    
    for(i=0;line[i]!=0;i++){
        if(line[i] == '\n' || line[i] == '\r')
            line[i] = '\0';
    }
}

void showHelp(void){
    
    puts("For bug reports, questions and comments, please email: tuomas.hamala@gmail.com\n");
    puts("Required parameters:");
    puts("\t-likes [file]*\tGenotype likelihoods in Beagle format");
    puts("\t-vcf [file]*\tGenotype calls in VCF 4.1 format");
    puts("\t-vcfp [file]*\tSame as '-vcf', except prints out a native PaBScan file (suffix .pabscan)");
    puts("\t-in [file]*\tGenotype calls in native PaBScan format");
    puts("\t-pop1 [file]\tList of individuals (one per line) from focal population one");
    puts("\t-pop2 [file]\tList of individuals (one per line) from focal population two");
    puts("\t-pop3 [file]\tList of individuals (one per line) from the outgroup population");
    puts("\t-out [string]\tName for the output file (suffix .pbs)");
    puts("\n\t*one of these is required");
    puts("\nOptional parameters:");
    puts("\t-win [int]\tSNP based window size (default is 1)");
    puts("\t-step [int]\tSNP based step size (default is 1)");
    puts("\t-ms [file]\tSimulated neutral data in ms format");
    puts("\t-msp [file]\tSame as '-ms', except prints the null-distributions to file (suffix .nulldist)");
    puts("\t-mc [int]\tNumber of permutation cycles for Monte Carlo testing (cannot be used with '-ms')");
    puts("\t-div [int]\tDivergence measure: [0] Hudsons's Fst [1] Weir & Cockerham's Fst [2] Dxy (default is 0)");
    puts("\t-min [int]\tMinimum number of individuals required per population (default is 1)");
    puts("\t-maf [double]\tMinimum minor allele frequency required per site (default is 0.01)");
    puts("\nUsage example:");
    puts("./pabscan -likes example.likes -pop1 list1.txt -pop2 list2.txt -pop3 list3.txt -out output -win 50 -step 51 -ms example.ms -div 2 -min 5 -maf 0.05\n");
}
