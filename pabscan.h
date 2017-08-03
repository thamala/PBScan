//Version 2017.07.27

typedef struct{
    int position, **genolist;
    double **likelist;
    char chromo[21];
}Geno_s;

typedef struct{
    int start, stop;
    double pbs, p;
    char chromo[21];
}Pbs_s;

typedef struct{
    double n1, n2, p1, p2;
}Freq_s;

typedef struct{
    double hw, hb, dxy;
}Div_s;

void openFiles(int argc, char *argv[]);
char **readPop(FILE *pop_file, int *n);
Geno_s *readVcf(FILE *geno_file, char **pop1, char **pop2, char **pop3, char name[], int indiv_n, int *snp_i, int *poplist, int print);
Geno_s *readNative(FILE *geno_file, char **pop1, char **pop2, char **pop3, int indiv_n, int *snp_i, int *poplist);
Geno_s *readBeagle(FILE *geno_file, char **pop1, char **pop2, char **pop3, int indiv_n, int *snp_i, int *poplist);
double **processMs(FILE *ms_file, char name[], int pop1_n, int pop2_n, int pop3_n, int *sims, int div, int msp, double maf);
Pbs_s **estimatePBS(Geno_s *sites, int poplist[], double **ms, int sims, int snp_n, int indiv_n, int win, int step, int perm, int div, int min, int like, double maf);
Div_s estimateDiv(Geno_s sites, int poplist[], double maf, int min, int indiv_n, int div, int like, int pop1, int pop2, int perm);
Freq_s freqCalls(Geno_s sites, int poplist[], double maf, int min, int indiv_n, int pop1, int pop2, int perm);
Freq_s freqLikes(Geno_s sites, int poplist[], double maf, int min, int indiv_n, int pop1, int pop2, int perm);
void printOutput(FILE *out_file, Pbs_s **estimates, int win, int perm, int sims);
void lineTerminator(char *line);
void showHelp(void);
