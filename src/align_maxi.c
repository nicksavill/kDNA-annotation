// gcc -O3 -o ../kDNA_annotation/align_maxi align_maxi.c -fopenmp -lyaml
// run mRNA_convert.py before running this program

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <yaml.h>

#define MIN(a, b) (a) < (b)? (a): (b)
#define NMRNA 30

short **shortmatrix(int sizex, int sizey)
{
  int i;
  short **m = (short**)calloc(sizex, sizeof(short*));
  m[0] = (short*)calloc(sizex * sizey, sizeof(short));
  for (i = 1; i < sizex; i++)
        m[i] = m[i-1] + sizey;
  return m;
}

char* rev(char *dest, char *src, int n) {
    int i;
    for (i = 0; i < n; i++)
        dest[n-1-i] = src[i];
    dest[n] = '\0';
    return dest;
}

char* comp(char *dest, char *src, int n) {
    int i;
    for (i = 0; i <= n; i++) {
        if      (src[i] == 'A') dest[i] = 'T';
        else if (src[i] == 'T') dest[i] = 'A';
        else if (src[i] == 'G') dest[i] = 'C';
        else if (src[i] == 'C') dest[i] = 'G';
        else                    dest[i] = src[i];
    }
    return dest;
}

int match(char a, char b) {
    return (a == b || a == 'G' && b == 'A' || a == 'T' && b == 'C');
}

void record(char *mRNA_seq, char *frag_seq, FILE *fp, short score, char *maxicircle_id, char *mRNA_name, char direction, int min_gRNA_length) {
    char alignA[score+1], alignB[score+1];
    char alignAp[score+1], alignBp[score+1];
    int i, j, length, mm[score], nmm = 1;
    int mm_allowed = 3;

    strncpy(alignA, mRNA_seq, score);
    strncpy(alignB, frag_seq, score);
    alignA[score] = '\0';
    alignB[score] = '\0';

    // set up positions of mismatches along this record
    mm[0] = -1;
    for (i = 0; i < score; i++)
        if (!match(alignA[i], alignB[i])) mm[nmm++] = i;
    mm[nmm++] = i;

#   pragma omp critical (output)
    {
        for (i = 0; i <= mm_allowed; i++)
            for (j = 0; j < nmm-1-i; j++) {
                length = mm[i+j+1]-(mm[j]+1);
                if (length >= min_gRNA_length) {
                    strncpy(alignAp, alignA+mm[j]+1, length);
                    strncpy(alignBp, alignB+mm[j]+1, length);
                    alignAp[length] = '\0';
                    alignBp[length] = '\0';

                    fprintf(fp, "%s\n", maxicircle_id);
                    fprintf(fp, "%s %s %s %c\n", mRNA_name, alignAp, alignBp, direction);
                    fflush(fp);
                }
            }
    }
}

void align(char *maxicircle_seq, int n, int nmRNA, ssize_t nbases[], char *mRNA_seq[], char mRNA_name[][100],
            FILE *fp, char *maxicircle_id, int min_gRNA_length, char direction) {
    int i, j, k, m, score, lenfrag, row, col, len_diag;
    char *fragment;
    short **c;

    lenfrag = 120;

#   pragma omp parallel for \
    firstprivate(maxicircle_seq, n, nbases, mRNA_seq, mRNA_name, fp, maxicircle_id, direction, lenfrag) \
    private(fragment, c, i, j, m, score, row, col, len_diag)

    for (k = 0; k < n-lenfrag/2; k += lenfrag/2) {
        fragment = maxicircle_seq+k;
        if (lenfrag > n-k) lenfrag = n-k;

        for (m = 0; m < nmRNA; m++) {
            c = shortmatrix(nbases[m]+1, lenfrag+1);

            // construct scoring matrix
            for (i = nbases[m]-1; i >= 0; i--)
                for (j = lenfrag-1; j >= 0; j--) {
                    score = c[i+1][j+1];
                    if (mRNA_seq[m][i] == 'N' || fragment[j] == 'N')
                        continue;

                    if (match(mRNA_seq[m][i], fragment[j]))
                        // if match always add 1
                        c[i][j] = score+1;
                    else if (i < nbases[m]-1 && j < lenfrag-1) {
                        // we may have an internal mismatch
                        if (match(mRNA_seq[m][i+1], fragment[j+1]))
                            // if last was a match score 1
                            c[i][j] = score+1;
                        else
                            // else last was a mismatch, set that score to 0
                            c[i+1][j+1] = 0;
                    }

                }

            // search for all alignments of lengths >= min_gRNA_length
            for (row = 0; row < nbases[m]; row++) {
                len_diag = MIN(nbases[m]-row, lenfrag);
                if (len_diag >= min_gRNA_length) {
                    i = 0;
                    while (i < len_diag) {
                        score = c[row+i][i];
                        if (score == 0)
                            i++;
                        else {
                            if (min_gRNA_length <= score)
                                record(mRNA_seq[m]+row+i, fragment+i, fp, score, maxicircle_id, mRNA_name[m]+1, direction, min_gRNA_length);
                            i += score;
                        }
                    }
                }
            }
            for (col = 1; col < lenfrag; col++) {
                len_diag = MIN(nbases[m], lenfrag-col);
                if (len_diag >= min_gRNA_length) {
                    i = 0;
                    while (i < len_diag) {
                        score = c[i][col+i];
                        if (score == 0)
                            i++;
                        else {
                            if (min_gRNA_length <= score)
                                record(mRNA_seq[m]+i, fragment+col+i, fp, score, maxicircle_id, mRNA_name[m]+1, direction, min_gRNA_length);
                            i += score;
                        }
                    }
                } 
            }
            free(c[0]);
            free(c);
        }
    }        
}

char* make_filename(char* name1, char* name2, char* name3) {
    size_t l1 = strlen(name1);
    size_t l2 = strlen(name2);
    size_t l3 = strlen(name3);
    char *file = (char*) malloc((l1+l2+l3+3)*sizeof(char));
    file[0] = '\0';
    file = strcat(file, name1);
    file = strcat(file, "/");
    file = strcat(file, name2);
    file = strcat(file, "/");
    file = strcat(file, name3);
    return file;
}

int main(int argc, char *argv[]) {
    int i, count, nmRNA, min_gRNA_length;
    char *mRNA_seq[] = {[0 ... NMRNA] = NULL}, mRNA_name[NMRNA][100];
    char maxicircle_seq[2<<17] = "", maxicircle_id[2<<10], *maxicircle_seq_rc;
    char *buf = malloc(sizeof(char));
    char *project_dir;
    char *work_dir;
    char *edited_mRNA_file;
    char *maxicircle_file;
    char *alignment_file;
    char *min_gRNA_length_str;
    ssize_t n = 1, nbases[NMRNA] = {0};
    FILE *fh;

    if (argc == 1)
        fh = fopen("config.yaml", "r");
    else
        fh = fopen(argv[1], "r");

    yaml_parser_t parser;
    yaml_token_t  token;

    // Initialize parser
    if (!yaml_parser_initialize(&parser)) {
        fputs("Failed to initialize parser!\n", stderr);
        exit(0);
    }
    if (fh == NULL) {
        fputs("Failed to open file!\n", stderr);
        exit(0);
    }

    // Set input file
    yaml_parser_set_input_file(&parser, fh);
 
    int state = 0;
    char** datap = NULL;
    char* tk;
    do {
        yaml_parser_scan(&parser, &token);
        switch(token.type) {
            case YAML_KEY_TOKEN:
                state = 0; 
                break;
            case YAML_VALUE_TOKEN:   
                state = 1; 
                break;
            case YAML_SCALAR_TOKEN:
                tk = token.data.scalar.value;
                if (state == 0) {
                    // we have a key token, setup datap to point to a variable
                    // if this token is not needed data is NULL 
                    if (!strcmp(tk, "project"))
                        datap = &project_dir;
                    else if (!strcmp(tk, "working directory"))
                        datap = &work_dir;
                    else if (!strcmp(tk, "edited mRNA fasta file"))
                        datap = &edited_mRNA_file;
                    else if (!strcmp(tk, "maxicircle clean fasta file"))
                        datap = &maxicircle_file;
                    else if (!strcmp(tk, "maxicircle alignments file"))
                        datap = &alignment_file;
                    else if (!strcmp(tk, "minimum gRNA length"))
                        datap = &min_gRNA_length_str;
                    else
                        datap = NULL;
                } else {
                    // state=1, we have a value token
                    // if datap is NULL ignore otherwise get its value
                    if (datap != NULL) {
                        *datap = strdup(tk);
                    }
                }
                break;
            default: 
                break;
        }
        if (token.type != YAML_STREAM_END_TOKEN)
            yaml_token_delete(&token);
    } while(token.type != YAML_STREAM_END_TOKEN);
    yaml_token_delete(&token);
    yaml_parser_delete(&parser);
    fclose(fh);

    if (project_dir == NULL) {
        printf("project not found in config file\n");
        exit(0);
    }
    if (work_dir == NULL) {
        printf("working directory not found in config file\n");
        exit(0);
    }
    if (edited_mRNA_file == NULL) {
        printf("edited mRNA fasta file not found in config file\n");
        exit(0);
    }
    if (maxicircle_file == NULL) {
        printf("maxicircle clean fasta file not found in config file\n");
        exit(0);
    }
    if (alignment_file == NULL) {
        printf("minimum gRNA length not found in config file\n");
        exit(0);
    }
    if (min_gRNA_length_str == NULL) {
        printf("not found\n");
        exit(0);
    }

    min_gRNA_length = atoi(min_gRNA_length_str);

    char *edited_mRNA_filename = make_filename(project_dir, work_dir, edited_mRNA_file);
    char *maxicircle_filename = make_filename(project_dir, work_dir, maxicircle_file);
    char *alignment_filename = make_filename(project_dir, work_dir, alignment_file);

    FILE *fpa = fopen(edited_mRNA_filename, "r");
    if (!fpa) {
        fprintf(stderr, "maxi_align: File does not exist: %s\n", edited_mRNA_filename);
        exit(0);
    }
    FILE *fpb = fopen(maxicircle_filename, "r");
    if (!fpb) {
        fprintf(stderr, "maxi_align: File does not exist: %s\n", maxicircle_filename);
        exit(0);
    }
    FILE *fpc  = fopen(alignment_filename, "w");
    if (!fpc) {
        fprintf(stderr, "maxi_align: Cannot open file: %s\n", alignment_filename);
        exit(0);
    }

    nmRNA = 0;
    while(fscanf(fpa, "%s\n", mRNA_name[nmRNA]) != EOF) {
        nbases[nmRNA] = getline(&mRNA_seq[nmRNA], &nbases[nmRNA], fpa)-1;
        mRNA_seq[nmRNA][nbases[nmRNA]] = '\0';
        nmRNA++;
    }
    fclose(fpa);
    printf("Number of mRNAs = %d\n", nmRNA);

    while ((n = getline(&buf, &n, fpb)) != -1) {
        buf[--n] = '\0';

        if (buf[0] == '>')
            strcpy(maxicircle_id, buf);
        else
            strcat(maxicircle_seq, buf);
    }   
    fclose(fpb);

    n = strlen(maxicircle_seq);
    maxicircle_seq_rc = malloc((n+1)*sizeof(char));

    align(maxicircle_seq, n, nmRNA, nbases, mRNA_seq, mRNA_name, fpc, maxicircle_id, min_gRNA_length, 't');

    comp(maxicircle_seq_rc, maxicircle_seq, n);
    rev(maxicircle_seq, maxicircle_seq_rc, n);
    align(maxicircle_seq, n, nmRNA, nbases, mRNA_seq, mRNA_name, fpc, maxicircle_id, min_gRNA_length, 'c');

    free(maxicircle_seq_rc);

    return 0;
}
