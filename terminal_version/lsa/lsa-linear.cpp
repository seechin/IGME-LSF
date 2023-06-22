#include    <stdio.h>
#include    <math.h>
#include    <string.h>
#include    <unistd.h>
#include    <stdlib.h>
#include    <sys/types.h>
#include    <sys/stat.h>
#include    <fcntl.h>
#include    <ctype.h>
#include    <pthread.h>
#include    "String.cpp"
#include    "matrix.cpp"

template <class Type> class Array { // return: number of pages
  #define ArrayPageSize 10000
  public: int n, nmax; Type * a;
    Array(){ n = nmax = 0; a = NULL; }
    int insert(Type e){
        if (!a || n >= nmax){ nmax += ArrayPageSize;
            Type * an = (Type*) malloc(sizeof(Type) * nmax);
            for (int i=0; i<n; i++) an[i] = a[i];
            free(a); a = an;
        }
        a[n++] = e;
        return nmax / ArrayPageSize;
    }
    void init(int nl){
        if (a) free(a); nmax = nl; n = nl;
        a = nmax==0? NULL : (Type*) malloc(sizeof(Type) * nmax);
    }
    void dispose(){ if (a) free(a); n = nmax = 0; a = NULL; }
};

const char * software_name = "lsa-linear";
const char * szHelpBrief = "%s [[-s] file] [-x colx_list -y coly]\n";
const char * szHelp = "%s (c) Cao Siqin 2015.5.23-2022.1.31\n\
    least square analysis: multi variable fitting\n\
Parameters:\n\
   [-s]     a file of training set, can be stdin/con\n\
            can be: -s-tab, -s-csv. -s = -s-tab\n\
            can be -S... to force show equation and details\n\
    -f      a file of the test set, can be ::input\n\
            can be: -f-tab, -f-csv. -f = -f-tab\n\
            can be -F... to show independent variables\n\
    -x      col of the variable list, default: -x 0\n\
            col 0 means a constant\n\
    -y      specify the y col, default: -y 1\n\
  -%%...     output format of data\n\
Example:\n\
  cat data | awk '{printf(\"%%f,%%f,%%f\\n\",$1,$1*$1,$2)}' |\n\
    lsa-linear -scsv -x 0 1 2 -y 3 -f ::input # y=A+Bx+Cx^2 fitting\n\
";

// -single, -double  set the output precision\n

#define MAX_FIT_TERMS 1000

double mathpower(double x, int n){
    switch (n){
      case 0: return 1;
      case 1: return x;
      case 2: return x*x;
      case 3: return x*x*x;
      default: return pow(x, n);
    }
}
void print_double(double x, char * term=" ", FILE * dev=stdout){
    fprintf(dev, fabs(x)<1e-5?"%11.4e%s":"%11f%s", x, term);
}

const char * string_nonspace_begin(const char * src){
    const char * ret = src; while (*ret && *ret==' ') ret ++;
    return ret;
}

int main(int argc, char *argv[]){
    int error = 0; bool analysis_csv = false; bool analysis_f_csv = false; bool show_testing_detail = false; bool show_equation = false; bool show_extra = false;
    bool single_repc = false;
    char * szfile = NULL; const char * szffile = NULL; const char * output_format = "%.14g"; int i_output_format_option = -1;
    int colx = -1; int coly = 0; int output = 0;
    int nlist[MAX_FIT_TERMS]; int nnlist = 1; memset(nlist, 0, sizeof(nlist)); nlist[0] = 1;
    if (argc > 1){
      for (int i = 1; i < argc; i ++){
        if (argv[i][0] == '-'){
            if (StringNS::string(argv[i]) == "-?"){
                printf(szHelpBrief, software_name);
                error = 1;
            } else if (StringNS::string(argv[i]) == "--?" || StringNS::string(argv[i]) == "-h" || StringNS::string(argv[i]) == "-help" || StringNS::string(argv[i]) == "--help"){
                printf(szHelp, software_name);
                printf(szHelpBrief, software_name);
                error = 1;
            } else if (StringNS::string(argv[i]) == "-x" || StringNS::string(argv[i]) == "--x"){
                nnlist = 0;
                while (i+1<argc && argv[i+1][0]!='-'){
                    i++; StringNS::string nsl[MAX_FIT_TERMS];
                    int nw = StringNS::analysis_line(argv[i], nsl, sizeof(nsl)/sizeof(nsl[0]), true); if (nw > MAX_FIT_TERMS-nnlist) nw = MAX_FIT_TERMS-nnlist;
                    if (nw>0 && is_string_number(nsl[0])){
                        for (int i=0; i<nw; i++) nlist[i+nnlist] = atoi(nsl[i].text) - 1;
                        nnlist += nw;
                    }
                    //printf("nlist[%d]:", nnlist); for (int i=0; i<nnlist; i++) printf(" %d", nlist[i]); printf("\n");
                }
            } else if ((StringNS::string(argv[i], 2) == "-y" && argv[i][2]>='0'&&argv[i][3]<='9') || (StringNS::string(argv[i], 3) == "--y" && argv[i][3]>='0'&&argv[i][3]<='9')){
                char * col_tx = (StringNS::string(argv[i], 2) == "-y")? &argv[i][2] : &argv[i][3];
                coly = atoi(col_tx) - 1;
            } else if (StringNS::string(argv[i]) == "-y" || StringNS::string(argv[i]) == "--y"){
                if (i+1<argc && argv[i+1][0] != '-'){ i++; coly = atoi(argv[i]) - 1; }
            } else if (StringNS::string(argv[i]) == "-s" || StringNS::string(argv[i]) == "--s" || StringNS::string(argv[i]) == "-s-tab" || StringNS::string(argv[i]) == "--s-tab" || StringNS::string(argv[i]) == "-stab" || StringNS::string(argv[i]) == "--stab"){
                analysis_csv = false;
                if (argv[i][1]=='S' || argv[i][2]=='S'){
                    show_equation = true;
                    show_extra = true;
                } else {
                    show_equation = false;
                    show_extra = false;
                }
                if (i+1<argc && argv[i+1][0] != '-'){ i++;
                    if (!szfile) szfile = argv[i]; else { fprintf(stderr, "%s [%d]: warning : -s %s already specified.\n", software_name, i, szfile); }
                }
            } else if (StringNS::string(argv[i]) == "-s-csv" || StringNS::string(argv[i]) == "--s-csv" || StringNS::string(argv[i]) == "-scsv" || StringNS::string(argv[i]) == "--scsv"){
                analysis_csv = true;
                if (argv[i][1]=='S' || argv[i][2]=='S') show_equation = true; else show_equation = false;
                if (i+1<argc && argv[i+1][0] != '-'){ i++;
                    if (!szfile) szfile = argv[i]; else { fprintf(stderr, "%s [%d]: warning : -s %s already specified.\n", software_name, i, szfile); }
                }
            } else if (StringNS::string(argv[i]) == "-f" || StringNS::string(argv[i]) == "--f" || StringNS::string(argv[i]) == "-f-tab" || StringNS::string(argv[i]) == "--f-tab" || StringNS::string(argv[i]) == "-ftab" || StringNS::string(argv[i]) == "--ftab"){
                analysis_f_csv = false;
                if (argv[i][1]=='F' || argv[i][2]=='F') show_testing_detail = true; else show_testing_detail = false;
                if (!szffile) szffile = "::input"; if (i+1<argc && argv[i+1][0] != '-'){ i++; szffile = argv[i]; }
            } else if (StringNS::string(argv[i]) == "-f-csv" || StringNS::string(argv[i]) == "--f-csv" || StringNS::string(argv[i]) == "-fcsv" || StringNS::string(argv[i]) == "--fcsv"){
                analysis_f_csv = true;
                if (argv[i][1]=='F' || argv[i][2]=='F') show_testing_detail = true; else show_testing_detail = false;
                if (!szffile) szffile = "::input"; if (i+1<argc && argv[i+1][0] != '-'){ i++; szffile = argv[i]; }
            } else if (StringNS::string(argv[i]) == "-single" || StringNS::string(argv[i]) == "--single"){
                single_repc = true;
            } else if (StringNS::string(argv[i]) == "-double" || StringNS::string(argv[i]) == "--double"){
                single_repc = false;
            } else if (StringNS::string(argv[i], 2) == "-%"){
                output_format = &argv[i][1];
                i_output_format_option = i;
            } else if (StringNS::string(argv[i], 3) == "--%"){
                output_format = &argv[i][2];
                i_output_format_option = i;
            } else {
                fprintf(stderr, "%s [%d]: error : unrecognizable option %s\n", software_name, i, argv[i]); error = -1;
            }
       } else {
        if (!szfile) szfile = argv[i]; else { fprintf(stderr, "%s [%d]: warning : file %s already specified.\n", software_name, i, szfile); }
       }
      }
    } else {
        printf(szHelpBrief, software_name);
        error = 1;
    }
    if (error) return error;

    for (int i=0; output_format[i]; i++){
        if (output_format[i]=='s' || output_format[i]=='S' || output_format[i]=='n' || output_format[i]=='N'){
            fprintf(stderr, "%s [%d] : possible fault detected: %s\n", software_name, i_output_format_option+1, i_output_format_option>=0?argv[i_output_format_option]:output_format);
            return -1;
        } else if (output_format[i]=='f' || output_format[i]=='F' || output_format[i]=='g' || output_format[i]=='G' || output_format[i]=='e' || output_format[i]=='E'){ break;
        } else if (output_format[i]=='p' || output_format[i]=='P'){ break;
        } else if (output_format[i]=='d' || output_format[i]=='D' || output_format[i]=='i' || output_format[i]=='I' || output_format[i]=='u' || output_format[i]=='U' || output_format[i]=='x' || output_format[i]=='X' || output_format[i]=='o' || output_format[i]=='O' || output_format[i]=='c' || output_format[i]=='C' || output_format[i]=='a' || output_format[i]=='A'){
            fprintf(stderr, "%s [%d] : warning : incorrect output format: %s\n", software_name, i_output_format_option+1, i_output_format_option>=0?argv[i_output_format_option]:output_format);
            break;
        }
    }

    for (int i=0; i<nnlist; i++) if (nlist[i]<0) nlist[i] = -1; for (int i=0; i<nnlist; i++){ for (int j=i+1; j<nnlist; j++){ if (nlist[i] == nlist[j]){ for (int k=j; k<nnlist-1; k++) nlist[k] = nlist[k+1]; nnlist --; break; } } }
    //if (!szfile){ fprintf(stderr, "%s : no file to process.\n", software_name); return -1; }
    FILE * file = stdin;
    if (szfile){
        if (StringNS::string(szfile)=="stdin" || StringNS::string(szfile)=="con"){
            file = stdin;
        } else {
            file = fopen(szfile, "r"); if (!file){ fprintf(stderr, "%s : error : can't open %s\n", software_name, szfile); return -1; }
        }
    }

    Array <double> seq[MAX_FIT_TERMS], seqy; for (int i=0; i<MAX_FIT_TERMS; i++) seq[i].init(0); seqy.init(0);
    StringNS::string sl[1000]; char input[4096];
    int maxln = 0; for (int i=0; i<nnlist; i++) if (nlist[i]>maxln) maxln = nlist[i];
    while (fgets(input, sizeof(input), file)){
        int nw = 0;
        if(analysis_csv) nw = StringNS::analysis_csv_line(input, sl, sizeof(sl)/sizeof(sl[0]), true); else nw = StringNS::analysis_line(input, sl, sizeof(sl)/sizeof(sl[0]), true);
        if (nw<1) continue;
        if (sl[0].text[0] == '#'){
            continue;
        } else if (sl[0].text[0] == ';' || sl[0].text[0] == '@'){
            continue;
        } else if (nw > maxln && nw > coly){
            seqy.insert(atof(sl[coly].text));
            for (int icol=0; icol<nnlist; icol++) seq[icol].insert(nlist[icol]<0? 1 : atof(sl[nlist[icol]].text));
        }
    }
    fclose(file);
    //for (int i=0; i<seq[0].n; i++){ printf("line %d : %12f :", i, seqy.a[i]); for (int ic=0; ic<nnlist; ic++) printf(" %12f", seq[ic].a[i]); printf("\n"); } return 0;

    MatrixNS::Matrix A, Z, M, Mi, Mit;
    A.init(nnlist); Z.init(nnlist); M.init(nnlist); Mi.init(nnlist); Mit.init(nnlist);

    Z=0; M=0; for (int i=0; i<seq[0].n; i++){
        for (int m=0; m<nnlist; m++){
            Z.a[m][0] += seqy.a[i] * seq[m].a[i];
            for (int n=0; n<nnlist; n++){
                M.a[m][n] += seq[m].a[i] * seq[n].a[i];
            }
        }
    }
    if (seq[0].n>0) for (int m=0; m<nnlist; m++){ Z.a[m][0] /= seq[0].n;
        for (int n=0; n<nnlist; n++) M.a[m][n] /= seq[0].n;
    }

    if (seq[0].n<=M.n){
        fprintf(stderr, "%s : error : data (%d) less than required (%d)\n", software_name, seq[0].n, M.n);
    } else if (!M.inverse(Mi, Mit)){
        fprintf(stderr, "%s : error : matrix irreversible\n", software_name);
    } else {
        A = 0; Mi.product(Z, A);
        bool show_fit_result = false;
        if (szffile){
            if (StringNS::string(szffile) == "::input" || StringNS::string(szffile) == ":input"){
                show_fit_result = true;
            } else {
                FILE * ffile = fopen(szffile, "r");
                if (ffile){
                    for (int i=0; i<nnlist; i++) seq[i].n = 0; seqy.n = 0;
                    while (fgets(input, sizeof(input), ffile)){
                        int nw = 0;
                        if(analysis_f_csv) nw = StringNS::analysis_csv_line(input, sl, sizeof(sl)/sizeof(sl[0]), true); else nw = StringNS::analysis_line(input, sl, sizeof(sl)/sizeof(sl[0]), true);
                        if (nw<1) continue;
                        if (sl[0].text[0] == '#'){
                            continue;
                        } else if (sl[0].text[0] == ';' || sl[0].text[0] == '@'){
                            continue;
                        } else if (nw > maxln && nw > coly){
                            seqy.insert(atof(sl[coly].text));
                            for (int icol=0; icol<nnlist; icol++) seq[icol].insert(nlist[icol]<0? 1 : atof(sl[nlist[icol]].text));
                        }
                    }
                    fclose(ffile);
                    show_fit_result = true;
                } else { fprintf(stderr, "%s : warning : cannot open -f %s\n", software_name, szffile); }
            }
        }
        if (!show_fit_result || show_equation){
            /*double SSReg = 0; double SST = 0;
            double y_average = 0; for (int i=0; i<seq[0].n; i++) y_average += seqy.a[i]/seq[0].n;
            for (int i=0; i<seq[0].n; i++){
                double y = seqy.a[i];
                double x = 0; for (int icol=0; icol<nnlist; icol++) x += A.a[icol][0] * seq[icol].a[i];
                SSReg += (x - y_average)*(x - y_average);
                SST += (y - y_average)*(y - y_average);
            }
            double r = sqrt(fabs(SSReg/SST));*/

            double sxx, sx, syy, sy, sxy; sxx = sx = syy = sy = sxy = 0; int count = 0;
            for (int i=0; i<seq[0].n; i++){
                double y = seqy.a[i];
                double x = 0; for (int icol=0; icol<nnlist; icol++) x += A.a[icol][0] * seq[icol].a[i];
                sxx += x*x; sx += x; syy += y*y; sy += y; sxy += x*y;
                count ++;
            }
            double r = (sxy - sx*sy/count) / sqrt((syy - sy*sy/count) * (sxx - sx*sx/count));
            double rmse = sqrt(fabs(sxx+syy-2*sxy)/count);


            FILE * o = stdout; if (show_fit_result) o = stderr;

          // show a human understandable equation
            fprintf(o, "#col(%d) =", coly+1); char buffer_tmp[512]; memset(buffer_tmp, 0, sizeof(buffer_tmp));
            bool first_term = true;
            for (int i=0; i<nnlist; i++){
              if (fabs(A.a[i][0]) > (single_repc?1e-7:1e-12)){
                //
                if (fabs(fabs(A.a[i][0])-1) < (single_repc?1e-7:1e-12)){ buffer_tmp[0] = 0;
                } else {
                    snprintf(buffer_tmp, sizeof(buffer_tmp), output_format, fabs(A.a[i][0]));
                }
                if (nlist[i]>=0) sprintf(&buffer_tmp[strlen(buffer_tmp)], buffer_tmp[0]? " col(%d)" : "col(%d)", nlist[i]+1);
                fprintf(o, " %s%s", A.a[i][0]<0? (first_term? "-" : "- ") : (first_term? "" : "+ "), string_nonspace_begin(buffer_tmp));
                first_term = false;



                /*
                snprintf(buffer_tmp, sizeof(buffer_tmp), output_format, fabs(A.a[i][0]));
                if (nlist[i]>=0) sprintf(&buffer_tmp[strlen(buffer_tmp)], buffer_tmp[0]? " col(%d)" : "col(%d)", nlist[i]+1);
                if (buffer_tmp[0]!=' '){
                    fprintf(o, " %s%s", A.a[i][0]<0? "-" : (first_term? "" : "+"), buffer_tmp);
                } else {
                    int k = 0; for (k=0; buffer_tmp[k]&&buffer_tmp[k]==' '; k++);
                    if (k>0) k--; buffer_tmp[k] = A.a[i][0]<0? '-' : (first_term? ' ' : '+');
                    fprintf(o, " %s", buffer_tmp);
                }
                //fprintf(o, " %s%s", A.a[i][0]<0? (first_term? "-" : "- ") : (first_term? "" : "+ "), buffer_tmp);
                first_term = false;
                */
              }
            }
            snprintf(buffer_tmp, sizeof(buffer_tmp), output_format, r); fprintf(o, " , R= %s", string_nonspace_begin(buffer_tmp));
            if (show_extra){
                fprintf(o, "\n### ");
                double s = count>2?sqrt(fabs((syy - sxy*sxy/sxx)/(count-2))) : -1;
                //snprintf(buffer_tmp, sizeof(buffer_tmp), output_format, s); fprintf(o, " s= %s", string_nonspace_begin(buffer_tmp));
                double SE_a = s * sqrt(1.0/count + (sx/count)*(sx/count)/sxx);
                snprintf(buffer_tmp, sizeof(buffer_tmp), output_format, SE_a); fprintf(o, " SE(intercept)= %s", string_nonspace_begin(buffer_tmp));
                double SE_b = s / sqrt(sxx);
                snprintf(buffer_tmp, sizeof(buffer_tmp), output_format, SE_b); fprintf(o, " SE(slope)= %s", string_nonspace_begin(buffer_tmp));
            }
            fprintf(o, "\n");
          // at full mode: show the list of fitting parameters and more information
            if (show_fit_result){
                fprintf(o, "  #a[]  =");
                for (int i=0; i<nnlist; i++){
                    fprintf(o, " "); fprintf(o, output_format, A.a[i][0]);
                }
                fprintf(o, "\n  #R    = "); fprintf(o, output_format, r);
                fprintf(o, "\n  #RMSE = "); fprintf(o, output_format, rmse);
                fprintf(o, "\n");
            }
        }
        if (show_fit_result){
            double sxyr = 0; double sxy = 0; double syy = 0; double sy = 0; int n = 0;
            for (int i=0; i<seq[0].n; i++){
                double y = seqy.a[i];
                double x = 0; for (int icol=0; icol<nnlist; icol++) x += A.a[icol][0] * seq[icol].a[i];
                printf(output_format, y);
                printf(" ");
                printf(output_format, x);
                if (show_testing_detail){
                    printf(" <-");
                    for (int icol=0; icol<nnlist; icol++){
                        printf(" "); printf(output_format, seq[icol].a[i]);
                    }
                }
                printf("\n");
                sxy += (x-y) * (x-y); sxyr += x-y; n ++; sy += y; syy += y*y;
            }
            /*if (show_testing_detail && n>0){
                fprintf(stderr, "# stdev ");
                fprintf(stderr, output_format, sqrt(fabs(sxy/n - sxyr/n*sxyr/n)));
                fprintf(stderr, " (");
                fprintf(stderr, output_format, sy/n);
                fprintf(stderr, "+");
                fprintf(stderr, output_format, sqrt(fabs(syy/n - sy/n*sy/n)));
                fprintf(stderr, ")\n");
            }*/
        }
    }



    for (int i=0; i<MAX_FIT_TERMS; i++) seq[i].dispose(); seqy.dispose();
    A.dispose(); Z.dispose(); M.dispose(); Mi.dispose(); Mit.dispose();
    return 0;
}
