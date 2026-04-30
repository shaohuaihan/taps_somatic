#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>

#define EPS 1e-10

// =====================
// context
// =====================
enum {CPG=0, CHG=1, CHH=2, OTHER=3};

// =====================
// 参数结构
// =====================
typedef struct {

    char tumor_bam[256];
    char normal_bam[256];
    char ref[256];

    double min_vaf;
    int min_depth;
    double llr_thres;
    double dispersion;
    double contamination;

    double strand_art;
    double orient_art;

    double cal_a;
    double cal_b;

} params_t;

// =====================
void init_params(params_t *p) {

    memset(p,0,sizeof(params_t));

    p->min_vaf = 0.01;
    p->min_depth = 10;
    p->llr_thres = 6.0;
    p->dispersion = 50;
    p->contamination = 0.01;

    p->strand_art = 0.9;
    p->orient_art = 0.95;

    p->cal_a = 1.0;
    p->cal_b = 0.0;
}

// =====================
// CLI解析
// =====================
void parse_args(int argc,char *argv[],params_t *p){

    static struct option long_opts[] = {
        {"tumor", required_argument, 0, 't'},
        {"normal", required_argument, 0, 'n'},
        {"ref", required_argument, 0, 'r'},
        {"min-vaf", required_argument, 0, 1},
        {"min-depth", required_argument, 0, 2},
        {"llr", required_argument, 0, 3},
        {"dispersion", required_argument, 0, 4},
        {"contamination", required_argument, 0, 5},
        {0,0,0,0}
    };

    int opt;
    while((opt=getopt_long(argc,argv,"t:n:r:",long_opts,NULL))!=-1){

        switch(opt){

            case 't': strcpy(p->tumor_bam,optarg); break;
            case 'n': strcpy(p->normal_bam,optarg); break;
            case 'r': strcpy(p->ref,optarg); break;

            case 1: p->min_vaf=atof(optarg); break;
            case 2: p->min_depth=atoi(optarg); break;
            case 3: p->llr_thres=atof(optarg); break;
            case 4: p->dispersion=atof(optarg); break;
            case 5: p->contamination=atof(optarg); break;

            default:
                fprintf(stderr,"Usage: -t tumor -n normal -r ref\n");
                exit(1);
        }
    }

    if(strlen(p->tumor_bam)==0 ||
       strlen(p->normal_bam)==0 ||
       strlen(p->ref)==0){

        fprintf(stderr,"Error: missing required arguments\n");
        exit(1);
    }
}

// =====================
// 数学
// =====================
double ll_binomial(int k,int n,double p){
    return k*log(p+EPS)+(n-k)*log(1-p+EPS);
}

double ll_beta_binom(int k,int n,double a,double b){
    return lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1)
        +lgamma(k+a)+lgamma(n-k+b)
        -lgamma(n+a+b)
        -(lgamma(a)+lgamma(b)-lgamma(a+b));
}

double calibrate(double llr,double a,double b){
    return 1.0/(1.0+exp(-a*llr+b));
}

// =====================
// context
// =====================
int get_ctx(char *s){

    if(s[1]!='C') return OTHER;
    if(s[2]=='G') return CPG;
    if(s[3]=='G') return CHG;
    return CHH;
}

// =====================
// stats
// =====================
typedef struct{
    int depth,alt;
    int fwd_alt,rev_alt;
    int f1r2,f2r1;
    double pos_sum;
    int alt_count;
    double alt_w,depth_w;
} stats_t;

double base_err(int Q){return pow(10,-Q/10.0);}
double map_err(int MQ){return pow(10,-MQ/10.0);}

// =====================
void collect(const bam_pileup1_t *plp,int n,char ref,stats_t *s){

    for(int i=0;i<n;i++){

        const bam_pileup1_t *p=&plp[i];

        if(p->is_del||p->is_refskip) continue;

        uint8_t *seq=bam_get_seq(p->b);
        uint8_t *qual=bam_get_qual(p->b);

        char base=seq_nt16_str[bam_seqi(seq,p->qpos)];

        int Q=qual[p->qpos];
        int MQ=p->b->core.qual;

        double w=(1-base_err(Q))*(1-map_err(MQ));

        s->depth++;
        s->depth_w+=1;

        if(base!=ref){

            s->alt++;
            s->alt_w+=w;

            if(bam_is_rev(p->b)) s->rev_alt++;
            else s->fwd_alt++;

            if(p->b->core.flag&BAM_FREAD1) s->f1r2++;
            else s->f2r1++;

            double rel=p->qpos/(double)p->b->core.l_qseq;
            s->pos_sum+=rel;
            s->alt_count++;
        }
    }
}

// =====================
// 学习context error
// =====================
void learn_error(char *bamfile,char *reffile,double ctx_err[4]){

    samFile *fp=sam_open(bamfile,"r");
    bam_hdr_t *hdr=sam_hdr_read(fp);
    bam1_t *b=bam_init1();
    faidx_t *fai=fai_load(reffile);

    long total[4]={0};
    long ct[4]={0};

    while(sam_read1(fp,hdr,b)>=0){

        uint8_t *seq=bam_get_seq(b);
        int len=b->core.l_qseq;
        int pos=b->core.pos;

        for(int i=0;i<len;i++){

            char base=seq_nt16_str[bam_seqi(seq,i)];

            int ref_pos=pos+i;

            int l;
            char *ref=faidx_fetch_seq(fai,
                hdr->target_name[b->core.tid],
                ref_pos-1,ref_pos+2,&l);

            if(!ref||l<4) continue;

            int ctx=get_ctx(ref);

            total[ctx]++;

            if(ref[1]=='C' && base=='T')
                ct[ctx]++;

            free(ref);
        }
    }

    for(int i=0;i<4;i++)
        ctx_err[i]=ct[i]/(double)(total[i]+1);

    fprintf(stderr,
        "[LEARN] CpG=%.4f CHG=%.4f CHH=%.4f\n",
        ctx_err[0],ctx_err[1],ctx_err[2]);
}

// =====================
// 主程序
// =====================
int main(int argc,char *argv[]){

    params_t P;
    init_params(&P);
    parse_args(argc,argv,&P);

    double ctx_err[4]={0};

    learn_error(P.normal_bam,P.ref,ctx_err);

    samFile *tfp=sam_open(P.tumor_bam,"r");
    samFile *nfp=sam_open(P.normal_bam,"r");

    if(!tfp||!nfp){
        fprintf(stderr,"Error: cannot open BAM\n");
        return 1;
    }

    bam_hdr_t *hdr=sam_hdr_read(tfp);
    faidx_t *fai=fai_load(P.ref);

    if(!fai){
        fprintf(stderr,"Error: cannot load reference\n");
        return 1;
    }

    void *data[2]={tfp,nfp};
    bam_mplp_t mplp=bam_mplp_init(2,(bam_plp_auto_f)bam_read1,data);

    int tid,pos,n_plp;
    const bam_pileup1_t *plp[2];

    printf("##fileformat=VCFv4.2\n");
    printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

    while(bam_mplp_auto(mplp,&tid,&pos,&n_plp,plp)>0){

        stats_t t={0},n={0};

        int l;
        char *seq=faidx_fetch_seq(fai,
            hdr->target_name[tid],
            pos-1,pos+2,&l);

        if(!seq||l<4) continue;

        char refb=seq[1];

        collect(plp[0],n_plp,refb,&t);
        collect(plp[1],n_plp,refb,&n);

        if(t.depth<P.min_depth) continue;

        double vaf=t.alt_w/t.depth_w;
        if(vaf<P.min_vaf) continue;

        double adj_vaf=vaf-P.contamination*0.5;
        if(adj_vaf<0) adj_vaf=0;

        int ctx=get_ctx(seq);

        double p_err=ctx_err[ctx];

        double ll_strand=
            ll_binomial(t.fwd_alt,
                t.fwd_alt+t.rev_alt,
                P.strand_art);

        double ll_orient=
            ll_binomial(t.f1r2,
                t.f1r2+t.f2r1,
                P.orient_art);

        double mean_pos=t.alt_count>0?
            t.pos_sum/t.alt_count:0.5;

        double ll_pos=-5*fabs(mean_pos-0.5);

        double ll_var=
            ll_beta_binom(t.alt,t.depth,
                adj_vaf*P.dispersion,
                (1-adj_vaf)*P.dispersion);

        double ll_art=
            ll_beta_binom(t.alt,t.depth,
                p_err*P.dispersion,
                (1-p_err)*P.dispersion);

        double llr=
            ll_var-ll_art
            +ll_pos
            -ll_strand
            -ll_orient;

        double thres=P.llr_thres;

        if(ctx==CPG) thres+=1;
        if(ctx==CHH) thres-=1;

        if(llr<thres) continue;

        double post=calibrate(llr,P.cal_a,P.cal_b);

        printf("%s\t%d\t.\t%c\tT\t.\tPASS\tVAF=%.4f;LLR=%.2f;POST=%.3f;CTX=%d\n",
            hdr->target_name[tid],
            pos+1,
            refb,
            vaf,
            llr,
            post,
            ctx);

        free(seq);
    }

    bam_mplp_destroy(mplp);
    sam_close(tfp);
    sam_close(nfp);

    return 0;
}
