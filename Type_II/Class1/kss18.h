#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#define X_length 85

/*============================================================================*/
/* Field                                                                      */
/*============================================================================*/
typedef struct{
	mpz_t x0;
}Fp;
typedef struct{
	Fp x0;
	Fp x1;
	Fp x2;
}Fp3;
typedef struct{
	Fp3 x0;
	Fp3 x1;
	Fp3 x2;
}Fp9;
typedef struct{
	Fp9 x0;
	Fp9 x1;
}Fp18;
/*----------------------------------------------------------------------------*/
//Fp
void Fp_init(Fp *A);
void Fp_clear(Fp *A);
void Fp_printf(Fp *A,char *str);
void Fp_set(Fp *ANS,Fp *A);
void Fp_set_ui(Fp *ANS,unsigned long int UI);
void Fp_set_mpz(Fp *ANS,mpz_t A);
void Fp_set_neg(Fp *ANS,Fp *A);
void Fp_set_random(Fp *ANS,gmp_randstate_t state);
void Fp_mul(Fp *ANS,Fp *A,Fp *B);
void Fp_mul_ui(Fp *ANS,Fp *A,unsigned long int UI);
void Fp_mul_mpz(Fp *ANS,Fp *A,mpz_t B);
void Fp_add(Fp *ANS,Fp *A,Fp *B);
void Fp_add_ui(Fp *ANS,Fp *A,unsigned long int UI);
void Fp_add_mpz(Fp *ANS,Fp *A,mpz_t B);
void Fp_sub(Fp *ANS,Fp *A,Fp *B);
void Fp_sub_ui(Fp *ANS,Fp *A,unsigned long int UI);
void Fp_sub_mpz(Fp *ANS,Fp *A,mpz_t B);
void Fp_inv(Fp *ANS,Fp *A);
int  Fp_legendre(Fp *A);
int  Fp_isCNR(Fp *A);
void Fp_sqrt(Fp *ANS,Fp *A);
void Fp_pow(Fp *ANS,Fp *A,mpz_t scalar);
int  Fp_cmp(Fp *A,Fp *B);
int  Fp_cmp_ui(Fp *A,unsigned long int UI);
int  Fp_cmp_mpz(Fp *A,mpz_t B);
int  Fp_cmp_zero(Fp *A);
int  Fp_cmp_one(Fp *A);
/*----------------------------------------------------------------------------*/
//Fp3
void Fp3_init(Fp3 *A);
void Fp3_clear(Fp3 *A);
void Fp3_printf(Fp3 *A,char *str);
void Fp3_set(Fp3 *ANS,Fp3 *A);
void Fp3_set_ui(Fp3 *ANS,unsigned long int UI);
void Fp3_set_mpz(Fp3 *ANS,mpz_t A);
void Fp3_set_neg(Fp3 *ANS,Fp3 *A);
void Fp3_set_random(Fp3 *ANS,gmp_randstate_t state);
void Fp3_mul(Fp3 *ANS,Fp3 *A,Fp3 *B);
void Fp3_mul_ui(Fp3 *ANS,Fp3 *A,unsigned long int UI);
void Fp3_mul_mpz(Fp3 *ANS,Fp3 *A,mpz_t B);
void Fp3_mul_basis(Fp3 *ANS,Fp3 *A);
void Fp3_inv_basis(Fp3 *ANS,Fp3 *A);//TODO
void Fp3_sqr(Fp3 *ANS,Fp3 *A);
void Fp3_add(Fp3 *ANS,Fp3 *A,Fp3 *B);
void Fp3_add_ui(Fp3 *ANS,Fp3 *A,unsigned long int UI);
void Fp3_add_mpz(Fp3 *ANS,Fp3 *A,mpz_t B);
void Fp3_sub(Fp3 *ANS,Fp3 *A,Fp3 *B);
void Fp3_sub_ui(Fp3 *ANS,Fp3 *A,unsigned long int UI);
void Fp3_sub_mpz(Fp3 *ANS,Fp3 *A,mpz_t B);
void Fp3_inv(Fp3 *ANS,Fp3 *A);
int  Fp3_legendre(Fp3 *A);
int  Fp3_isCNR(Fp3 *A);
void Fp3_sqrt(Fp3 *ANS,Fp3 *A);
void Fp3_pow(Fp3 *ANS,Fp3 *A,mpz_t scalar);
int  Fp3_cmp(Fp3 *A,Fp3 *B);
int  Fp3_cmp_ui(Fp3 *A,unsigned long int UI);
int  Fp3_cmp_mpz(Fp3 *A,mpz_t B);
int  Fp3_cmp_zero(Fp3 *A);
int  Fp3_cmp_one(Fp3 *A);
void Fp3_frobenius_map_p1(Fp3 *ANS,Fp3 *A);
void Fp3_frobenius_map_p2(Fp3 *ANS,Fp3 *A);
void Fp3_frobenius_map_p3(Fp3 *ANS,Fp3 *A);
/*----------------------------------------------------------------------------*/
//Fp9
void Fp9_init(Fp9 *A);
void Fp9_clear(Fp9 *A);
void Fp9_printf(Fp9 *A,char *str);
void Fp9_set(Fp9 *ANS,Fp9 *A);
void Fp9_set_ui(Fp9 *ANS,unsigned long int UI);
void Fp9_set_mpz(Fp9 *ANS,mpz_t A);
void Fp9_set_neg(Fp9 *ANS,Fp9 *A);
void Fp9_set_random(Fp9 *ANS,gmp_randstate_t state);
void Fp9_mul(Fp9 *ANS,Fp9 *A,Fp9 *B);
void Fp9_mul_ui(Fp9 *ANS,Fp9 *A,unsigned long int UI);
void Fp9_mul_mpz(Fp9 *ANS,Fp9 *A,mpz_t B);
void Fp9_mul_basis(Fp9 *ANS,Fp9 *A);
void Fp9_sqr(Fp9 *ANS,Fp9 *A);
void Fp9_add(Fp9 *ANS,Fp9 *A,Fp9 *B);
void Fp9_add_ui(Fp9 *ANS,Fp9 *A,unsigned long int UI);
void Fp9_add_mpz(Fp9 *ANS,Fp9 *A,mpz_t B);
void Fp9_sub(Fp9 *ANS,Fp9 *A,Fp9 *B);
void Fp9_sub_ui(Fp9 *ANS,Fp9 *A,unsigned long int UI);
void Fp9_sub_mpz(Fp9 *ANS,Fp9 *A,mpz_t B);
void Fp9_inv(Fp9 *ANS,Fp9 *A);
int  Fp9_legendre(Fp9 *A);
int  Fp9_isCNR(Fp9 *A);
void Fp9_sqrt(Fp9 *ANS,Fp9 *A);
void Fp9_pow(Fp9 *ANS,Fp9 *A,mpz_t scalar);
int  Fp9_cmp(Fp9 *A,Fp9 *B);
int  Fp9_cmp_ui(Fp9 *A,unsigned long int UI);
int  Fp9_cmp_mpz(Fp9 *A,mpz_t B);
int  Fp9_cmp_zero(Fp9 *A);
int  Fp9_cmp_one(Fp9 *A);
/*----------------------------------------------------------------------------*/
//Fp18
void Fp18_init(Fp18 *A);
void Fp18_clear(Fp18 *A);
void Fp18_printf(Fp18 *A,char *str);
void Fp18_set(Fp18 *ANS,Fp18 *A);
void Fp18_set_ui(Fp18 *ANS,unsigned long int UI);
void Fp18_set_mpz(Fp18 *ANS,mpz_t A);
void Fp18_set_neg(Fp18 *ANS,Fp18 *A);
void Fp18_set_random(Fp18 *ANS,gmp_randstate_t state);
void Fp18_mul(Fp18 *ANS,Fp18 *A,Fp18 *B);
void Fp18_mul_ui(Fp18 *ANS,Fp18 *A,unsigned long int UI);
void Fp18_mul_mpz(Fp18 *ANS,Fp18 *A,mpz_t B);
void Fp18_sqr(Fp18 *ANS,Fp18 *A);
void Fp18_sqr_cyclotomic(Fp18 *ANS,Fp18 *A);
void Fp18_add(Fp18 *ANS,Fp18 *A,Fp18 *B);
void Fp18_add_ui(Fp18 *ANS,Fp18 *A,unsigned long int UI);
void Fp18_add_mpz(Fp18 *ANS,Fp18 *A,mpz_t B);
void Fp18_sub(Fp18 *ANS,Fp18 *A,Fp18 *B);
void Fp18_sub_ui(Fp18 *ANS,Fp18 *A,unsigned long int UI);
void Fp18_sub_mpz(Fp18 *ANS,Fp18 *A,mpz_t B);
void Fp18_inv(Fp18 *ANS,Fp18 *A);
int  Fp18_legendre(Fp18 *A);
int  Fp18_isCNR(Fp18 *A);
void Fp18_sqrt(Fp18 *ANS,Fp18 *A);
void Fp18_pow(Fp18 *ANS,Fp18 *A,mpz_t scalar);
int  Fp18_cmp(Fp18 *A,Fp18 *B);
int  Fp18_cmp_ui(Fp18 *A,unsigned long int UI);
int  Fp18_cmp_mpz(Fp18 *A,mpz_t B);
int  Fp18_cmp_zero(Fp18 *A);
int  Fp18_cmp_one(Fp18 *A);
void Fp18_frobenius_map_p1(Fp18 *ANS,Fp18 *A);
void Fp18_frobenius_map_p2(Fp18 *ANS,Fp18 *A);
void Fp18_frobenius_map_p3(Fp18 *ANS,Fp18 *A);
void Fp18_frobenius_map_p4(Fp18 *ANS,Fp18 *A);
void Fp18_frobenius_map_p5(Fp18 *ANS,Fp18 *A);
void Fp18_frobenius_map_p6(Fp18 *ANS,Fp18 *A);
void Fp18_frobenius_map_p9(Fp18 *ANS,Fp18 *A);

/*============================================================================*/
/* Elliptic Curve                                                             */
/*============================================================================*/
typedef struct{
	Fp x,y;
	int infinity;
}EFp;

typedef struct{
    Fp3 x,y;
	int infinity;
}EFp3;

typedef struct{
    Fp9 x,y;
	int infinity;
}EFp9;

typedef struct{
    Fp18 x,y;
	int infinity;
}EFp18;
/*----------------------------------------------------------------------------*/
//EFp
void EFp_init(EFp *P);
void EFp_set(EFp *P,EFp *A);
void EFp_set_ui(EFp *ANS,unsigned long int UI);
void EFp_set_mpz(EFp *ANS,mpz_t A);
void EFp_set_neg(EFp *ANS,EFp *A);
void EFp_clear(EFp *P);
void EFp_printf(EFp *P,char *str);
void EFp_rational_point(EFp *P);
void EFp_ECD(EFp *ANS,EFp *P);
void EFp_ECA(EFp *ANS,EFp *P1,EFp *P2);
void EFp_SCM(EFp *ANS,EFp *P,mpz_t scalar);
//skew_frobenius_map
void EFp_skew_frobenius_map_p3(EFp *ANS,EFp *A);
/*----------------------------------------------------------------------------*/
//EFp3
void EFp3_init(EFp3 *P);
void EFp3_set(EFp3 *ANS,EFp3 *A);
void EFp3_set_ui(EFp3 *ANS,unsigned long int UI);
void EFp3_set_mpz(EFp3 *ANS,mpz_t A);
void EFp3_set_neg(EFp3 *ANS,EFp3 *A);
void EFp3_clear(EFp3 *P);
void EFp3_printf(EFp3 *P,char *str);
void EFp3_rational_point(EFp3 *P);
void EFp3_ECD(EFp3 *ANS,EFp3 *P);
void EFp3_ECA(EFp3 *ANS,EFp3 *P1,EFp3 *P2);
void EFp3_SCM(EFp3 *ANS,EFp3 *P,mpz_t scalar);
//skew_frobenius_map
void EFp3_skew_frobenius_map_p1(EFp3 *ANS,EFp3 *A);
void EFp3_skew_frobenius_map_p2(EFp3 *ANS,EFp3 *A);
void EFp3_skew_frobenius_map_p3(EFp3 *ANS,EFp3 *A);
/*----------------------------------------------------------------------------*/
//EFp9
void EFp9_init(EFp9 *P);
void EFp9_set(EFp9 *ANS,EFp9 *A);
void EFp9_set_ui(EFp9 *ANS,unsigned long int UI);
void EFp9_set_mpz(EFp9 *ANS,mpz_t A);
void EFp9_set_neg(EFp9 *ANS,EFp9 *A);
void EFp9_clear(EFp9 *P);
void EFp9_printf(EFp9 *P,char *str);
void EFp9_rational_point(EFp9 *P);
void EFp9_ECD(EFp9 *ANS,EFp9 *P);
void EFp9_ECA(EFp9 *ANS,EFp9 *P1,EFp9 *P2);
void EFp9_SCM(EFp9 *ANS,EFp9 *P,mpz_t scalar);
/*----------------------------------------------------------------------------*/
//EFp18
void EFp18_init(EFp18 *P);
void EFp18_set(EFp18 *ANS,EFp18 *A);
void EFp18_set_ui(EFp18 *ANS,unsigned long int UI);
void EFp18_set_mpz(EFp18 *ANS,mpz_t A);
void EFp18_set_neg(EFp18 *ANS,EFp18 *A);
void EFp18_clear(EFp18 *P);
void EFp18_printf(EFp18 *P,char *str);
void EFp18_rational_point(EFp18 *P);
void EFp18_generate_G1(EFp18 *P);
void EFp18_generate_G2(EFp18 *Q);
void EFp18_ECD(EFp18 *ANS,EFp18 *P);
void EFp18_ECA(EFp18 *ANS,EFp18 *P1,EFp18 *P2);
void EFp18_SCM(EFp18 *ANS,EFp18 *P,mpz_t scalar);
/*============================================================================*/
/* Pairing functions                                                          */
/*============================================================================*/
enum f_state{
	f_p1,f_p2,f_p3,f_p4,f_p5,f_p6,f_p7,f_p8,f_p9,f_p10,f_p11,f_p12,f_13,f_14,f_15,f_16,f_17,f_18
};
enum x_state{
	x_1,x_2,x_3,x_4,x_5,x_6,x_7
};
int X_binary[X_length+1];
mpz_t X,prime,order,trace;
mpz_t EFp_total,EFp3_total,EFp9_total,EFp18_total;
Fp curve_b;
Fp3 frobenius_constant[18][6];
Fp frobenius_constant_p3[6],frobenius_constant_p6[6];
Fp3 Tau1_2,Tau1_2_inv;
Fp3 ONE;
Fp seven_inv;
mpz_t epsilon1,epsilon2;

Fp TMP1_FP,TMP2_FP,TMP3_FP,TMP4_FP,TMP5_FP,TMP6_FP,TMP7_FP,TMP8_FP,TMP9_FP,TMP10_FP,TMP11_FP,TMP12_FP,TMP13_FP,TMP14_FP;
Fp3 TMP1_FP3,TMP2_FP3,TMP3_FP3,TMP4_FP3,TMP5_FP3,TMP6_FP3,TMP7_FP3;
Fp9 TMP1_FP9,TMP2_FP9,TMP3_FP9,TMP4_FP9;
Fp18 TMP1_FP18,TMP2_FP18,TMP3_FP18;

EFp TMP1_EFP,TMP2_EFP;
EFp3 TMP1_EFP3,TMP2_EFP3;
EFp9 TMP1_EFP9,TMP2_EFP9;
EFp18 TMP1_EFP18,TMP2_EFP18;
/*----------------------------------------------------------------------------*/
//twist
void EFp18_to_EFp3(EFp3 *ANS,EFp18 *A);
void EFp3_to_EFp18(EFp18 *ANS,EFp3 *A);
void EFp18_to_EFp(EFp *ANS,EFp18 *A);
void EFp_to_EFp18(EFp18 *ANS,EFp *A);
/*----------------------------------------------------------------------------*/
//Pseudo 12-sparse
void Pseudo_12_sparse_mapping(EFp *P,EFp3 *Q,Fp *L);
void Pseudo_12_sparse_mul(Fp18 *ANS,Fp18 *A,Fp18 *B);
void ff_ltt(Fp18 *f,EFp3 *T,EFp *P,Fp *L);
void f_ltq(Fp18 *f,EFp3 *T,EFp3 *Q,EFp *P,Fp *L);
/*----------------------------------------------------------------------------*/
//miller
void Miller_algo_for_plain_ate(Fp18 *ANS,EFp18 *Q,EFp18 *P);
void Miller_algo_for_opt_ate(Fp18 *ANS,EFp18 *Q,EFp18 *P);
/*----------------------------------------------------------------------------*/
//final exp
void Fp18_pow_X(Fp18 *ANS,Fp18 *A);
void Final_exp_easy(Fp18 *ANS,Fp18 *f);
void Final_exp_hard_plain(Fp18 *ANS,Fp18 *f);
void Final_exp_hard_optimal(Fp18 *ANS,Fp18 *f);
/*----------------------------------------------------------------------------*/
//pairing
void Plain_ate_pairing(Fp18 *ANS,EFp18 *P,EFp18 *Q);
void Opt_ate_pairing(Fp18 *ANS,EFp18 *P,EFp18 *Q);
/*----------------------------------------------------------------------------*/
//JSF
void Joint_sparse_form(int **binary,mpz_t S[2],int *loop_length);
//G1 SCM
void EFp18_G1_SCM_plain(EFp18 *ANS,EFp18 *P,mpz_t scalar);
void EFp18_G1_SCM_2split(EFp18 *ANS,EFp18 *P,mpz_t scalar);
void EFp18_G1_SCM_2split_JSF(EFp18 *ANS,EFp18 *P,mpz_t scalar);
/*----------------------------------------------------------------------------*/
//G2 SCM
void EFp18_G2_SCM_plain(EFp18 *ANS,EFp18 *Q,mpz_t scalar);
void EFp18_G2_SCM_2split(EFp18 *ANS,EFp18 *Q,mpz_t scalar);
void EFp18_G2_SCM_2split_JSF(EFp18 *ANS,EFp18 *Q,mpz_t scalar);
void EFp18_G2_SCM_6split(EFp18 *ANS,EFp18 *Q,mpz_t scalar);
/*----------------------------------------------------------------------------*/
//G3 EXP
void Fp18_G3_EXP_plain(Fp18 *ANS,Fp18 *A,mpz_t scalar);
void Fp18_G3_EXP_2split(Fp18 *ANS,Fp18 *A,mpz_t scalar);
void Fp18_G3_EXP_2split_JSF(Fp18 *ANS,Fp18 *A,mpz_t scalar);
void Fp18_G3_EXP_6split(Fp18 *ANS,Fp18 *A,mpz_t scalar);
/*----------------------------------------------------------------------------*/
//init/set/clear
void KSS18_init();
void KSS18_clear();
void KSS18_print_parameters();
void init_parameters();
void generate_X();
int  generate_prime();
int  generate_order();
void generate_trace();
void set_basis();
void weil();
void get_epsilon();
void set_frobenius_constant();
void set_curve_coefficient();
/*----------------------------------------------------------------------------*/
struct timeval tv_start,tv_end;
float MILLER_TATE,MILLER_PLAINATE,MILLER_OPTATE,MILLER_XATE;
float FINALEXP_PLAIN_EASY,FINALEXP_PLAIN_HARD,FINALEXP_OPT_EASY,FINALEXP_OPT_HARD;
float G1SCM_PLAIN,G1SCM_2SPLIT,G1SCM_2SPLIT_JSF;
float G2SCM_PLAIN,G2SCM_2SPLIT,G2SCM_2SPLIT_JSF,G2SCM_6SPLIT;
float G3EXP_PLAIN,G3EXP_2SPLIT,G3EXP_2SPLIT_JSF,G3EXP_6SPLIT;
struct mpz_Cost{
    unsigned long int mpz_mul;
    unsigned long int mpz_mul_ui;
    unsigned long int mpz_sqr;
    unsigned long int mpz_add;
    unsigned long int mpz_add_ui;
    unsigned long int mpz_invert;
};
struct Fp_Cost{
    unsigned long int Fp_mul;
    unsigned long int Fp_mul_mpz;
    unsigned long int Fp_mul_ui;
    unsigned long int Fp_sqr;
    unsigned long int Fp_basis;
    unsigned long int Fp_add;
    unsigned long int Fp_add_mpz;
    unsigned long int Fp_add_ui;
    unsigned long int Fp_inv;
    unsigned long int Fp_neg;
};
struct mpz_Cost mpz_cost;
struct Fp_Cost Fp_cost,Fp_cost_ave;

/*----------------------------------------------------------------------------*/
//time
float timedifference_msec(struct timeval t0, struct timeval t1);
float timedifference_usec(struct timeval t0, struct timeval t1);
/*----------------------------------------------------------------------------*/
//cost
void Init_mpz_Cost(struct mpz_Cost *cost);
void Print_mpz_Cost(struct mpz_Cost *cost,char *str);
void Init_Fp_Cost(struct Fp_Cost *cost);
void Print_Fp_Cost(struct Fp_Cost *cost,char *str);
/*----------------------------------------------------------------------------*/
//test
void test_Field();
void test_Frobenius_map();
void test_skew_frobenius_map();
void test_rational_point();
void test_twist();
void test_plain_ate_pairing();
void test_opt_ate_pairing();
void test_G1_SCM();
void test_G2_SCM();
void test_G3_EXP();
void computation_time();
void computation_cost();
void finite_field_cost();
void ecc_cost();
