#include "kss18.h"

/*============================================================================*/
/* main                                                                       */
/*============================================================================*/
int main(void){
    KSS18_init();
    KSS18_print_parameters();
    
    //test_Field();
    //test_Frobenius_map();
    //test_skew_frobenius_map();
    //test_rational_point();
    //test_twist();
    //test_plain_ate_pairing();
    test_opt_ate_pairing();
    //test_G1_SCM();
    //test_G2_SCM();
    //test_G3_EXP();
    //computation_time();
    //computation_cost();
    //finite_field_cost();
    //ecc_cost();
      
    KSS18_clear();
    return 0;
}

/*============================================================================*/
/* Field                                                                      */
/*============================================================================*/
/*----------------------------------------------------------------------------*/
//Fp
void Fp_init(Fp *A){
    mpz_init(A->x0);
}

void Fp_clear(Fp *A){
    mpz_clear(A->x0);
}

void Fp_printf(Fp *A,char *str){
    gmp_printf("%s%Zd",str,A->x0);
}

void Fp_set(Fp *ANS,Fp *A){
    mpz_set(ANS->x0,A->x0);
}

void Fp_set_ui(Fp *ANS,unsigned long int UI){
    mpz_set_ui(ANS->x0,UI);
}

void Fp_set_mpz(Fp *ANS,mpz_t A){
    mpz_set(ANS->x0,A);
}

//------------------------------------------------------//
//time

void Fp_set_neg(Fp *ANS,Fp *A){
    mpz_sub(ANS->x0,prime,A->x0);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_set_random(Fp *ANS,gmp_randstate_t state){
    mpz_urandomm(ANS->x0,state,prime);
}

void Fp_mul(Fp *ANS,Fp *A,Fp *B){
    mpz_mul(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_mul_ui(Fp *ANS,Fp *A,unsigned long int UI){
    mpz_mul_ui(ANS->x0,A->x0,UI);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_mul_mpz(Fp *ANS,Fp *A,mpz_t B){
    mpz_mul(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_add(Fp *ANS,Fp *A,Fp *B){
    mpz_add(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_add_ui(Fp *ANS,Fp *A,unsigned long int UI){
    mpz_add_ui(ANS->x0,A->x0,UI);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_add_mpz(Fp *ANS,Fp *A,mpz_t B){
    mpz_add(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_sub(Fp *ANS,Fp *A,Fp *B){
    mpz_sub(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_sub_ui(Fp *ANS,Fp *A,unsigned long int UI){
    mpz_sub_ui(ANS->x0,A->x0,UI);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_sub_mpz(Fp *ANS,Fp *A,mpz_t B){
    mpz_sub(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_inv(Fp *ANS,Fp *A){
    mpz_invert(ANS->x0,A->x0,prime);
}
//------------------------------------------------------//
//cost
/*
void Fp_set_neg(Fp *ANS,Fp *A){
    mpz_cost.mpz_add++;
    Fp_cost.Fp_neg++;
    mpz_sub(ANS->x0,prime,A->x0);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_set_random(Fp *ANS,gmp_randstate_t state){
    mpz_urandomm(ANS->x0,state,prime);
}

void Fp_mul(Fp *ANS,Fp *A,Fp *B){
    if(mpz_cmp(A->x0,B->x0)==0){
        mpz_cost.mpz_sqr++;
        Fp_cost.Fp_sqr++;
    }else{
        mpz_cost.mpz_mul++;
        Fp_cost.Fp_mul++;
    }
    mpz_mul(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_mul_ui(Fp *ANS,Fp *A,unsigned long int UI){
    Fp_cost.Fp_mul_ui++;
    mpz_cost.mpz_mul_ui++;
    mpz_mul_ui(ANS->x0,A->x0,UI);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_mul_mpz(Fp *ANS,Fp *A,mpz_t B){
    if(mpz_cmp(A->x0,B)==0){
        mpz_cost.mpz_sqr++;
        Fp_cost.Fp_sqr++;
    }else{
        mpz_cost.mpz_mul++;
        Fp_cost.Fp_mul_mpz++;
    }
    mpz_mul(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_add(Fp *ANS,Fp *A,Fp *B){
    Fp_cost.Fp_add++;
    mpz_cost.mpz_add++;
    mpz_add(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_add_ui(Fp *ANS,Fp *A,unsigned long int UI){
    Fp_cost.Fp_add_ui++;
    mpz_cost.mpz_add_ui++;
    mpz_add_ui(ANS->x0,A->x0,UI);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_add_mpz(Fp *ANS,Fp *A,mpz_t B){
    Fp_cost.Fp_add_mpz++;
    mpz_cost.mpz_add++;
    mpz_add(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_sub(Fp *ANS,Fp *A,Fp *B){
    Fp_cost.Fp_add++;
    mpz_cost.mpz_add++;
    mpz_sub(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_sub_ui(Fp *ANS,Fp *A,unsigned long int UI){
    Fp_cost.Fp_add_ui++;
    mpz_cost.mpz_add_ui++;
    mpz_sub_ui(ANS->x0,A->x0,UI);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_sub_mpz(Fp *ANS,Fp *A,mpz_t B){
    Fp_cost.Fp_add_mpz++;
    mpz_cost.mpz_add++;
    mpz_sub(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,prime);
}

void Fp_inv(Fp *ANS,Fp *A){
    Fp_cost.Fp_inv++;
    mpz_cost.mpz_invert++;
    mpz_invert(ANS->x0,A->x0,prime);
}
*/
//------------------------------------------------------//

int  Fp_legendre(Fp *A){
    return mpz_legendre(A->x0,prime);
}

int  Fp_isCNR(Fp *A){
    Fp tmp;
    Fp_init(&tmp);
    mpz_t exp;
    mpz_init(exp);
    
    mpz_sub_ui(exp,prime,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp_pow(&tmp,A,exp);
    
    if(Fp_cmp_one(&tmp)==0){
        mpz_clear(exp);
        Fp_clear(&tmp);
        return 1;
    }else{
        mpz_clear(exp);
        Fp_clear(&tmp);
        return -1;
    }
}

void Fp_sqrt(Fp *ANS,Fp *A){
    Fp x,y,t,k,n,tmp;
    Fp_init(&x);
    Fp_init(&y);
    Fp_init(&t);
    Fp_init(&k);
    Fp_init(&n);
    Fp_init(&tmp);
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    Fp_set_random(&n,state);
    
    while(Fp_legendre(&n)!=-1){
        Fp_set_random(&n,state);
    }
    mpz_sub_ui(q,prime,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    Fp_pow(&y,&n,q);
    mpz_set_ui(z,e);
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp_pow(&x,A,exp);
    Fp_mul(&tmp,&x,&x);
    Fp_mul(&k,&tmp,A);
    Fp_mul(&x,&x,A);
    while(mpz_cmp_ui(k.x0,1)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        Fp_pow(&tmp,&k,exp);
        while(mpz_cmp_ui(tmp.x0,1)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            Fp_pow(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        Fp_pow(&t,&y,result);
        Fp_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        Fp_mul(&x,&x,&t);
        Fp_mul(&k,&k,&y);
    }
    Fp_set(ANS,&x);
    
    mpz_clear(exp);
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(result);
    Fp_clear(&x);
    Fp_clear(&y);
    Fp_clear(&t);
    Fp_clear(&k);
    Fp_clear(&n);
    Fp_clear(&tmp);
}

void Fp_pow(Fp *ANS,Fp *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    Fp tmp;
    Fp_init(&tmp);
    
    Fp_set(&tmp,A);
    
    for(i=1; binary[i]!='\0'; i++){
        Fp_mul(&tmp,&tmp,&tmp);
        if(binary[i]=='1'){
            Fp_mul(&tmp,A,&tmp);
        }
    }
    Fp_set(ANS,&tmp);
    
    Fp_clear(&tmp);
}


int  Fp_cmp(Fp *A,Fp *B){
    if(mpz_cmp(A->x0,B->x0)==0){
        return 0;   
    }
    return 1;
}

int  Fp_cmp_ui(Fp *A,unsigned long int UI){
    if(mpz_cmp_ui(A->x0,UI)==0){
        return 0;
    }
    return 1;
}

int  Fp_cmp_mpz(Fp *A,mpz_t B){
    if(mpz_cmp(A->x0,B)==0){
        return 0;
    }
    return 1;
}

int  Fp_cmp_zero(Fp *A){
    if(mpz_cmp_ui(A->x0,0)==0){
        return 0;
    }
    return 1;
}

int  Fp_cmp_one(Fp *A){
    if(mpz_cmp_ui(A->x0,1)==0){
        return 0;
    }
    return 1;
}
/*----------------------------------------------------------------------------*/
//Fp3
void Fp3_init(Fp3 *A){
    Fp_init(&A->x0);
    Fp_init(&A->x1);
    Fp_init(&A->x2);
}

void Fp3_clear(Fp3 *A){
    Fp_clear(&A->x0);
    Fp_clear(&A->x1);
    Fp_clear(&A->x2);
}

void Fp3_printf(Fp3 *A,char *str){
    gmp_printf("%s(",str);
    Fp_printf(&A->x0,"");
    gmp_printf(",");
    Fp_printf(&A->x1,"");
    gmp_printf(",");
    Fp_printf(&A->x2,"");
    gmp_printf(")");
}

void Fp3_set(Fp3 *ANS,Fp3 *A){
    Fp_set(&ANS->x0,&A->x0);
    Fp_set(&ANS->x1,&A->x1);
    Fp_set(&ANS->x2,&A->x2);
}

void Fp3_set_ui(Fp3 *ANS,unsigned long int UI){
    Fp_set_ui(&ANS->x0,UI);
    Fp_set_ui(&ANS->x1,UI);
    Fp_set_ui(&ANS->x2,UI);
}

void Fp3_set_mpz(Fp3 *ANS,mpz_t A){
    Fp_set_mpz(&ANS->x0,A);
    Fp_set_mpz(&ANS->x1,A);
    Fp_set_mpz(&ANS->x2,A);
}

void Fp3_set_neg(Fp3 *ANS,Fp3 *A){
    Fp_set_neg(&ANS->x0,&A->x0);
    Fp_set_neg(&ANS->x1,&A->x1);
    Fp_set_neg(&ANS->x2,&A->x2);
}

void Fp3_set_random(Fp3 *ANS,gmp_randstate_t state){
    Fp_set_random(&ANS->x0,state);
    Fp_set_random(&ANS->x1,state);
    Fp_set_random(&ANS->x2,state);
}

void Fp3_mul(Fp3 *ANS,Fp3 *A,Fp3 *B){
    //AB = (c0, c1, c2)
    //   = (a0 + a1α + a2α^2)(b0 + b1α + b2α^2)
    //   = a0b0 + 2(a1b2 + a2b1) + (a0a1 + a1b0 + 2a2b2)α + (a0b2 + a2b0 + a1b1)α^2
    Fp_mul(&TMP2_FP,&A->x0,&B->x0);
    Fp_mul(&TMP3_FP,&A->x1,&B->x1);
    Fp_mul(&TMP4_FP,&A->x2,&B->x2);
    
    Fp_add(&TMP5_FP,&A->x1,&A->x2);
    Fp_add(&TMP1_FP,&B->x1,&B->x2);
    Fp_mul(&TMP5_FP,&TMP5_FP,&TMP1_FP);
    
    Fp_add(&TMP6_FP,&A->x0,&A->x1);
    Fp_add(&TMP1_FP,&B->x0,&B->x1);
    Fp_mul(&TMP6_FP,&TMP6_FP,&TMP1_FP);
    
    Fp_add(&TMP7_FP,&A->x0,&A->x2);
    Fp_add(&TMP1_FP,&B->x0,&B->x2);
    Fp_mul(&TMP7_FP,&TMP7_FP,&TMP1_FP);
    
    Fp_sub(&ANS->x0,&TMP5_FP,&TMP3_FP);
    Fp_sub(&ANS->x0,&ANS->x0,&TMP4_FP);
    Fp_add(&ANS->x0,&ANS->x0,&ANS->x0);
    Fp_add(&ANS->x0,&ANS->x0,&TMP2_FP);
    
    Fp_add(&ANS->x1,&TMP4_FP,&TMP4_FP);
    Fp_add(&ANS->x1,&ANS->x1,&TMP6_FP);
    Fp_sub(&ANS->x1,&ANS->x1,&TMP2_FP);
    Fp_sub(&ANS->x1,&ANS->x1,&TMP3_FP);
    
    Fp_add(&ANS->x2,&TMP3_FP,&TMP7_FP);
    Fp_sub(&ANS->x2,&ANS->x2,&TMP2_FP);
    Fp_sub(&ANS->x2,&ANS->x2,&TMP4_FP);
}

void Fp3_mul_ui(Fp3 *ANS,Fp3 *A,unsigned long int UI){
    Fp_mul_ui(&ANS->x0,&A->x0,UI);
    Fp_mul_ui(&ANS->x1,&A->x1,UI);
    Fp_mul_ui(&ANS->x2,&A->x2,UI);
}

void Fp3_mul_mpz(Fp3 *ANS,Fp3 *A,mpz_t B){
    Fp_mul_mpz(&ANS->x0,&A->x0,B);
    Fp_mul_mpz(&ANS->x1,&A->x1,B);
    Fp_mul_mpz(&ANS->x2,&A->x2,B);
}

void Fp3_mul_basis(Fp3 *ANS,Fp3 *A){
    Fp_set(&TMP1_FP,&A->x0);
    Fp_set(&TMP2_FP,&A->x1);
    Fp_set(&TMP3_FP,&A->x2);
    Fp_add(&ANS->x0,&TMP3_FP,&TMP3_FP);
    Fp_set(&ANS->x1,&TMP1_FP);
    Fp_set(&ANS->x2,&TMP2_FP);
}

void Fp3_inv_basis(Fp3 *ANS,Fp3 *A){
    Fp_mul_mpz(&TMP1_FP,&A->x0,Alpha_inv.x2.x0);
    Fp_mul_mpz(&TMP2_FP,&A->x1,Alpha_inv.x2.x0);
    Fp_mul_mpz(&TMP3_FP,&A->x2,Alpha_inv.x2.x0);
    Fp_add(&ANS->x0,&TMP2_FP,&TMP2_FP);
    Fp_add(&ANS->x1,&TMP3_FP,&TMP3_FP);
    Fp_set(&ANS->x2,&TMP1_FP);
}

void Fp3_sqr(Fp3 *ANS,Fp3 *A){
    //Chung-Hasan
    Fp_mul(&TMP1_FP,&A->x0,&A->x0);
    Fp_mul(&TMP4_FP,&A->x2,&A->x2);
    Fp_add(&TMP5_FP,&A->x1,&A->x1);
    Fp_mul(&TMP2_FP,&TMP5_FP,&A->x2);
    Fp_mul(&TMP3_FP,&A->x0,&TMP5_FP);
    Fp_add(&TMP5_FP,&A->x0,&A->x1);
    Fp_add(&TMP5_FP,&TMP5_FP,&A->x2);
    
    //x0
    Fp_add(&ANS->x0,&TMP2_FP,&TMP2_FP);
    Fp_add(&ANS->x0,&ANS->x0,&TMP1_FP);
    //x1
    Fp_add(&ANS->x1,&TMP4_FP,&TMP4_FP);
    Fp_add(&ANS->x1,&ANS->x1,&TMP3_FP);
    //x2
    Fp_mul(&ANS->x2,&TMP5_FP,&TMP5_FP);
    Fp_add(&TMP5_FP,&TMP1_FP,&TMP4_FP);
    Fp_add(&TMP5_FP,&TMP5_FP,&TMP2_FP);
    Fp_add(&TMP5_FP,&TMP5_FP,&TMP3_FP);
    Fp_sub(&ANS->x2,&ANS->x2,&TMP5_FP);
}

void Fp3_add(Fp3 *ANS,Fp3 *A,Fp3 *B){
    Fp_add(&ANS->x0,&A->x0,&B->x0);
    Fp_add(&ANS->x1,&A->x1,&B->x1);
    Fp_add(&ANS->x2,&A->x2,&B->x2);
}

void Fp3_add_ui(Fp3 *ANS,Fp3 *A,unsigned long int UI){
    Fp_add_ui(&ANS->x0,&A->x0,UI);
    Fp_add_ui(&ANS->x1,&A->x1,UI);
    Fp_add_ui(&ANS->x2,&A->x2,UI);
}

void Fp3_add_mpz(Fp3 *ANS,Fp3 *A,mpz_t B){
    Fp_add_mpz(&ANS->x0,&A->x0,B);
    Fp_add_mpz(&ANS->x1,&A->x1,B);
    Fp_add_mpz(&ANS->x2,&A->x2,B);
}

void Fp3_sub(Fp3 *ANS,Fp3 *A,Fp3 *B){
    Fp_sub(&ANS->x0,&A->x0,&B->x0);
    Fp_sub(&ANS->x1,&A->x1,&B->x1);
    Fp_sub(&ANS->x2,&A->x2,&B->x2);
}

void Fp3_sub_ui(Fp3 *ANS,Fp3 *A,unsigned long int UI){
    Fp_sub_ui(&ANS->x0,&A->x0,UI);
    Fp_sub_ui(&ANS->x1,&A->x1,UI);
    Fp_sub_ui(&ANS->x2,&A->x2,UI);
}

void Fp3_sub_mpz(Fp3 *ANS,Fp3 *A,mpz_t B){
    Fp_sub_mpz(&ANS->x0,&A->x0,B);
    Fp_sub_mpz(&ANS->x1,&A->x1,B);
    Fp_sub_mpz(&ANS->x2,&A->x2,B);
}

void Fp3_inv(Fp3 *ANS,Fp3 *A){
    Fp_mul(&TMP1_FP,&A->x0,&A->x0);
    Fp_mul(&TMP4_FP,&A->x1,&A->x2);
    Fp_add(&TMP4_FP,&TMP4_FP,&TMP4_FP);
    Fp_sub(&TMP1_FP,&TMP1_FP,&TMP4_FP);
    
    Fp_mul(&TMP2_FP,&A->x2,&A->x2);
    Fp_add(&TMP2_FP,&TMP2_FP,&TMP2_FP);
    Fp_mul(&TMP4_FP,&A->x0,&A->x1);
    Fp_sub(&TMP2_FP,&TMP2_FP,&TMP4_FP);
    
    Fp_mul(&TMP3_FP,&A->x1,&A->x1);
    Fp_mul(&TMP4_FP,&A->x0,&A->x2);
    Fp_sub(&TMP3_FP,&TMP3_FP,&TMP4_FP);
    
    Fp_mul(&TMP5_FP,&A->x1,&TMP3_FP);
    Fp_mul(&TMP4_FP,&A->x2,&TMP2_FP);
    Fp_add(&TMP5_FP,&TMP5_FP,&TMP4_FP);
    Fp_add(&TMP5_FP,&TMP5_FP,&TMP5_FP);
    Fp_mul(&TMP4_FP,&A->x0,&TMP1_FP);
    Fp_add(&TMP5_FP,&TMP5_FP,&TMP4_FP);
    
    Fp_inv(&TMP5_FP,&TMP5_FP);
    
    Fp_mul(&ANS->x0,&TMP5_FP,&TMP1_FP);
    Fp_mul(&ANS->x1,&TMP5_FP,&TMP2_FP);
    Fp_mul(&ANS->x2,&TMP5_FP,&TMP3_FP);
}

int  Fp3_legendre(Fp3 *A){
    Fp_mul(&TMP2_FP,&A->x0,&A->x0);
    Fp_mul(&TMP2_FP,&TMP2_FP,&A->x0);
    
    Fp_mul(&TMP1_FP,&A->x1,&A->x1);
    Fp_mul(&TMP1_FP,&TMP1_FP,&A->x1);
    Fp_add(&TMP1_FP,&TMP1_FP,&TMP1_FP);
    Fp_add(&TMP2_FP,&TMP2_FP,&TMP1_FP);
    
    Fp_mul(&TMP1_FP,&A->x2,&A->x2);
    Fp_mul(&TMP1_FP,&TMP1_FP,&A->x2);
    Fp_add(&TMP1_FP,&TMP1_FP,&TMP1_FP);
    Fp_add(&TMP1_FP,&TMP1_FP,&TMP1_FP);
    Fp_add(&TMP2_FP,&TMP2_FP,&TMP1_FP);
    
    Fp_mul(&TMP1_FP,&A->x0,&A->x1);
    Fp_mul(&TMP1_FP,&TMP1_FP,&A->x2);
    Fp_add(&TMP1_FP,&TMP1_FP,&TMP1_FP);
    Fp_sub(&TMP2_FP,&TMP2_FP,&TMP1_FP);
    Fp_sub(&TMP2_FP,&TMP2_FP,&TMP1_FP);
    Fp_sub(&TMP2_FP,&TMP2_FP,&TMP1_FP);
    
    return Fp_legendre(&TMP2_FP);
}

int  Fp3_isCNR(Fp3 *A){
    Fp3 tmp;
    Fp3_init(&tmp);
    mpz_t exp;
    mpz_init(exp);
    
    mpz_pow_ui(exp,prime,3);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp3_pow(&tmp,A,exp);
    
    if(Fp3_cmp_one(&tmp)==0){
        mpz_clear(exp);
        Fp3_clear(&tmp);
        return 1;
    }else{
        mpz_clear(exp);
        Fp3_clear(&tmp);
        return -1;
    }
}

void Fp3_sqrt(Fp3 *ANS,Fp3 *A){//TODO
    Fp3 tmp_A;
    Fp3_init(&tmp_A);
    Fp s,tmp;
    Fp_init(&tmp);
    Fp_init(&s);
    mpz_t p_plus1_over2;
    mpz_init(p_plus1_over2);
    mpz_add_ui(p_plus1_over2,prime,1);
    mpz_tdiv_q_ui(p_plus1_over2,p_plus1_over2,2);
    
    //s=A^{p^2+p+1}
    Fp_mul(&s,&A->x0,&A->x0);
    Fp_mul(&s,&s,&A->x0);
    
    Fp_mul(&tmp,&A->x1,&A->x1);
    Fp_mul(&tmp,&tmp,&A->x1);
    Fp_add(&tmp,&tmp,&tmp);
    Fp_add(&s,&s,&tmp);
    
    Fp_mul(&tmp,&A->x2,&A->x2);
    Fp_mul(&tmp,&tmp,&A->x2);
    Fp_add(&tmp,&tmp,&tmp);
    Fp_add(&tmp,&tmp,&tmp);
    Fp_add(&s,&s,&tmp);
    
    Fp_mul(&tmp,&A->x0,&A->x1);
    Fp_mul(&tmp,&tmp,&A->x2);
    Fp_add(&tmp,&tmp,&tmp);
    Fp_sub(&s,&s,&tmp);
    Fp_sub(&s,&s,&tmp);
    Fp_sub(&s,&s,&tmp);
    
    //sqrt{s}^{-1}
    Fp_sqrt(&s,&s);
    Fp_inv(&s,&s);
    
    //A^{(p^2+p+2)/2}
    Fp_set(&tmp_A.x0,&A->x0);
    Fp_mul(&tmp_A.x1,&A->x1,&frobenius_constant[f_p1][1]);
    Fp_mul(&tmp_A.x2,&A->x2,&frobenius_constant[f_p1][2]);
    Fp3_pow(&tmp_A,&tmp_A,p_plus1_over2);
    Fp3_mul(&tmp_A,&tmp_A,A);
    
    Fp_mul(&ANS->x0,&tmp_A.x0,&s);
    Fp_mul(&ANS->x1,&tmp_A.x1,&s);
    Fp_mul(&ANS->x2,&tmp_A.x2,&s);
    
    Fp_clear(&tmp);
    Fp_clear(&s);
    Fp3_clear(&tmp_A);
    mpz_clear(p_plus1_over2);
}

void Fp3_pow(Fp3 *ANS,Fp3 *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    Fp3 tmp;
    Fp3_init(&tmp);
    Fp3_set(&tmp,A);
    
    for(i=1; binary[i]!='\0'; i++){
        Fp3_sqr(&tmp,&tmp);
        if(binary[i]=='1'){
            Fp3_mul(&tmp,A,&tmp);
        }
    }
    
    Fp3_set(ANS,&tmp);
    Fp3_clear(&tmp);
}

int  Fp3_cmp(Fp3 *A,Fp3 *B){
    if(Fp_cmp(&A->x0,&B->x0)==0 && Fp_cmp(&A->x1,&B->x1)==0 && Fp_cmp(&A->x2,&B->x2)==0){
        return 0;   
    }
    return 1;
}

int  Fp3_cmp_ui(Fp3 *A,unsigned long int UI){
    if(Fp_cmp_ui(&A->x0,UI)==0 && Fp_cmp_ui(&A->x1,UI)==0 && Fp_cmp_ui(&A->x2,UI)==0){
        return 0;
    }
    return 1;
}

int  Fp3_cmp_mpz(Fp3 *A,mpz_t B){
    if(Fp_cmp_mpz(&A->x0,B)==0 && Fp_cmp_mpz(&A->x1,B)==0 && Fp_cmp_mpz(&A->x2,B)==0){
        return 0;
    }
    return 1;
}

int  Fp3_cmp_zero(Fp3 *A){
    if(Fp_cmp_zero(&A->x0)==0 && Fp_cmp_zero(&A->x1)==0 && Fp_cmp_zero(&A->x2)==0){
        return 0;
    }
    return 1;
}

int  Fp3_cmp_one(Fp3 *A){
    if(Fp_cmp_one(&A->x0)==0 && Fp_cmp_zero(&A->x1)==0 && Fp_cmp_zero(&A->x2)==0){
        return 0;
    }
    return 1;
}

/*----------------------------------------------------------------------------*/
//Fp9
void Fp9_init(Fp9 *A){
    Fp3_init(&A->x0);
    Fp3_init(&A->x1);
    Fp3_init(&A->x2);
}

void Fp9_clear(Fp9 *A){
    Fp3_clear(&A->x0);
    Fp3_clear(&A->x1);
    Fp3_clear(&A->x2);
}

void Fp9_printf(Fp9 *A,char *str){
    gmp_printf("%s(",str);
    Fp3_printf(&A->x0,"");
    gmp_printf(",");
    Fp3_printf(&A->x1,"");
    gmp_printf(",");
    Fp3_printf(&A->x2,"");
    gmp_printf(")");
}

void Fp9_set(Fp9 *ANS,Fp9 *A){
    Fp3_set(&ANS->x0,&A->x0);
    Fp3_set(&ANS->x1,&A->x1);
    Fp3_set(&ANS->x2,&A->x2);
}

void Fp9_set_ui(Fp9 *ANS,unsigned long int UI){
    Fp3_set_ui(&ANS->x0,UI);
    Fp3_set_ui(&ANS->x1,UI);
    Fp3_set_ui(&ANS->x2,UI);
}

void Fp9_set_mpz(Fp9 *ANS,mpz_t A){
    Fp3_set_mpz(&ANS->x0,A);
    Fp3_set_mpz(&ANS->x1,A);
    Fp3_set_mpz(&ANS->x2,A);
}

void Fp9_set_neg(Fp9 *ANS,Fp9 *A){
    Fp3_set_neg(&ANS->x0,&A->x0);
    Fp3_set_neg(&ANS->x1,&A->x1);
    Fp3_set_neg(&ANS->x2,&A->x2);
}

void Fp9_set_random(Fp9 *ANS,gmp_randstate_t state){
    Fp3_set_random(&ANS->x0,state);
    Fp3_set_random(&ANS->x1,state);
    Fp3_set_random(&ANS->x2,state);
}

void Fp9_mul(Fp9 *ANS,Fp9 *A,Fp9 *B){
    Fp3_mul(&TMP1_FP3,&A->x0,&B->x0);
    Fp3_mul(&TMP2_FP3,&A->x1,&B->x1);
    Fp3_mul(&TMP3_FP3,&A->x2,&B->x2);
    
    Fp3_add(&TMP5_FP3,&A->x0,&A->x1);
    Fp3_add(&TMP4_FP3,&B->x0,&B->x1);
    Fp3_mul(&TMP5_FP3,&TMP5_FP3,&TMP4_FP3);
    
    Fp3_add(&TMP6_FP3,&A->x1,&A->x2);
    Fp3_add(&TMP4_FP3,&B->x1,&B->x2);
    Fp3_mul(&TMP6_FP3,&TMP6_FP3,&TMP4_FP3);
    
    Fp3_add(&TMP7_FP3,&B->x0,&B->x2);
    Fp3_add(&TMP4_FP3,&A->x0,&A->x2);
    Fp3_mul(&TMP7_FP3,&TMP7_FP3,&TMP4_FP3);
    //x0
    Fp3_sub(&TMP6_FP3,&TMP6_FP3,&TMP2_FP3);
    Fp3_sub(&TMP6_FP3,&TMP6_FP3,&TMP3_FP3);
    Fp3_mul_basis(&TMP4_FP3,&TMP6_FP3);
    Fp3_add(&ANS->x0,&TMP1_FP3,&TMP4_FP3);
    //x1
    Fp3_sub(&TMP5_FP3,&TMP5_FP3,&TMP1_FP3);
    Fp3_sub(&TMP5_FP3,&TMP5_FP3,&TMP2_FP3);
    Fp3_mul_basis(&TMP4_FP3,&TMP3_FP3);
    Fp3_add(&ANS->x1,&TMP4_FP3,&TMP5_FP3);
    //x2
    Fp3_sub(&TMP7_FP3,&TMP7_FP3,&TMP1_FP3);
    Fp3_sub(&TMP7_FP3,&TMP7_FP3,&TMP3_FP3);
    Fp3_add(&ANS->x2,&TMP2_FP3,&TMP7_FP3);
}

void Fp9_mul_ui(Fp9 *ANS,Fp9 *A,unsigned long int UI){
    Fp3_mul_ui(&ANS->x0,&A->x0,UI);
    Fp3_mul_ui(&ANS->x1,&A->x1,UI);
    Fp3_mul_ui(&ANS->x2,&A->x2,UI);
}

void Fp9_mul_mpz(Fp9 *ANS,Fp9 *A,mpz_t B){
    Fp3_mul_mpz(&ANS->x0,&A->x0,B);
    Fp3_mul_mpz(&ANS->x1,&A->x1,B);
    Fp3_mul_mpz(&ANS->x2,&A->x2,B);
}

void Fp9_mul_basis(Fp9 *ANS,Fp9 *A){
    Fp3_set(&TMP1_FP3,&A->x0);
    Fp3_set(&TMP2_FP3,&A->x1);
    Fp3_set(&TMP3_FP3,&A->x2);
    Fp3_mul_basis(&ANS->x0,&TMP3_FP3);
    Fp3_set(&ANS->x1,&TMP1_FP3);
    Fp3_set(&ANS->x2,&TMP2_FP3);
}

void Fp9_sqr(Fp9 *ANS,Fp9 *A){
    Fp3_sqr(&TMP1_FP3,&A->x0);
    Fp3_sqr(&TMP4_FP3,&A->x2);
    Fp3_add(&TMP5_FP3,&A->x1,&A->x1);
    Fp3_mul(&TMP2_FP3,&TMP5_FP3,&A->x2);
    Fp3_mul(&TMP3_FP3,&A->x0,&TMP5_FP3);
    Fp3_add(&TMP5_FP3,&A->x0,&A->x1);
    Fp3_add(&TMP5_FP3,&TMP5_FP3,&A->x2);
    
    Fp3_mul_basis(&ANS->x0,&TMP2_FP3);
    Fp3_add(&ANS->x0,&ANS->x0,&TMP1_FP3);
    
    Fp3_mul_basis(&ANS->x1,&TMP4_FP3);
    Fp3_add(&ANS->x1,&ANS->x1,&TMP3_FP3);
    
    Fp3_sqr(&ANS->x2,&TMP5_FP3);
    Fp3_add(&TMP5_FP3,&TMP1_FP3,&TMP4_FP3);
    Fp3_add(&TMP5_FP3,&TMP5_FP3,&TMP2_FP3);
    Fp3_add(&TMP5_FP3,&TMP5_FP3,&TMP3_FP3);
    Fp3_sub(&ANS->x2,&ANS->x2,&TMP5_FP3);
}

void Fp9_add(Fp9 *ANS,Fp9 *A,Fp9 *B){
    Fp3_add(&ANS->x0,&A->x0,&B->x0);
    Fp3_add(&ANS->x1,&A->x1,&B->x1);
    Fp3_add(&ANS->x2,&A->x2,&B->x2);
}

void Fp9_add_ui(Fp9 *ANS,Fp9 *A,unsigned long int UI){
    Fp3_add_ui(&ANS->x0,&A->x0,UI);
    Fp3_add_ui(&ANS->x1,&A->x1,UI);
    Fp3_add_ui(&ANS->x2,&A->x2,UI);
}

void Fp9_add_mpz(Fp9 *ANS,Fp9 *A,mpz_t B){
    Fp3_add_mpz(&ANS->x0,&A->x0,B);
    Fp3_add_mpz(&ANS->x1,&A->x1,B);
    Fp3_add_mpz(&ANS->x2,&A->x2,B);
}

void Fp9_sub(Fp9 *ANS,Fp9 *A,Fp9 *B){
    Fp3_sub(&ANS->x0,&A->x0,&B->x0);
    Fp3_sub(&ANS->x1,&A->x1,&B->x1);
    Fp3_sub(&ANS->x2,&A->x2,&B->x2);
}

void Fp9_sub_ui(Fp9 *ANS,Fp9 *A,unsigned long int UI){
    Fp3_sub_ui(&ANS->x0,&A->x0,UI);
    Fp3_sub_ui(&ANS->x1,&A->x1,UI);
    Fp3_sub_ui(&ANS->x2,&A->x2,UI);
}

void Fp9_sub_mpz(Fp9 *ANS,Fp9 *A,mpz_t B){
    Fp3_sub_mpz(&ANS->x0,&A->x0,B);
    Fp3_sub_mpz(&ANS->x1,&A->x1,B);
    Fp3_sub_mpz(&ANS->x2,&A->x2,B);
}

void Fp9_inv(Fp9 *ANS,Fp9 *A){
    Fp3_sqr(&TMP1_FP3,&A->x0);
    Fp3_mul(&TMP4_FP3,&A->x1,&A->x2);
    Fp3_mul_basis(&TMP4_FP3,&TMP4_FP3);
    Fp3_sub(&TMP1_FP3,&TMP1_FP3,&TMP4_FP3);
    
    Fp3_sqr(&TMP2_FP3,&A->x2);
    Fp3_mul_basis(&TMP2_FP3,&TMP2_FP3);
    Fp3_mul(&TMP4_FP3,&A->x0,&A->x1);
    Fp3_sub(&TMP2_FP3,&TMP2_FP3,&TMP4_FP3);
    
    Fp3_sqr(&TMP3_FP3,&A->x1);
    Fp3_mul(&TMP4_FP3,&A->x0,&A->x2);
    Fp3_sub(&TMP3_FP3,&TMP3_FP3,&TMP4_FP3);
    
    Fp3_mul(&TMP5_FP3,&A->x1,&TMP3_FP3);
    Fp3_mul(&TMP4_FP3,&A->x2,&TMP2_FP3);
    Fp3_add(&TMP5_FP3,&TMP5_FP3,&TMP4_FP3);
    Fp3_mul_basis(&TMP5_FP3,&TMP5_FP3);
    Fp3_mul(&TMP4_FP3,&A->x0,&TMP1_FP3);
    Fp3_add(&TMP5_FP3,&TMP5_FP3,&TMP4_FP3);
    
    Fp3_inv(&TMP5_FP3,&TMP5_FP3);
    
    Fp3_mul(&ANS->x0,&TMP5_FP3,&TMP1_FP3);
    Fp3_mul(&ANS->x1,&TMP5_FP3,&TMP2_FP3);
    Fp3_mul(&ANS->x2,&TMP5_FP3,&TMP3_FP3);
}

int  Fp9_legendre(Fp9 *A){
    Fp3_sqr(&TMP2_FP3,&A->x0);
    Fp3_mul(&TMP2_FP3,&TMP2_FP3,&A->x0);
    
    Fp3_sqr(&TMP1_FP3,&A->x1);
    Fp3_mul(&TMP1_FP3,&TMP1_FP3,&A->x1);
    Fp3_mul_basis(&TMP1_FP3,&TMP1_FP3);
    Fp3_add(&TMP2_FP3,&TMP2_FP3,&TMP1_FP3);
    
    Fp3_sqr(&TMP1_FP3,&A->x2);
    Fp3_mul(&TMP1_FP3,&TMP1_FP3,&A->x2);
    Fp3_mul_basis(&TMP1_FP3,&TMP1_FP3);
    Fp3_mul_basis(&TMP1_FP3,&TMP1_FP3);
    Fp3_add(&TMP2_FP3,&TMP2_FP3,&TMP1_FP3);
    
    Fp3_mul(&TMP1_FP3,&A->x0,&A->x1);
    Fp3_mul(&TMP1_FP3,&TMP1_FP3,&A->x2);
    Fp3_mul_basis(&TMP1_FP3,&TMP1_FP3);
    Fp3_sub(&TMP2_FP3,&TMP2_FP3,&TMP1_FP3);
    Fp3_sub(&TMP2_FP3,&TMP2_FP3,&TMP1_FP3);
    Fp3_sub(&TMP2_FP3,&TMP2_FP3,&TMP1_FP3);
    
    return Fp3_legendre(&TMP2_FP3);
}

int  Fp9_isCNR(Fp9 *A){
    Fp9 tmp;
    Fp9_init(&tmp);
    mpz_t exp;
    mpz_init(exp);
    
    mpz_pow_ui(exp,prime,9);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp9_pow(&tmp,A,exp);
    
    if(Fp9_cmp_one(&tmp)==0){
        mpz_clear(exp);
        Fp9_clear(&tmp);
        return 1;
    }else{
        mpz_clear(exp);
        Fp9_clear(&tmp);
        return -1;
    }
}

void Fp9_sqrt(Fp9 *ANS,Fp9 *A){//TODO
    Fp9 tmp_A;
    Fp9_init(&tmp_A);
    Fp3 tmp,s;
    Fp3_init(&tmp);
    Fp3_init(&s);
    mpz_t p3_plus1_over2;
    mpz_init(p3_plus1_over2);
    mpz_pow_ui(p3_plus1_over2,prime,3);
    mpz_add_ui(p3_plus1_over2,p3_plus1_over2,1);
    mpz_tdiv_q_ui(p3_plus1_over2,p3_plus1_over2,2);
    
    //s=A^{p^6+p^3+1}
    Fp3_sqr(&s,&A->x0);
    Fp3_mul(&s,&s,&A->x0);
    
    Fp3_sqr(&tmp,&A->x1);
    Fp3_mul(&tmp,&tmp,&A->x1);
    Fp3_mul_basis(&tmp,&tmp);
    Fp3_add(&s,&s,&tmp);
    
    Fp3_sqr(&tmp,&A->x2);
    Fp3_mul(&tmp,&tmp,&A->x2);
    Fp3_mul_basis(&tmp,&tmp);
    Fp3_mul_basis(&tmp,&tmp);
    Fp3_add(&s,&s,&tmp);
    
    Fp3_mul(&tmp,&A->x0,&A->x1);
    Fp3_mul(&tmp,&tmp,&A->x2);
    Fp3_mul_basis(&tmp,&tmp);
    Fp3_sub(&s,&s,&tmp);
    Fp3_sub(&s,&s,&tmp);
    Fp3_sub(&s,&s,&tmp);
    
    //sqrt{s}^{-1}
    Fp3_sqrt(&s,&s);
    Fp3_inv(&s,&s);
    
    //A^{(p^6+p^3+2)/2}
    Fp3_set(&tmp_A.x0,&A->x0);
    Fp3_mul_mpz(&tmp_A.x1,&A->x1,frobenius_constant[f_p3][3].x0);
    Fp3_mul_mpz(&tmp_A.x2,&A->x2,frobenius_constant[f_p3][6].x0);
    Fp9_pow(&tmp_A,&tmp_A,p3_plus1_over2);
    Fp9_mul(&tmp_A,&tmp_A,A);
    
    Fp3_mul(&ANS->x0,&tmp_A.x0,&s);
    Fp3_mul(&ANS->x1,&tmp_A.x1,&s);
    Fp3_mul(&ANS->x2,&tmp_A.x2,&s);
    
    Fp3_clear(&tmp);
    Fp3_clear(&s);
    Fp9_clear(&tmp_A);
    mpz_clear(p3_plus1_over2);
}

void Fp9_pow(Fp9 *ANS,Fp9 *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    Fp9 tmp;
    Fp9_init(&tmp);
    Fp9_set(&tmp,A);
    
    for(i=1; binary[i]!='\0'; i++){
        Fp9_sqr(&tmp,&tmp);
        if(binary[i]=='1'){
            Fp9_mul(&tmp,A,&tmp);
        }
    }
    
    Fp9_set(ANS,&tmp);
    Fp9_clear(&tmp);
}

int  Fp9_cmp(Fp9 *A,Fp9 *B){
    if(Fp3_cmp(&A->x0,&B->x0)==0 && Fp3_cmp(&A->x1,&B->x1)==0 && Fp3_cmp(&A->x2,&B->x2)==0){
        return 0;   
    }
    return 1;
}

int  Fp9_cmp_ui(Fp9 *A,unsigned long int UI){
    if(Fp3_cmp_ui(&A->x0,UI)==0 && Fp3_cmp_ui(&A->x1,UI)==0 && Fp3_cmp_ui(&A->x2,UI)==0){
        return 0;
    }
    return 1;
}

int  Fp9_cmp_mpz(Fp9 *A,mpz_t B){
    if(Fp3_cmp_mpz(&A->x0,B)==0 && Fp3_cmp_mpz(&A->x1,B)==0 && Fp3_cmp_mpz(&A->x2,B)==0){
        return 0;
    }
    return 1;
}

int  Fp9_cmp_zero(Fp9 *A){
    if(Fp3_cmp_zero(&A->x0)==0 && Fp3_cmp_zero(&A->x1)==0 && Fp3_cmp_zero(&A->x2)==0){
        return 0;
    }
    return 1;
}

int  Fp9_cmp_one(Fp9 *A){
    if(Fp3_cmp_one(&A->x0)==0 && Fp3_cmp_zero(&A->x1)==0 && Fp3_cmp_zero(&A->x2)==0){
        return 0;
    }
    return 1;
}
/*----------------------------------------------------------------------------*/
//Fp18
void Fp18_init(Fp18 *A){
    Fp9_init(&A->x0);
    Fp9_init(&A->x1);
}

void Fp18_clear(Fp18 *A){
    Fp9_clear(&A->x0);
    Fp9_clear(&A->x1);
}

void Fp18_printf(Fp18 *A,char *str){
    gmp_printf("%s(",str);
    Fp9_printf(&A->x0,"");
    gmp_printf(",");
    Fp9_printf(&A->x1,"");
    gmp_printf(")");
}

void Fp18_set(Fp18 *ANS,Fp18 *A){
    Fp9_set(&ANS->x0,&A->x0);
    Fp9_set(&ANS->x1,&A->x1);
}

void Fp18_set_ui(Fp18 *ANS,unsigned long int UI){
    Fp9_set_ui(&ANS->x0,UI);
    Fp9_set_ui(&ANS->x1,UI);
}

void Fp18_set_mpz(Fp18 *ANS,mpz_t A){
    Fp9_set_mpz(&ANS->x0,A);
    Fp9_set_mpz(&ANS->x1,A);    
}

void Fp18_set_neg(Fp18 *ANS,Fp18 *A){
    Fp9_set_neg(&ANS->x0,&A->x0);
    Fp9_set_neg(&ANS->x1,&A->x1);
}

void Fp18_set_random(Fp18 *ANS,gmp_randstate_t state){
    Fp9_set_random(&ANS->x0,state);
    Fp9_set_random(&ANS->x1,state);
}

void Fp18_mul(Fp18 *ANS,Fp18 *A,Fp18 *B){
    Fp9_mul(&TMP2_FP9,&A->x1,&B->x1);
	Fp9_add(&TMP1_FP9,&A->x0,&A->x1);
	Fp9_add(&ANS->x1,&B->x0,&B->x1);
	Fp9_mul(&ANS->x1,&TMP1_FP9,&ANS->x1);
	Fp9_mul(&TMP1_FP9,&A->x0,&B->x0);
	
	Fp9_mul_basis(&ANS->x0,&TMP2_FP9);
	Fp9_add(&ANS->x0,&ANS->x0,&TMP1_FP9);
	
	Fp9_sub(&ANS->x1,&ANS->x1,&TMP1_FP9);
	Fp9_sub(&ANS->x1,&ANS->x1,&TMP2_FP9);
}

void Fp18_mul_ui(Fp18 *ANS,Fp18 *A,unsigned long int UI){
    Fp9_mul_ui(&ANS->x0,&A->x0,UI);
    Fp9_mul_ui(&ANS->x1,&A->x1,UI);
}

void Fp18_mul_mpz(Fp18 *ANS,Fp18 *A,mpz_t B){
    Fp9_mul_mpz(&ANS->x0,&A->x0,B);
    Fp9_mul_mpz(&ANS->x1,&A->x1,B);
}

void Fp18_sqr(Fp18 *ANS,Fp18 *A){
    Fp9_add(&TMP1_FP9,&A->x0,&A->x1);
	Fp9_mul_basis(&TMP2_FP9,&A->x1);
	Fp9_add(&TMP2_FP9,&TMP2_FP9,&A->x0);
	Fp9_mul(&TMP3_FP9,&A->x0,&A->x1);
	
	Fp9_mul(&ANS->x0,&TMP1_FP9,&TMP2_FP9);
	Fp9_sub(&ANS->x0,&ANS->x0,&TMP3_FP9);
	Fp9_mul_basis(&TMP1_FP9,&TMP3_FP9);
	Fp9_sub(&ANS->x0,&ANS->x0,&TMP1_FP9);
	
	Fp9_add(&ANS->x1,&TMP3_FP9,&TMP3_FP9);
}

void Fp18_sqr_cyclotomic(Fp18 *ANS,Fp18 *A){
    Fp9_add(&TMP1_FP9,&A->x0,&A->x1);
    Fp9_sqr(&TMP1_FP9,&TMP1_FP9);
    Fp9_sqr(&TMP2_FP9,&A->x1);
    Fp9_mul_basis(&ANS->x1,&TMP2_FP9);
    Fp9_add(&ANS->x0,&ANS->x1,&ANS->x1);
    Fp_add_ui(&ANS->x0.x0.x0,&ANS->x0.x0.x0,1);
    
    Fp9_sub(&ANS->x1,&TMP1_FP9,&ANS->x1);
    Fp9_sub(&ANS->x1,&ANS->x1,&TMP2_FP9);
    Fp_sub_ui(&ANS->x1.x0.x0,&ANS->x1.x0.x0,1);
}

void Fp18_add(Fp18 *ANS,Fp18 *A,Fp18 *B){
    Fp9_add(&ANS->x0,&A->x0,&B->x0);
    Fp9_add(&ANS->x1,&A->x1,&B->x1);
}

void Fp18_add_ui(Fp18 *ANS,Fp18 *A,unsigned long int UI){
    Fp9_add_ui(&ANS->x0,&A->x0,UI);
    Fp9_add_ui(&ANS->x1,&A->x1,UI);
}

void Fp18_add_mpz(Fp18 *ANS,Fp18 *A,mpz_t B){
    Fp9_add_mpz(&ANS->x0,&ANS->x0,B);
    Fp9_add_mpz(&ANS->x1,&ANS->x1,B);
}

void Fp18_sub(Fp18 *ANS,Fp18 *A,Fp18 *B){
    Fp9_sub(&ANS->x0,&A->x0,&B->x0);
    Fp9_sub(&ANS->x1,&A->x1,&B->x1);
}

void Fp18_sub_ui(Fp18 *ANS,Fp18 *A,unsigned long int UI){
    Fp9_sub_ui(&ANS->x0,&ANS->x0,UI);
    Fp9_sub_ui(&ANS->x1,&ANS->x1,UI);
}

void Fp18_sub_mpz(Fp18 *ANS,Fp18 *A,mpz_t B){
    Fp9_sub_mpz(&ANS->x0,&ANS->x0,B);
    Fp9_sub_mpz(&ANS->x1,&ANS->x1,B);
}

void Fp18_inv(Fp18 *ANS,Fp18 *A){
    Fp9_set(&TMP1_FP9,&A->x0);
    Fp9_set_neg(&TMP2_FP9,&A->x1);
    
    Fp9_mul(&TMP3_FP9,&TMP1_FP9,&A->x0);
    Fp9_mul(&TMP4_FP9,&TMP2_FP9,&A->x1);
    Fp9_mul_basis(&TMP4_FP9,&TMP4_FP9);
    Fp9_add(&TMP3_FP9,&TMP3_FP9,&TMP4_FP9);
    Fp9_inv(&TMP3_FP9,&TMP3_FP9);
    Fp9_mul(&ANS->x0,&TMP1_FP9,&TMP3_FP9);
    Fp9_mul(&ANS->x1,&TMP2_FP9,&TMP3_FP9);
}

int  Fp18_legendre(Fp18 *A){
    Fp9_sqr(&TMP1_FP9,&A->x0);
    Fp9_sqr(&TMP2_FP9,&A->x1);
    Fp9_mul_basis(&TMP2_FP9,&TMP2_FP9);
    Fp9_sub(&TMP1_FP9,&TMP1_FP9,&TMP2_FP9);
    
    return Fp9_legendre(&TMP1_FP9);
}

void Fp18_sqrt(Fp18 *ANS,Fp18 *A){
    Fp18 x,y,t,k,n,tmp;
    Fp18_init(&x);
    Fp18_init(&y);
    Fp18_init(&t);
    Fp18_init(&k);
    Fp18_init(&n);
    Fp18_init(&tmp);
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    Fp18_set_random(&n,state);
    while(Fp18_legendre(&n)!=-1){
        Fp18_set_random(&n,state);
    }
    mpz_pow_ui(q,prime,18);
    mpz_sub_ui(q,q,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    Fp18_pow(&y,&n,q);
    mpz_set_ui(z,e);    
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp18_pow(&x,A,exp);
    Fp18_mul(&tmp,&x,&x);
    Fp18_mul(&k,&tmp,A);
    Fp18_mul(&x,&x,A);
    while(Fp18_cmp_one(&k)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        Fp18_pow(&tmp,&k,exp);
        while(Fp18_cmp_one(&tmp)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            Fp18_pow(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        Fp18_pow(&t,&y,result);
        Fp18_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        Fp18_mul(&x,&x,&t); 
        Fp18_mul(&k,&k,&y);
    }
    Fp18_set(ANS,&x);
    
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(exp);
    mpz_clear(result);
    Fp18_clear(&x);
    Fp18_clear(&y);
    Fp18_clear(&t);
    Fp18_clear(&k);
    Fp18_clear(&n);
    Fp18_clear(&tmp);
}

void Fp18_pow(Fp18 *ANS,Fp18 *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    Fp18 tmp;
    Fp18_init(&tmp);
    Fp18_set(&tmp,A);
    
    for(i=1; binary[i]!='\0'; i++){
        Fp18_sqr(&tmp,&tmp);
        if(binary[i]=='1'){
            Fp18_mul(&tmp,A,&tmp);
        }
    }
    
    Fp18_set(ANS,&tmp);
    Fp18_clear(&tmp);
}

int  Fp18_cmp(Fp18 *A,Fp18 *B){
    if(Fp9_cmp(&A->x0,&B->x0)==0 && Fp9_cmp(&A->x1,&B->x1)==0){
        return 0;   
    }
    return 1;
}

int  Fp18_cmp_ui(Fp18 *A,unsigned long int UI){
    if(Fp9_cmp_ui(&A->x0,UI)==0 && Fp9_cmp_ui(&A->x1,UI)==0){
        return 0;
    }
    return 1;
}

int  Fp18_cmp_mpz(Fp18 *A,mpz_t B){
    if(Fp9_cmp_mpz(&A->x0,B)==0 && Fp9_cmp_mpz(&A->x1,B)==0){
        return 0;
    }
    return 1;
}

int  Fp18_cmp_zero(Fp18 *A){
    if(Fp9_cmp_zero(&A->x0)==0 && Fp9_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}

int  Fp18_cmp_one(Fp18 *A){
    if(Fp9_cmp_one(&A->x0)==0 && Fp9_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}

void Fp18_frobenius_map_p1(Fp18 *ANS,Fp18 *A){
    //1
    Fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    Fp_mul(&ANS->x0.x0.x1,&A->x0.x0.x1,&frobenius_constant[f_p1][1]);
    Fp_mul(&ANS->x0.x0.x2,&A->x0.x0.x2,&frobenius_constant[f_p1][2]);
    
    //beta
    Fp_mul(&ANS->x0.x1.x0,&A->x0.x1.x0,&frobenius_constant[f_p1][3]);
    Fp_mul(&ANS->x0.x1.x1,&A->x0.x1.x1,&frobenius_constant[f_p1][4]);
    Fp_mul(&ANS->x0.x1.x2,&A->x0.x1.x2,&frobenius_constant[f_p1][5]);
    Fp3_set(&TMP1_FP3,&ANS->x0.x1);
    Fp_add(&ANS->x0.x1.x0,&TMP1_FP3.x1,&TMP1_FP3.x1);
    Fp_add(&ANS->x0.x1.x1,&TMP1_FP3.x2,&TMP1_FP3.x2);
    Fp_set(&ANS->x0.x1.x2,&TMP1_FP3.x0);
    
    //beta^2
    Fp_mul(&ANS->x0.x2.x0,&A->x0.x2.x0,&frobenius_constant[f_p1][6]);
    Fp_mul(&ANS->x0.x2.x1,&A->x0.x2.x1,&frobenius_constant[f_p1][7]);
    Fp_mul(&ANS->x0.x2.x2,&A->x0.x2.x2,&frobenius_constant[f_p1][8]);
    Fp3_set(&TMP1_FP3,&ANS->x0.x2);
    Fp_add(&ANS->x0.x2.x0,&TMP1_FP3.x2,&TMP1_FP3.x2);
    Fp_set(&ANS->x0.x2.x1,&TMP1_FP3.x0);
    Fp_set(&ANS->x0.x2.x2,&TMP1_FP3.x1);
    
    //gamma
    Fp_mul(&ANS->x1.x0.x0,&A->x1.x0.x0,&frobenius_constant[f_p1][9]);
    Fp_mul(&ANS->x1.x0.x1,&A->x1.x0.x1,&frobenius_constant[f_p1][10]);
    Fp_mul(&ANS->x1.x0.x2,&A->x1.x0.x2,&frobenius_constant[f_p1][11]);
    Fp3_set(&TMP1_FP3,&ANS->x1.x0);
    Fp_add(&ANS->x1.x0.x0,&TMP1_FP3.x2,&TMP1_FP3.x2);
    Fp_set(&ANS->x1.x0.x1,&TMP1_FP3.x0);
    Fp_set(&ANS->x1.x0.x2,&TMP1_FP3.x1);
    
    //betagamma
    Fp_mul(&ANS->x1.x1.x0,&A->x1.x1.x0,&frobenius_constant[f_p1][12]);
    Fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
    Fp_mul(&ANS->x1.x1.x2,&A->x1.x1.x2,&frobenius_constant[f_p1][14]);
    
    //beta^2gamma
    Fp_mul(&ANS->x1.x2.x0,&A->x1.x2.x0,&frobenius_constant[f_p1][15]);
    Fp_mul(&ANS->x1.x2.x1,&A->x1.x2.x1,&frobenius_constant[f_p1][16]);
    Fp_mul(&ANS->x1.x2.x2,&A->x1.x2.x2,&frobenius_constant[f_p1][17]);
    Fp3_set(&TMP1_FP3,&ANS->x1.x2);
    Fp_add(&ANS->x1.x2.x0,&TMP1_FP3.x1,&TMP1_FP3.x1);
    Fp_add(&ANS->x1.x2.x1,&TMP1_FP3.x2,&TMP1_FP3.x2);
    Fp_set(&ANS->x1.x2.x2,&TMP1_FP3.x0);
}

void Fp18_frobenius_map_p2(Fp18 *ANS,Fp18 *A){
    //1
    Fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    Fp_mul(&ANS->x0.x0.x1,&A->x0.x0.x1,&frobenius_constant[f_p2][1]);
    Fp_mul(&ANS->x0.x0.x2,&A->x0.x0.x2,&frobenius_constant[f_p2][2]);
    
    //beta
    Fp_mul(&ANS->x0.x1.x0,&A->x0.x1.x0,&frobenius_constant[f_p2][3]);
    Fp_mul(&ANS->x0.x1.x1,&A->x0.x1.x1,&frobenius_constant[f_p2][4]);
    Fp_mul(&ANS->x0.x1.x2,&A->x0.x1.x2,&frobenius_constant[f_p2][5]);
    Fp3_set(&TMP1_FP3,&ANS->x0.x1);
    Fp_add(&ANS->x0.x1.x0,&TMP1_FP3.x2,&TMP1_FP3.x2);
    Fp_set(&ANS->x0.x1.x1,&TMP1_FP3.x0);
    Fp_set(&ANS->x0.x1.x2,&TMP1_FP3.x1);
    
    //beta^2
    Fp_mul(&ANS->x0.x2.x0,&A->x0.x2.x0,&frobenius_constant[f_p2][6]);
    Fp_mul(&ANS->x0.x2.x1,&A->x0.x2.x1,&frobenius_constant[f_p2][7]);
    Fp_mul(&ANS->x0.x2.x2,&A->x0.x2.x2,&frobenius_constant[f_p2][8]);
    Fp3_set(&TMP1_FP3,&ANS->x0.x2);
    Fp_add(&ANS->x0.x2.x0,&TMP1_FP3.x1,&TMP1_FP3.x1);
    Fp_add(&ANS->x0.x2.x1,&TMP1_FP3.x2,&TMP1_FP3.x2);
    Fp_set(&ANS->x0.x2.x2,&TMP1_FP3.x0);
    
    //gamma
    Fp_mul(&ANS->x1.x0.x0,&A->x1.x0.x0,&frobenius_constant[f_p2][9]);
    Fp_mul(&ANS->x1.x0.x1,&A->x1.x0.x1,&frobenius_constant[f_p2][10]);
    Fp_mul(&ANS->x1.x0.x2,&A->x1.x0.x2,&frobenius_constant[f_p2][11]);
    Fp3_set(&TMP1_FP3,&ANS->x1.x0);
    Fp_add(&ANS->x1.x0.x0,&TMP1_FP3.x1,&TMP1_FP3.x1);
    Fp_add(&ANS->x1.x0.x1,&TMP1_FP3.x2,&TMP1_FP3.x2);
    Fp_set(&ANS->x1.x0.x2,&TMP1_FP3.x0);
    
    //betagamma
    Fp_mul(&ANS->x1.x1.x0,&A->x1.x1.x0,&frobenius_constant[f_p2][12]);
    Fp_set(&ANS->x1.x1.x1,&A->x1.x1.x1);
    Fp_mul(&ANS->x1.x1.x2,&A->x1.x1.x2,&frobenius_constant[f_p2][14]);
    
    //beta^2gamma
    Fp_mul(&ANS->x1.x2.x0,&A->x1.x2.x0,&frobenius_constant[f_p2][15]);
    Fp_mul(&ANS->x1.x2.x1,&A->x1.x2.x1,&frobenius_constant[f_p2][16]);
    Fp_mul(&ANS->x1.x2.x2,&A->x1.x2.x2,&frobenius_constant[f_p2][17]);
    Fp3_set(&TMP1_FP3,&ANS->x1.x2);
    Fp_add(&ANS->x1.x2.x0,&TMP1_FP3.x2,&TMP1_FP3.x2);
    Fp_set(&ANS->x1.x2.x1,&TMP1_FP3.x0);
    Fp_set(&ANS->x1.x2.x2,&TMP1_FP3.x1);
}

void Fp18_frobenius_map_p3(Fp18 *ANS,Fp18 *A){
    Fp3_set(&ANS->x0.x0,&A->x0.x0);
    Fp3_mul_mpz(&ANS->x0.x1,&A->x0.x1,frobenius_constant[f_p3][3].x0);
    Fp3_mul_mpz(&ANS->x0.x2,&A->x0.x2,frobenius_constant[f_p3][6].x0);
    Fp3_mul_mpz(&ANS->x1.x0,&A->x1.x0,frobenius_constant[f_p3][9].x0);
    Fp3_set_neg(&ANS->x1.x1,&A->x1.x1);
    Fp3_mul_mpz(&ANS->x1.x2,&A->x1.x2,frobenius_constant[f_p3][15].x0);
}

void Fp18_frobenius_map_p4(Fp18 *ANS,Fp18 *A){
    //1
    Fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    Fp_mul(&ANS->x0.x0.x1,&A->x0.x0.x1,&frobenius_constant[f_p4][1]);
    Fp_mul(&ANS->x0.x0.x2,&A->x0.x0.x2,&frobenius_constant[f_p4][2]);
    
    //beta
    Fp_mul(&ANS->x0.x1.x0,&A->x0.x1.x0,&frobenius_constant[f_p4][3]);
    Fp_mul(&ANS->x0.x1.x1,&A->x0.x1.x1,&frobenius_constant[f_p4][4]);
    Fp_mul(&ANS->x0.x1.x2,&A->x0.x1.x2,&frobenius_constant[f_p4][5]);
    Fp3_set(&TMP1_FP3,&ANS->x0.x1);
    Fp_add(&ANS->x0.x1.x0,&TMP1_FP3.x1,&TMP1_FP3.x1);
    Fp_add(&ANS->x0.x1.x1,&TMP1_FP3.x2,&TMP1_FP3.x2);
    Fp_set(&ANS->x0.x1.x2,&TMP1_FP3.x0);
    
    //beta^2
    Fp_mul(&ANS->x0.x2.x0,&A->x0.x2.x0,&frobenius_constant[f_p4][6]);
    Fp_mul(&ANS->x0.x2.x1,&A->x0.x2.x1,&frobenius_constant[f_p4][7]);
    Fp_mul(&ANS->x0.x2.x2,&A->x0.x2.x2,&frobenius_constant[f_p4][8]);
    Fp3_set(&TMP1_FP3,&ANS->x0.x2);
    Fp_add(&ANS->x0.x2.x0,&TMP1_FP3.x2,&TMP1_FP3.x2);
    Fp_set(&ANS->x0.x2.x1,&TMP1_FP3.x0);
    Fp_set(&ANS->x0.x2.x2,&TMP1_FP3.x1);
    
    //gamma
    Fp_mul(&ANS->x1.x0.x0,&A->x1.x0.x0,&frobenius_constant[f_p4][9]);
    Fp_mul(&ANS->x1.x0.x1,&A->x1.x0.x1,&frobenius_constant[f_p4][10]);
    Fp_mul(&ANS->x1.x0.x2,&A->x1.x0.x2,&frobenius_constant[f_p4][11]);
    Fp3_set(&TMP1_FP3,&ANS->x1.x0);
    Fp_add(&ANS->x1.x0.x0,&TMP1_FP3.x2,&TMP1_FP3.x2);
    Fp_set(&ANS->x1.x0.x1,&TMP1_FP3.x0);
    Fp_set(&ANS->x1.x0.x2,&TMP1_FP3.x1);
    
    //betagamma
    Fp_mul(&ANS->x1.x1.x0,&A->x1.x1.x0,&frobenius_constant[f_p4][12]);
    Fp_set(&ANS->x1.x1.x1,&A->x1.x1.x1);
    Fp_mul(&ANS->x1.x1.x2,&A->x1.x1.x2,&frobenius_constant[f_p4][14]);
    
    //beta^2gamma
    Fp_mul(&ANS->x1.x2.x0,&A->x1.x2.x0,&frobenius_constant[f_p4][15]);
    Fp_mul(&ANS->x1.x2.x1,&A->x1.x2.x1,&frobenius_constant[f_p4][16]);
    Fp_mul(&ANS->x1.x2.x2,&A->x1.x2.x2,&frobenius_constant[f_p4][17]);
    Fp3_set(&TMP1_FP3,&ANS->x1.x2);
    Fp_add(&ANS->x1.x2.x0,&TMP1_FP3.x1,&TMP1_FP3.x1);
    Fp_add(&ANS->x1.x2.x1,&TMP1_FP3.x2,&TMP1_FP3.x2);
    Fp_set(&ANS->x1.x2.x2,&TMP1_FP3.x0);
}

void Fp18_frobenius_map_p5(Fp18 *ANS,Fp18 *A){
    //1
    Fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    Fp_mul(&ANS->x0.x0.x1,&A->x0.x0.x1,&frobenius_constant[f_p5][1]);
    Fp_mul(&ANS->x0.x0.x2,&A->x0.x0.x2,&frobenius_constant[f_p5][2]);
    
    //beta
    Fp_mul(&ANS->x0.x1.x0,&A->x0.x1.x0,&frobenius_constant[f_p5][3]);
    Fp_mul(&ANS->x0.x1.x1,&A->x0.x1.x1,&frobenius_constant[f_p5][4]);
    Fp_mul(&ANS->x0.x1.x2,&A->x0.x1.x2,&frobenius_constant[f_p5][5]);
    Fp3_set(&TMP1_FP3,&ANS->x0.x1);
    Fp_add(&ANS->x0.x1.x0,&TMP1_FP3.x2,&TMP1_FP3.x2);
    Fp_set(&ANS->x0.x1.x1,&TMP1_FP3.x0);
    Fp_set(&ANS->x0.x1.x2,&TMP1_FP3.x1);
    
    //beta^2
    Fp_mul(&ANS->x0.x2.x0,&A->x0.x2.x0,&frobenius_constant[f_p5][6]);
    Fp_mul(&ANS->x0.x2.x1,&A->x0.x2.x1,&frobenius_constant[f_p5][7]);
    Fp_mul(&ANS->x0.x2.x2,&A->x0.x2.x2,&frobenius_constant[f_p5][8]);
    Fp3_set(&TMP1_FP3,&ANS->x0.x2);
    Fp_add(&ANS->x0.x2.x0,&TMP1_FP3.x1,&TMP1_FP3.x1);
    Fp_add(&ANS->x0.x2.x1,&TMP1_FP3.x2,&TMP1_FP3.x2);
    Fp_set(&ANS->x0.x2.x2,&TMP1_FP3.x0);
    
    //gamma
    Fp_mul(&ANS->x1.x0.x0,&A->x1.x0.x0,&frobenius_constant[f_p5][9]);
    Fp_mul(&ANS->x1.x0.x1,&A->x1.x0.x1,&frobenius_constant[f_p5][10]);
    Fp_mul(&ANS->x1.x0.x2,&A->x1.x0.x2,&frobenius_constant[f_p5][11]);
    Fp3_set(&TMP1_FP3,&ANS->x1.x0);
    Fp_add(&ANS->x1.x0.x0,&TMP1_FP3.x1,&TMP1_FP3.x1);
    Fp_add(&ANS->x1.x0.x1,&TMP1_FP3.x2,&TMP1_FP3.x2);
    Fp_set(&ANS->x1.x0.x2,&TMP1_FP3.x0);
    
    //betagamma
    Fp_mul(&ANS->x1.x1.x0,&A->x1.x1.x0,&frobenius_constant[f_p5][12]);
    Fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
    Fp_mul(&ANS->x1.x1.x2,&A->x1.x1.x2,&frobenius_constant[f_p5][14]);
    
    //beta^2gamma
    Fp_mul(&ANS->x1.x2.x0,&A->x1.x2.x0,&frobenius_constant[f_p5][15]);
    Fp_mul(&ANS->x1.x2.x1,&A->x1.x2.x1,&frobenius_constant[f_p5][16]);
    Fp_mul(&ANS->x1.x2.x2,&A->x1.x2.x2,&frobenius_constant[f_p5][17]);
    Fp3_set(&TMP1_FP3,&ANS->x1.x2);
    Fp_add(&ANS->x1.x2.x0,&TMP1_FP3.x2,&TMP1_FP3.x2);
    Fp_set(&ANS->x1.x2.x1,&TMP1_FP3.x0);
    Fp_set(&ANS->x1.x2.x2,&TMP1_FP3.x1);
}

void Fp18_frobenius_map_p6(Fp18 *ANS,Fp18 *A){
    Fp3_set(&ANS->x0.x0,&A->x0.x0);
    Fp3_mul_mpz(&ANS->x0.x1,&A->x0.x1,frobenius_constant[f_p6][3].x0);
    Fp3_mul_mpz(&ANS->x0.x2,&A->x0.x2,frobenius_constant[f_p6][6].x0);
    Fp3_mul_mpz(&ANS->x1.x0,&A->x1.x0,frobenius_constant[f_p6][9].x0);
    Fp3_set(&ANS->x1.x1,&A->x1.x1);
    Fp3_mul_mpz(&ANS->x1.x2,&A->x1.x2,frobenius_constant[f_p6][15].x0);
}

void Fp18_frobenius_map_p9(Fp18 *ANS,Fp18 *A){
    Fp9_set(&ANS->x0,&A->x0);
    Fp9_set_neg(&ANS->x1,&A->x1);
}

/*============================================================================*/
/* Elliptic curve                                                             */
/*============================================================================*/
/*----------------------------------------------------------------------------*/
//EFp
void EFp_init(EFp *P){
    Fp_init(&P->x);
    Fp_init(&P->y);
    P->infinity=0;
}

void EFp_set(EFp *ANS,EFp *A){
    Fp_set(&ANS->x,&A->x);
    Fp_set(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp_set_ui(EFp *ANS,unsigned long int UI){
    Fp_set_ui(&ANS->x,UI);
    Fp_set_ui(&ANS->y,UI);
    ANS->infinity=0;
}

void EFp_set_mpz(EFp *ANS,mpz_t A){
    Fp_set_mpz(&ANS->x,A);
    Fp_set_mpz(&ANS->y,A);
    ANS->infinity=0;
}

void EFp_set_neg(EFp *ANS,EFp *A){
    Fp_set(&ANS->x,&A->x);
    Fp_set_neg(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp_clear(EFp *P){
    Fp_clear(&P->x);
    Fp_clear(&P->y);
}

void EFp_printf(EFp *P,char *str){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp_printf(&P->x,"");
        printf(",");
        Fp_printf(&P->y,"");
        printf(")");
    }else{
        printf("0");
    }
}

void EFp_rational_point(EFp *P){
    Fp tmp1,tmp2,tmp_x;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    Fp_init(&tmp_x);
    gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
	
    while(1){
        Fp_set_random(&P->x,state);
        Fp_mul(&tmp1,&P->x,&P->x);
        Fp_mul(&tmp2,&tmp1,&P->x);
        Fp_add(&tmp_x,&tmp2,&curve_b);
        if(Fp_legendre(&tmp_x)==1){
            Fp_sqrt(&P->y,&tmp_x);
            break;
        }
    }
    
    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
    Fp_clear(&tmp_x);
}

void EFp_ECD(EFp *ANS,EFp *P){
    if(Fp_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFp_set(&TMP1_EFP,P);
        
    Fp_add(&TMP1_FP,&TMP1_EFP.y,&TMP1_EFP.y);
    Fp_inv(&TMP1_FP,&TMP1_FP);
    Fp_mul(&TMP2_FP,&TMP1_EFP.x,&TMP1_EFP.x);
    Fp_add(&TMP3_FP,&TMP2_FP,&TMP2_FP);
    Fp_add(&TMP2_FP,&TMP2_FP,&TMP3_FP);
    Fp_mul(&TMP3_FP,&TMP1_FP,&TMP2_FP);
    Fp_mul(&TMP1_FP,&TMP3_FP,&TMP3_FP);
    Fp_add(&TMP2_FP,&TMP1_EFP.x,&TMP1_EFP.x);
    Fp_sub(&ANS->x,&TMP1_FP,&TMP2_FP);
    Fp_sub(&TMP1_FP,&TMP1_EFP.x,&ANS->x);
    Fp_mul(&TMP2_FP,&TMP3_FP,&TMP1_FP);
    Fp_sub(&ANS->y,&TMP2_FP,&TMP1_EFP.y);
}

void EFp_ECA(EFp *ANS,EFp *P1,EFp *P2){
    if(P1->infinity==1){
        EFp_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp_set(ANS,P1);
        return;
    }else if(Fp_cmp(&P1->x,&P2->x)==0){
        if(Fp_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp_ECD(ANS,P1);
            return;
        }
    }
    
    EFp_set(&TMP1_EFP,P1);
    EFp_set(&TMP2_EFP,P2);
    
    Fp_sub(&TMP1_FP,&TMP2_EFP.x,&TMP1_EFP.x);
    Fp_inv(&TMP1_FP,&TMP1_FP);
    Fp_sub(&TMP2_FP,&TMP2_EFP.y,&TMP1_EFP.y);
    Fp_mul(&TMP3_FP,&TMP1_FP,&TMP2_FP);
    Fp_mul(&TMP1_FP,&TMP3_FP,&TMP3_FP);
    Fp_sub(&TMP2_FP,&TMP1_FP,&TMP1_EFP.x);
    Fp_sub(&ANS->x,&TMP2_FP,&TMP2_EFP.x);
    Fp_sub(&TMP1_FP,&TMP1_EFP.x,&ANS->x);
    Fp_mul(&TMP2_FP,&TMP3_FP,&TMP1_FP);
    Fp_sub(&ANS->y,&TMP2_FP,&TMP1_EFP.y);
}

void EFp_SCM(EFp *ANS,EFp *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp_set(ANS,P);
        return;
    }
    
    EFp Tmp_P,Next_P;
    EFp_init(&Tmp_P);
    EFp_set(&Tmp_P,P);
    EFp_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    
    EFp_set(&Next_P,&Tmp_P);
    for(i=1; binary[i]!='\0'; i++){
        EFp_ECD(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp_ECA(&Next_P,&Next_P,&Tmp_P);
        }
    }
    
    EFp_set(ANS,&Next_P);
    
    EFp_clear(&Next_P);
    EFp_clear(&Tmp_P);
}

//skew frobenius map
void EFp_skew_frobenius_map_p3(EFp *ANS,EFp *A){
    Fp_mul_mpz(&ANS->x,&A->x,epsilon2);
    Fp_set_neg(&ANS->y,&A->y);
}

/*----------------------------------------------------------------------------*/
//EFp3
void EFp3_init(EFp3 *P){
    Fp3_init(&P->x);
    Fp3_init(&P->y);
    P->infinity=0;
}

void EFp3_set(EFp3 *ANS,EFp3 *A){
    Fp3_set(&ANS->x,&A->x);
    Fp3_set(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp3_set_ui(EFp3 *ANS,unsigned long int UI){
    Fp3_set_ui(&ANS->x,UI);
    Fp3_set_ui(&ANS->y,UI);
    ANS->infinity=0;
}

void EFp3_set_mpz(EFp3 *ANS,mpz_t A){
    Fp3_set_mpz(&ANS->x,A);
    Fp3_set_mpz(&ANS->y,A);
    ANS->infinity=0;
}

void EFp3_set_neg(EFp3 *ANS,EFp3 *A){
    Fp3_set(&ANS->x,&A->x);
    Fp3_set_neg(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp3_clear(EFp3 *P){
    Fp3_clear(&P->x);
    Fp3_clear(&P->y);
}

void EFp3_printf(EFp3 *P,char *str){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp3_printf(&P->x,"");
        printf(",");
        Fp3_printf(&P->y,"");
        printf(")");
    }else{
        printf("0");
    }
}

void EFp3_rational_point(EFp3 *P){
    Fp3 tmp1,tmp2;
    Fp3_init(&tmp1);
    Fp3_init(&tmp2);
    gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    while(1){
        Fp3_set_random(&P->x,state);
        Fp3_sqr(&tmp1,&P->x);
        Fp3_mul(&tmp2,&tmp1,&P->x);
        Fp_add(&tmp2.x0,&tmp2.x0,&curve_b);
        if(Fp3_legendre(&tmp2)==1){
            Fp3_sqrt(&P->y,&tmp2);
            break;
        }
    }
    
    Fp3_clear(&tmp1);
    Fp3_clear(&tmp2);
}

void EFp3_ECD(EFp3 *ANS,EFp3 *P){
    if(Fp3_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFp3_set(&TMP1_EFP3,P);
    
    Fp3_add(&TMP1_FP3,&TMP1_EFP3.y,&TMP1_EFP3.y);
    
    Fp3_inv(&TMP1_FP3,&TMP1_FP3);
    Fp3_sqr(&TMP2_FP3,&TMP1_EFP3.x);
    Fp3_add(&TMP3_FP3,&TMP2_FP3,&TMP2_FP3);
    Fp3_add(&TMP2_FP3,&TMP2_FP3,&TMP3_FP3);
    Fp3_mul(&TMP3_FP3,&TMP1_FP3,&TMP2_FP3);
    
    Fp3_sqr(&TMP1_FP3,&TMP3_FP3);
    Fp3_add(&TMP2_FP3,&TMP1_EFP3.x,&TMP1_EFP3.x);
    Fp3_sub(&ANS->x,&TMP1_FP3,&TMP2_FP3);
    
    Fp3_sub(&TMP1_FP3,&TMP1_EFP3.x,&ANS->x);
    Fp3_mul(&TMP2_FP3,&TMP3_FP3,&TMP1_FP3);
    Fp3_sub(&ANS->y,&TMP2_FP3,&TMP1_EFP3.y);
}

void EFp3_ECA(EFp3 *ANS,EFp3 *P1,EFp3 *P2){
    if(P1->infinity==1){
        EFp3_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp3_set(ANS,P1);
        return;
    }else if(Fp3_cmp(&P1->x,&P2->x)==0){
        if(Fp3_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp3_ECD(ANS,P1);
            return;
        }
    }
    
    EFp3_set(&TMP1_EFP3,P1);
    EFp3_set(&TMP2_EFP3,P2);
    
    Fp3_sub(&TMP1_FP3,&TMP2_EFP3.x,&TMP1_EFP3.x);
    Fp3_inv(&TMP1_FP3,&TMP1_FP3);
    Fp3_sub(&TMP2_FP3,&TMP2_EFP3.y,&TMP1_EFP3.y);
    Fp3_mul(&TMP3_FP3,&TMP1_FP3,&TMP2_FP3);
    Fp3_sqr(&TMP1_FP3,&TMP3_FP3);
    Fp3_sub(&TMP2_FP3,&TMP1_FP3,&TMP1_EFP3.x);
    Fp3_sub(&ANS->x,&TMP2_FP3,&TMP2_EFP3.x);
    Fp3_sub(&TMP1_FP3,&TMP1_EFP3.x,&ANS->x);
    Fp3_mul(&TMP2_FP3,&TMP3_FP3,&TMP1_FP3);
    Fp3_sub(&ANS->y,&TMP2_FP3,&TMP1_EFP3.y);
}

void EFp3_SCM(EFp3 *ANS,EFp3 *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp3_set(ANS,P);
        return;
    }
    
    EFp3 Tmp_P,Next_P;
    EFp3_init(&Tmp_P);
    EFp3_set(&Tmp_P,P);
    EFp3_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    
    EFp3_set(&Next_P,&Tmp_P);
    for(i=1; binary[i]!='\0'; i++){
        EFp3_ECD(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp3_ECA(&Next_P,&Next_P,&Tmp_P);
        }
    }
    EFp3_set(ANS,&Next_P);
    
    EFp3_clear(&Next_P);
    EFp3_clear(&Tmp_P);
}

//skew_frobenius_map
void EFp3_skew_frobenius_map_p1(EFp3 *ANS,EFp3 *A){
    //beta
    Fp3_inv_basis(&ANS->x,&A->x);
    Fp_mul(&ANS->x.x0,&ANS->x.x0,&frobenius_constant[f_p1][6]);
    Fp_mul(&ANS->x.x1,&ANS->x.x1,&frobenius_constant[f_p1][7]);
    Fp_mul(&ANS->x.x2,&ANS->x.x2,&frobenius_constant[f_p1][8]);
    Fp3_set(&TMP1_FP3,&ANS->x);
    Fp_add(&ANS->x.x0,&TMP1_FP3.x2,&TMP1_FP3.x2);
    Fp_set(&ANS->x.x1,&TMP1_FP3.x0);
    Fp_set(&ANS->x.x2,&TMP1_FP3.x1);
    Fp3_mul_basis(&ANS->x,&ANS->x);
    
    //betagamma
    Fp3_inv_basis(&ANS->y,&A->y);
    Fp_mul(&ANS->y.x0,&ANS->y.x0,&frobenius_constant[f_p1][12]);
    Fp_set_neg(&ANS->y.x1,&ANS->y.x1);
    Fp_mul(&ANS->y.x2,&ANS->y.x2,&frobenius_constant[f_p1][14]);
    Fp3_mul_basis(&ANS->y,&ANS->y);
}

void EFp3_skew_frobenius_map_p2(EFp3 *ANS,EFp3 *A){
    Fp3_inv_basis(&ANS->x,&A->x);
    Fp_mul(&ANS->x.x0,&ANS->x.x0,&frobenius_constant[f_p2][6]);
    Fp_mul(&ANS->x.x1,&ANS->x.x1,&frobenius_constant[f_p2][7]);
    Fp_mul(&ANS->x.x2,&ANS->x.x2,&frobenius_constant[f_p2][8]);
    Fp3_set(&TMP1_FP3,&ANS->x);
    Fp_add(&ANS->x.x0,&TMP1_FP3.x1,&TMP1_FP3.x1);
    Fp_add(&ANS->x.x1,&TMP1_FP3.x2,&TMP1_FP3.x2);
    Fp_set(&ANS->x.x2,&TMP1_FP3.x0);
    Fp3_mul_basis(&ANS->x,&ANS->x);
    
    Fp3_inv_basis(&ANS->y,&A->y);
    Fp_mul_mpz(&ANS->y.x0,&ANS->y.x0,frobenius_constant[f_p2][12].x0);
    Fp_mul_mpz(&ANS->y.x1,&ANS->y.x1,frobenius_constant[f_p2][13].x0);
    Fp_mul_mpz(&ANS->y.x2,&ANS->y.x2,frobenius_constant[f_p2][14].x0);
    Fp3_mul_basis(&ANS->y,&ANS->y);
}

void EFp3_skew_frobenius_map_p3(EFp3 *ANS,EFp3 *A){
    Fp3_mul_mpz(&ANS->x,&A->x,frobenius_constant[f_p3][6].x0);
    Fp3_set_neg(&ANS->y,&A->y);
}

/*----------------------------------------------------------------------------*/
//EFp9
void EFp9_init(EFp9 *P){
    Fp9_init(&P->x);
    Fp9_init(&P->y);
    P->infinity=0;
}

void EFp9_set(EFp9 *ANS,EFp9 *A){
    Fp9_set(&ANS->x,&A->x);
    Fp9_set(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp9_set_ui(EFp9 *ANS,unsigned long int UI){
    Fp9_set_ui(&ANS->x,UI);
    Fp9_set_ui(&ANS->y,UI);
    ANS->infinity=0;
}

void EFp9_set_mpz(EFp9 *ANS,mpz_t A){
    Fp9_set_mpz(&ANS->x,A);
    Fp9_set_mpz(&ANS->y,A);
    ANS->infinity=0;
}

void EFp9_set_neg(EFp9 *ANS,EFp9 *A){
    Fp9_set(&ANS->x,&A->x);
    Fp9_set_neg(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp9_clear(EFp9 *P){
    Fp9_clear(&P->x);
    Fp9_clear(&P->y);
}

void EFp9_printf(EFp9 *P,char *str){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp9_printf(&P->x,"");
        printf(",");
        Fp9_printf(&P->y,"");
        printf(")");
    }else{
        printf("0");
    }
}

void EFp9_rational_point(EFp9 *P){
    Fp9 tmp1,tmp2;
    Fp9_init(&tmp1);
    Fp9_init(&tmp2);
    gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    while(1){
        Fp9_set_random(&P->x,state);
        Fp9_sqr(&tmp1,&P->x);
        Fp9_mul(&tmp2,&tmp1,&P->x);
        Fp_add(&tmp2.x0.x0,&tmp2.x0.x0,&curve_b);
        if(Fp9_legendre(&tmp2)==1){
            Fp9_sqrt(&P->y,&tmp2);
            break;
        }
    }
    
    Fp9_clear(&tmp1);
    Fp9_clear(&tmp2);
}

void EFp9_ECD(EFp9 *ANS,EFp9 *P){
    if(Fp9_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFp9_set(&TMP1_EFP9,P);
    
    Fp9_add(&TMP1_FP9,&TMP1_EFP9.y,&TMP1_EFP9.y);
    
    Fp9_inv(&TMP1_FP9,&TMP1_FP9);
    Fp9_sqr(&TMP2_FP9,&TMP1_EFP9.x);
    Fp9_add(&TMP3_FP9,&TMP2_FP9,&TMP2_FP9);
    Fp9_add(&TMP2_FP9,&TMP2_FP9,&TMP3_FP9);
    Fp9_mul(&TMP3_FP9,&TMP1_FP9,&TMP2_FP9);
    
    Fp9_sqr(&TMP1_FP9,&TMP3_FP9);
    Fp9_add(&TMP2_FP9,&TMP1_EFP9.x,&TMP1_EFP9.x);
    Fp9_sub(&ANS->x,&TMP1_FP9,&TMP2_FP9);
    
    Fp9_sub(&TMP1_FP9,&TMP1_EFP9.x,&ANS->x);
    Fp9_mul(&TMP2_FP9,&TMP3_FP9,&TMP1_FP9);
    Fp9_sub(&ANS->y,&TMP2_FP9,&TMP1_EFP9.y);
}

void EFp9_ECA(EFp9 *ANS,EFp9 *P1,EFp9 *P2){
    if(P1->infinity==1){
        EFp9_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp9_set(ANS,P1);
        return;
    }else if(Fp9_cmp(&P1->x,&P2->x)==0){
        if(Fp9_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp9_ECD(ANS,P1);
            return;
        }
    }
    
    EFp9_set(&TMP1_EFP9,P1);
    EFp9_set(&TMP2_EFP9,P2);
    
    Fp9_sub(&TMP1_FP9,&TMP2_EFP9.x,&TMP1_EFP9.x);
    Fp9_inv(&TMP1_FP9,&TMP1_FP9);
    Fp9_sub(&TMP2_FP9,&TMP2_EFP9.y,&TMP1_EFP9.y);
    Fp9_mul(&TMP3_FP9,&TMP1_FP9,&TMP2_FP9);
    Fp9_sqr(&TMP1_FP9,&TMP3_FP9);
    Fp9_sub(&TMP2_FP9,&TMP1_FP9,&TMP1_EFP9.x);
    Fp9_sub(&ANS->x,&TMP2_FP9,&TMP2_EFP9.x);
    Fp9_sub(&TMP1_FP9,&TMP1_EFP9.x,&ANS->x);
    Fp9_mul(&TMP2_FP9,&TMP3_FP9,&TMP1_FP9);
    Fp9_sub(&ANS->y,&TMP2_FP9,&TMP1_EFP9.y);
}

void EFp9_SCM(EFp9 *ANS,EFp9 *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp9_set(ANS,P);
        return;
    }
    
    EFp9 Tmp_P,Next_P;
    EFp9_init(&Tmp_P);
    EFp9_set(&Tmp_P,P);
    EFp9_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    
    EFp9_set(&Next_P,&Tmp_P);
    for(i=1; binary[i]!='\0'; i++){
        EFp9_ECD(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp9_ECA(&Next_P,&Next_P,&Tmp_P);
        }
    }
    EFp9_set(ANS,&Next_P);
    
    EFp9_clear(&Next_P);
    EFp9_clear(&Tmp_P);
}
/*----------------------------------------------------------------------------*/
//EFp18
void EFp18_init(EFp18 *P){
    Fp18_init(&P->x);
    Fp18_init(&P->y);
    P->infinity=0;
}

void EFp18_set(EFp18 *ANS,EFp18 *A){
    Fp18_set(&ANS->x,&A->x);
    Fp18_set(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp18_set_ui(EFp18 *ANS,unsigned long int UI){
    Fp18_set_ui(&ANS->x,UI);
    Fp18_set_ui(&ANS->y,UI);
    ANS->infinity=0;
}

void EFp18_set_mpz(EFp18 *ANS,mpz_t A){
    Fp18_set_mpz(&ANS->x,A);
    Fp18_set_mpz(&ANS->y,A);
    ANS->infinity=0;
}

void EFp18_set_neg(EFp18 *ANS,EFp18 *A){
    Fp18_set(&ANS->x,&A->x);
    Fp18_set_neg(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void EFp18_clear(EFp18 *P){
    Fp18_clear(&P->x);
    Fp18_clear(&P->y);
}

void EFp18_printf(EFp18 *P,char *str){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        Fp18_printf(&P->x,"");
        printf(",");
        Fp18_printf(&P->y,"");
        printf(")");
    }else{
        printf("0");
    }
}

void EFp18_rational_point(EFp18 *P){
    Fp18 tmp1,tmp2;
    Fp18_init(&tmp1);
    Fp18_init(&tmp2);
    gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    while(1){
        Fp18_set_random(&P->x,state);
        Fp18_sqr(&tmp1,&P->x);
        Fp18_mul(&tmp2,&tmp1,&P->x);
        Fp_add(&tmp2.x0.x0.x0,&tmp2.x0.x0.x0,&curve_b);
        if(Fp18_legendre(&tmp2)==1){
            Fp18_sqrt(&P->y,&tmp2);
            break;
        }
    }
    
    Fp18_clear(&tmp1);
    Fp18_clear(&tmp2);
}

void EFp18_generate_G1(EFp18 *P){
    EFp tmp_P;
	EFp_init(&tmp_P);
	mpz_t exp;
	mpz_init(exp);
	
	EFp_rational_point(&tmp_P);
	EFp18_set_ui(P,0);
	mpz_tdiv_q(exp,EFp_total,order);
	EFp_SCM(&tmp_P,&tmp_P,exp);
	Fp_set(&P->x.x0.x0.x0,&tmp_P.x);
	Fp_set(&P->y.x0.x0.x0,&tmp_P.y);
	P->infinity=tmp_P.infinity;
	
	EFp_clear(&tmp_P);
	mpz_clear(exp);
}

void EFp18_generate_G2(EFp18 *Q){
    EFp18 random_P,P,frobenius_P;
    EFp18_init(&random_P);
    EFp18_init(&P);
    EFp18_init(&frobenius_P);
    mpz_t exp;
    mpz_init(exp);
    
    EFp18_rational_point(&random_P);
    mpz_pow_ui(exp,order,2);
    mpz_tdiv_q(exp,EFp18_total,exp);
    EFp18_SCM(&P,&random_P,exp);
    Fp18_frobenius_map_p1(&frobenius_P.x,&P.x);
    Fp18_frobenius_map_p1(&frobenius_P.y,&P.y);
    EFp18_set_neg(&P,&P);
    EFp18_ECA(Q,&P,&frobenius_P);
    
    mpz_clear(exp);
    EFp18_clear(&random_P);
    EFp18_clear(&P);
    EFp18_clear(&frobenius_P);
}

void EFp18_ECD(EFp18 *ANS,EFp18 *P){
    if(Fp18_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    EFp18_set(&TMP1_EFP18,P);
    
    Fp18_add(&TMP1_FP18,&TMP1_EFP18.y,&TMP1_EFP18.y);
    
    Fp18_inv(&TMP1_FP18,&TMP1_FP18);
    Fp18_sqr(&TMP2_FP18,&TMP1_EFP18.x);
    Fp18_add(&TMP3_FP18,&TMP2_FP18,&TMP2_FP18);
    Fp18_add(&TMP2_FP18,&TMP2_FP18,&TMP3_FP18);
    Fp18_mul(&TMP3_FP18,&TMP1_FP18,&TMP2_FP18);
    
    Fp18_sqr(&TMP1_FP18,&TMP3_FP18);
    Fp18_add(&TMP2_FP18,&TMP1_EFP18.x,&TMP1_EFP18.x);
    Fp18_sub(&ANS->x,&TMP1_FP18,&TMP2_FP18);
    
    Fp18_sub(&TMP1_FP18,&TMP1_EFP18.x,&ANS->x);
    Fp18_mul(&TMP2_FP18,&TMP3_FP18,&TMP1_FP18);
    Fp18_sub(&ANS->y,&TMP2_FP18,&TMP1_EFP18.y);
}

void EFp18_ECA(EFp18 *ANS,EFp18 *P1,EFp18 *P2){
    if(P1->infinity==1){
        EFp18_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        EFp18_set(ANS,P1);
        return;
    }else if(Fp18_cmp(&P1->x,&P2->x)==0){
        if(Fp18_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            EFp18_ECD(ANS,P1);
            return;
        }
    }
    
    EFp18_set(&TMP1_EFP18,P1);
    EFp18_set(&TMP2_EFP18,P2);
    
    Fp18_sub(&TMP1_FP18,&TMP2_EFP18.x,&TMP1_EFP18.x);
    Fp18_inv(&TMP1_FP18,&TMP1_FP18);
    Fp18_sub(&TMP2_FP18,&TMP2_EFP18.y,&TMP1_EFP18.y);
    Fp18_mul(&TMP3_FP18,&TMP1_FP18,&TMP2_FP18);
    Fp18_sqr(&TMP1_FP18,&TMP3_FP18);
    Fp18_sub(&TMP2_FP18,&TMP1_FP18,&TMP1_EFP18.x);
    Fp18_sub(&ANS->x,&TMP2_FP18,&TMP2_EFP18.x);
    Fp18_sub(&TMP1_FP18,&TMP1_EFP18.x,&ANS->x);
    Fp18_mul(&TMP2_FP18,&TMP3_FP18,&TMP1_FP18);
    Fp18_sub(&ANS->y,&TMP2_FP18,&TMP1_EFP18.y);
}

void EFp18_SCM(EFp18 *ANS,EFp18 *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        EFp18_set(ANS,P);
        return;
    }
    
    EFp18 Tmp_P,Next_P;
    EFp18_init(&Tmp_P);
    EFp18_set(&Tmp_P,P);
    EFp18_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length];
    mpz_get_str(binary,2,scalar);
    
    EFp18_set(&Next_P,&Tmp_P);
    for(i=1; binary[i]!='\0'; i++){
        EFp18_ECD(&Next_P,&Next_P);
        if(binary[i]=='1'){
            EFp18_ECA(&Next_P,&Next_P,&Tmp_P);
        }
    }
    EFp18_set(ANS,&Next_P);
    
    EFp18_clear(&Next_P);
    EFp18_clear(&Tmp_P);
}

/*============================================================================*/
/* Pairing functions                                                          */
/*============================================================================*/
//twist
void EFp18_to_EFp3(EFp3 *ANS,EFp18 *A){
    Fp3_set_ui(&ANS->x,0);
	Fp3_set(&ANS->x,&A->x.x0.x2);
	Fp3_mul_basis(&ANS->x,&ANS->x);
	Fp3_set_ui(&ANS->y,0);
	Fp3_set(&ANS->y,&A->y.x1.x1);
	Fp3_mul_basis(&ANS->y,&ANS->y);
	ANS->infinity=A->infinity;
}

void EFp3_to_EFp18(EFp18 *ANS,EFp3 *A){
    Fp18_set_ui(&ANS->x,0);
	Fp3_set(&ANS->x.x0.x2,&A->x);
	Fp3_inv_basis(&ANS->x.x0.x2,&ANS->x.x0.x2);
	Fp18_set_ui(&ANS->y,0);
	Fp3_set(&ANS->y.x1.x1,&A->y);
	Fp3_inv_basis(&ANS->y.x1.x1,&ANS->y.x1.x1);
	ANS->infinity=A->infinity;
}

void EFp18_to_EFp(EFp *ANS,EFp18 *A){
    Fp_set(&ANS->x,&A->x.x0.x0.x0);
    Fp_set(&ANS->y,&A->y.x0.x0.x0);
    ANS->infinity=A->infinity;
}

void EFp_to_EFp18(EFp18 *ANS,EFp *A){
    Fp18_set_ui(&ANS->x,0);
    Fp18_set_ui(&ANS->y,0);
    Fp_set(&ANS->x.x0.x0.x0,&A->x);
    Fp_set(&ANS->y.x0.x0.x0,&A->y);
    ANS->infinity=A->infinity;
}

/*----------------------------------------------------------------------------*/
//Pseudo 12-sparse
void Pseudo_12_sparse_mapping(EFp *P,EFp3 *Q,Fp *L){
    EFp_set(&TMP1_EFP,P);
	EFp3_set(&TMP1_EFP3,Q);
	
	Fp_mul(&TMP1_FP,&TMP1_EFP.x,&TMP1_EFP.y);
	Fp_inv(&TMP1_FP,&TMP1_FP);
	Fp_mul(&TMP2_FP,&TMP1_EFP.x,&TMP1_EFP.x);
	Fp_mul(&TMP2_FP,&TMP2_FP,&TMP1_FP);
	Fp_mul(&TMP3_FP,&TMP1_EFP.y,&TMP1_FP);
	Fp_mul(&TMP4_FP,&TMP2_FP,&TMP2_FP);
	
	Fp3_mul_mpz(&Q->x,&TMP1_EFP3.x,TMP4_FP.x0);
	Fp_mul(&TMP5_FP,&TMP2_FP,&TMP4_FP);
	Fp3_mul_mpz(&Q->y,&TMP1_EFP3.y,TMP5_FP.x0);
	
	Fp_mul(&P->x,&TMP4_FP,&TMP1_EFP.x);
	Fp_set(&P->y,&P->x);
	
	Fp_mul(L,&TMP3_FP,&TMP1_EFP.y);
	Fp_mul(L,L,L);
	Fp_mul(L,L,&TMP3_FP);
}

void Pseudo_12_sparse_mul(Fp18 *ANS,Fp18 *A,Fp18 *B){
    //A= f0 + f1γ^2 + f2γ^4 + f3γ　+ f4γ^3 + f5γ^5
	//B= 1                        +  aγ^3 +  bγ^5
	// x0.x0  x0.x1  x0.x2  x1.x0   x1.x1   x1.x2
	Fp3_mul(&TMP1_FP3,&A->x0.x0,&B->x1.x1);
	Fp3_mul(&TMP2_FP3,&A->x0.x1,&B->x1.x2);
	Fp3_add(&TMP3_FP3,&A->x0.x0,&A->x0.x1);
	Fp3_add(&TMP4_FP3,&B->x1.x1,&B->x1.x2);
	Fp3_mul(&TMP3_FP3,&TMP3_FP3,&TMP4_FP3);
	Fp3_sub(&TMP3_FP3,&TMP3_FP3,&TMP1_FP3);
	Fp3_sub(&TMP3_FP3,&TMP3_FP3,&TMP2_FP3);
	
	Fp3_add(&TMP1_FP18.x1.x2,&TMP3_FP3,&A->x1.x2);
	Fp3_add(&TMP1_FP18.x1.x1,&TMP1_FP3,&A->x1.x1);
	Fp3_mul(&TMP3_FP3,&A->x0.x2,&B->x1.x2);
	Fp3_mul_basis(&TMP3_FP3,&TMP3_FP3);
	Fp3_add(&TMP1_FP18.x1.x1,&TMP1_FP18.x1.x1,&TMP3_FP3);
	Fp3_mul(&TMP1_FP3,&A->x0.x2,&B->x1.x1);
	Fp3_add(&TMP1_FP3,&TMP1_FP3,&TMP2_FP3);
	Fp3_mul_basis(&TMP1_FP3,&TMP1_FP3);
	Fp3_add(&TMP1_FP18.x1.x0,&TMP1_FP3,&A->x1.x0);
	
	Fp3_mul(&TMP1_FP3,&A->x1.x0,&B->x1.x1);
	Fp3_mul(&TMP2_FP3,&A->x1.x1,&B->x1.x2);
	Fp3_add(&TMP3_FP3,&A->x1.x0,&A->x1.x1);
	Fp3_mul(&TMP3_FP3,&TMP3_FP3,&TMP4_FP3);
	Fp3_sub(&TMP3_FP3,&TMP3_FP3,&TMP1_FP3);
	Fp3_sub(&TMP3_FP3,&TMP3_FP3,&TMP2_FP3);
	
	Fp3_mul_basis(&TMP3_FP3,&TMP3_FP3);
	Fp3_add(&TMP1_FP18.x0.x0,&TMP3_FP3,&A->x0.x0);
	
	Fp3_mul(&TMP3_FP3,&A->x1.x2,&B->x1.x1);
	Fp3_add(&TMP3_FP3,&TMP2_FP3,&TMP3_FP3);
	Fp3_mul_basis(&TMP3_FP3,&TMP3_FP3);
	Fp3_add(&TMP1_FP18.x0.x1,&TMP3_FP3,&A->x0.x1);
	Fp3_mul(&TMP4_FP3,&A->x1.x2,&B->x1.x2);
	Fp3_mul_basis(&TMP4_FP3,&TMP4_FP3);
	
	Fp3_add(&TMP1_FP3,&TMP1_FP3,&TMP4_FP3);
	Fp3_add(&TMP1_FP18.x0.x2,&TMP1_FP3,&A->x0.x2);
	
	Fp18_set(ANS,&TMP1_FP18);
}

void ff_ltt(Fp18 *f,EFp3 *T,EFp *P,Fp *L){
	EFp3_set(&TMP1_EFP3,T);
	
	Fp3_add(&TMP1_FP3,&TMP1_EFP3.y,&TMP1_EFP3.y);
	Fp3_inv(&TMP1_FP3,&TMP1_FP3);
	Fp3_sqr(&TMP2_FP3,&TMP1_EFP3.x);
	Fp3_add(&TMP3_FP3,&TMP2_FP3,&TMP2_FP3);
	Fp3_add(&TMP2_FP3,&TMP3_FP3,&TMP2_FP3);
	Fp3_mul(&TMP3_FP3,&TMP1_FP3,&TMP2_FP3);
	Fp3_add(&TMP4_FP3,&TMP1_EFP3.x,&TMP1_EFP3.x);
	Fp3_sqr(&T->x,&TMP3_FP3);
	Fp3_sub(&T->x,&T->x,&TMP4_FP3);
	Fp3_mul(&TMP5_FP3,&TMP3_FP3,&TMP1_EFP3.x);
	Fp3_sub(&TMP5_FP3,&TMP5_FP3,&TMP1_EFP3.y);
	Fp3_mul(&T->y,&TMP3_FP3,&T->x);
	Fp3_sub(&T->y,&TMP5_FP3,&T->y);
	
	Fp_set_ui(&TMP2_FP18.x0.x0.x0,1);
	Fp3_set_neg(&TMP2_FP18.x1.x2,&TMP3_FP3);
	Fp3_inv_basis(&TMP2_FP18.x1.x2,&TMP2_FP18.x1.x2);
	Fp3_mul_mpz(&TMP2_FP18.x1.x1,&TMP5_FP3,L->x0);
	Fp3_inv_basis(&TMP2_FP18.x1.x1,&TMP2_FP18.x1.x1);
    
	Fp18_sqr(f,f);
	Pseudo_12_sparse_mul(f,f,&TMP2_FP18);
}

void f_ltq(Fp18 *f,EFp3 *T,EFp3 *Q,EFp *P,Fp *L){
	EFp3_set(&TMP1_EFP3,T);
    
	Fp3_sub(&TMP1_FP3,&Q->x,&TMP1_EFP3.x);
	Fp3_inv(&TMP1_FP3,&TMP1_FP3);
	Fp3_sub(&TMP2_FP3,&Q->y,&TMP1_EFP3.y);
	Fp3_mul(&TMP3_FP3,&TMP1_FP3,&TMP2_FP3);
	Fp3_add(&TMP4_FP3,&TMP1_EFP3.x,&Q->x);
	Fp3_sqr(&T->x,&TMP3_FP3);
	Fp3_sub(&T->x,&T->x,&TMP4_FP3);
	Fp3_mul(&TMP5_FP3,&TMP3_FP3,&TMP1_EFP3.x);
	Fp3_sub(&TMP5_FP3,&TMP5_FP3,&TMP1_EFP3.y);
	Fp3_mul(&T->y,&TMP3_FP3,&T->x);
	Fp3_sub(&T->y,&TMP5_FP3,&T->y);
	
	Fp_set_ui(&TMP2_FP18.x0.x0.x0,1);
	Fp3_set_neg(&TMP2_FP18.x1.x2,&TMP3_FP3);
	Fp3_inv_basis(&TMP2_FP18.x1.x2,&TMP2_FP18.x1.x2);
	Fp3_mul_mpz(&TMP2_FP18.x1.x1,&TMP5_FP3,L->x0);
	Fp3_inv_basis(&TMP2_FP18.x1.x1,&TMP2_FP18.x1.x1);
	
	Pseudo_12_sparse_mul(f,f,&TMP2_FP18);
}

/*----------------------------------------------------------------------------*/
//miller
void Miller_algo_for_plain_ate(Fp18 *ANS,EFp18 *Q,EFp18 *P){
    gettimeofday(&tv_start,NULL);
    
    EFp3 T;
	EFp3_init(&T);
	EFp3 mapped_Q;
	EFp3_init(&mapped_Q);
	EFp mapped_P;
	EFp_init(&mapped_P);
	Fp18 f;
	Fp18_init(&f);
	Fp L;
	Fp_init(&L);
	//loop time
	mpz_t loop;
	mpz_init(loop);
	mpz_sub_ui(loop,trace,1);
	int i,length;
	length=(int)mpz_sizeinbase(loop,2);
	char binary[length];
	mpz_get_str(binary,2,loop);
	
	//set mapped_P & mapped_Q
	EFp18_to_EFp(&mapped_P,P);
	EFp18_to_EFp3(&mapped_Q,Q);
	Pseudo_12_sparse_mapping(&mapped_P,&mapped_Q,&L);
	//set T & f
	EFp3_set(&T,&mapped_Q);
	Fp_set_ui(&f.x0.x0.x0,1);
	
	//miller
	for(i=1; binary[i]!='\0'; i++){
		ff_ltt(&f,&T,&mapped_P,&L);
		if(binary[i]=='1'){
			f_ltq(&f,&T,&mapped_Q,&mapped_P,&L);
		}
	}
	
	Fp18_set(ANS,&f);
	
	Fp18_clear(&f);
	EFp3_clear(&T);
	EFp3_clear(&mapped_Q);
	EFp_clear(&mapped_P);
	Fp_clear(&L);
	mpz_clear(loop);
    
    gettimeofday(&tv_end,NULL);
    MILLER_PLAINATE=timedifference_msec(tv_start,tv_end);
}

void Miller_algo_for_opt_ate(Fp18 *ANS,EFp18 *Q,EFp18 *P){
    gettimeofday(&tv_start,NULL);
    
    EFp18 Buf;
    EFp18_init(&Buf);
    EFp3 T;
    EFp3_init(&T);
    EFp3 mapped_Q,mapped_Q_neg,mapped_Q1,mapped_Q2_neg,P3_mapped_Q,tmp;
    EFp3_init(&mapped_Q);
    EFp3_init(&mapped_Q_neg);
    EFp3_init(&mapped_Q1);
    EFp3_init(&mapped_Q2_neg);
    EFp3_init(&P3_mapped_Q);
    EFp3_init(&tmp);
    EFp mapped_P;
    EFp_init(&mapped_P);
    Fp18 f,tmp_f;
    Fp18_init(&f);
    Fp18_init(&tmp_f);
    Fp L;
    Fp_init(&L);
    int i;
    
    //set
    EFp18_to_EFp(&mapped_P,P);//set mapped_P
    EFp18_to_EFp3(&mapped_Q,Q);//set mapped_Q
    Pseudo_12_sparse_mapping(&mapped_P,&mapped_Q,&L);
    EFp3_set_neg(&mapped_Q_neg,&mapped_Q);//set mapped_Q_neg
    EFp3_set(&T,&mapped_Q);     //set T
    Fp_set_ui(&f.x0.x0.x0,1);
    //miller
    for(i=X_length-1; i>=0; i--){
        switch(X_binary[i]){
            case 0:
                ff_ltt(&f,&T,&mapped_P,&L);
                break;
            case 1:
                ff_ltt(&f,&T,&mapped_P,&L);
                f_ltq(&f,&T,&mapped_Q,&mapped_P,&L);
                break;
            case -1:
                ff_ltt(&f,&T,&mapped_P,&L);
                f_ltq(&f,&T,&mapped_Q_neg,&mapped_P,&L);
                break;
            default:
                break;
        }
        
    }
    
    //l
    Fp18_set_ui(&tmp_f,0);
    Fp_set_ui(&tmp_f.x0.x0.x0,1);
    EFp3_skew_frobenius_map_p1(&tmp,&mapped_Q);
    EFp3_ECD(&P3_mapped_Q,&tmp);
    EFp3_ECA(&P3_mapped_Q,&P3_mapped_Q,&tmp);
    f_ltq(&tmp_f,&T,&P3_mapped_Q,&mapped_P,&L);
    Fp18_mul(&f,&f,&tmp_f);
    
    
    //f3,Q(P)^p
    EFp3_set(&T,&mapped_Q);
    Fp18_set_ui(&tmp_f,0);
    Fp_set_ui(&tmp_f.x0.x0.x0,1);
    ff_ltt(&tmp_f,&T,&mapped_P,&L);
    f_ltq(&tmp_f,&T,&mapped_Q,&mapped_P,&L);
    Fp18_frobenius_map_p1(&tmp_f,&tmp_f);
    
    Fp18_mul(ANS,&f,&tmp_f);
    
    EFp18_clear(&Buf);
    Fp18_clear(&f);
    Fp18_clear(&tmp_f);
    EFp3_clear(&T);
    EFp3_clear(&mapped_Q);
    EFp3_clear(&mapped_Q_neg);
    EFp3_clear(&mapped_Q1);
    EFp3_clear(&mapped_Q2_neg);
    EFp3_clear(&P3_mapped_Q);
    EFp3_clear(&tmp);
    EFp_clear(&mapped_P);
    Fp_clear(&L);
    
    gettimeofday(&tv_end,NULL);
    MILLER_OPTATE=timedifference_msec(tv_start,tv_end);
}

/*----------------------------------------------------------------------------*/
//final exp
void Fp18_pow_X(Fp18 *ANS,Fp18 *A){
    int i;
	Fp18_frobenius_map_p9(&TMP2_FP18,A);
	Fp18_set(&TMP1_FP18,A);
	for(i=X_length-1; i>=0; i--){
		switch(X_binary[i]){
			case 0:
				Fp18_sqr_cyclotomic(&TMP1_FP18,&TMP1_FP18);
				break;
			case 1:
				Fp18_sqr_cyclotomic(&TMP1_FP18,&TMP1_FP18);
				Fp18_mul(&TMP1_FP18,&TMP1_FP18,A);
				break;
			case -1:
				Fp18_sqr_cyclotomic(&TMP1_FP18,&TMP1_FP18);
				Fp18_mul(&TMP1_FP18,&TMP1_FP18,&TMP2_FP18);
				break;
			default:
				break;
		}
	}
	Fp18_set(ANS,&TMP1_FP18);
}

void Final_exp_easy(Fp18 *ANS,Fp18 *f){
    gettimeofday(&tv_start,NULL);
    
    Fp18 t0,t1;
    Fp18_init(&t0);
    Fp18_init(&t1);
    
	//p^9-1
	Fp18_frobenius_map_p9(&t0,f);
	Fp18_inv(&t1,f);
	Fp18_mul(f,&t0,&t1);
	//p^3+1
	Fp18_frobenius_map_p3(&t0,f);
	Fp18_mul(ANS,&t0,f);
	
	Fp18_clear(&t0);
    Fp18_clear(&t1);
	
	gettimeofday(&tv_end,NULL);
    FINALEXP_OPT_EASY=timedifference_msec(tv_start,tv_end);
}

void Final_exp_hard_plain(Fp18 *ANS,Fp18 *f){
    gettimeofday(&tv_start,NULL);
    
    mpz_t exp,buf;
    mpz_init(exp);
    mpz_init(buf);
	
	mpz_pow_ui(exp,prime,6);
	mpz_pow_ui(buf,prime,3);
	mpz_sub(exp,exp,buf);
	mpz_add_ui(exp,exp,1);
	mpz_tdiv_q(exp,exp,order);
	Fp18_pow(ANS,f,exp);
	
	mpz_clear(exp);
	mpz_clear(buf);
	
	gettimeofday(&tv_end,NULL);
    FINALEXP_PLAIN_HARD=timedifference_msec(tv_start,tv_end); 
}

void Final_exp_hard_optimal(Fp18 *ANS,Fp18 *f){
    gettimeofday(&tv_start,NULL);
    
    int i;
    Fp18 A,B,C,t0,t1,t2,t3,t4,t5,t6,t7,t8,tmp_f[7],tmp;
    Fp18_init(&A);
    Fp18_init(&B);
    Fp18_init(&C);
    Fp18_init(&t0);
    Fp18_init(&t1);
    Fp18_init(&t2);
    Fp18_init(&t3);
    Fp18_init(&t4);
    Fp18_init(&t5);
    Fp18_init(&t6);
    Fp18_init(&t7);
    Fp18_init(&t8);
    Fp18_init(&tmp);
    for(i=0; i<7; i++){
        Fp18_init(&tmp_f[i]);
	}
    
	//(p^6-p^3+1)/r
    Fp18_pow_X(&tmp_f[x_1],f);
	Fp18_pow_X(&tmp_f[x_2],&tmp_f[x_1]);
	Fp18_pow_X(&tmp_f[x_3],&tmp_f[x_2]);
	Fp18_pow_X(&tmp_f[x_4],&tmp_f[x_3]);
	Fp18_pow_X(&tmp_f[x_5],&tmp_f[x_4]);
	Fp18_pow_X(&tmp_f[x_6],&tmp_f[x_5]);
	Fp18_pow_X(&tmp_f[x_7],&tmp_f[x_6]);
	
	Fp18_frobenius_map_p9(&A,f);//y23
	Fp18_frobenius_map_p1(&A,&A);
	Fp18_sqr_cyclotomic(&t0,&A);
	
	Fp18_frobenius_map_p9(&B,f);//y21
	Fp18_frobenius_map_p4(&B,&B);
	Fp18_mul(&t0,&t0,&B);
	
	Fp18_frobenius_map_p9(&C,&tmp_f[x_1]);//y22
	Fp18_frobenius_map_p1(&C,&C);
	Fp18_mul(&t1,&t0,&C);
	
	Fp18_frobenius_map_p3(&B,&tmp_f[x_1]);//y19
	Fp18_mul(&t0,&t1,&B);
	
	Fp18_frobenius_map_p9(&tmp,&tmp_f[x_2]);//y14
	Fp18_frobenius_map_p1(&tmp,&tmp);
	Fp18_frobenius_map_p9(&B,&tmp_f[x_3]);
	Fp18_frobenius_map_p2(&B,&B);
	Fp18_mul(&B,&B,&tmp);
	Fp18_mul(&t1,&t1,&B);
	
	Fp18_frobenius_map_p2(&A,f);//y4
	
	Fp18_frobenius_map_p9(&B,&tmp_f[x_1]);//y20
	Fp18_frobenius_map_p4(&B,&B);
	Fp18_mul(&t6,&A,&B);
	Fp18_mul(&t0,&t0,&B);
	
	Fp18_frobenius_map_p9(&tmp,&tmp_f[x_5]);//y2
	Fp18_frobenius_map_p4(&tmp,&tmp);
	Fp18_frobenius_map_p5(&A,f);
	Fp18_mul(&A,&A,&tmp);
	Fp18_mul(&t4,&A,&B);
	
	Fp18_set(&B,&tmp_f[x_1]);//y17
	Fp18_mul(&t2,&t0,&B);
	Fp18_mul(&t0,&t0,&t1);
	
	Fp18_frobenius_map_p3(&A,&tmp_f[x_2]);//y18
	Fp18_mul(&t1,&A,&C);
	
	Fp18_frobenius_map_p9(&B,&tmp_f[x_4]);//y9
	Fp18_frobenius_map_p2(&B,&B);
	Fp18_mul(&t3,&A,&B);
	Fp18_mul(&t2,&t1,&t2);
	
	Fp18_frobenius_map_p9(&B,&tmp_f[x_4]);//y8
	Fp18_frobenius_map_p4(&B,&B);
	Fp18_mul(&t5,&t1,&B);
	
	Fp18_frobenius_map_p9(&B,&tmp_f[x_2]);//y16
	Fp18_frobenius_map_p2(&B,&B);
	Fp18_mul(&t1,&t2,&B);
	
	Fp18_frobenius_map_p3(&B,&tmp_f[x_4]);//y7
	Fp18_mul(&t8,&t2,&B);
	
	Fp18_set(&B,&tmp_f[x_2]);//y15
	Fp18_mul(&t2,&t1,&B);
	
	Fp18_frobenius_map_p9(&B,&tmp_f[x_4]);//y11
	Fp18_frobenius_map_p1(&B,&B);
	Fp18_mul(&t1,&t1,&B);
	Fp18_mul(&t0,&t2,&t0);
	
	Fp18_frobenius_map_p3(&B,&tmp_f[x_5]);//y6
	Fp18_mul(&t7,&t2,&B);
	Fp18_sqr_cyclotomic(&t0,&t0);
	
	Fp18_frobenius_map_p9(&B,&tmp_f[x_2]);//y13
	Fp18_frobenius_map_p4(&B,&B);
	Fp18_mul(&t2,&t0,&B);
	
	Fp18_frobenius_map_p9(&tmp,&tmp_f[x_3]);//y12
	Fp18_frobenius_map_p1(&tmp,&tmp);
	Fp18_frobenius_map_p3(&B,&tmp_f[x_3]);
	Fp18_mul(&B,&B,&tmp);
	Fp18_mul(&t0,&t2,&B);
	Fp18_mul(&t2,&t2,&t8);
	Fp18_mul(&t1,&t0,&t1);
	Fp18_mul(&t0,&t0,&t7);
	Fp18_mul(&t3,&t1,&t3);
    Fp18_mul(&t1,&t1,&t6);
    
    Fp18_frobenius_map_p9(&tmp,&tmp_f[x_3]);//y10
    Fp18_frobenius_map_p4(&tmp,&tmp);
    Fp18_mul(&B,&tmp,&tmp_f[x_3]);
    Fp18_mul(&t6,&t3,&B);
    
    Fp18_frobenius_map_p3(&A,&tmp_f[x_6]);//y1
    Fp18_mul(&t3,&A,&B);
    Fp18_mul(&t2,&t6,&t2);
    
    Fp18_frobenius_map_p5(&B,&tmp_f[x_3]);//y3	
	Fp18_frobenius_map_p9(&tmp,&tmp_f[x_6]);
    Fp18_frobenius_map_p2(&tmp,&tmp);
    Fp18_mul(&B,&B,&tmp);
    Fp18_frobenius_map_p9(&tmp,&tmp_f[x_5]);
    Fp18_frobenius_map_p1(&tmp,&tmp);
    Fp18_mul(&B,&B,&tmp);
    Fp18_mul(&B,&B,&tmp_f[x_5]);
    Fp18_mul(&t6,&t6,&B);
    Fp18_mul(&t2,&t2,&t5);
    
    Fp18_frobenius_map_p5(&B,&tmp_f[x_4]);//y0
    Fp18_frobenius_map_p9(&tmp,&tmp_f[x_7]);
    Fp18_frobenius_map_p2(&tmp,&tmp);
    Fp18_mul(&B,&B,&tmp);
    Fp18_mul(&B,&B,&tmp_f[x_6]);
    Fp18_mul(&t5,&t5,&B);
    Fp18_sqr_cyclotomic(&t2,&t2);
    
    Fp18_frobenius_map_p5(&B,&tmp_f[x_2]);//y5
    Fp18_frobenius_map_p9(&tmp,&tmp_f[x_5]);
    Fp18_frobenius_map_p2(&tmp,&tmp);
    Fp18_mul(&B,&B,&tmp);
    Fp18_mul(&B,&B,&tmp_f[x_4]);
    Fp18_mul(&t2,&t2,&B);
    Fp18_sqr_cyclotomic(&t0,&t0);
    Fp18_mul(&t0,&t0,&t6);
    Fp18_mul(&t1,&t2,&t1);
    Fp18_mul(&t2,&t2,&t5);
    Fp18_sqr_cyclotomic(&t1,&t1);
    Fp18_mul(&t1,&t1,&t4);
    Fp18_mul(&t1,&t1,&t0);
    Fp18_mul(&t0,&t0,&t3);
    Fp18_mul(&t0,&t0,&t1);
    Fp18_mul(&t1,&t1,&t2);
    Fp18_sqr_cyclotomic(&t0,&t0);
    Fp18_mul(ANS,&t0,&t1);
	
	Fp18_clear(&A);
    Fp18_clear(&B);
    Fp18_clear(&C);
    Fp18_clear(&t0);
    Fp18_clear(&t1);
    Fp18_clear(&t2);
    Fp18_clear(&t3);
    Fp18_clear(&t4);
    Fp18_clear(&t5);
    Fp18_clear(&t6);
    Fp18_clear(&t7);
    Fp18_clear(&t8);
    Fp18_clear(&tmp);
    for(i=0; i<7; i++){
        Fp18_clear(&tmp_f[i]);
	}
    
    gettimeofday(&tv_end,NULL);
    FINALEXP_OPT_HARD=timedifference_msec(tv_start,tv_end);
}

/*----------------------------------------------------------------------------*/
//pairing
void Plain_ate_pairing(Fp18 *ANS,EFp18 *P,EFp18 *Q){
    Miller_algo_for_plain_ate(ANS,P,Q);
	Final_exp_easy(ANS,ANS);
	Final_exp_hard_plain(ANS,ANS);
}

void Opt_ate_pairing(Fp18 *ANS,EFp18 *P,EFp18 *Q){
    Miller_algo_for_opt_ate(ANS,P,Q);
	Final_exp_easy(ANS,ANS);
	Final_exp_hard_optimal(ANS,ANS);
}
/*============================================================================*/
/* G1 SCM                                                                     */
/*============================================================================*/
void JSF(int **binary,mpz_t scalar[2],int *loop_length){
	int i,j;
	unsigned long int u;
	mpz_t mod_2,mod_4,mod_8;
	mpz_init(mod_2);
	mpz_init(mod_4);
	mpz_init(mod_8);
	
	mpz_t k[2];
	mpz_init(k[0]);
	mpz_init(k[1]);
	//set
	j=0;
	mpz_set(k[0],scalar[0]);
	mpz_set(k[1],scalar[1]);
	
	while(mpz_cmp_ui(k[0],0)>0 || mpz_cmp_ui(k[1],0)>0){
		for(i=0; i<2; i++){
			mpz_mod_ui(mod_2,k[i],2);
			if(mpz_cmp_ui(mod_2,0)==0){
				u=0;
			}else{
				mpz_mod_ui(mod_4,k[i],4);
				u=mpz_get_ui(mod_4);
				if(u==3){
					u=-1;
				}
				mpz_mod_ui(mod_8,k[i],8);
				mpz_mod_ui(mod_4,k[1-i],4);
				if((mpz_cmp_ui(mod_8,3)==0 || mpz_cmp_ui(mod_8,5)==0) && mpz_cmp_ui(mod_4,2)==0){
					u=-u;
				}
			}
			binary[i][j]=u;
		}
		for(i=0; i<2; i++){
			u=binary[i][j];
			switch (u){
				case 1:
					mpz_sub_ui(k[i],k[i],1);
					break;
				case -1:
					mpz_add_ui(k[i],k[i],1);
					break;
				default:
					break;
			}
			mpz_tdiv_q_ui(k[i],k[i],2);
		}
		j=j+1;
	}
	*loop_length=j-1;
	
	mpz_clear(mod_2);
	mpz_clear(mod_4);
	mpz_clear(mod_8);
	mpz_clear(k[0]);
	mpz_clear(k[1]);	
}

void EFp18_G1_SCM_plain(EFp18 *ANS,EFp18 *P,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    EFp tmp_P;
	EFp_init(&tmp_P);
	
	EFp18_to_EFp(&tmp_P,P);
	EFp_SCM(&tmp_P,&tmp_P,scalar);
	EFp_to_EFp18(ANS,&tmp_P);
	
	EFp_clear(&tmp_P);
    
    gettimeofday(&tv_end,NULL);
    G1SCM_PLAIN=timedifference_msec(tv_start,tv_end);
}

void EFp18_G1_SCM_2split(EFp18 *ANS,EFp18 *P,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    //s=s0+s1[x^4]
	int i,length_s[2],loop_length;
	EFp next_tmp_P,tmp_P,tmp_P_3x,tmp;
	EFp_init(&next_tmp_P);
	EFp_init(&tmp_P);
	EFp_init(&tmp_P_3x);
	EFp_init(&tmp);
	mpz_t s[2],buf;
	mpz_init(buf);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	//table
	EFp table[4];
	for(i=0; i<4; i++){
		EFp_init(&table[i]);
	}
	
	//set
	EFp18_to_EFp(&tmp_P,P);				//tmp_P
	mpz_set_ui(buf,18);
	EFp_skew_frobenius_map_p3(&tmp_P_3x,&tmp_P);//tmp_P_3x
	EFp_SCM(&tmp,&tmp_P,buf);
	EFp_ECA(&tmp_P_3x,&tmp_P_3x,&tmp);
	EFp_set_neg(&tmp_P_3x,&tmp_P_3x);
	
	//set table
	table[0].infinity=1;						//00
	EFp_set(&table[1],&tmp_P);			//01
	EFp_set(&table[2],&tmp_P_3x);			//10
	EFp_ECA(&table[3],&tmp_P,&tmp_P_3x);	//11
	
	//s0,s1
	mpz_set(buf,X);
	mpz_pow_ui(buf,buf,3);
	mpz_tdiv_qr(s[1],s[0],scalar,buf);
	
	//binary
	loop_length=0;
	for(i=0; i<2; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	
	//set binary
	char binary_s[2][loop_length+1];
	char str[5],*e;
	int binary[loop_length+1];
	for(i=0; i<2; i++){
		if(length_s[i]==loop_length){
			mpz_get_str(binary_s[i],2,s[i]);
		}else{
			char binary_buf[loop_length+1];
			mpz_get_str(binary_buf,2,s[i]);
			memset(binary_s[i],'0',sizeof(binary_s[i]));
			memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
		}
	}
	for(i=0; i<loop_length; i++){
		sprintf(str,"%c%c",binary_s[1][i],binary_s[0][i]);
		binary[i]=strtol(str,&e,2);
	}
	EFp_set(&next_tmp_P,&table[binary[0]]);
	
	//SCM
	for(i=1; i<loop_length; i++){
		EFp_ECD(&next_tmp_P,&next_tmp_P);
		EFp_ECA(&next_tmp_P,&next_tmp_P,&table[binary[i]]);
	}
	
	EFp_to_EFp18(ANS,&next_tmp_P);
	
	mpz_clear(buf);
	EFp_clear(&next_tmp_P);
	EFp_clear(&tmp_P);
	EFp_clear(&tmp_P_3x);
	EFp_clear(&tmp);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
	for(i=0; i<4; i++){
		EFp_clear(&table[i]);
	}
    gettimeofday(&tv_end,NULL);
    G1SCM_2SPLIT=timedifference_msec(tv_start,tv_end);
}

void EFp18_G1_SCM_2split_JSF(EFp18 *ANS,EFp18 *P,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    //s=s0+s1[x^4]
	int i,length_s[2],loop_length;
	EFp next_tmp_P,tmp_P,tmp_P_neg,tmp_P_3x,tmp_P_3x_neg,tmp;
	EFp_init(&next_tmp_P);
	EFp_init(&tmp_P);
	EFp_init(&tmp_P_neg);
	EFp_init(&tmp_P_3x);
	EFp_init(&tmp_P_3x_neg);
	EFp_init(&tmp);
	mpz_t s[2],buf;
	mpz_init(buf);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	//table
	EFp table[9];
	for(i=0; i<9; i++){
		EFp_init(&table[i]);
	}
	
	//set
	EFp18_to_EFp(&tmp_P,P);					//tmp_P
	EFp_set_neg(&tmp_P_neg,&tmp_P);			//tmp_P_neg
	mpz_set_ui(buf,18);
	EFp_skew_frobenius_map_p3(&tmp_P_3x,&tmp_P);//tmp_P_3x
	EFp_SCM(&tmp,&tmp_P,buf);
	EFp_ECA(&tmp_P_3x,&tmp_P_3x,&tmp);
	EFp_set(&tmp_P_3x_neg,&tmp_P_3x);
	EFp_set_neg(&tmp_P_3x,&tmp_P_3x);
	
	//set table
	table[0].infinity=1;						//00
	EFp_set(&table[1],&tmp_P);				//01
	EFp_set(&table[2],&tmp_P_3x);				//10
	EFp_ECA(&table[3],&tmp_P_3x,&tmp_P);		//11
	EFp_set(&table[4],&tmp_P_neg);			//0-1
	EFp_set(&table[5],&tmp_P_3x_neg);			//-10
	EFp_ECA(&table[6],&tmp_P_3x_neg,&tmp_P_neg);	//-1-1
	EFp_ECA(&table[7],&tmp_P_3x,&tmp_P_neg);		//1-1
	EFp_ECA(&table[8],&tmp_P_3x_neg,&tmp_P);		//-11
	
	//s0,s1
	mpz_set(buf,X);
	mpz_pow_ui(buf,buf,3);
	mpz_tdiv_qr(s[1],s[0],scalar,buf);
	
	//get loop_length
	loop_length=0;
	for(i=0; i<2; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	//JSF
	int JSF_length;
	int JSF_binary[2][loop_length+1];
	for(i=0; i<loop_length; i++){
		JSF_binary[0][i]=0;
		JSF_binary[1][i]=0;
	}
	int *JSF_pointer[2];
	JSF_pointer[0]=JSF_binary[0];
	JSF_pointer[1]=JSF_binary[1];
	JSF(JSF_pointer,s,&JSF_length);
	int binary[JSF_length+1];
	
	for(i=JSF_length; i>=0; i--){
		if(JSF_binary[1][i]==0 && JSF_binary[0][i]==0) 		binary[i]=0;
		else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==1) 	binary[i]=1;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==0) 	binary[i]=2;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==1)	binary[i]=3;
		else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==-1)	binary[i]=4;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==0)	binary[i]=5;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==-1)	binary[i]=6;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==-1)	binary[i]=7;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==1)	binary[i]=8;
	}
    EFp_set(&next_tmp_P,&table[binary[JSF_length]]);
	//SCM
	for(i=JSF_length-1; i>=0; i--){
		EFp_ECD(&next_tmp_P,&next_tmp_P);
		EFp_ECA(&next_tmp_P,&next_tmp_P,&table[binary[i]]);
	}
	EFp_to_EFp18(ANS,&next_tmp_P);
	
	mpz_clear(buf);
	EFp_clear(&next_tmp_P);
	EFp_clear(&tmp_P);
	EFp_clear(&tmp_P_neg);
	EFp_clear(&tmp_P_3x);
	EFp_clear(&tmp_P_3x_neg);
	EFp_clear(&tmp);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
	for(i=0; i<9; i++){
		EFp_clear(&table[i]);
	}
    
    gettimeofday(&tv_end,NULL);
    G1SCM_2SPLIT_JSF=timedifference_msec(tv_start,tv_end);
}
/*============================================================================*/
/* G2 SCM                                                                     */
/*============================================================================*/
void EFp18_G2_SCM_plain(EFp18 *ANS,EFp18 *Q,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    EFp3 tmp_Q;
	EFp3_init(&tmp_Q);
	
	EFp18_to_EFp3(&tmp_Q,Q);
	EFp3_SCM(&tmp_Q,&tmp_Q,scalar);
	EFp3_to_EFp18(ANS,&tmp_Q);
	
	EFp3_clear(&tmp_Q);
    
    gettimeofday(&tv_end,NULL);
    G2SCM_PLAIN=timedifference_msec(tv_start,tv_end);
}

void EFp18_G2_SCM_2split(EFp18 *ANS,EFp18 *Q,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    int i,length_s[2],loop_length;
	EFp3 next_twisted_Q,twisted_Q,skew_Q,tmp;
	EFp3_init(&next_twisted_Q);
	EFp3_init(&twisted_Q);
	EFp3_init(&skew_Q);
	EFp3_init(&tmp);
	mpz_t s[2],buf;
	mpz_init(buf);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	//table
	EFp3 table[4];
	for(i=0; i<4; i++){
		EFp3_init(&table[i]);
	}
	
	//set
	EFp18_to_EFp3(&twisted_Q,Q);				//twisted_Q
	mpz_set_ui(buf,18);
	EFp3_skew_frobenius_map_p3(&skew_Q,&twisted_Q);//skew_Q
	EFp3_SCM(&tmp,&twisted_Q,buf);
	EFp3_ECA(&skew_Q,&skew_Q,&tmp);
	EFp3_set_neg(&skew_Q,&skew_Q);
	
	//set table
	table[0].infinity=1;						//00
	EFp3_set(&table[1],&twisted_Q);			//01
	EFp3_set(&table[2],&skew_Q);				//10
	EFp3_ECA(&table[3],&twisted_Q,&skew_Q);		//11
	
	//s0,s1
	mpz_set(buf,X);
	mpz_pow_ui(buf,buf,3);
	mpz_tdiv_qr(s[1],s[0],scalar,buf);
	
	//binary
	loop_length=0;
	for(i=0; i<2; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	//set binary
	char binary_s[2][loop_length+1];
	char str[5],*e;
	int binary[loop_length+1];
	for(i=0; i<2; i++){
		if(length_s[i]==loop_length){
			mpz_get_str(binary_s[i],2,s[i]);
		}else{
			char binary_buf[loop_length+1];
			mpz_get_str(binary_buf,2,s[i]);
			memset(binary_s[i],'0',sizeof(binary_s[i]));
			memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
		}
	}
	for(i=0; i<loop_length; i++){
		sprintf(str,"%c%c",binary_s[1][i],binary_s[0][i]);
		binary[i]=(int)strtol(str,&e,2);
	}
	EFp3_set(&next_twisted_Q,&table[binary[0]]);
	
	//SCM
	for(i=1; i<loop_length; i++){
		EFp3_ECD(&next_twisted_Q,&next_twisted_Q);
		EFp3_ECA(&next_twisted_Q,&next_twisted_Q,&table[binary[i]]);
	}
	
	EFp3_to_EFp18(ANS,&next_twisted_Q);
	
	mpz_clear(buf);
	EFp3_clear(&next_twisted_Q);
	EFp3_clear(&twisted_Q);
	EFp3_clear(&skew_Q);
	EFp3_clear(&tmp);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
	for(i=0; i<4; i++){
		EFp3_clear(&table[i]);
	}
    
    gettimeofday(&tv_end,NULL);
    G2SCM_2SPLIT=timedifference_msec(tv_start,tv_end);
}

void EFp18_G2_SCM_2split_JSF(EFp18 *ANS,EFp18 *Q,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    int i,length_s[2],loop_length;
	EFp3 next_tmp_Q,tmp_Q,tmp_Q_neg,skew_Q,skew_Q_neg,tmp;
	EFp3_init(&next_tmp_Q);
	EFp3_init(&tmp_Q);
	EFp3_init(&tmp_Q_neg);
	EFp3_init(&skew_Q);
	EFp3_init(&skew_Q_neg);
	EFp3_init(&tmp);
	mpz_t s[2],buf;
	mpz_init(buf);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	//table
	EFp3 table[9];
	for(i=0; i<9; i++){
		EFp3_init(&table[i]);
	}
	
	//set
	EFp18_to_EFp3(&tmp_Q,Q);					//tmp_Q
	EFp3_set_neg(&tmp_Q_neg,&tmp_Q);			//tmp_Q_neg
	mpz_set_ui(buf,18);
	EFp3_skew_frobenius_map_p3(&skew_Q,&tmp_Q);//skew_Q
	EFp3_SCM(&tmp,&tmp_Q,buf);
	EFp3_ECA(&skew_Q,&skew_Q,&tmp);
	EFp3_set(&skew_Q_neg,&skew_Q);
	EFp3_set_neg(&skew_Q,&skew_Q);
	
	//set table
	table[0].infinity=1;						//00
	EFp3_set(&table[1],&tmp_Q);				//01
	EFp3_set(&table[2],&skew_Q);				//10
	EFp3_ECA(&table[3],&skew_Q,&tmp_Q);		//11
	EFp3_set(&table[4],&tmp_Q_neg);			//0-1
	EFp3_set(&table[5],&skew_Q_neg);			//-10
	EFp3_ECA(&table[6],&skew_Q_neg,&tmp_Q_neg);	//-1-1
	EFp3_ECA(&table[7],&skew_Q,&tmp_Q_neg);		//1-1
	EFp3_ECA(&table[8],&skew_Q_neg,&tmp_Q);		//-11
	
	//s0,s1
	mpz_set(buf,X);
	mpz_pow_ui(buf,buf,3);
	mpz_tdiv_qr(s[1],s[0],scalar,buf);
	
	//get loop_length
	loop_length=0;
	for(i=0; i<2; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	//JSF
	int JSF_length;
	int JSF_binary[2][loop_length+1];
	for(i=0; i<loop_length; i++){
		JSF_binary[0][i]=0;
		JSF_binary[1][i]=0;
	}
	int *JSF_pointer[2];
	JSF_pointer[0]=JSF_binary[0];
	JSF_pointer[1]=JSF_binary[1];
	JSF(JSF_pointer,s,&JSF_length);
	int binary[JSF_length+1];
	for(i=JSF_length; i>=0; i--){
		if(JSF_binary[1][i]==0 && JSF_binary[0][i]==0) 		binary[i]=0;
		else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==1) 	binary[i]=1;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==0) 	binary[i]=2;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==1)	binary[i]=3;
		else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==-1)	binary[i]=4;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==0)	binary[i]=5;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==-1)	binary[i]=6;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==-1)	binary[i]=7;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==1)	binary[i]=8;
	}
	EFp3_set(&next_tmp_Q,&table[binary[JSF_length]]);
	//SCM
	for(i=JSF_length-1; i>=0; i--){
		EFp3_ECD(&next_tmp_Q,&next_tmp_Q);
		EFp3_ECA(&next_tmp_Q,&next_tmp_Q,&table[binary[i]]);
	}
	EFp3_to_EFp18(ANS,&next_tmp_Q);
	
	mpz_clear(buf);
	EFp3_clear(&next_tmp_Q);
	EFp3_clear(&tmp_Q);
	EFp3_clear(&tmp_Q_neg);
	EFp3_clear(&skew_Q);
	EFp3_clear(&skew_Q_neg);
	EFp3_clear(&tmp);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
	for(i=0; i<9; i++){
		EFp3_clear(&table[i]);
	}
    
    gettimeofday(&tv_end,NULL);
    G2SCM_2SPLIT_JSF=timedifference_msec(tv_start,tv_end);
}

void EFp18_G2_SCM_6split(EFp18 *ANS,EFp18 *Q,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    //s=s0+s1[x]+s2[x^2]+s3[x^3]+s4[x^4]+s5[x^5]
	int i,length_s[6],loop_length;
	EFp3 next_twisted_Q,twisted_Q_p3,tmp_Q_p3,tmp_Q,tmp,twisted_Q_x[6];
	EFp3_init(&twisted_Q_p3);
	EFp3_init(&tmp_Q_p3);
	EFp3_init(&tmp_Q);
	EFp3_init(&next_twisted_Q);
	EFp3_init(&tmp);
	for(i=0; i<6; i++){
		EFp3_init(&twisted_Q_x[i]);
	}
	EFp3 table0to2[8],table3to5[8];
	for(i=0; i<8; i++){
		EFp3_init(&table0to2[i]);
		EFp3_init(&table3to5[i]);
	}
	mpz_t A,B,C,D,s[6],x_1,x_2,x_3,exp;
	mpz_init(A);
	mpz_init(B);
	mpz_init(C);
	mpz_init(D);
	mpz_init(x_1);
	mpz_init(x_2);
	mpz_init(x_3);
	mpz_init(exp);
	for(i=0; i<6; i++){
		mpz_init(s[i]);
	}
	mpz_set(x_1,X);
	mpz_mul(x_2,x_1,x_1);
	mpz_mul(x_3,x_2,x_1);
	
	EFp18_to_EFp3(&twisted_Q_x[0],Q);
	EFp3_skew_frobenius_map_p3(&twisted_Q_p3,&twisted_Q_x[0]);
	//twisted_Q_x1=[p(p^3-3)]Q
	EFp3_ECD(&tmp,&twisted_Q_x[0]);
	EFp3_ECA(&tmp_Q,&tmp,&twisted_Q_x[0]);
	EFp3_set_neg(&tmp,&tmp_Q);
	EFp3_ECA(&tmp,&twisted_Q_p3,&tmp);
	EFp3_skew_frobenius_map_p1(&twisted_Q_x[1],&tmp);	
		
	//twisted_Q_x1=[p^2(-5p^3+8)]Q=[p^2(5(-p^3+1)+3]Q
	EFp3_set_neg(&tmp_Q_p3,&twisted_Q_p3);
	EFp3_ECA(&tmp,&tmp_Q_p3,&twisted_Q_x[0]);
	mpz_set_ui(exp,5);
	EFp3_SCM(&tmp,&tmp,exp);
	EFp3_ECA(&tmp,&tmp,&tmp_Q);
	EFp3_skew_frobenius_map_p2(&twisted_Q_x[2],&tmp);
		
	//twisted_Q_x3=[-p^3-18]Q=[-(p^3+6*(3p))]Q
	mpz_set_ui(exp,6);
	EFp3_SCM(&tmp,&tmp_Q,exp);
	EFp3_ECA(&tmp,&tmp,&twisted_Q_p3);
	EFp3_set_neg(&twisted_Q_x[3],&tmp);
	
	//twisted_Q_x4=[p(-16p^3+55)]Q=[p(16(-p^3+3)+2*3+1)]Q
	mpz_set_ui(exp,16);
	EFp3_ECA(&tmp,&tmp_Q_p3,&tmp_Q);
	EFp3_SCM(&tmp,&tmp,exp);
	EFp3_ECD(&tmp_Q,&tmp_Q);
	EFp3_ECA(&twisted_Q_x[4],&tmp,&tmp_Q);
	EFp3_ECA(&twisted_Q_x[4],&twisted_Q_x[4],&twisted_Q_x[0]);
	EFp3_skew_frobenius_map_p1(&twisted_Q_x[4],&twisted_Q_x[4]);
	
	//twisted_Q_x5=[p^2(87p^3-149)]Q=[p^2(-5*16(-p^3+3)+7(p^3+2*6+1)]Q
	mpz_set_ui(exp,5);
	EFp3_SCM(&tmp,&tmp,exp);
	EFp3_set_neg(&tmp,&tmp);
	mpz_set_ui(exp,7);
	EFp3_ECA(&tmp_Q_p3,&twisted_Q_p3,&twisted_Q_x[0]);
	EFp3_ECD(&tmp_Q,&tmp_Q);
	EFp3_ECA(&tmp_Q_p3,&tmp_Q_p3,&tmp_Q);
	EFp3_SCM(&tmp_Q_p3,&tmp_Q_p3,exp);
	EFp3_ECA(&tmp,&tmp,&tmp_Q_p3);
	EFp3_skew_frobenius_map_p2(&twisted_Q_x[5],&tmp);
	
	
	//set table
	table0to2[0].infinity=1;							        //000
	EFp3_set(&table0to2[1],&twisted_Q_x[0]);				//001
	EFp3_set(&table0to2[2],&twisted_Q_x[1]);				//010
	EFp3_ECA(&table0to2[3],&twisted_Q_x[0],&twisted_Q_x[1]);//011
	EFp3_set(&table0to2[4],&twisted_Q_x[2]);				//100
	EFp3_ECA(&table0to2[5],&twisted_Q_x[2],&twisted_Q_x[0]);//101
	EFp3_ECA(&table0to2[6],&twisted_Q_x[2],&twisted_Q_x[1]);//110
	EFp3_ECA(&table0to2[7],&table0to2[6],&twisted_Q_x[0]);	//111
	
	table3to5[0].infinity=1;							        //000
	EFp3_set(&table3to5[1],&twisted_Q_x[3]);				//001
	EFp3_set(&table3to5[2],&twisted_Q_x[4]);				//010
	EFp3_ECA(&table3to5[3],&twisted_Q_x[3],&twisted_Q_x[4]);//011
	EFp3_set(&table3to5[4],&twisted_Q_x[5]);				//100
	EFp3_ECA(&table3to5[5],&twisted_Q_x[5],&twisted_Q_x[3]);//101
	EFp3_ECA(&table3to5[6],&twisted_Q_x[5],&twisted_Q_x[4]);//110
	EFp3_ECA(&table3to5[7],&table3to5[6],&twisted_Q_x[3]);	//111
	
	//s0,s1,s2,s3,s4,s5
	mpz_tdiv_qr(B,A,scalar,x_3);
	mpz_tdiv_qr(s[2],C,A,x_2);
	mpz_tdiv_qr(s[1],s[0],C,x_1);
	mpz_tdiv_qr(s[5],C,B,x_2);
	mpz_tdiv_qr(s[4],s[3],C,x_1);
	
	//binary
	loop_length=0;
	for(i=0; i<6; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	
	//set binary
	char binary_s[6][loop_length+1];
	char str[4],*e;
	int binary0to2[loop_length],binary3to5[loop_length];
	for(i=0; i<6; i++){
		if(length_s[i]==loop_length){
			mpz_get_str(binary_s[i],2,s[i]);
		}else{
			char binary_buf[loop_length+1];
			mpz_get_str(binary_buf,2,s[i]);
			memset(binary_s[i],'0',sizeof(binary_s[i]));
			memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
		}
	}
	
	for(i=0; i<loop_length; i++){
		sprintf(str,"%c%c%c",binary_s[2][i],binary_s[1][i],binary_s[0][i]);
		binary0to2[i]=strtol(str,&e,2);
		sprintf(str,"%c%c%c",binary_s[5][i],binary_s[4][i],binary_s[3][i]);
		binary3to5[i]=strtol(str,&e,2);
	}
	
	EFp3_ECA(&next_twisted_Q,&table0to2[binary0to2[0]],&table3to5[binary3to5[0]]);
	
	//miller
	for(i=1; i<loop_length; i++){
		EFp3_ECD(&next_twisted_Q,&next_twisted_Q);
		EFp3_ECA(&next_twisted_Q,&next_twisted_Q,&table0to2[binary0to2[i]]);
		EFp3_ECA(&next_twisted_Q,&next_twisted_Q,&table3to5[binary3to5[i]]);
	}
	EFp3_to_EFp18(ANS,&next_twisted_Q);
	ANS->infinity=next_twisted_Q.infinity;
	
	mpz_clear(A);
	mpz_clear(B);
	mpz_clear(C);
	mpz_clear(D);
	mpz_clear(x_1);
	mpz_clear(x_2);
	mpz_clear(x_3);
	mpz_clear(exp);
	EFp3_clear(&twisted_Q_p3);
	EFp3_clear(&tmp_Q_p3);
	EFp3_clear(&tmp_Q);
	EFp3_clear(&next_twisted_Q);
	EFp3_clear(&tmp);
	for(i=0; i<6; i++){
		EFp3_clear(&twisted_Q_x[i]);
	}
	for(i=0; i<6; i++){
		mpz_clear(s[i]);
	}
	for(i=0; i<8; i++){
		EFp3_clear(&table0to2[i]);
		EFp3_clear(&table3to5[i]);
	}
    
    gettimeofday(&tv_end,NULL);
    G2SCM_6SPLIT=timedifference_msec(tv_start,tv_end);
}

/*============================================================================*/
/* G3 EXP                                                                     */
/*============================================================================*/
void Fp18_G3_EXP_plain(Fp18 *ANS,Fp18 *A,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    int i,length;
	length=(int)mpz_sizeinbase(scalar,2);
	char binary[length];
	mpz_get_str(binary,2,scalar);
	Fp18 buf;
	Fp18_init(&buf);
	Fp18_set(&buf,A);
	
	for(i=1; binary[i]!='\0'; i++){
		Fp18_sqr_cyclotomic(&buf,&buf);
		if(binary[i]=='1'){
			Fp18_mul(&buf,A,&buf);
		}
	}
	
	Fp18_set(ANS,&buf);
	Fp18_clear(&buf);
    
    gettimeofday(&tv_end,NULL);
    G3EXP_PLAIN=timedifference_msec(tv_start,tv_end);
}

void Fp18_G3_EXP_2split(Fp18 *ANS,Fp18 *A,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    int i,length_s[2],loop_length;
	Fp18 Buf,next_f,f,frobenius_f,tmp;
	Fp18_init(&Buf);
	Fp18_init(&next_f);
	Fp18_init(&f);
	Fp18_init(&frobenius_f);
	Fp18_init(&tmp);
	mpz_t s[2],buf;
	mpz_init(buf);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	//table
	Fp18 table[4];
	for(i=0; i<4; i++){
		Fp18_init(&table[i]);
	}
	
	//set
	Fp18_set(&f,A);						//f
	mpz_set_ui(buf,18);
	Fp18_frobenius_map_p3(&frobenius_f,&f);//skew_Q
	Fp18_pow(&tmp,&f,buf);
	Fp18_mul(&frobenius_f,&frobenius_f,&tmp);
	Fp18_frobenius_map_p9(&frobenius_f,&frobenius_f);
	
	//set table
	Fp_set_ui(&table[0].x0.x0.x0,1);			//00
	Fp18_set(&table[1],&f);					//01
	Fp18_set(&table[2],&frobenius_f);			//10
	Fp18_mul(&table[3],&f,&frobenius_f);		//11
	
	//s0,s1
	mpz_set(buf,X);
	mpz_pow_ui(buf,buf,3);
	mpz_tdiv_qr(s[1],s[0],scalar,buf);
	
	//binary
	loop_length=0;
	for(i=0; i<2; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	//set binary
	char binary_s[2][loop_length+1];
	char str[5],*e;
	int binary[loop_length+1];
	for(i=0; i<2; i++){
		if(length_s[i]==loop_length){
			mpz_get_str(binary_s[i],2,s[i]);
		}else{
			char binary_buf[loop_length+1];
			mpz_get_str(binary_buf,2,s[i]);
			memset(binary_s[i],'0',sizeof(binary_s[i]));
			memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
		}
	}
	for(i=0; i<loop_length; i++){
		sprintf(str,"%c%c",binary_s[1][i],binary_s[0][i]);
		binary[i]=(int)strtol(str,&e,2);
	}
	Fp18_set(&next_f,&table[binary[0]]);
	
	//EXP
	for(i=1; i<loop_length; i++){
		Fp18_sqr_cyclotomic(&next_f,&next_f);
		Fp18_mul(&next_f,&next_f,&table[binary[i]]);
	}
	
	Fp18_set(ANS,&next_f);
	
	mpz_clear(buf);
	Fp18_clear(&Buf);
	Fp18_clear(&next_f);
	Fp18_clear(&f);
	Fp18_clear(&frobenius_f);
	Fp18_clear(&tmp);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
	for(i=0; i<4; i++){
		Fp18_clear(&table[i]);
	}
    
    gettimeofday(&tv_end,NULL);
    G3EXP_2SPLIT=timedifference_msec(tv_start,tv_end);
}

void Fp18_G3_EXP_2split_JSF(Fp18 *ANS,Fp18 *A,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
    int i,length_s[2],loop_length;
	Fp18 next_f,f,f_inv,frobenius_f,frobenius_f_inv,tmp;
	Fp18_init(&next_f);
	Fp18_init(&f);
	Fp18_init(&f_inv);
	Fp18_init(&frobenius_f);
	Fp18_init(&frobenius_f_inv);
	Fp18_init(&tmp);
	mpz_t s[2],buf;
	mpz_init(buf);
	for(i=0; i<2; i++){
		mpz_init(s[i]);
	}
	//table
	Fp18 table[9];
	for(i=0; i<9; i++){
		Fp18_init(&table[i]);
	}
	
	//set
	Fp18_set(&f,A);							//f
	Fp18_frobenius_map_p9(&f_inv,&f);						//f_inv
	mpz_set_ui(buf,18);
	Fp18_frobenius_map_p3(&frobenius_f,&f);//skew_Q
	Fp18_pow(&tmp,&f,buf);
	Fp18_mul(&frobenius_f,&frobenius_f,&tmp);
	Fp18_set(&frobenius_f_inv,&frobenius_f);
	Fp18_frobenius_map_p9(&frobenius_f,&frobenius_f);						//f_inv
	
	//set table
	Fp_set_ui(&table[0].x0.x0.x0,1);			//00
	Fp18_set(&table[1],&f);					//01
	Fp18_set(&table[2],&frobenius_f);			//10
	Fp18_mul(&table[3],&frobenius_f,&f);		//11
	Fp18_set(&table[4],&f_inv);				//0-1
	Fp18_set(&table[5],&frobenius_f_inv);		//-10
	Fp18_mul(&table[6],&frobenius_f_inv,&f_inv);	//-1-1
	Fp18_mul(&table[7],&frobenius_f,&f_inv);	//1-1
	Fp18_mul(&table[8],&frobenius_f_inv,&f);	//-11
	
	//s0,s1
	mpz_set(buf,X);
	mpz_pow_ui(buf,buf,3);
	mpz_tdiv_qr(s[1],s[0],scalar,buf);
	
	//get loop_length
	loop_length=0;
	for(i=0; i<2; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	//JSF
	int JSF_length;
	int JSF_binary[2][loop_length+1];
	for(i=0; i<loop_length; i++){
		JSF_binary[0][i]=0;
		JSF_binary[1][i]=0;
	}
	int *JSF_pointer[2];
	JSF_pointer[0]=JSF_binary[0];
	JSF_pointer[1]=JSF_binary[1];
	JSF(JSF_pointer,s,&JSF_length);
	int binary[JSF_length+1];
	for(i=JSF_length; i>=0; i--){
		if(JSF_binary[1][i]==0 && JSF_binary[0][i]==0) 		binary[i]=0;
		else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==1) 	binary[i]=1;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==0) 	binary[i]=2;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==1)	binary[i]=3;
		else if(JSF_binary[1][i]==0 && JSF_binary[0][i]==-1)	binary[i]=4;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==0)	binary[i]=5;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==-1)	binary[i]=6;
		else if(JSF_binary[1][i]==1 && JSF_binary[0][i]==-1)	binary[i]=7;
		else if(JSF_binary[1][i]==-1 && JSF_binary[0][i]==1)	binary[i]=8;
	}
	Fp18_set(&next_f,&table[binary[JSF_length]]);
	//SCM
	for(i=JSF_length-1; i>=0; i--){
		Fp18_sqr_cyclotomic(&next_f,&next_f);
		Fp18_mul(&next_f,&next_f,&table[binary[i]]);
	}
	Fp18_set(ANS,&next_f);
	
	mpz_clear(buf);
	Fp18_clear(&next_f);
	Fp18_clear(&f);
	Fp18_clear(&f_inv);
	Fp18_clear(&frobenius_f);
	Fp18_clear(&frobenius_f_inv);
	Fp18_clear(&tmp);
	for(i=0; i<2; i++){
		mpz_clear(s[i]);
	}
	for(i=0; i<9; i++){
		Fp18_clear(&table[i]);
	}
    
    gettimeofday(&tv_end,NULL);
    G3EXP_2SPLIT_JSF=timedifference_msec(tv_start,tv_end);
}

void Fp18_G3_EXP_6split(Fp18 *ANS,Fp18 *A,mpz_t scalar){
    gettimeofday(&tv_start,NULL);
    
	int i,length_s[6],loop_length;
	Fp18 next_f,f_p3,tmp_f_p3,tmp_f,tmp,f_x[6];
	Fp18_init(&f_p3);
	Fp18_init(&tmp_f_p3);
	Fp18_init(&tmp_f);
	Fp18_init(&next_f);
	Fp18_init(&tmp);
	for(i=0; i<6; i++){
		Fp18_init(&f_x[i]);
	}
	Fp18 table0to2[8],table3to5[8];
	for(i=0; i<8; i++){
		Fp18_init(&table0to2[i]);
		Fp18_init(&table3to5[i]);
	}
	mpz_t a,b,c,d,s[6],x_1,x_2,x_3,exp;
	mpz_init(a);
	mpz_init(b);
	mpz_init(c);
	mpz_init(d);
	mpz_init(x_1);
	mpz_init(x_2);
	mpz_init(x_3);
	mpz_init(exp);
	for(i=0; i<6; i++){
		mpz_init(s[i]);
	}
	mpz_set(x_1,X);
	mpz_mul(x_2,x_1,x_1);
	mpz_mul(x_3,x_2,x_1);
	
	Fp18_set(&f_x[0],A);
	Fp18_frobenius_map_p3(&f_p3,&f_x[0]);
	//f_x1=[p(p^3-3)]Q
	Fp18_sqr_cyclotomic(&tmp,&f_x[0]);
	Fp18_mul(&tmp_f,&tmp,&f_x[0]);
	Fp18_frobenius_map_p9(&tmp,&tmp_f);
	Fp18_mul(&tmp,&f_p3,&tmp);
	Fp18_frobenius_map_p1(&f_x[1],&tmp);	
		
	//f_x1=[p^2(-5p^3+8)]Q=[p^2(5(-p^3+1)+3]Q
	Fp18_frobenius_map_p9(&tmp_f_p3,&f_p3);
	Fp18_mul(&tmp,&tmp_f_p3,&f_x[0]);
	mpz_set_ui(exp,5);
	Fp18_pow(&tmp,&tmp,exp);
	Fp18_mul(&tmp,&tmp,&tmp_f);
	Fp18_frobenius_map_p2(&f_x[2],&tmp);
		
	//f_x3=[-p^3-18]Q=[-(p^3+6*(3p))]Q
	mpz_set_ui(exp,6);
	Fp18_pow(&tmp,&tmp_f,exp);
	Fp18_mul(&tmp,&tmp,&f_p3);
	Fp18_frobenius_map_p9(&f_x[3],&tmp);
	
	//f_x4=[p(-16p^3+55)]Q=[p(16(-p^3+3)+2*3+1)]Q
	mpz_set_ui(exp,16);
	Fp18_mul(&tmp,&tmp_f_p3,&tmp_f);
	Fp18_pow(&tmp,&tmp,exp);
	Fp18_sqr_cyclotomic(&tmp_f,&tmp_f);
	Fp18_mul(&f_x[4],&tmp,&tmp_f);
	Fp18_mul(&f_x[4],&f_x[4],&f_x[0]);
	Fp18_frobenius_map_p1(&f_x[4],&f_x[4]);
	
	//f_x5=[p^2(87p^3-149)]Q=[p^2(-5*16(-p^3+3)+7(p^3+2*6+1)]Q
	mpz_set_ui(exp,5);
	Fp18_pow(&tmp,&tmp,exp);
	Fp18_frobenius_map_p9(&tmp,&tmp);
	mpz_set_ui(exp,7);
	Fp18_mul(&tmp_f_p3,&f_p3,&f_x[0]);
	Fp18_sqr_cyclotomic(&tmp_f,&tmp_f);
	Fp18_mul(&tmp_f_p3,&tmp_f_p3,&tmp_f);
	Fp18_pow(&tmp_f_p3,&tmp_f_p3,exp);
	Fp18_mul(&tmp,&tmp,&tmp_f_p3);
	Fp18_frobenius_map_p2(&f_x[5],&tmp);
	
	
	//set table
	Fp_set_ui(&table0to2[0].x0.x0.x0,1);                    //000
	Fp18_set(&table0to2[1],&f_x[0]);				//001
	Fp18_set(&table0to2[2],&f_x[1]);				//010
	Fp18_mul(&table0to2[3],&f_x[0],&f_x[1]);//011
	Fp18_set(&table0to2[4],&f_x[2]);				//100
	Fp18_mul(&table0to2[5],&f_x[2],&f_x[0]);//101
	Fp18_mul(&table0to2[6],&f_x[2],&f_x[1]);//110
	Fp18_mul(&table0to2[7],&table0to2[6],&f_x[0]);	//111
	
	Fp_set_ui(&table3to5[0].x0.x0.x0,1);                    //000
	Fp18_set(&table3to5[1],&f_x[3]);				//001
	Fp18_set(&table3to5[2],&f_x[4]);				//010
	Fp18_mul(&table3to5[3],&f_x[3],&f_x[4]);//011
	Fp18_set(&table3to5[4],&f_x[5]);				//100
	Fp18_mul(&table3to5[5],&f_x[5],&f_x[3]);//101
	Fp18_mul(&table3to5[6],&f_x[5],&f_x[4]);//110
	Fp18_mul(&table3to5[7],&table3to5[6],&f_x[3]);	//111
	
	//s0,s1,s2,s3,s4,s5
	mpz_tdiv_qr(b,a,scalar,x_3);
	mpz_tdiv_qr(s[2],c,a,x_2);
	mpz_tdiv_qr(s[1],s[0],c,x_1);
	mpz_tdiv_qr(s[5],c,b,x_2);
	mpz_tdiv_qr(s[4],s[3],c,x_1);
	
	//binary
	loop_length=0;
	for(i=0; i<6; i++){
		length_s[i]=(int)mpz_sizeinbase(s[i],2);
		if(loop_length<length_s[i]){
			loop_length=length_s[i];
		}
	}
	
	//set binary
	char binary_s[6][loop_length+1];
	char str[4],*e;
	int binary0to2[loop_length],binary3to5[loop_length];
	for(i=0; i<6; i++){
		if(length_s[i]==loop_length){
			mpz_get_str(binary_s[i],2,s[i]);
		}else{
			char binary_buf[loop_length+1];
			mpz_get_str(binary_buf,2,s[i]);
			memset(binary_s[i],'0',sizeof(binary_s[i]));
			memmove(binary_s[i]+loop_length-length_s[i],binary_buf,sizeof(binary_buf));
		}
	}
	
	for(i=0; i<loop_length; i++){
		sprintf(str,"%c%c%c",binary_s[2][i],binary_s[1][i],binary_s[0][i]);
		binary0to2[i]=strtol(str,&e,2);
		sprintf(str,"%c%c%c",binary_s[5][i],binary_s[4][i],binary_s[3][i]);
		binary3to5[i]=strtol(str,&e,2);
	}
	
	Fp18_mul(&next_f,&table0to2[binary0to2[0]],&table3to5[binary3to5[0]]);
	
	//miller
	for(i=1; i<loop_length; i++){
		Fp18_sqr_cyclotomic(&next_f,&next_f);
		Fp18_mul(&next_f,&next_f,&table0to2[binary0to2[i]]);
		Fp18_mul(&next_f,&next_f,&table3to5[binary3to5[i]]);
	}
	Fp18_set(ANS,&next_f);
	
	mpz_clear(a);
	mpz_clear(b);
	mpz_clear(c);
	mpz_clear(d);
	mpz_clear(x_1);
	mpz_clear(x_2);
	mpz_clear(x_3);
	mpz_clear(exp);
	Fp18_clear(&f_p3);
	Fp18_clear(&tmp_f_p3);
	Fp18_clear(&tmp_f);
	Fp18_clear(&next_f);
	Fp18_clear(&tmp);
	for(i=0; i<6; i++){
		Fp18_clear(&f_x[i]);
	}
	for(i=0; i<6; i++){
		mpz_clear(s[i]);
	}
	for(i=0; i<8; i++){
		Fp18_clear(&table0to2[i]);
		Fp18_clear(&table3to5[i]);
	}
    
    gettimeofday(&tv_end,NULL);
    G3EXP_6SPLIT=timedifference_msec(tv_start,tv_end);
}

/*============================================================================*/
/* KSS18 init                                                                 */
/*============================================================================*/
//init/set/clear
void KSS18_init(){
    init_parameters();
    generate_X();
    
    if(generate_prime()==1 && generate_order()==1){
        generate_trace();
        weil();
        get_epsilon();
        set_basis();
        set_frobenius_constant();
        set_curve_coefficient();
    }else{
        KSS18_clear();
        printf("error : prime\nexit\n");
        exit(1);
    }
}

void init_parameters(){
    int i,j;
    
    mpz_init(X);
    mpz_init(prime);
    mpz_init(order);
    mpz_init(trace);
    
    mpz_init(EFp_total);
    mpz_init(EFp3_total);
    mpz_init(EFp9_total);
    mpz_init(EFp18_total);
    Fp_init(&curve_b);
    
    for(i=0; i<X_length+1; i++){
        X_binary[i]=0;
    }
    
    mpz_init(epsilon1);
    mpz_init(epsilon2);
    
    for(i=0; i<18; i++){
        for(j=0; j<18; j++){
            Fp_init(&frobenius_constant[i][j]);
        }
        for(j=0; j<2; j++){
            Fp3_init(&skew_frobenius_constant[i][j]);
        }
    }
    
    Fp3_init(&Alpha);
    Fp3_init(&Alpha_inv);
    
    Fp_init(&TMP1_FP);
    Fp_init(&TMP2_FP);
    Fp_init(&TMP3_FP);
    Fp_init(&TMP4_FP);
    Fp_init(&TMP5_FP);
    Fp_init(&TMP6_FP);
    Fp_init(&TMP7_FP);
    
    Fp3_init(&TMP1_FP3);
    Fp3_init(&TMP2_FP3);
    Fp3_init(&TMP3_FP3);
    Fp3_init(&TMP4_FP3);
    Fp3_init(&TMP5_FP3);
    Fp3_init(&TMP6_FP3);
    Fp3_init(&TMP7_FP3);
    
    Fp9_init(&TMP1_FP9);
    Fp9_init(&TMP2_FP9);
    Fp9_init(&TMP3_FP9);
    Fp9_init(&TMP4_FP9);
    
    Fp18_init(&TMP1_FP18);
    Fp18_init(&TMP2_FP18);
    Fp18_init(&TMP3_FP18);
    
    EFp_init(&TMP1_EFP);
    EFp_init(&TMP2_EFP);
    EFp3_init(&TMP1_EFP3);
    EFp3_init(&TMP2_EFP3);
    EFp9_init(&TMP1_EFP9);
    EFp9_init(&TMP2_EFP9);
    EFp18_init(&TMP1_EFP18);
    EFp18_init(&TMP2_EFP18);
}

void KSS18_print_parameters(){
    printf("====================================================================================\n");
    printf("KSS18\n\n");
    gmp_printf("parameters\n");
    gmp_printf("X     (%dbit length) : %Zd \n",(int)mpz_sizeinbase(X,2),X);
    gmp_printf("prime (%dbit length) : %Zd \n",(int)mpz_sizeinbase(prime,2),prime);
    gmp_printf("order (%dbit length) : %Zd \n",(int)mpz_sizeinbase(order,2),order);
    gmp_printf("trace (%dbit length) : %Zd \n",(int)mpz_sizeinbase(trace,2),trace);
    
    gmp_printf("\nKSS-18 curve\n");
    gmp_printf("E:y^2=x^3+2^5\n");
    gmp_printf("twisted curve\n");
    gmp_printf("E':y^2=x^3+2^5*alpha\n");
    
    gmp_printf("\nmodulo polynomial\n");
    gmp_printf("Fp3 = Fp[alpha]/(alpha^3-2)\n");
    gmp_printf("Fp9 = Fp3[beta]/(beta^3-alpha)\n");
    gmp_printf("Fp18= Fp9[gamma]/(gamma^2-beta)\n");
}

void KSS18_clear(){
    int i,j;
    
    mpz_clear(X);
    mpz_clear(prime);
    mpz_clear(order);
    mpz_clear(trace);
    
    mpz_clear(EFp_total);
    mpz_clear(EFp3_total);
    mpz_clear(EFp9_total);
    mpz_clear(EFp18_total);
    Fp_clear(&curve_b);
    
    mpz_clear(epsilon1);
    mpz_clear(epsilon2);
    
    for(i=0; i<18; i++){
        for(j=0; j<18; j++){
            Fp_clear(&frobenius_constant[i][j]);
        }
        for(j=0; j<2; j++){
            Fp3_init(&skew_frobenius_constant[i][j]);
        }
    }
    
    Fp3_clear(&Alpha);
    Fp3_clear(&Alpha_inv);
    
    Fp_clear(&TMP1_FP);
    Fp_clear(&TMP2_FP);
    Fp_clear(&TMP3_FP);
    Fp_clear(&TMP4_FP);
    Fp_clear(&TMP5_FP);
    Fp_clear(&TMP6_FP);
    Fp_clear(&TMP7_FP);
    
    Fp3_clear(&TMP1_FP3);
    Fp3_clear(&TMP2_FP3);
    Fp3_clear(&TMP3_FP3);
    Fp3_clear(&TMP4_FP3);
    Fp3_clear(&TMP5_FP3);
    Fp3_clear(&TMP6_FP3);
    Fp3_clear(&TMP7_FP3);
    
    Fp9_clear(&TMP1_FP9);
    Fp9_clear(&TMP2_FP9);
    Fp9_clear(&TMP3_FP9);
    Fp9_clear(&TMP4_FP9);
    
    Fp18_clear(&TMP1_FP18);
    Fp18_clear(&TMP2_FP18);
    Fp18_clear(&TMP3_FP18);
    
    EFp_clear(&TMP1_EFP);
    EFp_clear(&TMP2_EFP);
    EFp3_clear(&TMP1_EFP3);
    EFp3_clear(&TMP2_EFP3);
    EFp9_clear(&TMP1_EFP9);
    EFp9_clear(&TMP2_EFP9);
    EFp18_clear(&TMP1_EFP18);
    EFp18_clear(&TMP2_EFP18);
}

void generate_X(){
    int i;
    mpz_t buf;
    mpz_init(buf);
    
    //X_binary
    X_binary[85]=1;
    X_binary[74]=-1;
    X_binary[71]=-1;
    X_binary[45]=1;
    X_binary[1]=-1;
    
    //KSS18.X
    mpz_set_ui(X,0);
    for(i=X_length; i>=0; i--){
        if(X_binary[i]==1){
            mpz_ui_pow_ui(buf,2,i);
            mpz_add(X,X,buf);
        }else if(X_binary[i]==-1){
            mpz_ui_pow_ui(buf,2,i);
            mpz_sub(X,X,buf);
        }
    }
    
    mpz_clear(buf);
}

int generate_prime(){
    //p=1/21(x^8+5x^7+7x^6+37x^5+188x^4+259x^3+343x^2+1763x+2401)
    mpz_t buf,result;
    mpz_init(buf);
    mpz_init(result);
    
    mpz_pow_ui(buf,X,8);
    mpz_set(result,buf);
    mpz_pow_ui(buf,X,7);
    mpz_mul_ui(buf,buf,5);
    mpz_add(result,result,buf);
    mpz_pow_ui(buf,X,6);
    mpz_mul_ui(buf,buf,7);
    mpz_add(result,result,buf);
    mpz_pow_ui(buf,X,5);
    mpz_mul_ui(buf,buf,37);
    mpz_add(result,result,buf);
    mpz_pow_ui(buf,X,4);
    mpz_mul_ui(buf,buf,188);
    mpz_add(result,result,buf);
    mpz_pow_ui(buf,X,3);
    mpz_mul_ui(buf,buf,259);
    mpz_add(result,result,buf);
    mpz_pow_ui(buf,X,2);
    mpz_mul_ui(buf,buf,343);
    mpz_add(result,result,buf);
    mpz_mul_ui(buf,X,1763);
    mpz_add(result,result,buf);
    mpz_add_ui(result,result,2401);
    mpz_tdiv_q_ui(result,result,21);
    
    //isprime
    if(mpz_probab_prime_p(result,25)==0){
        mpz_clear(buf);
        mpz_clear(result);
        return 0;
    }else{
        mpz_set(prime,result);
        mpz_clear(buf);
        mpz_clear(result);
        return 1;
    }
}

int generate_order(){
    //r=(x^6+37x^3+343)/343
    mpz_t buf,result;
    mpz_init(buf);
    mpz_init(result);
    
    mpz_pow_ui(buf,X,6);
    mpz_set(result,buf);
    mpz_pow_ui(buf,X,3);
    mpz_mul_ui(buf,buf,37);
    mpz_add(result,result,buf);
    mpz_add_ui(result,result,343);
    mpz_tdiv_q_ui(result,result,343);
    
    //isprime
    if(mpz_probab_prime_p(result,25)==0){
        mpz_clear(buf);
        mpz_clear(result);
        return 0;
    }else{
        mpz_set(order,result);
        mpz_clear(buf);
        mpz_clear(result);
        return 1;
    }
}

void generate_trace(){
    //t=1/7(x^4+16x+7)
    mpz_t buf,result;
    mpz_init(buf);
    mpz_init(result);
    
    mpz_pow_ui(result,X,4);
    mpz_mul_ui(buf,X,16);
    mpz_add(result,result,buf);
    mpz_add_ui(result,result,7);
    mpz_tdiv_q_ui(trace,result,7);
    
    mpz_clear(buf);
    mpz_clear(result);
}

void weil(){
    mpz_t t1,t3,t9,t18,p1,p3,p9,buf;
    mpz_init(t1);
    mpz_init(t3);
    mpz_init(t9);
    mpz_init(t18);
    mpz_init(p1);
    mpz_init(p3);
    mpz_init(p9);
    mpz_init(buf);
    
    mpz_set(p1,prime);
    mpz_set(t1,trace);
    
    //t3←α^3+β^3
    mpz_pow_ui(t3,t1,3);
    mpz_mul(buf,t1,p1);
    mpz_mul_ui(buf,buf,3);
    mpz_sub(t3,t3,buf);
    mpz_pow_ui(p3,p1,3);
    
    //t9←α^9+β^9
    mpz_pow_ui(t9,t3,3);
    mpz_mul(buf,t3,p3);
    mpz_mul_ui(buf,buf,3);
    mpz_sub(t9,t9,buf);
    mpz_pow_ui(p9,p3,3);
    
    //α^18+β^18
    mpz_pow_ui(t18,t9,2);
    mpz_mul_ui(buf,p9,2);
    mpz_sub(t18,t18,buf);
    
    //#EFp
    mpz_sub(buf,p1,t1);
    mpz_add_ui(EFp_total,buf,1);
    //#EFp3
    mpz_sub(buf,p3,t3);
    mpz_add_ui(EFp3_total,buf,1);
    //#EFp9
    mpz_sub(buf,p9,t9);
    mpz_add_ui(EFp9_total,buf,1);
    //#EFp18
    mpz_pow_ui(buf,p9,2);
    mpz_sub(buf,buf,t18);
    mpz_add_ui(EFp18_total,buf,1);
    
    mpz_clear(t1);
    mpz_clear(t3);
    mpz_clear(t9);
    mpz_clear(t18);
    mpz_clear(p1);
    mpz_clear(p3);
    mpz_clear(p9);
    mpz_clear(buf);
}

void get_epsilon(){
    Fp inv,buf,result1,result2;
    Fp_init(&inv);
    Fp_init(&buf);
    Fp_init(&result1);
    Fp_init(&result2);
    
    Fp_set_ui(&buf,2);
    Fp_inv(&inv,&buf);
    mpz_sub_ui(buf.x0,prime,3);
    
    Fp_sqrt(&buf,&buf);
    Fp_sub_ui(&buf,&buf,1);
    Fp_mul(&result1,&buf,&inv);
    Fp_mul(&result2,&result1,&result1);
    
    mpz_set(epsilon1,result1.x0);
    mpz_set(epsilon2,result2.x0);
    
    Fp_clear(&inv);
    Fp_clear(&buf);
    Fp_clear(&result1);
    Fp_clear(&result2);
}

void set_basis(){
    Fp_set_ui(&Alpha.x1,1);
    Fp3_inv(&Alpha_inv,&Alpha);
}

void set_frobenius_constant(){
    int i,j;
    unsigned long int s;
    Fp tmp_alpha1,tmp_alpha2,tmp[6];
    Fp3 tmp_frobenius_constant[6],tmp_beta,tmp_gamma,tmp1,tmp2;
    mpz_t tmp_p[18],exp;
    //init
    Fp_init(&tmp_alpha1);
    Fp_init(&tmp_alpha2);
    for(i=0; i<6; i++){
        Fp_init(&tmp[i]);
    }
    for(i=0; i<6; i++){
        Fp3_init(&tmp_frobenius_constant[i]);
    }
    Fp3_init(&tmp_beta);
    Fp3_init(&tmp_gamma);
    Fp3_init(&tmp1);
    Fp3_init(&tmp2);
    mpz_init(exp);
    for(i=0; i<18; i++){
        mpz_init(tmp_p[i]);
    }
    for(i=0; i<18; i++){
	    s=(unsigned long int)i+1;
	    mpz_pow_ui(tmp_p[i],prime,s);
	}
    
    for(i=0; i<18; i++){
        //alpha
        mpz_sub_ui(exp,tmp_p[i],1);
	    mpz_tdiv_q_ui(exp,exp,3);
	    Fp_set_ui(&tmp_alpha1,2);
	    Fp_pow(&tmp_alpha1,&tmp_alpha1,exp);
	    Fp_mul(&tmp_alpha2,&tmp_alpha1,&tmp_alpha1);
	    //beta
	    mpz_sub_ui(exp,tmp_p[i],1);
	    mpz_tdiv_q_ui(exp,exp,3);
	    Fp3_pow(&tmp_beta,&Alpha,exp);
	    //gamma
	    mpz_tdiv_q_ui(exp,exp,2);
	    Fp3_pow(&tmp_gamma,&Alpha,exp);
	    Fp_set_ui(&tmp_frobenius_constant[0].x0,1);
	    Fp3_set(&tmp_frobenius_constant[1],&tmp_beta);
	    Fp3_sqr(&tmp_frobenius_constant[2],&tmp_beta);
	    Fp3_set(&tmp_frobenius_constant[3],&tmp_gamma);
	    Fp3_mul(&tmp_frobenius_constant[4],&tmp_gamma,&tmp_frobenius_constant[1]);
	    Fp3_mul(&tmp_frobenius_constant[5],&tmp_gamma,&tmp_frobenius_constant[2]);
	    
	    for(j=0; j<6; j++){
	        if(Fp_cmp_zero(&tmp_frobenius_constant[j].x0)!=0){
	            Fp_set(&tmp[j],&tmp_frobenius_constant[j].x0);
	        }else if(Fp_cmp_zero(&tmp_frobenius_constant[j].x1)!=0){
	            Fp_set(&tmp[j],&tmp_frobenius_constant[j].x1);
	        }else if(Fp_cmp_zero(&tmp_frobenius_constant[j].x2)!=0){
	            Fp_set(&tmp[j],&tmp_frobenius_constant[j].x2);
	        }
	    }
	    //f_pn
	    for(j=0; j<6; j++){
	        Fp_set(&frobenius_constant[i][3*j],&tmp[j]);
	        Fp_mul(&frobenius_constant[i][3*j+1],&tmp_alpha1,&tmp[j]);
	        Fp_mul(&frobenius_constant[i][3*j+2],&tmp_alpha2,&tmp[j]);
	    }
	    
	    Fp_inv(&tmp_alpha1,&tmp_alpha1);
	    Fp3_mul_mpz(&skew_frobenius_constant[i][0],&tmp_frobenius_constant[2],tmp_alpha1.x0);
	    Fp3_mul_mpz(&skew_frobenius_constant[i][1],&tmp_frobenius_constant[4],tmp_alpha1.x0);
	}
	
    Fp_clear(&tmp_alpha1);
    Fp_clear(&tmp_alpha2);
    Fp3_clear(&tmp_beta);
    Fp3_clear(&tmp_gamma);
    mpz_clear(exp);
    for(i=0; i<18; i++){
        mpz_clear(tmp_p[i]);
    }
    for(i=0; i<6; i++){
        Fp3_clear(&tmp_frobenius_constant[i]);
    }
}

void set_curve_coefficient(){
    Fp_set_ui(&curve_b,32);
}

/*============================================================================*/
/* Test functions                                                             */
/*============================================================================*/
/*----------------------------------------------------------------------------*/
//time
float timedifference_msec(struct timeval tv_start, struct timeval tv_end){
    return (tv_end.tv_sec - tv_start.tv_sec) * 1000.0f + (tv_end.tv_usec - tv_start.tv_usec) / 1000.0f;
}

float timedifference_usec(struct timeval tv_start, struct timeval tv_end){
    return (tv_end.tv_sec - tv_start.tv_sec) + (tv_end.tv_usec - tv_start.tv_usec);
}

/*----------------------------------------------------------------------------*/
//cost
void Init_mpz_Cost(struct mpz_Cost *cost){
    cost->mpz_mul=0;
    cost->mpz_mul_ui=0;
    cost->mpz_sqr=0;
    cost->mpz_add=0;
    cost->mpz_add_ui=0;
    cost->mpz_invert=0;
}

void Print_mpz_Cost(struct mpz_Cost *cost,char *str){
    printf("%s",str);
    printf("mpz_mul,mpz_mul_ui,mpz_sqr,mpz_add,mpz_add_ui,mpz_invert\n");
    printf("%ld,",cost->mpz_mul);
    printf("%ld,",cost->mpz_mul_ui);
    printf("%ld,",cost->mpz_sqr);
    printf("%ld,",cost->mpz_add);
    printf("%ld,",cost->mpz_add_ui);
    printf("%ld",cost->mpz_invert);
    printf("\n");
}

void Init_Fp_Cost(struct Fp_Cost *cost){
    cost->Fp_mul=0;
    cost->Fp_mul_mpz=0;
    cost->Fp_mul_ui=0;
    cost->Fp_sqr=0;
    cost->Fp_basis=0;
    cost->Fp_add=0;
    cost->Fp_add_mpz=0;
    cost->Fp_add_ui=0;
    cost->Fp_inv=0;
    cost->Fp_neg=0;
}

void Print_Fp_Cost(struct Fp_Cost *cost,char *str){
    printf("%s",str);
    //printf("Fp_mul,Fp_mul_mpz,Fp_mul_ui,Fp_sqr,Fp_basis,Fp_add,Fp_add_mpz,Fp_add_ui,Fp_inv,Fp_neg\n");
    printf("Fp_mul,Fp_mul_mpz,Fp_sqr,Fp_add,Fp_add_ui,Fp_inv,Fp_neg\n");
    printf("%ld,",cost->Fp_mul);
    printf("%ld,",cost->Fp_mul_mpz);
    //printf("%ld,",cost->Fp_mul_ui);
    printf("%ld,",cost->Fp_sqr);
    //printf("%ld,",cost->Fp_basis);
    printf("%ld,",cost->Fp_add);
    //printf("%ld,",cost->Fp_add_mpz);
    printf("%ld,",cost->Fp_add_ui);
    printf("%ld,",cost->Fp_inv);
    printf("%ld",cost->Fp_neg);
    printf("\n");
}

void test_Field(){
    printf("====================================================================================\n");
    Fp18 tmp_Fp18,test1,test2;
    Fp18_init(&tmp_Fp18);
    Fp18_init(&test1);
    Fp18_init(&test2);
    mpz_t exp;
    mpz_init(exp);
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    Fp18_set_random(&tmp_Fp18,state);
    Fp18_printf(&tmp_Fp18,"");
    printf("\n\n");
    
    printf("mul/sqr\n");
    Fp18_mul(&test1,&tmp_Fp18,&tmp_Fp18);
    Fp18_printf(&test1,"");
    printf("\n");
    
    Fp18_sqr(&test2,&tmp_Fp18);
    Fp18_printf(&test2,"");
    printf("\n");
    
    if(Fp18_cmp(&test1,&test2)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    printf("pow/inv\n");
    mpz_pow_ui(exp,prime,18);
    mpz_sub_ui(exp,exp,2);
    Fp18_pow(&test1,&tmp_Fp18,exp);
    Fp18_printf(&test1,"");
    printf("\n");
    
    Fp18_inv(&test2,&tmp_Fp18);
    Fp18_printf(&test2,"");
    printf("\n");
    
    if(Fp18_cmp(&test1,&test2)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(exp);
    Fp18_clear(&tmp_Fp18);
    Fp18_clear(&test1);
    Fp18_clear(&test2);
}

void test_Frobenius_map(){
    printf("====================================================================================\n");
    Fp18 tmp_Fp18,test1,test2;
    Fp18_init(&tmp_Fp18);
    Fp18_init(&test1);
    Fp18_init(&test2);
    mpz_t exp;
    mpz_init(exp);
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    Fp18_set_random(&tmp_Fp18,state);
    Fp18_printf(&tmp_Fp18,"");
    printf("\n\n");
    
    mpz_pow_ui(exp,prime,9);
    
    Fp18_pow(&test1,&tmp_Fp18,exp);
    Fp18_printf(&test1,"");
    printf("\n\n");
    
    Fp18_frobenius_map_p9(&test2,&tmp_Fp18);
    Fp18_printf(&test2,"");
    printf("\n\n");
    
    if(Fp18_cmp(&test1,&test2)==0){
        printf("success\n");
    }else{
        printf("failed\n");
    }
        
    mpz_clear(exp);
    Fp18_clear(&tmp_Fp18);
    Fp18_clear(&test1);
    Fp18_clear(&test2);
}

void test_skew_frobenius_map(){
    printf("====================================================================================\n");
    EFp18 Q,test1,test2;
    EFp18_init(&Q);
    EFp18_init(&test1);
    EFp18_init(&test2);
    EFp3 twisted_Q;
    EFp3_init(&twisted_Q);
    
    EFp18_generate_G2(&Q);
    EFp18_to_EFp3(&twisted_Q,&Q);
    
    Fp18_frobenius_map_p3(&test1.x,&Q.x);
    Fp18_frobenius_map_p3(&test1.y,&Q.y);
    EFp18_printf(&test1,"");
    printf("\n\n");
    
    EFp3_skew_frobenius_map_p3(&twisted_Q,&twisted_Q);
    EFp3_to_EFp18(&test2,&twisted_Q);
    EFp18_printf(&test2,"");
    printf("\n\n");
    
    if(Fp18_cmp(&test1.x,&test2.x)==0 && Fp18_cmp(&test1.y,&test2.y)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    EFp18_clear(&Q);
    EFp18_clear(&test1);
    EFp18_clear(&test2);
    EFp3_clear(&twisted_Q);
}

void test_rational_point(){
    printf("====================================================================================\n");
    printf("Rational point\n\n");
    EFp tmp_EFp;
    EFp3 tmp_EFp3;
    EFp9 tmp_EFp9;
    EFp18 tmp_EFp18;
    EFp_init(&tmp_EFp);
    EFp3_init(&tmp_EFp3);
    EFp9_init(&tmp_EFp9);
    EFp18_init(&tmp_EFp18);
    EFp18 test_G1,test_G2,random_P;
    EFp18_init(&test_G1);
    EFp18_init(&test_G2);
    EFp18_init(&random_P);
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test EFp\n");
    EFp_rational_point(&tmp_EFp);
    EFp_printf(&tmp_EFp,"");
    printf("\n");
    EFp_SCM(&tmp_EFp,&tmp_EFp,EFp_total);
    EFp_printf(&tmp_EFp,"");
    printf("\n");
    if(tmp_EFp.infinity==1){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test EFp3\n");
    EFp3_rational_point(&tmp_EFp3);
    EFp3_printf(&tmp_EFp3,"");
    printf("\n");
    EFp3_SCM(&tmp_EFp3,&tmp_EFp3,EFp3_total);
    EFp3_printf(&tmp_EFp3,"");
    printf("\n");
    if(tmp_EFp3.infinity==1){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test EFp9\n");
    EFp9_rational_point(&tmp_EFp9);
    EFp9_printf(&tmp_EFp9,"");
    printf("\n");
    EFp9_SCM(&tmp_EFp9,&tmp_EFp9,EFp9_total);
    EFp9_printf(&tmp_EFp9,"");
    printf("\n");
    if(tmp_EFp9.infinity==1){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test EFp18\n");
    EFp18_rational_point(&tmp_EFp18);
    EFp18_printf(&tmp_EFp18,"");
    printf("\n");
    EFp18_SCM(&tmp_EFp18,&tmp_EFp18,EFp18_total);
    EFp18_printf(&tmp_EFp18,"");
    printf("\n");
    if(tmp_EFp18.infinity==1){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test G1\n");
    EFp18_generate_G1(&test_G1);
    EFp18_printf(&test_G1,"");
    printf("\n");
    EFp18_SCM(&test_G1,&test_G1,order);
    EFp18_printf(&test_G1,"");
    printf("\n");
    if(test_G1.infinity==1){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test G2\n");
    EFp18_generate_G2(&test_G2);
    EFp18_printf(&test_G2,"");
    printf("\n");
    EFp18_SCM(&test_G2,&test_G2,order);
    EFp18_printf(&test_G2,"");
    printf("\n");
    
    if(test_G2.infinity==1){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    EFp_clear(&tmp_EFp);
    EFp3_clear(&tmp_EFp3);
    EFp9_clear(&tmp_EFp9);
    EFp18_clear(&tmp_EFp18);
    EFp18_clear(&test_G1);
    EFp18_clear(&test_G2);
    EFp18_clear(&random_P);
}

void test_twist(){
    printf("====================================================================================\n");
    EFp18 Q,test1,test2;
    EFp18_init(&Q);
    EFp18_init(&test1);
    EFp18_init(&test2);
    EFp3 twist_Q;
    EFp3_init(&twist_Q);
    
    
    EFp18_generate_G2(&Q);
    EFp18_printf(&Q,"Q\n");
    printf("\n\n");
    
    EFp18_to_EFp3(&twist_Q,&Q);
    EFp3_ECD(&twist_Q,&twist_Q);
    EFp3_to_EFp18(&test1,&twist_Q);
    EFp18_printf(&test1,"");
    printf("\n\n");
    
    EFp18_ECD(&test2,&Q);
    EFp18_printf(&test2,"");
    printf("\n\n");
    
    if(Fp18_cmp(&test1.x,&test2.x)==0 && Fp18_cmp(&test1.y,&test2.y)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    EFp18_clear(&Q);
    EFp18_clear(&test1);
    EFp18_clear(&test2);
    EFp3_clear(&twist_Q);
}

void test_plain_ate_pairing(){
    printf("====================================================================================\n");
    printf("Plain-ate pairing\n\n");
    EFp18 P,Q,s1P,s2P,s1Q,s2Q;
    EFp18_init(&P);
    EFp18_init(&Q);
    EFp18_init(&s1P);
    EFp18_init(&s2P);
    EFp18_init(&s1Q);
    EFp18_init(&s2Q);
    Fp18 Z,test1,test2,test3;
    Fp18_init(&Z);
    Fp18_init(&test1);
    Fp18_init(&test2);
    Fp18_init(&test3);
    mpz_t s12,s1,s2;
    mpz_init(s12);
    mpz_init(s1);
    mpz_init(s2);
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    mpz_urandomm(s1,state,order);
	mpz_urandomm(s2,state,order);
    mpz_mul(s12,s1,s2);
    mpz_mod(s12,s12,order);
    
    EFp18_generate_G1(&P);
    EFp18_printf(&P,"P\n");
    printf("\n\n");
    EFp18_generate_G2(&Q);
    EFp18_printf(&Q,"Q\n");
    printf("\n\n");
    EFp18_SCM(&s1P,&P,s1);
    EFp18_SCM(&s2P,&P,s2);
    EFp18_SCM(&s1Q,&Q,s1);
    EFp18_SCM(&s2Q,&Q,s2);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("plain_ate(Q,P)^s1*s2\n");
    Plain_ate_pairing(&Z,&Q,&P);
    Fp18_pow(&test1,&Z,s12);
    printf("Miller's Algo. (Plain-ate) : %.2f[ms]\n",MILLER_PLAINATE);
    printf("Final Exp. Plain (easy) : %.2f[ms]\n",FINALEXP_PLAIN_EASY);
    printf("Final Exp. Plain (hard) : %.2f[ms]\n",FINALEXP_PLAIN_HARD);
    Fp18_printf(&test1,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("plain_ate([s2]Q,[s1]P)\n");
    Plain_ate_pairing(&test2,&s2Q,&s1P);
    printf("Miller's Algo. (Plain-ate) : %.2f[ms]\n",MILLER_PLAINATE);
    printf("Final Exp. Plain (easy) : %.2f[ms]\n",FINALEXP_PLAIN_EASY);
    printf("Final Exp. Plain (hard) : %.2f[ms]\n",FINALEXP_PLAIN_HARD);
    Fp18_printf(&test2,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("plain_ate([s1]Q,[s2]P)\n");
    Plain_ate_pairing(&test3,&s1Q,&s2P);
    printf("Miller's Algo. (Plain-ate) : %.2f[ms]\n",MILLER_PLAINATE);
    printf("Final Exp. Plain (easy) : %.2f[ms]\n",FINALEXP_PLAIN_EASY);
    printf("Final Exp. Plain (hard) : %.2f[ms]\n",FINALEXP_PLAIN_HARD);
    Fp18_printf(&test3,"");
    printf("\n\n");
    
    printf("bilinear test\n");
    if(Fp18_cmp(&test1,&test2)==0 && Fp18_cmp(&test1,&test3)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(s12);
    mpz_clear(s1);
    mpz_clear(s2);
    EFp18_clear(&P);
    EFp18_clear(&Q);
    EFp18_clear(&s1P);
    EFp18_clear(&s2P);
    EFp18_clear(&s1Q);
    EFp18_clear(&s2Q);
    Fp18_clear(&Z);
    Fp18_clear(&test1);
    Fp18_clear(&test2);
    Fp18_clear(&test3);
}

void test_opt_ate_pairing(){
    printf("====================================================================================\n");
    printf("Opt-ate pairing\n\n");
    EFp18 P,Q,s1P,s2P,s1Q,s2Q;
    EFp18_init(&P);
    EFp18_init(&Q);
    EFp18_init(&s1P);
    EFp18_init(&s2P);
    EFp18_init(&s1Q);
    EFp18_init(&s2Q);
    Fp18 Z,test1,test2,test3;
    Fp18_init(&Z);
    Fp18_init(&test1);
    Fp18_init(&test2);
    Fp18_init(&test3);
    mpz_t s12,s1,s2;
    mpz_init(s12);
    mpz_init(s1);
    mpz_init(s2);
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    mpz_urandomm(s1,state,order);
	mpz_urandomm(s2,state,order);
    mpz_mul(s12,s1,s2);
    mpz_mod(s12,s12,order);
    
    EFp18_generate_G1(&P);
    EFp18_printf(&P,"P\n");
    printf("\n\n");
    EFp18_generate_G2(&Q);
    EFp18_printf(&Q,"Q\n");
    printf("\n\n");
    EFp18_SCM(&s1P,&P,s1);
    EFp18_SCM(&s2P,&P,s2);
    EFp18_SCM(&s1Q,&Q,s1);
    EFp18_SCM(&s2Q,&Q,s2);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("opt_ate(Q,P)^s1*s2\n");
    Opt_ate_pairing(&Z,&Q,&P);
    Fp18_pow(&test1,&Z,s12);
    printf("Miller's Algo. (Opt-ate) : %.2f[ms]\n",MILLER_OPTATE);
    printf("Final Exp. Opt (easy) : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. Opt (hard) : %.2f[ms]\n",FINALEXP_OPT_HARD);
    Fp18_printf(&test1,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("opt_ate([s2]Q,[s1]P)\n");
    Opt_ate_pairing(&test2,&s2Q,&s1P);
    printf("Miller's Algo. (Opt-ate) : %.2f[ms]\n",MILLER_OPTATE);
    printf("Final Exp. Opt (easy) : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. Opt (hard) : %.2f[ms]\n",FINALEXP_OPT_HARD);
    Fp18_printf(&test2,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("opt_ate([s1]Q,[s2]P)\n");
    Opt_ate_pairing(&test3,&s1Q,&s2P);
    printf("Miller's Algo. (Opt-ate) : %.2f[ms]\n",MILLER_OPTATE);
    printf("Final Exp. Opt (easy) : %.2f[ms]\n",FINALEXP_OPT_EASY);
    printf("Final Exp. Opt (hard) : %.2f[ms]\n",FINALEXP_OPT_HARD);
    Fp18_printf(&test3,"");
    printf("\n\n");
    
    printf("bilinear test\n");
    if(Fp18_cmp(&test1,&test2)==0 && Fp18_cmp(&test1,&test3)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(s12);
    mpz_clear(s1);
    mpz_clear(s2);
    EFp18_clear(&P);
    EFp18_clear(&Q);
    EFp18_clear(&s1P);
    EFp18_clear(&s2P);
    EFp18_clear(&s1Q);
    EFp18_clear(&s2Q);
    Fp18_clear(&Z);
    Fp18_clear(&test1);
    Fp18_clear(&test2);
    Fp18_clear(&test3);
}

void test_G1_SCM(){
    printf("====================================================================================\n");
    printf("G1 SCM\n\n");
    EFp18 P,test1,test2,test3;
    EFp18_init(&P);
    EFp18_init(&test1);
    EFp18_init(&test2);
    EFp18_init(&test3);
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    mpz_urandomm(scalar,state,order);
    
    printf("scalar : %d-bit\n",(int)mpz_sizeinbase(scalar,2));
    
    printf("generating rational point P in G1\n\n");
    EFp18_generate_G1(&P);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("plain\n");
    EFp18_G1_SCM_plain(&test1,&P,scalar);
    printf("G1 SCM (plain) : %.2f[ms]\n",G1SCM_PLAIN);
    EFp18_printf(&test1,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("2split\n");
    EFp18_G1_SCM_2split(&test2,&P,scalar);
    printf("G1 SCM (2split) : %.2f[ms]\n",G1SCM_2SPLIT);
    EFp18_printf(&test2,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("2split-JSF\n");
    EFp18_G1_SCM_2split_JSF(&test3,&P,scalar);
    printf("G1 SCM (2split-JSF) : %.2f[ms]\n",G1SCM_2SPLIT_JSF);
    EFp18_printf(&test3,"");
    printf("\n\n");
    
    if(Fp18_cmp(&test1.x,&test2.x)==0 && Fp18_cmp(&test1.y,&test2.y)==0
    && Fp18_cmp(&test1.x,&test3.x)==0 && Fp18_cmp(&test1.y,&test3.y)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(scalar);
    EFp18_clear(&P);
    EFp18_clear(&test1);
    EFp18_clear(&test2);
    EFp18_clear(&test3);
}

void test_G2_SCM(){
    printf("====================================================================================\n");
    printf("G2 SCM\n\n");
    EFp18 Q,test1,test2,test3,test4;
    EFp18_init(&Q);
    EFp18_init(&test1);
    EFp18_init(&test2);
    EFp18_init(&test3);
    EFp18_init(&test4);
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    mpz_urandomm(scalar,state,order);
    
    printf("scalar : %d-bit\n",(int)mpz_sizeinbase(scalar,2));
    
    printf("generating rational point Q in G2\n\n");
    EFp18_generate_G2(&Q);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("plain\n");
    EFp18_G2_SCM_plain(&test1,&Q,scalar);
    printf("G2 SCM (plain) : %.2f[ms]\n",G2SCM_PLAIN);
    EFp18_printf(&test1,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("2split\n");
    EFp18_G2_SCM_2split(&test2,&Q,scalar);
    printf("G2 SCM (2split) : %.2f[ms]\n",G2SCM_2SPLIT);
    EFp18_printf(&test2,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("2split-JSF\n");
    EFp18_G2_SCM_2split_JSF(&test3,&Q,scalar);
    printf("G2 SCM (2split-JSF) : %.2f[ms]\n",G2SCM_2SPLIT_JSF);
    EFp18_printf(&test3,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test4\n\n");
    printf("6split\n");
    EFp18_G2_SCM_6split(&test4,&Q,scalar);
    printf("G2 SCM (6split) : %.2f[ms]\n",G2SCM_6SPLIT);
    EFp18_printf(&test4,"");
    printf("\n\n");
    
    
    if(Fp18_cmp(&test1.x,&test2.x)==0 && Fp18_cmp(&test1.y,&test2.y)==0
    && Fp18_cmp(&test1.x,&test3.x)==0 && Fp18_cmp(&test1.y,&test3.y)==0
    && Fp18_cmp(&test1.x,&test4.x)==0 && Fp18_cmp(&test1.y,&test4.y)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(scalar);
    EFp18_clear(&Q);
    EFp18_clear(&test1);
    EFp18_clear(&test2);
    EFp18_clear(&test3);
    EFp18_clear(&test4);
}

void test_G3_EXP(){
    printf("====================================================================================\n");
    printf("G3 Exp.\n\n");
    EFp18 P,Q;
    EFp18_init(&P);
    EFp18_init(&Q);
    Fp18 Z,test1,test2,test3,test4;
    Fp18_init(&Z);
    Fp18_init(&test1);
    Fp18_init(&test2);
    Fp18_init(&test3);
    Fp18_init(&test4);
    mpz_t scalar;
    mpz_init(scalar);
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    mpz_urandomm(scalar,state,order);
    
    printf("scalar : %d-bit\n",(int)mpz_sizeinbase(scalar,2));
    
    printf("generating rational point P in G1\n\n");
    EFp18_generate_G1(&P);
    printf("generating rational point Q in G2\n\n");
    EFp18_generate_G2(&Q);
    printf("x-ate(Q,P)\n");
    Opt_ate_pairing(&Z,&Q,&P);
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test1\n\n");
    printf("plain\n");
    Fp18_G3_EXP_plain(&test1,&Z,scalar);
    printf("G3 SCM (plain) : %.2f[ms]\n",G3EXP_PLAIN);
    Fp18_printf(&test1,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test2\n\n");
    printf("2split\n");
    Fp18_G3_EXP_2split(&test2,&Z,scalar);
    printf("G3 SCM (2split) : %.2f[ms]\n",G3EXP_2SPLIT);
    Fp18_printf(&test2,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test3\n\n");
    printf("2split-JSF\n");
    Fp18_G3_EXP_2split_JSF(&test3,&Z,scalar);
    printf("G3 SCM (2split-JSF) : %.2f[ms]\n",G3EXP_2SPLIT_JSF);
    Fp18_printf(&test3,"");
    printf("\n\n");
    
    printf("------------------------------------------------------------------------------------\n");
    printf("test4\n\n");
    printf("6split\n");
    Fp18_G3_EXP_6split(&test4,&Z,scalar);
    printf("G3 SCM (6split) : %.2f[ms]\n",G3EXP_6SPLIT);
    Fp18_printf(&test4,"");
    printf("\n\n");
    
    if(Fp18_cmp(&test1,&test2)==0 && Fp18_cmp(&test1,&test3)==0 && Fp18_cmp(&test1,&test4)==0){
        printf("success\n\n");
    }else{
        printf("failed\n\n");
    }
    
    mpz_clear(scalar);
    EFp18_clear(&P);
    EFp18_clear(&Q);
    Fp18_clear(&Z);
    Fp18_clear(&test1);
    Fp18_clear(&test2);
    Fp18_clear(&test3);
    Fp18_clear(&test4);
}

void computation_time(){
    printf("====================================================================================\n");
    printf("Computation time\n\n"); 
    int count;
    float AVE_MILLER_PLAINATE=0,AVE_MILLER_OPTATE=0;
    float AVE_FINALEXP_PLAIN_EASY=0,AVE_FINALEXP_PLAIN_HARD=0,AVE_FINALEXP_OPT_EASY=0,AVE_FINALEXP_OPT_HARD=0;
    float AVE_G1SCM_PLAIN=0,AVE_G1SCM_2SPLIT=0,AVE_G1SCM_2SPLIT_JSF=0;
    float AVE_G2SCM_PLAIN=0,AVE_G2SCM_2SPLIT=0,AVE_G2SCM_2SPLIT_JSF=0,AVE_G2SCM_6SPLIT=0;
    float AVE_G3EXP_PLAIN=0,AVE_G3EXP_2SPLIT=0,AVE_G3EXP_2SPLIT_JSF=0,AVE_G3EXP_6SPLIT=0;
    EFp18 P,Q,tmp_P,tmp_Q;
    EFp18_init(&P);
    EFp18_init(&Q);
    EFp18_init(&tmp_P);
    EFp18_init(&tmp_Q);
    Fp18 Z,tmp_Z;
    Fp18_init(&Z);
    Fp18_init(&tmp_Z);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    mpz_t scalar;
    mpz_init(scalar);
    
    printf("generating rational point P in G1\n\n");
    EFp18_generate_G1(&P);
    EFp18_printf(&P,"P\n");
    printf("\n\n");
    
    printf("generating rational point Q in G2\n\n");
    EFp18_generate_G2(&Q);
    EFp18_printf(&Q,"Q\n");
    printf("\n\n");
    
    printf("Computing...\n");
    count=0;
    while(count<100){
        printf("%d\n",count);
        //plain
        Plain_ate_pairing(&Z,&Q,&P);
        AVE_MILLER_PLAINATE=AVE_MILLER_PLAINATE+MILLER_PLAINATE;
        AVE_FINALEXP_PLAIN_EASY=AVE_FINALEXP_PLAIN_EASY+FINALEXP_PLAIN_EASY;
        AVE_FINALEXP_PLAIN_HARD=AVE_FINALEXP_PLAIN_HARD+FINALEXP_PLAIN_HARD;
        
        //opt
        Opt_ate_pairing(&Z,&Q,&P);
        AVE_MILLER_OPTATE=AVE_MILLER_OPTATE+MILLER_OPTATE;
        AVE_FINALEXP_OPT_EASY=AVE_FINALEXP_OPT_EASY+FINALEXP_OPT_EASY;
        AVE_FINALEXP_OPT_HARD=AVE_FINALEXP_OPT_HARD+FINALEXP_OPT_HARD;
        
        mpz_urandomm(scalar,state,order);
        
        //G1_SCM
        EFp18_G1_SCM_plain(&tmp_P,&P,scalar);
        EFp18_G1_SCM_2split(&tmp_P,&P,scalar);
        EFp18_G1_SCM_2split_JSF(&tmp_P,&P,scalar);
        AVE_G1SCM_PLAIN=AVE_G1SCM_PLAIN+G1SCM_PLAIN;
        AVE_G1SCM_2SPLIT=AVE_G1SCM_2SPLIT+G1SCM_2SPLIT;
        AVE_G1SCM_2SPLIT_JSF=AVE_G1SCM_2SPLIT_JSF+G1SCM_2SPLIT_JSF;
        
        //G2_SCM
        EFp18_G2_SCM_plain(&tmp_Q,&Q,scalar);
        EFp18_G2_SCM_2split(&tmp_Q,&Q,scalar);
        EFp18_G2_SCM_2split_JSF(&tmp_Q,&Q,scalar);
        EFp18_G2_SCM_6split(&tmp_Q,&Q,scalar);
        AVE_G2SCM_PLAIN=AVE_G2SCM_PLAIN+G2SCM_PLAIN;
        AVE_G2SCM_2SPLIT=AVE_G2SCM_2SPLIT+G2SCM_2SPLIT;
        AVE_G2SCM_2SPLIT_JSF=AVE_G2SCM_2SPLIT_JSF+G2SCM_2SPLIT_JSF;
        AVE_G2SCM_6SPLIT=AVE_G2SCM_6SPLIT+G2SCM_6SPLIT;
        
        //G3_EXP
        Fp18_G3_EXP_plain(&tmp_Z,&Z,scalar);
        Fp18_G3_EXP_2split(&tmp_Z,&Z,scalar);
        Fp18_G3_EXP_2split_JSF(&tmp_Z,&Z,scalar);
        Fp18_G3_EXP_6split(&tmp_Z,&Z,scalar);
        AVE_G3EXP_PLAIN=AVE_G3EXP_PLAIN+G3EXP_PLAIN;
        AVE_G3EXP_2SPLIT=AVE_G3EXP_2SPLIT+G3EXP_2SPLIT;
        AVE_G3EXP_2SPLIT_JSF=AVE_G3EXP_2SPLIT_JSF+G3EXP_2SPLIT_JSF;
        AVE_G3EXP_6SPLIT=AVE_G3EXP_6SPLIT+G3EXP_6SPLIT;
        
        count++;
    }
    printf("\n");
    
    printf("Miller's Algo. (plain-ate) : %.2f[ms]\n",AVE_MILLER_PLAINATE/100);
    printf("Miller's Algo. (opt-ate)   : %.2f[ms]\n",AVE_MILLER_OPTATE/100);
    
    printf("\n");
    printf("Final Exp. Plain   (easy) : %.2f[ms]\n",AVE_FINALEXP_PLAIN_EASY/100);
    printf("Final Exp. Plain   (hard) : %.2f[ms]\n",AVE_FINALEXP_PLAIN_HARD/100);
    printf("Final Exp. Optimal (easy) : %.2f[ms]\n",AVE_FINALEXP_OPT_EASY/100);
    printf("Final Exp. Optimal (hard) : %.2f[ms]\n",AVE_FINALEXP_OPT_HARD/100);
    
    printf("G1 SCM (plain)            : %.2f[ms]\n",AVE_G1SCM_PLAIN/100);
    printf("G1 SCM (2split)           : %.2f[ms]\n",AVE_G1SCM_2SPLIT/100);
    printf("G1 SCM (2split-JSF)       : %.2f[ms]\n",AVE_G1SCM_2SPLIT_JSF/100);
    printf("\n");
    
    printf("G2 SCM (plain)            : %.2f[ms]\n",AVE_G2SCM_PLAIN/100);
    printf("G2 SCM (2split)           : %.2f[ms]\n",AVE_G2SCM_2SPLIT/100);
    printf("G2 SCM (2split-JSF)       : %.2f[ms]\n",AVE_G2SCM_2SPLIT_JSF/100);
    printf("G2 SCM (6split)           : %.2f[ms]\n",AVE_G2SCM_6SPLIT/100);
    printf("\n");
    
    printf("G3 EXP (plain)            : %.2f[ms]\n",AVE_G3EXP_PLAIN/100);
    printf("G3 EXP (2split)           : %.2f[ms]\n",AVE_G3EXP_2SPLIT/100);
    printf("G3 EXP (2split-JSF)       : %.2f[ms]\n",AVE_G3EXP_2SPLIT_JSF/100);
    printf("G3 EXP (6split)           : %.2f[ms]\n",AVE_G3EXP_6SPLIT/100);
    printf("\n");
    
    EFp18_clear(&P);
    EFp18_clear(&Q);
    EFp18_clear(&tmp_P);
    EFp18_clear(&tmp_Q);
    Fp18_clear(&Z);
    Fp18_clear(&tmp_Z);
    mpz_clear(scalar);
}

void computation_cost(){
    printf("====================================================================================\n");
    printf("Computation cost\n\n");
    int count;
    EFp18 P,Q,tmp_P,tmp_Q;
    EFp18_init(&P);
    EFp18_init(&Q);
    EFp18_init(&tmp_P);
    EFp18_init(&tmp_Q);
    Fp18 Z,tmp_Z;
    Fp18_init(&Z);
    Fp18_init(&tmp_Z);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    mpz_t scalar;
    mpz_init(scalar);
    
    printf("generating rational point P in G1\n\n");
    EFp18_generate_G1(&P);
    
    printf("generating rational point Q in G2\n\n");
    EFp18_generate_G2(&Q);
    
    //plain
    Init_Fp_Cost(&Fp_cost);
    Miller_algo_for_plain_ate(&Z,&P,&Q);
    Print_Fp_Cost(&Fp_cost,"Plain-ate\n");
    Init_Fp_Cost(&Fp_cost);
    Final_exp_easy(&Z,&Z);
    Print_Fp_Cost(&Fp_cost,"Final exp (easy)\n");
    Init_Fp_Cost(&Fp_cost);
    Final_exp_hard_plain(&Z,&Z);
    Print_Fp_Cost(&Fp_cost,"Final exp Plain (hard)\n");
    printf("\n");
    
    //opt
    Init_Fp_Cost(&Fp_cost);
    Miller_algo_for_opt_ate(&Z,&P,&Q);
    Print_Fp_Cost(&Fp_cost,"Opt-ate\n");
    Init_Fp_Cost(&Fp_cost);
    Final_exp_easy(&Z,&Z);
    Print_Fp_Cost(&Fp_cost,"Final exp (easy)\n");
    Init_Fp_Cost(&Fp_cost);
    Final_exp_hard_optimal(&Z,&Z);
    Print_Fp_Cost(&Fp_cost,"Final exp Plain (hard)\n");
    printf("\n");
    
    //G1 SCM
    Init_Fp_Cost(&Fp_cost_ave);
    count=0;
    mpz_urandomm(scalar,state,order);
    while(count<100){
        Init_Fp_Cost(&Fp_cost);
        EFp18_G1_SCM_2split_JSF(&P,&P,scalar);
        Fp_cost_ave.Fp_mul=Fp_cost_ave.Fp_mul+Fp_cost.Fp_mul;
        Fp_cost_ave.Fp_mul_mpz=Fp_cost_ave.Fp_mul_mpz+Fp_cost.Fp_mul_mpz;
        Fp_cost_ave.Fp_mul_ui=Fp_cost_ave.Fp_mul_ui+Fp_cost.Fp_mul_ui;
        Fp_cost_ave.Fp_sqr=Fp_cost_ave.Fp_sqr+Fp_cost.Fp_sqr;
        Fp_cost_ave.Fp_basis=Fp_cost_ave.Fp_basis+Fp_cost.Fp_basis;
        Fp_cost_ave.Fp_add=Fp_cost_ave.Fp_add+Fp_cost.Fp_add;
        Fp_cost_ave.Fp_add_mpz=Fp_cost_ave.Fp_add_mpz+Fp_cost.Fp_add_mpz;
        Fp_cost_ave.Fp_add_ui=Fp_cost_ave.Fp_add_ui+Fp_cost.Fp_add_ui;
        Fp_cost_ave.Fp_inv=Fp_cost_ave.Fp_inv+Fp_cost.Fp_inv;
        Fp_cost_ave.Fp_neg=Fp_cost_ave.Fp_neg+Fp_cost.Fp_neg;
        count++;
        mpz_urandomm(scalar,state,order);
    }
    Fp_cost_ave.Fp_mul=Fp_cost_ave.Fp_mul/count;
    Fp_cost_ave.Fp_mul_mpz=Fp_cost_ave.Fp_mul_mpz/count;
    Fp_cost_ave.Fp_mul_ui=Fp_cost_ave.Fp_mul_ui/count;
    Fp_cost_ave.Fp_sqr=Fp_cost_ave.Fp_sqr/count;
    Fp_cost_ave.Fp_basis=Fp_cost_ave.Fp_basis/count;
    Fp_cost_ave.Fp_add=Fp_cost_ave.Fp_add/count;
    Fp_cost_ave.Fp_add_mpz=Fp_cost_ave.Fp_add_mpz/count;
    Fp_cost_ave.Fp_add_ui=Fp_cost_ave.Fp_add_ui/count;
    Fp_cost_ave.Fp_inv=Fp_cost_ave.Fp_inv/count;
    Fp_cost_ave.Fp_neg=Fp_cost_ave.Fp_neg/count;
    
    Print_Fp_Cost(&Fp_cost_ave,"G1 SCM (2-split JSF)");
    printf("\n");
    
    //G2 SCM
    Init_Fp_Cost(&Fp_cost_ave);
    count=0;
    mpz_urandomm(scalar,state,order);
    while(count<100){
        Init_Fp_Cost(&Fp_cost);
        EFp18_G2_SCM_6split(&Q,&Q,scalar);
        Fp_cost_ave.Fp_mul=Fp_cost_ave.Fp_mul+Fp_cost.Fp_mul;
        Fp_cost_ave.Fp_mul_mpz=Fp_cost_ave.Fp_mul_mpz+Fp_cost.Fp_mul_mpz;
        Fp_cost_ave.Fp_mul_ui=Fp_cost_ave.Fp_mul_ui+Fp_cost.Fp_mul_ui;
        Fp_cost_ave.Fp_sqr=Fp_cost_ave.Fp_sqr+Fp_cost.Fp_sqr;
        Fp_cost_ave.Fp_basis=Fp_cost_ave.Fp_basis+Fp_cost.Fp_basis;
        Fp_cost_ave.Fp_add=Fp_cost_ave.Fp_add+Fp_cost.Fp_add;
        Fp_cost_ave.Fp_add_mpz=Fp_cost_ave.Fp_add_mpz+Fp_cost.Fp_add_mpz;
        Fp_cost_ave.Fp_add_ui=Fp_cost_ave.Fp_add_ui+Fp_cost.Fp_add_ui;
        Fp_cost_ave.Fp_inv=Fp_cost_ave.Fp_inv+Fp_cost.Fp_inv;
        Fp_cost_ave.Fp_neg=Fp_cost_ave.Fp_neg+Fp_cost.Fp_neg;
        count++;
        mpz_urandomm(scalar,state,order);
    }
    Fp_cost_ave.Fp_mul=Fp_cost_ave.Fp_mul/count;
    Fp_cost_ave.Fp_mul_mpz=Fp_cost_ave.Fp_mul_mpz/count;
    Fp_cost_ave.Fp_mul_ui=Fp_cost_ave.Fp_mul_ui/count;
    Fp_cost_ave.Fp_sqr=Fp_cost_ave.Fp_sqr/count;
    Fp_cost_ave.Fp_basis=Fp_cost_ave.Fp_basis/count;
    Fp_cost_ave.Fp_add=Fp_cost_ave.Fp_add/count;
    Fp_cost_ave.Fp_add_mpz=Fp_cost_ave.Fp_add_mpz/count;
    Fp_cost_ave.Fp_add_ui=Fp_cost_ave.Fp_add_ui/count;
    Fp_cost_ave.Fp_inv=Fp_cost_ave.Fp_inv/count;
    Fp_cost_ave.Fp_neg=Fp_cost_ave.Fp_neg/count;
    
    Print_Fp_Cost(&Fp_cost_ave,"G2 SCM (6-split)");
    printf("\n");
    
    //G3 EXP
    Init_Fp_Cost(&Fp_cost_ave);
    count=0;
    mpz_urandomm(scalar,state,order);
    while(count<100){
        Init_Fp_Cost(&Fp_cost);
        Fp18_G3_EXP_6split(&Z,&Z,scalar);
        Fp_cost_ave.Fp_mul=Fp_cost_ave.Fp_mul+Fp_cost.Fp_mul;
        Fp_cost_ave.Fp_mul_mpz=Fp_cost_ave.Fp_mul_mpz+Fp_cost.Fp_mul_mpz;
        Fp_cost_ave.Fp_mul_ui=Fp_cost_ave.Fp_mul_ui+Fp_cost.Fp_mul_ui;
        Fp_cost_ave.Fp_sqr=Fp_cost_ave.Fp_sqr+Fp_cost.Fp_sqr;
        Fp_cost_ave.Fp_basis=Fp_cost_ave.Fp_basis+Fp_cost.Fp_basis;
        Fp_cost_ave.Fp_add=Fp_cost_ave.Fp_add+Fp_cost.Fp_add;
        Fp_cost_ave.Fp_add_mpz=Fp_cost_ave.Fp_add_mpz+Fp_cost.Fp_add_mpz;
        Fp_cost_ave.Fp_add_ui=Fp_cost_ave.Fp_add_ui+Fp_cost.Fp_add_ui;
        Fp_cost_ave.Fp_inv=Fp_cost_ave.Fp_inv+Fp_cost.Fp_inv;
        Fp_cost_ave.Fp_neg=Fp_cost_ave.Fp_neg+Fp_cost.Fp_neg;
        count++;
        mpz_urandomm(scalar,state,order);
    }
    Fp_cost_ave.Fp_mul=Fp_cost_ave.Fp_mul/count;
    Fp_cost_ave.Fp_mul_mpz=Fp_cost_ave.Fp_mul_mpz/count;
    Fp_cost_ave.Fp_mul_ui=Fp_cost_ave.Fp_mul_ui/count;
    Fp_cost_ave.Fp_sqr=Fp_cost_ave.Fp_sqr/count;
    Fp_cost_ave.Fp_basis=Fp_cost_ave.Fp_basis/count;
    Fp_cost_ave.Fp_add=Fp_cost_ave.Fp_add/count;
    Fp_cost_ave.Fp_add_mpz=Fp_cost_ave.Fp_add_mpz/count;
    Fp_cost_ave.Fp_add_ui=Fp_cost_ave.Fp_add_ui/count;
    Fp_cost_ave.Fp_inv=Fp_cost_ave.Fp_inv/count;
    Fp_cost_ave.Fp_neg=Fp_cost_ave.Fp_neg/count;
    
    Print_Fp_Cost(&Fp_cost_ave,"G3 Exp. (6-split)");
    printf("\n");
    
    //G1_SCM
    /*printf("G1 SCM\n");
     Init_Fp_Cost(&Fp_cost);
     EFp18_G1_SCM_plain(&tmp_P,&P,scalar);
     Print_Fp_Cost(&Fp_cost,"plain\n");
     Init_Fp_Cost(&Fp_cost);
     EFp18_G1_SCM_2split(&tmp_P,&P,scalar);
     Print_Fp_Cost(&Fp_cost,"2split\n");
     Init_Fp_Cost(&Fp_cost);
     EFp18_G1_SCM_2split_JSF(&tmp_P,&P,scalar);
     Print_Fp_Cost(&Fp_cost,"2split-JSF\n");
     printf("\n");
     
     //G2_SCM
     printf("G2 SCM\n");
     Init_Fp_Cost(&Fp_cost);
     EFp18_G2_SCM_plain(&tmp_Q,&Q,scalar);
     Print_Fp_Cost(&Fp_cost,"plain\n");
     Init_Fp_Cost(&Fp_cost);
     EFp18_G2_SCM_2split(&tmp_Q,&Q,scalar);
     Print_Fp_Cost(&Fp_cost,"2split\n");
     Init_Fp_Cost(&Fp_cost);
     EFp18_G2_SCM_2split_JSF(&tmp_Q,&Q,scalar);
     Print_Fp_Cost(&Fp_cost,"2split-JSF\n");
     Init_Fp_Cost(&Fp_cost);
     EFp18_G2_SCM_6split(&tmp_Q,&Q,scalar);
     Print_Fp_Cost(&Fp_cost,"6split\n");
     printf("\n");
     
     //G3_EXP
     printf("G3 EXP\n");
     Init_Fp_Cost(&Fp_cost);
     Fp18_G3_EXP_plain(&tmp_Z,&Z,scalar);
     Print_Fp_Cost(&Fp_cost,"plain\n");
     Init_Fp_Cost(&Fp_cost);
     Fp18_G3_EXP_2split(&tmp_Z,&Z,scalar);
     Print_Fp_Cost(&Fp_cost,"2split\n");
     Init_Fp_Cost(&Fp_cost);
     Fp18_G3_EXP_2split_JSF(&tmp_Z,&Z,scalar);
     Print_Fp_Cost(&Fp_cost,"2split-JSF\n");
     Init_Fp_Cost(&Fp_cost);
     Fp18_G3_EXP_6split(&tmp_Z,&Z,scalar);
     Print_Fp_Cost(&Fp_cost,"6split\n");
     printf("\n");*/
    
    //Cyclotomic Squaring
    Init_Fp_Cost(&Fp_cost);
    Fp18_sqr_cyclotomic(&tmp_Z,&Z);
    Print_Fp_Cost(&Fp_cost,"Cyclotomic Squaring\n");
    
    EFp18_clear(&P);
    EFp18_clear(&Q);
    EFp18_clear(&tmp_P);
    EFp18_clear(&tmp_Q);
    Fp18_clear(&Z);
    Fp18_clear(&tmp_Z);
    mpz_clear(scalar);
}

void finite_field_cost(){
    printf("====================================================================================\n");
    Fp3 tmp1_Fp3,tmp2_Fp3;
    Fp3_init(&tmp1_Fp3);
    Fp3_init(&tmp2_Fp3);
    Fp9 tmp1_Fp9,tmp2_Fp9;
    Fp9_init(&tmp1_Fp9);
    Fp9_init(&tmp2_Fp9);
    Fp18 tmp1_Fp18,tmp2_Fp18;
    Fp18_init(&tmp1_Fp18);
    Fp18_init(&tmp2_Fp18);
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
	gmp_randseed_ui(state,(unsigned long)time(NULL));
	
	printf("Fp3\n");
    Fp3_set_random(&tmp1_Fp3,state);
    Fp3_set_random(&tmp2_Fp3,state);
    Init_Fp_Cost(&Fp_cost);
    Fp3_mul(&tmp1_Fp3,&tmp1_Fp3,&tmp2_Fp3);
    Print_Fp_Cost(&Fp_cost,"mul\n");
    Init_Fp_Cost(&Fp_cost);
    Fp3_sqr(&tmp1_Fp3,&tmp1_Fp3);
    Print_Fp_Cost(&Fp_cost,"sqr\n");
    Init_Fp_Cost(&Fp_cost);
    Fp3_inv(&tmp1_Fp3,&tmp2_Fp3);
    Print_Fp_Cost(&Fp_cost,"inv\n");
    Init_Fp_Cost(&Fp_cost);
    Fp3_mul_basis(&tmp1_Fp3,&tmp2_Fp3);
    Print_Fp_Cost(&Fp_cost,"mul by beta^3\n");
    Init_Fp_Cost(&Fp_cost);
    Fp3_inv_basis(&tmp1_Fp3,&tmp2_Fp3);
    Print_Fp_Cost(&Fp_cost,"mul by beta^{-3}\n");
    printf("\n");
	
	printf("Fp9\n");
    Fp9_set_random(&tmp1_Fp9,state);
    Fp9_set_random(&tmp2_Fp9,state);
    Init_Fp_Cost(&Fp_cost);
    Fp9_mul(&tmp1_Fp9,&tmp1_Fp9,&tmp2_Fp9);
    Print_Fp_Cost(&Fp_cost,"mul\n");
    Init_Fp_Cost(&Fp_cost);
    Fp9_sqr(&tmp1_Fp9,&tmp1_Fp9);
    Print_Fp_Cost(&Fp_cost,"sqr\n");
    Init_Fp_Cost(&Fp_cost);
    Fp9_inv(&tmp1_Fp9,&tmp2_Fp9);
    Print_Fp_Cost(&Fp_cost,"inv\n");
    printf("\n");
	
	printf("Fp18\n");
    Fp18_set_random(&tmp1_Fp18,state);
    Fp18_set_random(&tmp2_Fp18,state);
    Init_Fp_Cost(&Fp_cost);
    Fp18_mul(&tmp1_Fp18,&tmp1_Fp18,&tmp2_Fp18);
    Print_Fp_Cost(&Fp_cost,"mul\n");
    Init_Fp_Cost(&Fp_cost);
    Fp18_sqr(&tmp1_Fp18,&tmp1_Fp18);
    Print_Fp_Cost(&Fp_cost,"sqr\n");
    Init_Fp_Cost(&Fp_cost);
    Fp18_inv(&tmp1_Fp18,&tmp2_Fp18);
    Print_Fp_Cost(&Fp_cost,"inv\n");
    Init_Fp_Cost(&Fp_cost);
    Fp18_frobenius_map_p1(&tmp1_Fp18,&tmp2_Fp18);
    Print_Fp_Cost(&Fp_cost,"p-power\n");
    Init_Fp_Cost(&Fp_cost);
    Fp18_frobenius_map_p2(&tmp1_Fp18,&tmp2_Fp18);
    Print_Fp_Cost(&Fp_cost,"p2-power\n");
    Init_Fp_Cost(&Fp_cost);
    Fp18_frobenius_map_p3(&tmp1_Fp18,&tmp2_Fp18);
    Print_Fp_Cost(&Fp_cost,"p3-power\n");
    Init_Fp_Cost(&Fp_cost);
    Fp18_frobenius_map_p4(&tmp1_Fp18,&tmp2_Fp18);
    Print_Fp_Cost(&Fp_cost,"p4-power\n");
    Init_Fp_Cost(&Fp_cost);
    Fp18_frobenius_map_p5(&tmp1_Fp18,&tmp2_Fp18);
    Print_Fp_Cost(&Fp_cost,"p5-power\n");
    Init_Fp_Cost(&Fp_cost);
    Fp18_frobenius_map_p6(&tmp1_Fp18,&tmp2_Fp18);
    Print_Fp_Cost(&Fp_cost,"p6-power\n");
    Init_Fp_Cost(&Fp_cost);
    Fp18_frobenius_map_p9(&tmp1_Fp18,&tmp2_Fp18);
    Print_Fp_Cost(&Fp_cost,"p9-power\n");
    printf("\n");
    
    Fp3_clear(&tmp1_Fp3);
    Fp3_clear(&tmp2_Fp3);
    Fp9_clear(&tmp1_Fp9);
    Fp9_clear(&tmp2_Fp9);
    Fp18_clear(&tmp1_Fp18);
    Fp18_clear(&tmp2_Fp18);
}

void ecc_cost(){
    printf("====================================================================================\n");
    EFp3 tmp1_EFp3,tmp2_EFp3;
    EFp3_init(&tmp1_EFp3);
    EFp3_init(&tmp2_EFp3);
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    EFp3_rational_point(&tmp1_EFp3);
    sleep(1);
    EFp3_rational_point(&tmp2_EFp3);
    
    EFp3_printf(&tmp1_EFp3,"");
    printf("\n");
    EFp3_printf(&tmp2_EFp3,"");
    printf("\n");
    
    Init_Fp_Cost(&Fp_cost);
    EFp3_ECA(&tmp1_EFp3,&tmp1_EFp3,&tmp2_EFp3);
    Print_Fp_Cost(&Fp_cost,"ECA\n");
    Init_Fp_Cost(&Fp_cost);
    EFp3_ECD(&tmp1_EFp3,&tmp1_EFp3);
    Print_Fp_Cost(&Fp_cost,"ECD\n");
    Init_Fp_Cost(&Fp_cost);
    EFp3_set_neg(&tmp1_EFp3,&tmp1_EFp3);
    Print_Fp_Cost(&Fp_cost,"Neg.\n");
    
    Init_Fp_Cost(&Fp_cost);
    EFp3_skew_frobenius_map_p1(&tmp1_EFp3,&tmp1_EFp3);
    Print_Fp_Cost(&Fp_cost,"Skew frobenius p1\n");
    Init_Fp_Cost(&Fp_cost);
    EFp3_skew_frobenius_map_p2(&tmp1_EFp3,&tmp1_EFp3);
    Print_Fp_Cost(&Fp_cost,"Skew frobenius p2\n");
    Init_Fp_Cost(&Fp_cost);
    EFp3_skew_frobenius_map_p3(&tmp1_EFp3,&tmp1_EFp3);
    Print_Fp_Cost(&Fp_cost,"Skew frobenius p3\n");
    
    EFp3_clear(&tmp1_EFp3);
    EFp3_clear(&tmp2_EFp3);
}
