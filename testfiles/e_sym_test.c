#include "miller_header.h"

int main (void){
    // 𝑒([𝑏]𝑃, [𝑎]𝑄) = 𝑒([𝑎]𝑃, [𝑏]𝑄) = 𝑒(𝑃, 𝑄)^𝑎𝑏
    struct efp12 P, Q;
    efp12_init(&P);efp12_init(&Q);
    struct efp12 tempP, tempQ;
    efp12_init(&tempP); efp12_init(&tempQ);
    struct fp12 f1, f2, f3;
    fp12_init(&f1);fp12_init(&f2);fp12_init(&f3);

    mpz_t a, b;mpz_inits(a, b, NULL);
    mpz_set_str(a, "100", 10);
    mpz_set_str(b, "200", 10);

    mpz_t z;mpz_t p;mpz_t r;mpz_t t;
    mpz_inits(z, p, r, t, NULL);
    gen_params(z, p, r, t);

    G_random(&P);
    G_random(&Q);

    // １回目
    printf("\nf1start\n");
    efp12_scm(&tempP, P, b, p);
    efp12_scm(&tempQ, Q, a, p);
    symmetric_miller(&f1, tempP, tempQ);
    gmp_printf("pair1 : ");
    fp12_printf(f1);

    // ２回目
    printf("\nf2start\n");
    efp12_scm(&tempP, P, a, p);
    efp12_scm(&tempQ, Q, b, p);
    symmetric_miller(&f2, tempP, tempQ);
    gmp_printf("pair2 : ");
    fp12_printf(f2);

    // ３回目
    printf("\nf3start\n");
    mpz_t s;mpz_init(s);mpz_mul(s, a, b);
    gmp_printf("ab : %Zd\n", s);
    symmetric_miller(&f3, P, Q);
    fp12_pow(&f3, f3, s, p);
    gmp_printf("pair3 : ");
    fp12_printf(f3);

    if(fp12_cmp(f1, f2) == 0 && fp12_cmp(f2, f3) == 0){
        printf("\nSuccess!!\n");
    }
    else{
        printf("\nMiss...\n");
    }

    fp12_clear(&tempP.x);fp12_clear(&tempP.y);fp12_clear(&tempQ.x);fp12_clear(&tempQ.y);
    fp12_clear(&f1);fp12_clear(&f2);fp12_clear(&f3);
    mpz_clear(s);
}