#include "miller_header.h"
 
int main(void){
    // パラメータ定義
    gmp_randstate_t state;
    gmp_randinit_default(state);
    make_state(state);
    mpz_t z;mpz_t p;mpz_t r;mpz_t t;
    mpz_inits(z, p, r, t, NULL);
    gen_params(z, p, r, t);
    struct fp b;mpz_init(b.x0);
    mpz_set_str(b.x0, "2806781539090543763928146397551071025921865095800381583843579968964127551432039332258992094003963260740981125881345582810061579481053866112", 10);
 
 
    // ciにランダムな値を生成
    struct efp12 c1;efp12_init(&c1);
    efp12_random(&c1, b, p, state);
    printf("c1 : ");efp12_printf(c1);printf("\n");
 
    struct fp12 c2;fp12_init(&c2);
    fp12_random(&c2, p, state);
    printf("c2 : ");fp12_printf(c2);printf("\n");
 
    struct fp12 c3;fp12_init(&c3);
    fp12_random(&c3, p, state);
    printf("c3 : ");fp12_printf(c3);printf("\n");
 
    struct fp12 c4;fp12_init(&c4);
    fp12_random(&c4, p, state);
    printf("c4 : ");fp12_printf(c4);printf("\n");
 
    struct fp12 c5;fp12_init(&c5);
    fp12_random(&c5, p, state);
    printf("c5 : ");fp12_printf(c5);printf("\n");
 
    printf("\n");
 
 
    // gammaの計算
    mpz_t gamma;mpz_init(gamma);
    Gamma(gamma, c1, c2, c3, c4, p);
    printf("SHA-512 hash by gamma(mpz_t): ");
    gmp_printf("%Zd",gamma);;
    printf("\n");

    // 決定性の確認
    Gamma(gamma, c1, c2, c3, c4, p);
    printf("deterministic by gamma(mpz_t): ");
    gmp_printf("%Zd",gamma);;
    printf("\n");
 
 
    // functionの計算
    mpz_t function;mpz_init(function);
    Function(function, c5);
    printf("SHA-512 hash by function(mpz_t): ");
    gmp_printf("%Zd",function);;
    printf("\n");

    // 決定性の確認
    Function(function, c5);
    printf("deterministic by function(mpz_t): ");
    gmp_printf("%Zd",function);;
    printf("\n");
 
    printf("last check\n");
    printf("c1 : ");efp12_printf(c1);printf("\n");
 
    printf("c2 : ");fp12_printf(c2);printf("\n");
 
    printf("c3 : ");fp12_printf(c3);printf("\n");
 
    printf("c4 : ");fp12_printf(c4);printf("\n");
 
    printf("c5 : ");fp12_printf(c5);printf("\n");
 
 
    // メモリを解放
    mpz_clear(p);
    mpz_clear(gamma);
    efp12_clear(&c1);
    fp12_clear(&c2);
    fp12_clear(&c3);
    fp12_clear(&c4);
    fp12_clear(&c5);
 
    return 0;
}