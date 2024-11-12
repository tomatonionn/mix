#include "miller_header.h"

int main (void){
    // パラメータ定義
    printf("start\n");
    gmp_randstate_t state;
    gmp_randinit_default(state);
    make_state(state);

    mpz_t z;mpz_t p;mpz_t r;mpz_t t;
    mpz_inits(z, p, r, t, NULL);
    gen_params(z, p, r, t);
    print_params(z, p, r, t);
    printf("\n");

    mpz_t E;mpz_init(E);
    rank_number(E, 12, t, p);
    gmp_printf("有利点数:%Zd\n\n", E);

    struct fp b;mpz_init(b.x0);
    mpz_set_str(b.x0, "2806781539090543763928146397551071025921865095800381583843579968964127551432039332258992094003963260740981125881345582810061579481053866112", 10);

    struct efp12 R;efp12_init(&R);

    struct efp12 P;fp12_init(&P.x);fp12_init(&P.y);P.inf = 0;
    generate1(&P, b, z, r, p, state);
    printf("P : ");
    efp12_printf(P);

    struct efp12 Q;fp12_init(&Q.x);fp12_init(&Q.y);Q.inf = 0;
    generate2(&Q, b, E, r, p, state);
    printf("Q : ");
    efp12_printf(Q);
    
    struct fp12 S;fp12_init(&S);
    struct fp12 f;fp12_init(&f);


    // // 時間測定
    // struct timespec ts;
    // struct tm tm;
    // long double start_sec, start_nsec, end_sec, end_nsec, exe_time_buf0;
    // clock_gettime(CLOCK_REALTIME, &ts);
    // start_sec = ts.tv_sec;
    // start_nsec = ts.tv_nsec;
    // int loop = 1000;
    // for (int i = 0; i < loop; i ++){
    //     optimal_ate_miller(&f, z, P, Q, p);
    //     // final_exponentiation(&S, f, r, p);
    //     hard_final_exponentiation(&S, f, r, p, z);
    // }
    // clock_gettime(CLOCK_REALTIME, &ts);
    // end_sec = ts.tv_sec;
    // end_nsec = ts.tv_nsec;
    // exe_time_buf0 = ((end_sec - start_sec) + (end_nsec - start_nsec) / 1000000000.0)*1000.0/loop;
    // printf("Execution time!! : %Lf[ms]\n", exe_time_buf0);


    // // 乗算回数測定
    // int Maverage;
    // // ミラー
    // Mcounter = 0;
    // for (int i = 0; i < loop; i ++){
    //     optimal_ate_miller(&f, z, P, Q, p);
    // }
    // Maverage = Mcounter / loop;
    // printf("ミラーの掛け算の回数 : %d\n", Maverage);

    // // 最終べき
    // Mcounter = 0;
    // for (int i = 0; i < loop; i ++){
    //     // final_exponentiation(&S, f, r, p);
    //     hard_final_exponentiation(&S, f, r, p, z);
    // }
    // Maverage = Mcounter / loop;
    // printf("最終べき掛け算の回数 : %d\n", Maverage);


    optimal_ate_miller(&f, z, P, Q, p);
    printf("\nmiller : ");
    fp12_printf(f);

    // final_exponentiation(&S, f, r, p);
    // easy_final_exponentiation(&S, f, r, p);
    hard_final_exponentiation(&S, f, r, p, z);
    printf("\nPARE : ");
    fp12_printf(S);



    // struct fp12 f;fp12_init(&f);
    optimal_ate_miller(&f, z, P, Q, p);
    printf("\nmiller : ");
    fp12_printf(f);


    Mcounter = 0;
    // struct fp12 S;fp12_init(&S);
    final_exponentiation(&S, f, r, p);
    printf("\nPARE : ");
    fp12_printf(S);


    mpz_t s1, s2;mpz_inits(s1, s2, NULL);
    mpz_set_str(s1, "100", 10);
    mpz_set_str(s2, "200", 10);
    bilinearity(P, Q, s1, s2, z, r, p);

    printf("\ntest\n");
    struct fp12 W;fp12_init(&W);
    fp12_random(&W, p, state);
    fp12_printf(W);
}