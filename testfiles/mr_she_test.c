#include "../miller_header.h"
 
int main (void){

    mpz_t z;mpz_t p;mpz_t r;mpz_t t;
    mpz_inits(z, p, r, t, NULL);
    gen_params(z, p, r, t);

    // gmp_printf("z = %Zd\n", z);
    // gmp_printf("r = %Zd\n", r);

    struct fp12 M;fp12_init(&M);
    mpz_set_str(M.x0.x0.x0.x0, "1", 10);
    mpz_set_str(M.x0.x0.x1.x0, "2", 10);
    mpz_set_str(M.x0.x1.x0.x0, "3", 10);
    mpz_set_str(M.x0.x1.x1.x0, "4", 10);
    mpz_set_str(M.x0.x2.x0.x0, "5", 10);
    mpz_set_str(M.x0.x2.x1.x0, "6", 10);
    mpz_set_str(M.x1.x0.x0.x0, "7", 10);
    mpz_set_str(M.x1.x0.x1.x0, "8", 10);
    mpz_set_str(M.x1.x1.x0.x0, "9", 10);
    mpz_set_str(M.x1.x1.x1.x0, "10", 10);
    mpz_set_str(M.x1.x2.x0.x0, "11", 10);
    mpz_set_str(M.x1.x2.x1.x0, "12", 10);
    printf("Message : ");fp12_printf(M);printf("\n");
    printf("\n");

    mpz_t omega;mpz_init(omega);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    struct timeval tv;
    gettimeofday(&tv, NULL);
    usleep(1000);
    unsigned long seed = tv.tv_sec * 1000 + tv.tv_usec / 1000;
    gmp_randseed_ui(state, seed);
    mpz_urandomm(omega, state, r);
    gmp_printf("Key Word : %Zd\n", omega);
    printf("\n");
    

    // Key Generate
    printf("Public Key = (g, g_1, h_1, h_2, h_3, h_4, hk, f)\n");
    struct PubKey pk;pk_init(&pk);
    pk_init(&pk);
    struct SecKey sk;
    sk_init(&sk);
    KeyGen(&pk, &sk);
    printf("g : ");efp12_printf(pk.g);printf("\n");
    printf("g_1 : ");efp12_printf(pk.g1);printf("\n");
    printf("h_1 : ");efp12_printf(pk.h1);printf("\n");
    printf("h_2 : ");efp12_printf(pk.h2);printf("\n");
    printf("h_3 : ");efp12_printf(pk.h3);printf("\n");
    printf("h_4 : ");efp12_printf(pk.h4);printf("\n");
    printf("hk, f : SHA-512\n");
    printf("\n");

    printf("Secret Key = α\n");
    gmp_printf("α : %Zd\n", sk.alpha);
    printf("\n");


    // Homomorphic Key Generate
    printf("Homomorphic Key = (g^ω, (r_{ω,3}, h_{ω,3}), (r_{ω,4}, h_{ω,4}))\n");
    struct HomKey hk;hk_init(&hk);
    HomKeyGen(&hk, pk, sk, omega, p, r);
    printf("g^ω : ");efp12_printf(hk.g);printf("\n");
    gmp_printf("r_{ω,3} : %Zd\n", hk.rh3.r_w);
    printf("h_{ω,3} : ");efp12_printf(hk.rh3.h_w);printf("\n");
    gmp_printf("r_{ω,4} : %Zd\n", hk.rh4.r_w);
    printf("h_{ω,4} : ");efp12_printf(hk.rh4.h_w);printf("\n");
    printf("\n");

    // Encryption
    printf("Cypher Text = (c_1, c_2, c_3, c_4, τ)\n");
    struct Ciphertext ct;ct_init(&ct);
    Enc(&ct, pk, M, omega, p, r);
    printf("c_1 : ");efp12_printf(ct.c1);printf("\n");
    printf("c_2 : ");fp12_printf(ct.c2);printf("\n");
    printf("c_3 : ");fp12_printf(ct.c3);printf("\n");
    printf("c_4 : ");fp12_printf(ct.c4);printf("\n");
    gmp_printf("τ : %Zd\n", ct.tau);
    printf("\n");

    // Test
    int test = 0;
    test = Test(pk, hk, ct, p);
    printf("Test Result : ");
    if(test == 1){printf("Success\n");}
    else{printf("Failed\n");}
    printf("\n");

    // Decryption
    printf("Decryption\n");
    Dec(&M, pk, sk, omega, ct, p, r);
    printf("Decrypted Message : ");fp12_printf(M);printf("\n");
    printf("\n");

    // Evaluation
    struct fp12 M2;fp12_init(&M2);
    mpz_set_str(M2.x0.x0.x0.x0, "2", 10);
    mpz_set_str(M2.x0.x0.x1.x0, "3", 10);
    mpz_set_str(M2.x0.x1.x0.x0, "4", 10);
    mpz_set_str(M2.x0.x1.x1.x0, "5", 10);
    mpz_set_str(M2.x0.x2.x0.x0, "6", 10);
    mpz_set_str(M2.x0.x2.x1.x0, "7", 10);
    mpz_set_str(M2.x1.x0.x0.x0, "8", 10);
    mpz_set_str(M2.x1.x0.x1.x0, "9", 10);
    mpz_set_str(M2.x1.x1.x0.x0, "10", 10);
    mpz_set_str(M2.x1.x1.x1.x0, "11", 10);
    mpz_set_str(M2.x1.x2.x0.x0, "12", 10);
    mpz_set_str(M2.x1.x2.x1.x0, "1", 10);
    printf("Message2 : ");fp12_printf(M2);printf("\n");
    printf("\n");

    printf("Cypher Text2 = (c2_1, c2_2, c2_3, c2_4, τ2)\n");
    struct Ciphertext ct2;ct_init(&ct2);
    Enc(&ct2, pk, M2, omega, p, r);
    printf("c2_1 : ");efp12_printf(ct2.c1);printf("\n");
    printf("c2_2 : ");fp12_printf(ct2.c2);printf("\n");
    printf("c2_3 : ");fp12_printf(ct2.c3);printf("\n");
    printf("c2_4 : ");fp12_printf(ct2.c4);printf("\n");
    gmp_printf("τ2 : %Zd\n", ct.tau);
    printf("\n");

    struct Ciphertext ct_eval;ct_init(&ct_eval);
    Eval(&ct_eval, pk, hk, ct, ct2, p, r);
    printf("c_eval_1 : ");efp12_printf(ct_eval.c1);printf("\n");
    printf("c_eval_2 : ");fp12_printf(ct_eval.c2);printf("\n");
    printf("c_eval_3 : ");fp12_printf(ct_eval.c3);printf("\n");
    printf("c_eval_4 : ");fp12_printf(ct_eval.c4);printf("\n");
    gmp_printf("τ_eval : %Zd\n", ct_eval.tau);
    printf("\n");

    test = Test(pk, hk, ct_eval, p);
    printf("Test Result : ");
    if(test == 1){printf("Success\n");}
    else{printf("Failed\n");}
    printf("\n");

    struct fp12 M_eval;fp12_init(&M_eval);
    Dec(&M_eval, pk, sk, omega, ct_eval, p, r);
    printf("Decrypted Message : ");fp12_printf(M_eval);printf("\n");

    struct fp12 M3;fp12_init(&M3);
    fp12_mul(&M3, M, M2, p);
    printf("M3 : ");fp12_printf(M3);printf("\n");

    gmp_randclear(state);
}
