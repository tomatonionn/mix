#include "../miller_header.h"

void pk_init(struct PubKey *pk){
    efp12_init(&pk->g);
    efp12_init(&pk->g1);
    efp12_init(&pk->h1);
    efp12_init(&pk->h2);
    efp12_init(&pk->h3);
    efp12_init(&pk->h4);
}

void pk_clear(struct PubKey *pk){
    efp12_clear(&pk->g);
    efp12_clear(&pk->g1);
    efp12_clear(&pk->h1);
    efp12_clear(&pk->h2);
    efp12_clear(&pk->h3);
    efp12_clear(&pk->h4);
}

void sk_init(struct SecKey *sk){
    mpz_init(sk->alpha);
}

void sk_clear(struct SecKey *sk){
    mpz_clear(sk->alpha);
}

void hk_init(struct HomKey *hk){
    efp12_init(&hk->g);
    mpz_init(hk->rh3.r_w);
    efp12_init(&hk->rh3.h_w);
    mpz_init(hk->rh4.r_w);
    efp12_init(&hk->rh4.h_w);
}

void hk_clear(struct HomKey *hk){
    efp12_clear(&hk->g);
    mpz_clear(hk->rh3.r_w);
    efp12_clear(&hk->rh3.h_w);
    mpz_clear(hk->rh4.r_w);
    efp12_clear(&hk->rh4.h_w);
}

void ct_init(struct Ciphertext *ct){
    efp12_init(&ct->c1);
    fp12_init(&ct->c2);
    fp12_init(&ct->c3);
    fp12_init(&ct->c4);
    mpz_init(ct->tau);
}

void ct_clear(struct Ciphertext *ct){
    efp12_clear(&ct->c1);
    fp12_clear(&ct->c2);
    fp12_clear(&ct->c3);
    fp12_clear(&ct->c4);
    mpz_clear(ct->tau);
}

void symmetric_miller(struct fp12 *e, const struct efp12 P, const struct efp12 Q){
    
    mpz_t z;mpz_t p;mpz_t r;mpz_t t;
    mpz_inits(z, p, r, t, NULL);
    gen_params(z, p, r, t);
    
    struct efp12 P_tmp;efp12_init(&P_tmp);efp12_set(&P_tmp, P);
    struct efp12 Q_tmp;efp12_init(&Q_tmp);efp12_set(&Q_tmp, Q);
    struct fp12 e_tmp;fp12_init(&e_tmp);
    mpz_t r_prime;mpz_init(r_prime);

    // 1 : r' ← (p-1)^-1 (mod r)
    mpz_sub_ui(r_prime, p, 1);
    mpz_invert(r_prime, r_prime, r);

    // 2 : P' ← (p-πp)×r' ×P
    struct efp12 pi_p;efp12_init(&pi_p);
    fp12_Frobenius(&pi_p.x, P.x, p);
    fp12_Frobenius(&pi_p.y, P.y, p);
    fp12_neg(&pi_p.y, pi_p.y, p);
    efp12_scm(&P_tmp, P, p, p);
    efp12_eca(&P_tmp, P_tmp, pi_p, p);
    efp12_scm(&P_tmp, P_tmp, r_prime, p);
    // printf("P' : ");efp12_printf(P_tmp);printf("\n");
    
    // 3 : Q' ← (πp −1)×r' ×Q
    fp12_Frobenius(&pi_p.x, Q.x, p);
    fp12_Frobenius(&pi_p.y, Q.y, p);
    fp12_set(&Q_tmp.x, Q.x);
    fp12_neg(&Q_tmp.y, Q.y, p);
    efp12_eca(&Q_tmp, Q_tmp, pi_p, p);
    efp12_scm(&Q_tmp, Q_tmp, r_prime, p);
    // printf("Q' : ");efp12_printf(Q_tmp);printf("\n");

    // 4 : e ←e_asy(Q',P')
    optimal_ate_miller(&e_tmp, z, P_tmp, Q_tmp, p);
    final_exponentiation(&e_tmp, e_tmp, r, p);
    fp12_set(e, e_tmp);
    
    efp12_clear(&P_tmp);
    efp12_clear(&Q_tmp);
    fp12_clear(&e_tmp);
    mpz_clear(r_prime);
}

void G_random(struct efp12 *R){
    struct efp12 R_tmp;efp12_init(&R_tmp);
    struct efp12 P;efp12_init(&P);
    struct efp12 Q;efp12_init(&Q);

    gmp_randstate_t state;
    gmp_randinit_default(state);
    struct timeval tv;
    gettimeofday(&tv, NULL);
    usleep(100);
    unsigned long seed = tv.tv_sec * 1000 + tv.tv_usec / 1000;  // ミリ秒単位のシード
    gmp_randseed_ui(state, seed);

    mpz_t z;mpz_t p;mpz_t r;mpz_t t;
    mpz_inits(z, p, r, t, NULL);
    gen_params(z, p, r, t);

    mpz_t E;mpz_init(E);
    rank_number(E, 12, t, p);

    struct fp b;fp_init(&b);
    mpz_set_str(b.x0, "2806781539090543763928146397551071025921865095800381583843579968964127551432039332258992094003963260740981125881345582810061579481053866112", 10);

    generate1(&P, b, z, r, p, state);
    generate2(&Q, b, E, r, p, state);
    efp12_eca(&R_tmp, P, Q, p);
    efp12_set(R, R_tmp);

    efp12_clear(&R_tmp);
    efp12_clear(&P);
    efp12_clear(&Q);
    fp_clear(&b);
    mpz_clear(z);
    mpz_clear(p);
    mpz_clear(r);
    mpz_clear(t);
}

void KeyGen(struct PubKey *pk, struct SecKey *sk){
    
    // シード設定
    gmp_randstate_t state;
    gmp_randinit_default(state);
    struct timeval tv;
    gettimeofday(&tv, NULL);
    usleep(100);
    unsigned long seed = tv.tv_sec * 1000 + tv.tv_usec / 1000;
    gmp_randseed_ui(state, seed);

    // パラメータ設定
    mpz_t z;mpz_t p;mpz_t r;mpz_t t;
    mpz_inits(z, p, r, t, NULL);
    gen_params(z, p, r, t);
    struct fp b;mpz_init(b.x0);
    mpz_set_str(b.x0, "2806781539090543763928146397551071025921865095800381583843579968964127551432039332258992094003963260740981125881345582810061579481053866112", 10);

    // 公開鍵、秘密鍵生成
    struct PubKey pk_tmp;pk_init(&pk_tmp);
    struct SecKey sk_tmp;sk_init(&sk_tmp);

    // 出力
    G_random(&pk_tmp.g);       // g ← G
    G_random(&pk_tmp.h1);      // h1 ← G
    G_random(&pk_tmp.h2);      // h2 ← G
    G_random(&pk_tmp.h3);      // h3 ← G
    G_random(&pk_tmp.h4);      // h4 ← G
    mpz_urandomm(sk_tmp.alpha, state, r);       // α ← Zr
    efp12_scm(&pk_tmp.g1, pk_tmp.g, sk_tmp.alpha, p);   // g1 ← g^α
    // hk, f ← SHA-512

    // 代入
    efp12_set(&pk->g, pk_tmp.g);
    efp12_set(&pk->g1, pk_tmp.g1);
    efp12_set(&pk->h1, pk_tmp.h1);
    efp12_set(&pk->h2, pk_tmp.h2);
    efp12_set(&pk->h3, pk_tmp.h3);
    efp12_set(&pk->h4, pk_tmp.h4);
    mpz_set(sk->alpha, sk_tmp.alpha);

    // 解放
    pk_clear(&pk_tmp);
    sk_clear(&sk_tmp);
    fp_clear(&b);
    mpz_clears(z, p, r, t, NULL);
    gmp_randclear(state);
}

void HomKeyGen(struct HomKey *hk, struct PubKey pk, struct SecKey sk, mpz_t omega, mpz_t p, mpz_t r){
    
    // シード設定
    gmp_randstate_t state;
    gmp_randinit_default(state);
    struct timeval tv;
    gettimeofday(&tv, NULL);
    usleep(100);
    unsigned long seed = tv.tv_sec * 1000 + tv.tv_usec / 1000;
    gmp_randseed_ui(state, seed);

    mpz_t index;mpz_init(index);
    
    // 演算鍵（トラップドア）生成
    struct HomKey hk_tmp;hk_init(&hk_tmp);

    mpz_urandomm(hk_tmp.rh3.r_w, state, r);          // r_{ω,3} ← Zr
    mpz_urandomm(hk_tmp.rh4.r_w, state, r);          // r_{ω,4} ← Zr

    // h_{ω,3} ← (h_3 h^{-r_{ω,3}})^{1/(α-1))}
    efp12_scm(&hk_tmp.rh3.h_w, pk.g, hk_tmp.rh3.r_w, p);    // g^r_{ω,3}
    fp12_neg(&hk_tmp.rh3.h_w.y, hk_tmp.rh3.h_w.y, p);       // g^-r_{ω,3}
    efp12_eca(&hk_tmp.rh3.h_w, hk_tmp.rh3.h_w, pk.h3, p);   // h_3 + g^-r_{ω,3}
    mpz_sub(index, sk.alpha, omega);                        // α - ω
    mpz_invert(index, index, p);                            // index = 1 / α - ω
    efp12_scm(&hk_tmp.rh3.h_w, hk_tmp.rh3.h_w, index, p);   // h_{ω,3} = (h_3 h^{-r_{ω,3}})^{1/(α-1))}

    // h_{ω,4} ← (h_4 h^{-r_{ω,4}})^{1/(α-1))}
    efp12_scm(&hk_tmp.rh4.h_w, pk.g, hk_tmp.rh4.r_w, p);    // g^r_{ω,4}
    fp12_neg(&hk_tmp.rh4.h_w.y, hk_tmp.rh4.h_w.y, p);       // g^-r_{ω,4}
    efp12_eca(&hk_tmp.rh4.h_w, hk_tmp.rh4.h_w, pk.h4, p);   // h_4 + g^index
    mpz_sub(index, sk.alpha, omega);                        // α - ω
    mpz_invert(index, index, p);                            // index = 1 / α - ω
    efp12_scm(&hk_tmp.rh4.h_w, hk_tmp.rh4.h_w, index, p);   // h_{ω,4} = (h_4 + h^{-r_{ω,4}})^{1/(α-1))}

    efp12_scm(&hk_tmp.g, pk.g, omega, p);          // g^{ω}

    // 代入
    efp12_set(&hk->g, hk_tmp.g);
    mpz_set(hk->rh3.r_w, hk_tmp.rh3.r_w);
    efp12_set(&hk->rh3.h_w, hk_tmp.rh3.h_w);
    mpz_set(hk->rh4.r_w, hk_tmp.rh4.r_w);
    efp12_set(&hk->rh4.h_w, hk_tmp.rh4.h_w);

    // 解放
    hk_clear(&hk_tmp);
    mpz_clear(index);
    gmp_randclear(state);
}

void Enc(struct Ciphertext *ct, struct PubKey pk, struct fp12 M, mpz_t omega, mpz_t p, mpz_t r){

    // シード設定
    gmp_randstate_t state;
    gmp_randinit_default(state);
    struct timeval tv;
    gettimeofday(&tv, NULL);
    usleep(100);
    unsigned long seed = tv.tv_sec * 1000 + tv.tv_usec / 1000;  // ミリ秒単位のシード
    gmp_randseed_ui(state, seed);

    mpz_t index;mpz_init(index);
    struct efp12 efp12_tmp;efp12_init(&efp12_tmp);
    struct fp12 fp12_tmp;fp12_init(&fp12_tmp);
    
    // 暗号文生成
    struct Ciphertext ct_tmp;ct_init(&ct_tmp);

    // s ← Zr
    mpz_t s;mpz_init(s);
    mpz_urandomm(s, state, r);

    // c1 ← g1^s g^{-sω}
    efp12_scm(&ct_tmp.c1, pk.g1, s, p);                 // g1^s
    mpz_mul(index, s, omega);                           // sω
    efp12_scm(&efp12_tmp, pk.g, index, p);              // g^{sω}
    fp12_neg(&efp12_tmp.y, efp12_tmp.y, p);             // g^{-sω}
    efp12_eca(&ct_tmp.c1, ct_tmp.c1, efp12_tmp, p);     // c1 = g1^s + g^{-sω}

    // c2 ← e(g,g)^s
    symmetric_miller(&ct_tmp.c2, pk.g, pk.g);           // e(g,g)
    fp12_pow(&ct_tmp.c2, ct_tmp.c2, s, p);              // c2 = e(g,g)^s

    // c3 ← M・e(g,h1)^-s
    symmetric_miller(&ct_tmp.c3, pk.g, pk.h1);          // e(g,h1)
    fp12_pow(&fp12_tmp, ct_tmp.c3, s, p);               // e(g,h1)^s
    fp12_inv(&fp12_tmp, fp12_tmp, p);                   // e(g,h1)^-s
    fp12_mul(&ct_tmp.c3, ct_tmp.c3, M, p);              // c3 = M・e(g,h1)^-s

    // c4 ← e(g,h2)^s
    symmetric_miller(&ct_tmp.c4, pk.g, pk.h2);          // e(g,h2)
    fp12_pow(&ct_tmp.c4, ct_tmp.c4, s, p);              // c4 = e(g,h2)^s

    mpz_t delta;mpz_init(delta);
    Gamma(delta, ct_tmp.c1, ct_tmp.c2, ct_tmp.c3, ct_tmp.c4, p);        // δ ← Γ(c1,c2,c3,c4)
    // gmp_printf("delta : %Zd\n", delta);

    // c5 ← e(g,h3)^s e(g,h4)^sδ
    struct fp12 c5;fp12_init(&c5);
    symmetric_miller(&fp12_tmp, pk.g, pk.h3);           // e(g,h3)
    fp12_pow(&fp12_tmp, fp12_tmp, s, p);                // e(g,h3)^s
    symmetric_miller(&c5, pk.g, pk.h4);                 // e(g,h4)
    mpz_mul(index, s, delta);                           // sδ
    fp12_pow(&c5, c5, index, p);                        // e(g,h4)^sδ
    fp12_mul(&c5, fp12_tmp, c5, p);                     // c5 = e(g,h3)^s e(g,h4)^sδ

    Function(ct_tmp.tau, c5);                           // τ ← F(c5)

    // 代入
    efp12_set(&ct->c1, ct_tmp.c1);
    fp12_set(&ct->c2, ct_tmp.c2);
    fp12_set(&ct->c3, ct_tmp.c3);
    fp12_set(&ct->c4, ct_tmp.c4);
    mpz_set(ct->tau, ct_tmp.tau);

    // 解放
    mpz_clear(s);
    mpz_clear(index);
    efp12_clear(&efp12_tmp);
    fp12_clear(&fp12_tmp);
    mpz_clear(delta);
    fp12_clear(&c5);
    ct_clear(&ct_tmp);
    gmp_randclear(state);
}

int Test(struct PubKey pk, struct HomKey hk, struct Ciphertext ct, mpz_t p){
    
    mpz_t index;mpz_init(index);
    struct efp12 efp12_tmp;efp12_init(&efp12_tmp);
    struct fp12 fp12_tmp1;fp12_init(&fp12_tmp1);
    struct fp12 fp12_tmp2;fp12_init(&fp12_tmp2);
    
    // 暗号文検証
    int test = 0;


    // δ ← Γ(c1,c2,c3,c4)
    mpz_t delta;mpz_init(delta);
    Gamma(delta, ct.c1, ct.c2, ct.c3, ct.c4, p);

    // τ_ch ← f(e(c1,h_{ω,3) h_{ω,4}^δ) c2^{r_{ω,3}+r_{ω,4}δ})
    mpz_t tau_ch;mpz_init(tau_ch);
    efp12_scm(&efp12_tmp, hk.rh4.h_w, delta, p);            // h_{ω,4}^δ
    efp12_eca(&efp12_tmp, hk.rh3.h_w, efp12_tmp, p);        // h_{ω,3) h_{ω,4}^δ
    symmetric_miller(&fp12_tmp1, ct.c1, efp12_tmp);         // e(c1,h_{ω,3) h_{ω,4}^δ)
    mpz_mul(index, hk.rh4.r_w, delta);                      // r_{ω,4}δ
    mpz_add(index, index, hk.rh3.r_w);                      // r_{ω,3}+r_{ω,4}δ
    fp12_pow(&fp12_tmp2, ct.c2, index, p);                  // c2^{r_{ω,3}+r_{ω,4}δ}
    fp12_mul(&fp12_tmp1, fp12_tmp1, fp12_tmp2, p);          // e(c1,h_{ω,3) h_{ω,4}^δ) c2^{r_{ω,3}+r_{ω,4}δ}
    Function(tau_ch, fp12_tmp1);                            // τ_ch = f(e(c1,h_{ω,3) h_{ω,4}^δ) c2^{r_{ω,3}+r_{ω,4}δ})

    if(mpz_cmp(ct.tau, tau_ch) == 0){
        test = 1;
    }

    // 解放
    mpz_clear(index);
    efp12_clear(&efp12_tmp);
    fp12_clear(&fp12_tmp1);
    fp12_clear(&fp12_tmp2);
    mpz_clear(delta);
    mpz_clear(tau_ch);

    return test;
}

void Dec(struct fp12 *M, struct PubKey pk, struct SecKey sk, mpz_t omega, struct Ciphertext ct, mpz_t p, mpz_t r){
        
    gmp_randstate_t state;
    gmp_randinit_default(state);
    struct timeval tv;
    gettimeofday(&tv, NULL);
    usleep(100);
    unsigned long seed = tv.tv_sec * 1000 + tv.tv_usec / 1000;  // ミリ秒単位のシード
    gmp_randseed_ui(state, seed);

    mpz_t index;mpz_init(index);
    struct fp12 fp12_tmp;fp12_init(&fp12_tmp);
    struct efp12 efp12_tmp;efp12_init(&efp12_tmp);

    // 復号
    struct fp12 M_tmp;fp12_init(&M_tmp);

    mpz_t r_w1;mpz_init(r_w1);
    mpz_urandomm(r_w1, state, r);                          // r_{ω,1} ← Zr
    mpz_t r_w2;mpz_init(r_w2);
    mpz_urandomm(r_w2, state, r);                          // r_{ω,2} ← Zr
    mpz_t r_w3;mpz_init(r_w3);
    mpz_urandomm(r_w3, state, r);                          // r_{ω,3} ← Zr
    mpz_t r_w4;mpz_init(r_w4);
    mpz_urandomm(r_w4, state, r);                          // r_{ω,4} ← Zr


    // h_{ω,1} ← (h1 g^{-r_{ω,1}})^{1/(α-ω)}
    struct efp12 h_w1;efp12_init(&h_w1);
    mpz_neg(index, r_w1);                            // -r_{ω,1}
    mpz_mod(index, index, p);                            // index = -r_{ω,1} mod p
    efp12_scm(&h_w1, pk.g, index, p);                // g^{-r_{ω,1}}
    efp12_eca(&h_w1, pk.h1, h_w1, p);            // h1 g^{-r_{ω,1}}
    mpz_sub(index, sk.alpha, omega);                     // α - ω
    mpz_invert(index, index, p);                         // index = 1 / α - ω
    efp12_scm(&h_w1, h_w1, index, p);            // h_{ω,1} = (h1 g^{-r_{ω,1}})^{1/(α-ω)}

    // h_{ω,2} ← (h2 g^{-r_{ω,2}})^{1/(α-ω)}
    struct efp12 h_w2;efp12_init(&h_w2);
    mpz_neg(index, r_w2);                            // -r_{ω,2}
    mpz_mod(index, index, p);                            // index = -r_{ω,2} mod p
    efp12_scm(&h_w2, pk.g, index, p);                // g^{-r_{ω,2}}
    efp12_eca(&h_w2, pk.h2, h_w2, p);            // h2 g^{-r_{ω,2}}
    mpz_sub(index, sk.alpha, omega);                     // α - ω
    mpz_invert(index, index, p);                         // index = 1 / α - ω
    efp12_scm(&h_w2, h_w2, index, p);            // h_{ω,2} = (h2 g^{-r_{ω,2}})^{1/(α-ω)}

    // h_{ω,3} ← (h3 g^{-r_{ω,3}})^{1/(α-ω)}
    struct efp12 h_w3;efp12_init(&h_w3);
    mpz_neg(index, r_w3);                            // -r_{ω,3}
    mpz_mod(index, index, p);                            // index = -r_{ω,3} mod p
    efp12_scm(&h_w3, pk.g, index, p);                // g^{-r_{ω,3}}
    efp12_eca(&h_w3, pk.h3, h_w3, p);            // h3 g^{-r_{ω,3}}
    mpz_sub(index, sk.alpha, omega);                     // α - ω
    mpz_invert(index, index, p);                         // index = 1 / α - ω
    efp12_scm(&h_w3, h_w3, index, p);            // h_{ω,3} = (h3 g^{-r_{ω,3}})^{1/(α-ω)}

    // h_{ω,4} ← (h4 g^{-r_{ω,4}})^{1/(α-ω)}
    struct efp12 h_w4;efp12_init(&h_w4);
    mpz_neg(index, r_w4);                            // -r_{ω,4}
    mpz_mod(index, index, p);                            // index = -r_{ω,4} mod p
    efp12_scm(&h_w4, pk.g, index, p);                // g^{-r_{ω,4}}
    efp12_eca(&h_w4, pk.h4, h_w4, p);            // h4 g^{-r_{ω,4}}
    mpz_sub(index, sk.alpha, omega);                     // α - ω
    mpz_invert(index, index, p);                         // index = 1 / α - ω
    efp12_scm(&h_w4, h_w4, index, p);            // h_{ω,4} = (h4 g^{-r_{ω,4}})^{1/(α-ω)}

    mpz_t delta;mpz_init(delta);
    Gamma(delta, ct.c1, ct.c2, ct.c3, ct.c4, p);        // δ ← Γ(c1,c2,c3,c4)

    // c4' ← e(c1,h_{ω,2}) c2^{r_{ω,2}}
    struct fp12 c4_prime;fp12_init(&c4_prime);
    symmetric_miller(&c4_prime, ct.c1, h_w2);    // e(c1,h_{ω,2})
    fp12_pow(&fp12_tmp, ct.c2, r_w2, p);                  // c2^{r_{ω,2}}
    fp12_mul(&c4_prime, c4_prime, fp12_tmp, p);               // c4' = e(c1,h_{ω,2}) c2^{r_{ω,2}}

    // c5 ← e(c1,h_{ω,3} h_{ω,4}^δ) c2^{r_{ω,3}+r_{ω,4}δ}
    struct fp12 c5;fp12_init(&c5);
    efp12_scm(&efp12_tmp, h_w4, delta, p);                  // h_{ω,4}^δ
    efp12_eca(&efp12_tmp, h_w3, efp12_tmp, p);              // h_{ω,3) h_{ω,4}^δ
    symmetric_miller(&fp12_tmp, ct.c1, efp12_tmp);     // e(c1,h_{ω,3) h_{ω,4}^δ)
    mpz_mul(index, r_w4, delta);                            // r_{ω,4}δ
    mpz_add(index, index, r_w3);                            // r_{ω,3}+r_{ω,4}δ
    fp12_pow(&c5, ct.c2, index, p);                             // c2^{r_{ω,3}+r_{ω,4}δ}
    fp12_mul(&c5, fp12_tmp, c5, p);                             // c5' = e(c1,h_{ω,3} h_{ω,4}^δ) c2^{r_{ω,3}+r_{ω,4}δ}

    mpz_t tau_ch;mpz_init(tau_ch);
    Function(tau_ch, c5);                                       // τ_ch ← F(c5)
    if(fp12_cmp(ct.c4, c4_prime) != 0 || mpz_cmp(ct.tau, tau_ch) != 0){
        fp12_set(M, M_tmp);
        printf("Decryption failed\n");
    }

    else{
        // M ← c3 e(c1,h_{ω,1}) c2^{r_{ω,1}}
        fp12_pow(&fp12_tmp, ct.c2, r_w1, p);                // c2^{r_{ω,1}}
        symmetric_miller(&M_tmp, ct.c1, h_w1);         // e(c1,h_{ω,1})
        fp12_mul(&M_tmp, ct.c3, M_tmp, p);                              // c3 e(c1,h_{ω,1})
        fp12_mul(&M_tmp, M_tmp, fp12_tmp, p);                           // M = c3 e(c1,h_{ω,1}) c2^{r_{ω,1}}
        fp12_set(M, M_tmp);
    }

    mpz_clear(index);
    efp12_clear(&efp12_tmp);
    fp12_clear(&fp12_tmp);
    mpz_clear(r_w1);
    mpz_clear(r_w2);
    mpz_clear(r_w3);
    mpz_clear(r_w4);
    mpz_clear(delta);
    fp12_clear(&c4_prime);
    fp12_clear(&c5);
    mpz_clear(tau_ch);

}

void Eval(struct Ciphertext ct, struct PubKey pk, struct HomKey hk, struct Ciphertext ct1, struct Ciphertext ct2, mpz_t p, mpz_t r){
    
    mpz_t index;mpz_init(index);
    struct fp12 fp12_tmp;fp12_init(&fp12_tmp);
    struct efp12 efp12_tmp;efp12_init(&efp12_tmp);
    
    // δ_1 = Γ(c_{1,1}, c_{1,2}, c_{1,3}, c_{1,4})
    mpz_t delta1;mpz_init(delta1);
    Gamma(delta1, ct1.c1, ct1.c2, ct1.c3, ct1.c4, p);

    // c_{1,5} ← e(c_{1,1},h_{ω,3} h_{ω,4}^δ_1) c_{1,2}^{r_{ω,3}+r_{ω,4}δ_1}
    struct fp12 c1_5;fp12_init(&c1_5);
    efp12_scm(&efp12_tmp, hk.rh4.h_w, delta1, p);       // h_{ω,4}^δ_1
    efp12_eca(&efp12_tmp, efp12_tmp, hk.rh3.h_w, p);    // h_{ω,3} + h_{ω,4}^δ_1
    symmetric_miller(&c1_5, ct1.c1, efp12_tmp);             // e(c_{1,1},h_{ω,3} h_{ω,4}^δ_1)
    mpz_mul(index, hk.rh4.r_w, delta1);                 // r_{ω,4}δ_1
    mpz_add(index, index, hk.rh3.r_w);                  // r_{ω,3}+r_{ω,4}δ_1
    fp12_pow(&fp12_tmp, ct1.c2, index, p);                  // c_{1,2}^{r_{ω,3}+r_{ω,4}δ_1}
    fp12_mul(&c1_5, c1_5, fp12_tmp, p);                     // c_{1,5} = e(c_{1,1},h_{ω,3} h_{ω,4}^δ_1) c_{1,2}^{r_{ω,3}+r_{ω,4}δ_1}

    // δ_2 = Γ(c_{2,1}, c_{2,2}, c_{2,3}, c_{2,4})
    mpz_t delta2;mpz_init(delta2);
    Gamma(delta2, ct2.c1, ct2.c2, ct2.c3, ct2.c4, p);

    // c_{2,5} ← e(c_{2,1},h_{ω,3} h_{ω,4}^δ_2) c_{2,2}^{r_{ω,3}+r_{ω,4}δ_2}
    struct fp12 c2_5;fp12_init(&c2_5);
    efp12_scm(&efp12_tmp, hk.rh4.h_w, delta2, p);       // h_{ω,4}^δ_2
    efp12_eca(&efp12_tmp, efp12_tmp, hk.rh3.h_w, p);    // h_{ω,3} + h_{ω,4}^δ_2
    symmetric_miller(&c2_5, ct2.c1, efp12_tmp);             // e(c_{2,1},h_{ω,3} h_{ω,4}^δ_2)
    mpz_mul(index, hk.rh4.r_w, delta2);                 // r_{ω,4}δ_2
    mpz_add(index, index, hk.rh3.r_w);                  // r_{ω,3}+r_{ω,4}δ_2
    fp12_pow(&fp12_tmp, ct2.c2, index, p);                  // c_{2,2}^{r_{ω,3}+r_{ω,4}δ_2}
    fp12_mul(&c2_5, c2_5, fp12_tmp, p);                     // c_{2,5} = e(c_{2,1},h_{ω,3} h_{ω,4}^δ_2) c_{2,2}^{r_{ω,3}+r_{ω,4}δ_2}


    mpz_t tau_ch1;mpz_init(tau_ch1);
    Function(tau_ch1, c1_5);                                // τ_ch1 ← F(c_{1,5})
    mpz_t tau_ch2;mpz_init(tau_ch2);
    Function(tau_ch2, c2_5);                                // τ_ch2 ← F(c_{2,5})
    if(mpz_cmp(tau_ch1, ct1.tau) != 0 || mpz_cmp(tau_ch2, ct2.tau) != 0){
        printf("Evaluation failed\n");
    }

    else{

        gmp_randstate_t state;
        gmp_randinit_default(state);
        struct timeval tv;
        gettimeofday(&tv, NULL);
        usleep(100);
        unsigned long seed = tv.tv_sec * 1000 + tv.tv_usec / 1000;  // ミリ秒単位のシード
        gmp_randseed_ui(state, seed);

        // s ← Zr
        mpz_t s;mpz_init(s);
        mpz_urandomm(s, state, r);

        // c1 ← c_{1,1} c_{2,1} g_1^s g^{-sω}
        mpz_sub_ui(index, p, 1);                // -1
        efp12_scm(&ct.c1, pk.g1, index, p);     // g_1^{-1}
        efp12_eca(&ct.c1, ct1.c1, hk.g, p);     // g_1^{-1} g^ω
        mpz_neg(index, s);                      // -s
        mpz_mod(index, index, p);               // -s mod p
        efp12_scm(&ct.c1, ct.c1, index, p);     // g_1^s g^{-sω}
        efp12_eca(&ct.c1, ct.c1, ct1.c1, p);    // c_{1,1} g_1^s g^{-sω}
        efp12_eca(&ct.c1, ct.c1, ct2.c1, p);    // c_{1,1} c_{2,1} g_1^s g^{-sω}

        // c2 ← c_{1,2} c_{2,2} e(g,g)^s
        symmetric_miller(&ct.c2, pk.g, pk.g);    // e(g,g)
        fp12_pow(&ct.c2, ct.c2, s, p);            // e(g,g)^s
        fp12_mul(&ct.c2, ct.c2, ct1.c2, p);       // c_{1,2} e(g,g)^s
        fp12_mul(&ct.c2, ct.c2, ct2.c2, p);       // c_{1,2} c_{2,2} e(g,g)^s

        // c3 ← c_{1,3} c_{2,3} e(g,h_1)^{-s}
        symmetric_miller(&ct.c3, pk.g, pk.h1);    // e(g,h_1)
        mpz_neg(index, s);                        // -s
        mpz_mod(index, index, p);                 // -s mod p
        fp12_pow(&ct.c3, ct.c3, index, p);        // e(g,h_1)^{-s}
        fp12_mul(&ct.c3, ct.c3, ct1.c3, p);       // c_{1,3} e(g,h_1)^{-s}
        fp12_mul(&ct.c3, ct.c3, ct2.c3, p);       // c_{1,3} c_{2,3} e(g,h_1)^{-s}

        // c4 ← c_{1,4} c_{2,4} e(g,h_2)^s
        symmetric_miller(&ct.c4, pk.g, pk.h2);    // e(g,h_2)
        fp12_pow(&ct.c4, ct.c4, s, p);            // e(g,h_2)^s
        fp12_mul(&ct.c4, ct.c4, ct1.c4, p);       // c_{1,4} e(g,h_2)^s
        fp12_mul(&ct.c4, ct.c4, ct2.c4, p);       // c_{1,4} c_{2,4} e(g,h_2)^s

        // δ ← Γ(c1,c2,c3,c4)
        mpz_t delta;mpz_init(delta);
        Gamma(delta, ct.c1, ct.c2, ct.c3, ct.c4, p);
        
        // c5 ← e(c1,h_{ω,3} h_{ω,4}^δ) c2^{r_{ω,3}+r_{ω,4}δ}
        struct fp12 c5;fp12_init(&c5);
        efp12_scm(&efp12_tmp, hk.rh4.h_w, delta, p);        // h_{ω,4}^δ
        efp12_eca(&efp12_tmp, hk.rh3.h_w, efp12_tmp, p);    // h_{ω,3) h_{ω,4}^δ
        symmetric_miller(&fp12_tmp, ct.c1, efp12_tmp);          // e(c1,h_{ω,3} h_{ω,4}^δ)
        mpz_mul(index, hk.rh4.r_w, delta);                  // r_{ω,4}δ
        mpz_add(index, index, hk.rh3.r_w);                  // r_{ω,3}+r_{ω,4}δ
        fp12_pow(&c5, ct.c2, index, p);                         // c2^{r_{ω,3}+r_{ω,4}δ}
        fp12_mul(&c5, fp12_tmp, c5, p);                         // c5 = e(c1,h_{ω,3} h_{ω,4}^δ) c2^{r_{ω,3}+r_{ω,4}δ}

        Function(ct.tau, c5);                                   // τ ← F(c5)
    }

}
