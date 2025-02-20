#include <ELiPS/bls12.h>
#include <stdio.h>
#include <openssl/evp.h>
#include <openssl/sha.h>

typedef struct{
    efp12_t g;
    efp12_t g1;
    efp12_t h1;
    efp12_t h2;
    efp12_t h3;
    efp12_t h4;
    // hk in HK
    // f
}PubKey;

typedef struct{
    mpz_t alpha;
}SecKey;

typedef struct{
    mpz_t r_w;
    efp12_t h_w;
}Ran_Hash;

typedef struct{
    efp12_t g;
    Ran_Hash rh3;
    Ran_Hash rh4;
}HomKey;

typedef struct{
    efp12_t c1;
    fp12_t c2;
    fp12_t c3;
    fp12_t c4;
    mpz_t tau;
}Ciphertext;

typedef struct{
    fp12_t gg;
    fp12_t gh1;
    fp12_t gh2;
    fp12_t gh3;
    fp12_t gh4;
}PreValue;

void Gamma(mpz_t gamma, efp12_t c1, fp12_t c2, fp12_t c3, fp12_t c4){

    // (c1, c2, c3, c4) → fp12_sum → c_sum
    mp_limb_t limb_sum;
    fp12_t fp12_sum;fp12_init(&fp12_sum);

    // c1.x +c1.y + c2 + c3 + c4 → fp12_sum
    fp12_add(&fp12_sum, &c1.x, &c1.y);
    fp12_add(&fp12_sum, &fp12_sum, &c2);
    fp12_add(&fp12_sum, &fp12_sum, &c3);
    fp12_add(&fp12_sum, &fp12_sum, &c4);

    // fp12_sum.x0.x0.x0.x0 + fp12_sum.x0.x0.x1.x0 + ... → mpz_t c_sum
    mpn_add_n(&limb_sum, &fp12_sum.x0.x0.x0.x0[0], &fp12_sum.x0.x0.x1.x0[0], FPLIMB);
    mpn_add_n(&limb_sum, &limb_sum, &fp12_sum.x0.x1.x0.x0[0], FPLIMB);
    mpn_add_n(&limb_sum, &limb_sum, &fp12_sum.x0.x1.x1.x0[0], FPLIMB);
    mpn_add_n(&limb_sum, &limb_sum, &fp12_sum.x0.x2.x0.x0[0], FPLIMB);
    mpn_add_n(&limb_sum, &limb_sum, &fp12_sum.x0.x2.x1.x0[0], FPLIMB);
    mpn_add_n(&limb_sum, &limb_sum, &fp12_sum.x1.x0.x0.x0[0], FPLIMB);
    mpn_add_n(&limb_sum, &limb_sum, &fp12_sum.x1.x0.x1.x0[0], FPLIMB);
    mpn_add_n(&limb_sum, &limb_sum, &fp12_sum.x1.x1.x0.x0[0], FPLIMB);
    mpn_add_n(&limb_sum, &limb_sum, &fp12_sum.x1.x1.x1.x0[0], FPLIMB);
    mpn_add_n(&limb_sum, &limb_sum, &fp12_sum.x1.x2.x0.x0[0], FPLIMB);
    mpn_add_n(&limb_sum, &limb_sum, &fp12_sum.x1.x2.x1.x0[0], FPLIMB);


    // HASH(c_sum) → gamma
    // mpz_t型をバイト列に変換
    mpz_t mpz_sum; mpz_init(mpz_sum);
    mpz_set_ui(mpz_sum, limb_sum);
    size_t c_sum_len = (mpz_sizeinbase(mpz_sum, 2) + 7) / 8; // バイト数を計算
    unsigned char *c_sum_bytes = (unsigned char *)malloc(c_sum_len);
    mpz_export(c_sum_bytes, NULL, 1, 1, 0, 0, mpz_sum);

    // ハッシュ値を生成
    unsigned char char_gamma[EVP_MAX_MD_SIZE]; // 最大のハッシュサイズ
    EVP_MD_CTX *mdctx;
    const EVP_MD *md;
    unsigned int md_len;
    md = EVP_sha512();
    mdctx = EVP_MD_CTX_new();
    EVP_DigestInit_ex(mdctx, md, NULL);
    EVP_DigestUpdate(mdctx, c_sum_bytes, c_sum_len);
    EVP_DigestFinal_ex(mdctx, char_gamma, &md_len);
    EVP_MD_CTX_free(mdctx);

    // ハッシュ値をmpz_t型に変換
    mpz_import(gamma, SHA512_DIGEST_LENGTH, 1, 1, 0, 0, char_gamma);

    mpz_clear(mpz_sum);
}

void Function(mpz_t function, fp12_t c5){

    // c5.x0.x0.x0.x0 + c5.x0.x0.x0.x0 + ... → mpz_t limb_sum
    mp_limb_t limb_sum;
    mpn_add_n(&limb_sum, &c5.x0.x0.x0.x0[0], &c5.x0.x0.x1.x0[0], FPLIMB);
    mpn_add_n(&limb_sum, &limb_sum, &c5.x0.x1.x0.x0[0], FPLIMB);
    mpn_add_n(&limb_sum, &limb_sum, &c5.x0.x1.x1.x0[0], FPLIMB);
    mpn_add_n(&limb_sum, &limb_sum, &c5.x0.x2.x0.x0[0], FPLIMB);
    mpn_add_n(&limb_sum, &limb_sum, &c5.x0.x2.x1.x0[0], FPLIMB);
    mpn_add_n(&limb_sum, &limb_sum, &c5.x1.x0.x0.x0[0], FPLIMB);
    mpn_add_n(&limb_sum, &limb_sum, &c5.x1.x0.x1.x0[0], FPLIMB);
    mpn_add_n(&limb_sum, &limb_sum, &c5.x1.x1.x0.x0[0], FPLIMB);
    mpn_add_n(&limb_sum, &limb_sum, &c5.x1.x1.x1.x0[0], FPLIMB);
    mpn_add_n(&limb_sum, &limb_sum, &c5.x1.x2.x0.x0[0], FPLIMB);
    mpn_add_n(&limb_sum, &limb_sum, &c5.x1.x2.x1.x0[0], FPLIMB);


    // HASH(c5_sum) → function
    // mpz_t型をバイト列に変換
    mpz_t mpz_sum; mpz_init(mpz_sum);
    mpz_set_ui(mpz_sum, limb_sum);
    size_t c5_sum_len = (mpz_sizeinbase(mpz_sum, 2) + 7) / 8; // バイト数を計算
    unsigned char *c5_sum_bytes = (unsigned char *)malloc(c5_sum_len);
    mpz_export(c5_sum_bytes, NULL, 1, 1, 0, 0, mpz_sum);

    // ハッシュ値を生成
    unsigned char char_function[EVP_MAX_MD_SIZE]; // 最大のハッシュサイズ
    EVP_MD_CTX *mdctx;
    const EVP_MD *md;
    unsigned int md_len;
    md = EVP_sha512();
    mdctx = EVP_MD_CTX_new();
    EVP_DigestInit_ex(mdctx, md, NULL);
    EVP_DigestUpdate(mdctx, c5_sum_bytes, c5_sum_len);
    EVP_DigestFinal_ex(mdctx, char_function, &md_len);
    EVP_MD_CTX_free(mdctx);

    // ハッシュ値をmpz_t型に変換
    mpz_import(function, SHA512_DIGEST_LENGTH, 1, 1, 0, 0, char_function);

    mpz_clear(mpz_sum);
}

void pk_init(PubKey *pk){
    efp12_init(&pk->g);
    efp12_init(&pk->g1);
    efp12_init(&pk->h1);
    efp12_init(&pk->h2);
    efp12_init(&pk->h3);
    efp12_init(&pk->h4);
}

void sk_init(SecKey *sk){
    mpz_init(sk->alpha);
}

void hk_init(HomKey *hk){
    efp12_init(&hk->g);
    mpz_init(hk->rh3.r_w);
    efp12_init(&hk->rh3.h_w);
    mpz_init(hk->rh4.r_w);
    efp12_init(&hk->rh4.h_w);

}

void ct_init(Ciphertext *ct){
    efp12_init(&ct->c1);
    fp12_init(&ct->c2);
    fp12_init(&ct->c3);
    fp12_init(&ct->c4);
    mpz_init(ct->tau);
}

void pv_init(PreValue *pv){
    fp12_init(&pv->gg);
    fp12_init(&pv->gh1);
    fp12_init(&pv->gh2);
    fp12_init(&pv->gh3);
    fp12_init(&pv->gh4);
}

// モンゴメリ記法
void symmetric_paring(fp12_t *e, efp12_t P, efp12_t Q){
    efp12_t P_tmp;efp12_init(&P_tmp);
    efp12_t Q_tmp;efp12_init(&Q_tmp);
    fp12_t e_tmp;fp12_init(&e_tmp);
    mpz_t r_prime;mpz_init(r_prime);
 
    // 1 : r' ← (p-1)^-1 (mod r)
    mpz_sub_ui(r_prime, prime_z, 1);
    mpz_invert(r_prime, r_prime, order_z);
 
    // 2 : P' ← (p-πp)×r' ×P
    efp12_t pi_p;efp12_init(&pi_p);
    fp12_frobenius_map_p1(&pi_p.x, &P.x);
    fp12_frobenius_map_p1(&pi_p.y, &P.y);
    fp12_set_neg(&pi_p.y, &pi_p.y);
    efp12_scm(&P_tmp, &P, prime_z);
    efp12_eca(&P_tmp, &P_tmp, &pi_p);
    efp12_scm(&P_tmp, &P_tmp, r_prime);
    // efp12_println("P' = ", &P_tmp);
   
    // 3 : Q' ← (πp −1)×r' ×Q
    fp12_frobenius_map_p1(&pi_p.x, &Q.x);
    fp12_frobenius_map_p1(&pi_p.y, &Q.y);
    fp12_set(&Q_tmp.x, &Q.x);
    fp12_set_neg(&Q_tmp.y, &Q.y);
    efp12_eca(&Q_tmp, &Q_tmp, &pi_p);
    efp12_scm(&Q_tmp, &Q_tmp, r_prime);
    // efp12_println("Q' = ", &Q_tmp);
 
    efp_t P_twist;efp_init(&P_twist);
    efp2_t Q_twist;efp2_init(&Q_twist);
    efp12_to_efp(&P_twist, &P_tmp);
    fp_to_montgomery(&P_twist.x, &P_twist.x);
    fp_to_montgomery(&P_twist.y, &P_twist.y);
    efp12_to_efp2(&Q_twist, &Q_tmp);
    fp2_to_montgomery(&Q_twist.x, &Q_twist.x);
    fp2_to_montgomery(&Q_twist.y, &Q_twist.y);
 
    // 4 : e ←e_asy(Q',P')
    g1g2_to_g3_pairing(&e_tmp, &P_twist, &Q_twist);
    fp12_set(e, &e_tmp);
 
    mpz_clear(r_prime);
}

// 整数表現
void G_random(efp12_t *R, efp_t gen1, efp2_t gen2){
    fr_t s1;fr_init(&s1);
    fr_t s2;fr_init(&s2);
    fr_set_random(&s1, state);
    fr_set_random(&s2, state);

    efp_t R1;efp_init(&R1);
    g1_scm(&R1, &gen1, &s1);
    efp12_t R1_12;efp12_init(&R1_12);
    efp_to_efp12(&R1_12, &R1);
    fp12_mod_montgomery(&R1_12.x, &R1_12.x);
    fp12_mod_montgomery(&R1_12.y, &R1_12.y);

    efp2_t R2;efp2_init(&R2);
    g2_scm(&R2, &gen2, &s2);
    efp12_t R2_12;efp12_init(&R2_12);
    efp2_to_efp12(&R2_12, &R2);
    fp12_mod_montgomery(&R2_12.x, &R2_12.x);
    fp12_mod_montgomery(&R2_12.y, &R2_12.y);

    efp12_eca(R, &R1_12, &R2_12);
}

void KeyGen(PubKey *pk, SecKey *sk, efp_t gen1, efp2_t gen2){
    // 公開鍵、秘密鍵生成
    PubKey pk_tmp;pk_init(&pk_tmp);
    SecKey sk_tmp;sk_init(&sk_tmp);


    // g ← G
    fr_t s1;fr_init(&s1);
    fr_t s2;fr_init(&s2);
    fr_set_random(&s1, state);
    fr_set_random(&s2, state);
    efp_t R1;efp_init(&R1);
    g1_scm(&R1, &gen1, &s1);
    efp12_t R1_12;efp12_init(&R1_12);
    efp_to_efp12(&R1_12, &R1);
    fp12_mod_montgomery(&R1_12.x, &R1_12.x);
    fp12_mod_montgomery(&R1_12.y, &R1_12.y);
    efp2_t R2;efp2_init(&R2);
    g2_scm(&R2, &gen2, &s2);
    efp12_t R2_12;efp12_init(&R2_12);
    efp2_to_efp12(&R2_12, &R2);
    fp12_mod_montgomery(&R2_12.x, &R2_12.x);
    fp12_mod_montgomery(&R2_12.y, &R2_12.y);
    efp12_eca(&pk_tmp.g, &R1_12, &R2_12);


    fr_t alpha;fr_init(&alpha);
    fr_set_random(&alpha, state);
    fr_mul(&s1, &alpha, &s1);
    fr_mul(&s2, &alpha, &s2);
    g1_scm(&R1, &gen1, &s1);
    efp_to_efp12(&R1_12, &R1);
    fp12_mod_montgomery(&R1_12.x, &R1_12.x);
    fp12_mod_montgomery(&R1_12.y, &R1_12.y);
    g2_scm(&R2, &gen2, &s2);
    efp2_to_efp12(&R2_12, &R2);
    fp12_mod_montgomery(&R2_12.x, &R2_12.x);
    fp12_mod_montgomery(&R2_12.y, &R2_12.y);
    efp12_eca(&pk_tmp.g1, &R1_12, &R2_12);

    mpz_set_fr(sk_tmp.alpha, &alpha);

    G_random(&pk_tmp.h1, gen1, gen2);      // h1 ← G

    G_random(&pk_tmp.h2, gen1, gen2);      // h2 ← G   

    G_random(&pk_tmp.h3, gen1, gen2);      // h3 ← G

    G_random(&pk_tmp.h4, gen1, gen2);      // h4 ← G


    // mpz_urandomm(sk_tmp.alpha, state, order_z);       // α ← Zr

    // efp12_scm_lazy(&pk_tmp.g1, &pk_tmp.g, sk_tmp.alpha);   // g1 ← g^α

    // hk, f ← SHA-512

    // 代入
    efp12_set(&pk->g, &pk_tmp.g);
    efp12_set(&pk->g1, &pk_tmp.g1);
    efp12_set(&pk->h1, &pk_tmp.h1);
    efp12_set(&pk->h2, &pk_tmp.h2);
    efp12_set(&pk->h3, &pk_tmp.h3);
    efp12_set(&pk->h4, &pk_tmp.h4);
    mpz_set(sk->alpha, sk_tmp.alpha);

}

// 整数表現
void PreCal(PreValue *pv, PubKey pk){
    PreValue pv_tmp;
    
    // e(g, g)
    symmetric_paring(&pv_tmp.gg, pk.g, pk.g);
    fp12_mod_montgomery(&pv_tmp.gg, &pv_tmp.gg);

    // e(g, h1)
    symmetric_paring(&pv_tmp.gh1, pk.g, pk.h1);
    fp12_mod_montgomery(&pv_tmp.gh1, &pv_tmp.gh1);

    // e(g, h2)
    symmetric_paring(&pv_tmp.gh2, pk.g, pk.h2);
    fp12_mod_montgomery(&pv_tmp.gh2, &pv_tmp.gh2);

    // e(g, h3)
    symmetric_paring(&pv_tmp.gh3, pk.g, pk.h3);
    fp12_mod_montgomery(&pv_tmp.gh3, &pv_tmp.gh3);

    // e(g, h4)
    symmetric_paring(&pv_tmp.gh4, pk.g, pk.h4);
    fp12_mod_montgomery(&pv_tmp.gh4, &pv_tmp.gh4);

    fp12_set(&pv->gg, &pv_tmp.gg);
    fp12_set(&pv->gh1, &pv_tmp.gh1);
    fp12_set(&pv->gh2, &pv_tmp.gh2);
    fp12_set(&pv->gh3, &pv_tmp.gh3);
    fp12_set(&pv->gh4, &pv_tmp.gh4);
}

void HomKeyGen(HomKey *hk, PubKey pk, SecKey sk, mpz_t omega){
    mpz_t index;mpz_init(index);
    
    // 演算鍵（トラップドア）生成
    HomKey hk_tmp;hk_init(&hk_tmp);

    mpz_urandomm(hk_tmp.rh3.r_w, state, order_z);           // r_{ω,3} ← Zr

    mpz_urandomm(hk_tmp.rh4.r_w, state, order_z);           // r_{ω,4} ← Zr

    // h_{ω,3} ← (h_3 h^{-r_{ω,3}})^{1/(α-1))}
    mpz_neg(index, hk_tmp.rh3.r_w);                         // -r_{ω,3}
    mpz_mod(index, index, order_z);                         // index = -r_{ω,3} mod r
    efp12_scm(&hk_tmp.rh3.h_w, &pk.g, index);               // g^{-r_{ω,3}}
    efp12_eca(&hk_tmp.rh3.h_w, &hk_tmp.rh3.h_w, &pk.h3);    // h_3 + g^-r_{ω,3}
    mpz_sub(index, sk.alpha, omega);                        // α - ω
    mpz_mod(index, index, order_z);                         // index = α - ω mod r
    mpz_invert(index, index, order_z);                      // index = 1 / α - ω mod r
    efp12_scm(&hk_tmp.rh3.h_w, &hk_tmp.rh3.h_w, index);     // h_{ω,3} = (h_3 h^{-r_{ω,3}})^{1/(α-1))}

    // h_{ω,4} ← (h_4 h^{-r_{ω,4}})^{1/(α-1))}
    mpz_neg(index, hk_tmp.rh4.r_w);                         // -r_{ω,4}
    mpz_mod(index, index, order_z);                         // index = -r_{ω,4} mod r
    efp12_scm(&hk_tmp.rh4.h_w, &pk.g, index);               // g^{-r_{ω,4}}
    efp12_eca(&hk_tmp.rh4.h_w, &hk_tmp.rh4.h_w, &pk.h4);    // h_4 + g^{-r_{ω,4}}
    mpz_sub(index, sk.alpha, omega);                        // α - ω
    mpz_mod(index, index, order_z);                         // index = α - ω mod r
    mpz_invert(index, index, order_z);                      // index = 1 / α - ω
    efp12_scm(&hk_tmp.rh4.h_w, &hk_tmp.rh4.h_w, index);     // h_{ω,4} = (h_4 + h^{-r_{ω,4}})^{1/(α-1))}

    efp12_scm(&hk_tmp.g, &pk.g, omega);                     // g^{ω}

    // 代入
    efp12_set(&hk->g, &hk_tmp.g);
    mpz_set(hk->rh3.r_w, hk_tmp.rh3.r_w);
    efp12_set(&hk->rh3.h_w, &hk_tmp.rh3.h_w);
    mpz_set(hk->rh4.r_w, hk_tmp.rh4.r_w);
    efp12_set(&hk->rh4.h_w, &hk_tmp.rh4.h_w);

    // 解放
    mpz_clear(index);
}

void Enc(Ciphertext *ct, PubKey pk, fp12_t M, mpz_t omega, PreValue pv){
    mpz_t index;mpz_init(index);
    efp12_t efp12_tmp;efp12_init(&efp12_tmp);
    fp12_t fp12_tmp;fp12_init(&fp12_tmp);
    
    // 暗号文生成
    Ciphertext ct_tmp;ct_init(&ct_tmp);

    // s ← Zr
    mpz_t s;mpz_init(s);
    mpz_urandomm(s, state, order_z);

    // c1 ← g1^s g^{-sω}
    efp12_scm(&ct_tmp.c1, &pk.g1, s);                   // g1^s
    mpz_mul(index, s, omega);                           // sω
    mpz_neg(index, index);                              // -sω
    mpz_mod(index, index, order_z);                     // index = -sω mod r
    efp12_scm(&efp12_tmp, &pk.g, index);                // g^{-sω}
    efp12_eca(&ct_tmp.c1, &ct_tmp.c1, &efp12_tmp);      // c1 = g1^s + g^{-sω}

    // c2 ← e(g,g)^s
    fp12_pow(&ct_tmp.c2, &pv.gg, s);                    // c2 = e(g,g)^s

    // c3 ← M・e(g,h1)^-s
    fp12_pow(&ct_tmp.c3, &pv.gh1, s);                   // e(g,h1)^s
    fp12_inv(&ct_tmp.c3, &ct_tmp.c3);                   // e(g,h1)^-s
    fp12_mul(&ct_tmp.c3, &ct_tmp.c3, &M);               // c3 = M・e(g,h1)^-s

    // c4 ← e(g,h2)^s
    fp12_pow(&ct_tmp.c4, &pv.gh2, s);                   // c4 = e(g,h2)^s

    mpz_t delta;mpz_init(delta);
    Gamma(delta, ct_tmp.c1, ct_tmp.c2, ct_tmp.c3, ct_tmp.c4);        // δ ← Γ(c1,c2,c3,c4)

    // c5 ← e(g,h3)^s e(g,h4)^sδ
    fp12_t c5;fp12_init(&c5);
    fp12_pow(&fp12_tmp, &pv.gh3, s);                    // e(g,h3)^s
    mpz_mul(index, s, delta);                           // sδ
    fp12_pow(&c5, &pv.gh4, index);                      // e(g,h4)^sδ
    fp12_mul(&c5, &fp12_tmp, &c5);                      // c5 = e(g,h3)^s e(g,h4)^sδ

    Function(ct_tmp.tau, c5);                           // τ ← F(c5)

    // 代入
    efp12_set(&ct->c1, &ct_tmp.c1);
    fp12_set(&ct->c2, &ct_tmp.c2);
    fp12_set(&ct->c3, &ct_tmp.c3);
    fp12_set(&ct->c4, &ct_tmp.c4);
    mpz_set(ct->tau, ct_tmp.tau);

    // 解放
    mpz_clear(s);
    mpz_clear(index);
    mpz_clear(delta);
}

int Test(PubKey pk, HomKey hk, Ciphertext ct){
    mpz_t index;mpz_init(index);
    efp12_t efp12_tmp;efp12_init(&efp12_tmp);
    fp12_t fp12_tmp1;fp12_init(&fp12_tmp1);
    fp12_t fp12_tmp2;fp12_init(&fp12_tmp2);
    
    // 暗号文検証
    int test = 0;


    // δ ← Γ(c1,c2,c3,c4)
    mpz_t delta;mpz_init(delta);
    Gamma(delta, ct.c1, ct.c2, ct.c3, ct.c4);

    // τ_ch ← f(e(c1,h_{ω,3) h_{ω,4}^δ) c2^{r_{ω,3}+r_{ω,4}δ})
    mpz_t tau_ch;mpz_init(tau_ch);
    efp12_scm(&efp12_tmp, &hk.rh4.h_w, delta);          // h_{ω,4}^δ
    efp12_eca(&efp12_tmp, &hk.rh3.h_w, &efp12_tmp);     // h_{ω,3) h_{ω,4}^δ
    symmetric_paring(&fp12_tmp1, ct.c1, efp12_tmp);     // e(c1,h_{ω,3) h_{ω,4}^δ)
    fp12_mod_montgomery(&fp12_tmp1, &fp12_tmp1);
    mpz_mul(index, hk.rh4.r_w, delta);                  // r_{ω,4}δ
    mpz_add(index, index, hk.rh3.r_w);                  // r_{ω,3}+r_{ω,4}δ
    fp12_pow(&fp12_tmp2, &ct.c2, index);                // c2^{r_{ω,3}+r_{ω,4}δ}
    fp12_mul(&fp12_tmp1, &fp12_tmp1, &fp12_tmp2);       // e(c1,h_{ω,3) h_{ω,4}^δ) c2^{r_{ω,3}+r_{ω,4}δ}
    Function(tau_ch, fp12_tmp1);                        // τ_ch = f(e(c1,h_{ω,3) h_{ω,4}^δ) c2^{r_{ω,3}+r_{ω,4}δ})

    if(mpz_cmp(ct.tau, tau_ch) == 0){
        test = 1;
    }

    // 解放
    mpz_clear(index);

    return test;
}

void Dec(fp12_t *M, PubKey pk, SecKey sk, mpz_t omega, Ciphertext ct){
    mpz_t index;mpz_init(index);
    fp12_t fp12_tmp;fp12_init(&fp12_tmp);
    efp12_t efp12_tmp;efp12_init(&efp12_tmp);

    // 復号
    fp12_t M_tmp;fp12_init(&M_tmp);

    mpz_t r_w1;mpz_init(r_w1);
    mpz_urandomm(r_w1, state, order_z);         // r_{ω,1} ← Zr
    mpz_t r_w2;mpz_init(r_w2);
    mpz_urandomm(r_w2, state, order_z);         // r_{ω,2} ← Zr
    mpz_t r_w3;mpz_init(r_w3);
    mpz_urandomm(r_w3, state, order_z);         // r_{ω,3} ← Zr
    mpz_t r_w4;mpz_init(r_w4);
    mpz_urandomm(r_w4, state, order_z);         // r_{ω,4} ← Zr


    // h_{ω,1} ← (h1 g^{-r_{ω,1}})^{1/(α-ω)}
    efp12_t h_w1;efp12_init(&h_w1);
    mpz_neg(index, r_w1);                       // -r_{ω,1}
    mpz_mod(index, index, order_z);             // -r_{ω,1} mod r
    efp12_scm(&h_w1, &pk.g, index);             // g^{-r_{ω,1}}
    efp12_eca(&h_w1, &pk.h1, &h_w1);            // h1 g^{-r_{ω,1}}
    mpz_sub(index, sk.alpha, omega);            // α - ω
    mpz_mod(index, index, order_z);             // α - ω mod r
    mpz_invert(index, index, order_z);          // 1 / α - ω
    efp12_scm(&h_w1, &h_w1, index);             // (h1 g^{-r_{ω,1}})^{1/(α-ω)}

    // h_{ω,2} ← (h2 g^{-r_{ω,2}})^{1/(α-ω)}
    efp12_t h_w2;efp12_init(&h_w2);
    mpz_neg(index, r_w2);                       // r_{ω,2}
    mpz_mod(index, index, order_z);             // -r_{ω,2} mod r
    efp12_scm(&h_w2, &pk.g, index);             // g^{-r_{ω,2}}
    efp12_eca(&h_w2, &pk.h2, &h_w2);            // h2 g^{-r_{ω,2}}
    mpz_sub(index, sk.alpha, omega);            // α - ω
    mpz_mod(index, index, order_z);             // α - ω mod r
    mpz_invert(index, index, order_z);          // 1 / α - ω
    efp12_scm(&h_w2, &h_w2, index);             // (h2 g^{-r_{ω,2}})^{1/(α-ω)}

    // h_{ω,3} ← (h3 g^{-r_{ω,3}})^{1/(α-ω)}
    efp12_t h_w3;efp12_init(&h_w3);
    mpz_neg(index, r_w3);                       // r_{ω,3}
    mpz_mod(index, index, order_z);             // -r_{ω,3} mod r
    efp12_scm(&h_w3, &pk.g, index);             // g^{-r_{ω,3}}
    efp12_eca(&h_w3, &pk.h3, &h_w3);            // h3 g^{-r_{ω,3}}
    mpz_sub(index, sk.alpha, omega);            // α - ω
    mpz_mod(index, index, order_z);             // α - ω mod r
    mpz_invert(index, index, order_z);          // 1 / α - ω
    efp12_scm(&h_w3, &h_w3, index);             // (h3 g^{-r_{ω,3}})^{1/(α-ω)}

    // h_{ω,4} ← (h4 g^{-r_{ω,4}})^{1/(α-ω)}
    efp12_t h_w4;efp12_init(&h_w4);
    mpz_neg(index, r_w4);                       // r_{ω,4}
    mpz_mod(index, index, order_z);             // -r_{ω,4} mod r
    efp12_scm(&h_w4, &pk.g, index);             // g^{-r_{ω,4}}
    efp12_eca(&h_w4, &pk.h4, &h_w4);            // h4 g^{-r_{ω,4}}
    mpz_sub(index, sk.alpha, omega);            // α - ω
    mpz_mod(index, index, order_z);             // α - ω mod r
    mpz_invert(index, index, order_z);          // 1 / α - ω
    efp12_scm(&h_w4, &h_w4, index);             // (h4 g^{-r_{ω,4}})^{1/(α-ω)}

    mpz_t delta;mpz_init(delta);
    Gamma(delta, ct.c1, ct.c2, ct.c3, ct.c4);           // δ ← Γ(c1,c2,c3,c4)

    // c4' ← e(c1,h_{ω,2}) c2^{r_{ω,2}}
    fp12_t c4_prime;fp12_init(&c4_prime);
    symmetric_paring(&c4_prime, ct.c1, h_w2);           // e(c1,h_{ω,2})
    fp12_mod_montgomery(&c4_prime, &c4_prime);
    fp12_pow(&fp12_tmp, &ct.c2, r_w2);                  // c2^{r_{ω,2}}
    fp12_mul(&c4_prime, &c4_prime, &fp12_tmp);          // e(c1,h_{ω,2}) c2^{r_{ω,2}}

    // c5 ← e(c1,h_{ω,3} h_{ω,4}^δ) c2^{r_{ω,3}+r_{ω,4}δ}
    fp12_t c5;fp12_init(&c5);
    efp12_scm(&efp12_tmp, &h_w4, delta);                // h_{ω,4}^δ
    efp12_eca(&efp12_tmp, &h_w3, &efp12_tmp);           // h_{ω,3) h_{ω,4}^δ
    symmetric_paring(&fp12_tmp, ct.c1, efp12_tmp);      // e(c1,h_{ω,3) h_{ω,4}^δ)
    fp12_mod_montgomery(&fp12_tmp, &fp12_tmp);
    mpz_mul(index, r_w4, delta);                        // r_{ω,4}δ
    mpz_add(index, index, r_w3);                        // r_{ω,3}+r_{ω,4}δ
    fp12_pow(&c5, &ct.c2, index);                       // c2^{r_{ω,3}+r_{ω,4}δ}
    fp12_mul(&c5, &fp12_tmp, &c5);                      // e(c1,h_{ω,3} h_{ω,4}^δ) c2^{r_{ω,3}+r_{ω,4}δ}

    // τ_ch ← F(c5)
    mpz_t tau_ch;mpz_init(tau_ch);
    Function(tau_ch, c5);

    if(fp12_cmp(&ct.c4, &c4_prime) != 0 || mpz_cmp(ct.tau, tau_ch) != 0){
        fp12_set(M, &M_tmp);
        printf("DecryptionError\n");
        if(fp12_cmp(&ct.c4, &c4_prime) != 0){
            printf("c4 != c4'\n");
        }
        if(mpz_cmp(ct.tau, tau_ch) != 0){
            printf("tau != tau_ch\n");
        }
    }

    else{
        // M ← c3 e(c1,h_{ω,1}) c2^{r_{ω,1}}
        symmetric_paring(&fp12_tmp, ct.c1, h_w1);       // e(c1,h_{ω,1})
        fp12_mod_montgomery(&fp12_tmp, &fp12_tmp);
        fp12_mul(&M_tmp, &ct.c3, &fp12_tmp);            // c3 e(c1,h_{ω,1})
        fp12_pow(&fp12_tmp, &ct.c2, r_w1);              // c2^{r_{ω,1}}
        fp12_mul(&M_tmp, &M_tmp, &fp12_tmp);            // M = c3 e(c1,h_{ω,1}) c2^{r_{ω,1}}
        fp12_set(M, &M_tmp);
    }

    // 解放
    mpz_clear(index);
    mpz_clear(delta);
    mpz_clear(r_w1);
    mpz_clear(r_w2);
    mpz_clear(r_w3);
    mpz_clear(r_w4);
    mpz_clear(tau_ch);

}

void Eval(Ciphertext *ct, PubKey pk, HomKey hk, Ciphertext ct1, Ciphertext ct2, PreValue pv){
    
    mpz_t index;mpz_init(index);
    fp12_t fp12_tmp;fp12_init(&fp12_tmp);
    efp12_t efp12_tmp;efp12_init(&efp12_tmp);
    
    // δ_1 = Γ(c_{1,1}, c_{1,2}, c_{1,3}, c_{1,4})
    mpz_t delta1;mpz_init(delta1);
    Gamma(delta1, ct1.c1, ct1.c2, ct1.c3, ct1.c4);

    // c_{1,5} ← e(c_{1,1},h_{ω,3} h_{ω,4}^δ_1) c_{1,2}^{r_{ω,3}+r_{ω,4}δ_1}
    fp12_t c1_5;fp12_init(&c1_5);
    efp12_scm(&efp12_tmp, &hk.rh4.h_w, delta1);         // h_{ω,4}^δ_1
    efp12_eca(&efp12_tmp, &efp12_tmp, &hk.rh3.h_w);     // h_{ω,3} + h_{ω,4}^δ_1
    symmetric_paring(&c1_5, ct1.c1, efp12_tmp);         // e(c_{1,1},h_{ω,3} h_{ω,4}^δ_1)
    fp12_mod_montgomery(&c1_5, &c1_5);
    mpz_mul(index, hk.rh4.r_w, delta1);                 // r_{ω,4}δ_1
    mpz_add(index, index, hk.rh3.r_w);                  // r_{ω,3}+r_{ω,4}δ_1
    fp12_pow(&fp12_tmp, &ct1.c2, index);                // c_{1,2}^{r_{ω,3}+r_{ω,4}δ_1}
    fp12_mul(&c1_5, &c1_5, &fp12_tmp);                  // c_{1,5} = e(c_{1,1},h_{ω,3} h_{ω,4}^δ_1) c_{1,2}^{r_{ω,3}+r_{ω,4}δ_1}

    // δ_2 = Γ(c_{2,1}, c_{2,2}, c_{2,3}, c_{2,4})
    mpz_t delta2;mpz_init(delta2);
    Gamma(delta2, ct2.c1, ct2.c2, ct2.c3, ct2.c4);

    // c_{2,5} ← e(c_{2,1},h_{ω,3} h_{ω,4}^δ_2) c_{2,2}^{r_{ω,3}+r_{ω,4}δ_2}
    fp12_t c2_5;fp12_init(&c2_5);
    efp12_scm(&efp12_tmp, &hk.rh4.h_w, delta2);       // h_{ω,4}^δ_2
    efp12_eca(&efp12_tmp, &efp12_tmp, &hk.rh3.h_w);    // h_{ω,3} + h_{ω,4}^δ_2
    symmetric_paring(&c2_5, ct2.c1, efp12_tmp);         // e(c_{2,1},h_{ω,3} h_{ω,4}^δ_2)
    fp12_mod_montgomery(&c2_5, &c2_5);
    mpz_mul(index, hk.rh4.r_w, delta2);                 // r_{ω,4}δ_2
    mpz_add(index, index, hk.rh3.r_w);                  // r_{ω,3}+r_{ω,4}δ_2
    fp12_pow(&fp12_tmp, &ct2.c2, index);              // c_{2,2}^{r_{ω,3}+r_{ω,4}δ_2}
    fp12_mul(&c2_5, &c2_5, &fp12_tmp);                 // c_{2,5} = e(c_{2,1},h_{ω,3} h_{ω,4}^δ_2) c_{2,2}^{r_{ω,3}+r_{ω,4}δ_2}

    mpz_t tau_ch1;mpz_init(tau_ch1);
    Function(tau_ch1, c1_5);                            // τ_ch1 ← F(c_{1,5})
    mpz_t tau_ch2;mpz_init(tau_ch2);
    Function(tau_ch2, c2_5);                            // τ_ch2 ← F(c_{2,5})
    if(mpz_cmp(tau_ch1, ct1.tau) != 0 || mpz_cmp(tau_ch2, ct2.tau) != 0){
        printf("EvaluationError\n");
    }

    else{
        Ciphertext ct_tmp;ct_init(&ct_tmp);

        // s ← Zr
        mpz_t s;mpz_init(s);
        mpz_urandomm(s, state, order_z);

        // c1 ← c_{1,1} c_{2,1} g_1^s g^{-sω}
        efp12_set(&ct_tmp.c1, &pk.g1);
        fp12_set_neg(&ct_tmp.c1.y, &ct_tmp.c1.y);         // g_1^{-1}
        efp12_eca(&ct_tmp.c1, &ct_tmp.c1, &hk.g);      // g_1^{-1} g^ω
        mpz_neg(index, s);                              // -s
        mpz_mod(index, index, order_z);                       // -s mod r
        efp12_scm(&ct_tmp.c1, &ct_tmp.c1, index);     // g_1^s g^{-sω}
        efp12_eca(&ct_tmp.c1, &ct_tmp.c1, &ct1.c1);    // c_{1,1} g_1^s g^{-sω}
        efp12_eca(&ct_tmp.c1, &ct_tmp.c1, &ct2.c1);    // c_{1,1} c_{2,1} g_1^s g^{-sω}

        // c2 ← c_{1,2} c_{2,2} e(g,g)^s
        fp12_pow(&ct_tmp.c2, &pv.gg, s);          // e(g,g)^s
        fp12_mul(&ct_tmp.c2, &ct_tmp.c2, &ct1.c2);     // c_{1,2} e(g,g)^s
        fp12_mul(&ct_tmp.c2, &ct_tmp.c2, &ct2.c2);     // c_{1,2} c_{2,2} e(g,g)^s

        // c3 ← c_{1,3} c_{2,3} e(g,h_1)^{-s}
        mpz_neg(index, s);                              // -s
        mpz_mod(index, index, order_z);                       // -s mod r
        fp12_pow(&ct_tmp.c3, &pv.gh1, index);      // e(g,h_1)^{-s}
        fp12_mul(&ct_tmp.c3, &ct_tmp.c3, &ct1.c3);     // c_{1,3} e(g,h_1)^{-s}
        fp12_mul(&ct_tmp.c3, &ct_tmp.c3, &ct2.c3);     // c_{1,3} c_{2,3} e(g,h_1)^{-s}

        // c4 ← c_{1,4} c_{2,4} e(g,h_2)^s
        fp12_pow(&ct_tmp.c4, &pv.gh2, s);          // e(g,h_2)^s
        fp12_mul(&ct_tmp.c4, &ct_tmp.c4, &ct1.c4);     // c_{1,4} e(g,h_2)^s
        fp12_mul(&ct_tmp.c4, &ct_tmp.c4, &ct2.c4);     // c_{1,4} c_{2,4} e(g,h_2)^s

        // δ ← Γ(c1,c2,c3,c4)
        mpz_t delta;mpz_init(delta);
        Gamma(delta, ct_tmp.c1, ct_tmp.c2, ct_tmp.c3, ct_tmp.c4);
        
        // c5 ← e(c1,h_{ω,3} h_{ω,4}^δ) c2^{r_{ω,3}+r_{ω,4}δ}
        fp12_t c5;fp12_init(&c5);
        efp12_scm(&efp12_tmp, &hk.rh4.h_w, delta);            // h_{ω,4}^δ
        efp12_eca(&efp12_tmp, &hk.rh3.h_w, &efp12_tmp);        // h_{ω,3) h_{ω,4}^δ
        symmetric_paring(&fp12_tmp, ct_tmp.c1, efp12_tmp);      // e(c1,h_{ω,3} h_{ω,4}^δ)
        fp12_mod_montgomery(&fp12_tmp, &fp12_tmp);
        mpz_mul(index, hk.rh4.r_w, delta);                      // r_{ω,4}δ
        mpz_add(index, index, hk.rh3.r_w);                      // r_{ω,3}+r_{ω,4}δ
        fp12_pow(&c5, &ct_tmp.c2, index);                     // c2^{r_{ω,3}+r_{ω,4}δ}
        fp12_mul(&c5, &fp12_tmp, &c5);                         // c5 = e(c1,h_{ω,3} h_{ω,4}^δ) c2^{r_{ω,3}+r_{ω,4}δ}

        Function(ct_tmp.tau, c5);                               // τ ← F(c5)

        // 代入
        efp12_set(&ct->c1, &ct_tmp.c1);
        fp12_set(&ct->c2, &ct_tmp.c2);
        fp12_set(&ct->c3, &ct_tmp.c3);
        fp12_set(&ct->c4, &ct_tmp.c4);
        mpz_set(ct->tau, ct_tmp.tau);

        // 解放
        mpz_clear(s);
        mpz_clear(delta);
    }

    mpz_clear(index);
}


int main (){
    // bls parameters
    bls12_init();
    bls12_print_parameters();
    printf("====================================================================================\n");
    printf("\n");


    // 生成元
    efp_t gen1;efp_init(&gen1);
    efp2_t gen2;efp2_init(&gen2);
    g1_set_random(&gen1, state);
    g2_set_random(&gen2, state);


    // Message
    fp12_t M;fp12_init(&M);
    fp_set_ui(&M.x0.x0.x0, 1);
    fp_set_ui(&M.x0.x0.x1, 2);
    fp_set_ui(&M.x0.x1.x0, 3);
    fp_set_ui(&M.x0.x1.x1, 4);
    fp_set_ui(&M.x0.x2.x0, 5);
    fp_set_ui(&M.x0.x2.x1, 6);
    fp_set_ui(&M.x1.x0.x0, 7);
    fp_set_ui(&M.x1.x0.x1, 8);
    fp_set_ui(&M.x1.x1.x0, 9);
    fp_set_ui(&M.x1.x1.x1, 10);
    fp_set_ui(&M.x1.x2.x0, 11);
    fp_set_ui(&M.x1.x2.x1, 12);
    fp12_println("Message : ", &M);
    printf("\n");

    mpz_t omega;mpz_init(omega);
    // mpz_urandomm(omega, state, order_z);
    mpz_set_str(omega, "10", 10);
    gmp_printf("Key Word : %Zd\n", omega);
    printf("\n");
    

    // Key Generate
    printf("Public Key = (g, g_1, h_1, h_2, h_3, h_4, hk, f)\n");
    PubKey pk;pk_init(&pk);
    pk_init(&pk);
    SecKey sk;
    sk_init(&sk);
    KeyGen(&pk, &sk, gen1, gen2);
    efp12_println("g : ", &pk.g);
    efp12_println("g_1 : ", &pk.g1);
    efp12_println("h_1 : ", &pk.h1);
    efp12_println("h_2 : ", &pk.h2);
    efp12_println("h_3 : ", &pk.h3);
    efp12_println("h_4 : ", &pk.h4);
    printf("hk, f : SHA-512\n");
    printf("\n");

    printf("Secret Key = α\n");
    gmp_printf("α : %Zd\n", sk.alpha);
    printf("\n");


    // Pre Calculation
    PreValue pv;pv_init(&pv);
    PreCal(&pv, pk);


    // Homomorphic Key Generate
    printf("Homomorphic Key = (g^ω, (r_{ω,3}, h_{ω,3}), (r_{ω,4}, h_{ω,4}))\n");
    HomKey hk;hk_init(&hk);
    HomKeyGen(&hk, pk, sk, omega);
    efp12_println("g^ω : ", &hk.g);
    gmp_printf("r_{ω,3} : %Zd\n", hk.rh3.r_w);
    efp12_println("h_{ω,3} : ", &hk.rh3.h_w);
    gmp_printf("r_{ω,4} : %Zd\n", hk.rh4.r_w);
    efp12_println("h_{ω,4} : ", &hk.rh4.h_w);
    printf("\n");


    // Encryption
    printf("Cypher Text = (c_1, c_2, c_3, c_4, τ)\n");
    Ciphertext ct;ct_init(&ct);
    Enc(&ct, pk, M, omega, pv);
    efp12_println("c_1 : ", &ct.c1);
    fp12_println("c_2 : ", &ct.c2);
    fp12_println("c_3 : ", &ct.c3);
    fp12_println("c_4 : ", &ct.c4);
    gmp_printf("τ : %Zd\n", ct.tau);
    printf("\n");


    // Test
    int test = 0;
    test = Test(pk, hk, ct);
    printf("Test Result : ");
    if(test == 1){printf("Success\n");}
    else{printf("Failed\n");}
    printf("\n");


    // Decryption
    printf("Decryption\n");
    Dec(&M, pk, sk, omega, ct);
    fp12_println("Decrypted Message : ", &M);
    printf("\n");


    // Multiplication
    fp12_t M2;fp12_init(&M2);
    fp_set_ui(&M2.x0.x0.x0, 2);
    fp_set_ui(&M2.x0.x0.x1, 3);
    fp_set_ui(&M2.x0.x1.x0, 4);
    fp_set_ui(&M2.x0.x1.x1, 5);
    fp_set_ui(&M2.x0.x2.x0, 6);
    fp_set_ui(&M2.x0.x2.x1, 7);
    fp_set_ui(&M2.x1.x0.x0, 8);
    fp_set_ui(&M2.x1.x0.x1, 9);
    fp_set_ui(&M2.x1.x1.x0, 10);
    fp_set_ui(&M2.x1.x1.x1, 11);
    fp_set_ui(&M2.x1.x2.x0, 12);
    fp_set_ui(&M2.x1.x2.x1, 1);
    fp12_println("Message2 : ", &M2);
    printf("\n");

    printf("Cypher Text2 = (c2_1, c2_2, c2_3, c2_4, τ2)\n");
    Ciphertext ct2;ct_init(&ct2);
    Enc(&ct2, pk, M2, omega, pv);
    efp12_println("c2_1 : ", &ct2.c1);
    fp12_println("c2_2 : ", &ct2.c2);
    fp12_println("c2_3 : ", &ct2.c3);
    fp12_println("c2_4 : ", &ct2.c4);
    gmp_printf("τ2 : %Zd\n", ct.tau);
    printf("\n");

    Ciphertext ct_mul;ct_init(&ct_mul);
    Eval(&ct_mul, pk, hk, ct, ct2, pv);
    efp12_println("c_mul_1 : ", &ct_mul.c1);
    fp12_println("c_mul_2 : ", &ct_mul.c2);
    fp12_println("c_mul_3 : ", &ct_mul.c3);
    fp12_println("c_mul_4 : ", &ct_mul.c4);
    gmp_printf("τ_mul : %Zd\n", ct_mul.tau);
    printf("\n");

    test = Test(pk, hk, ct_mul);
    printf("Test Result : ");
    if(test == 1){printf("Success\n");}
    else{printf("Failed\n");}
    printf("\n");

    fp12_t M_mul;fp12_init(&M_mul);
    Dec(&M_mul, pk, sk, omega, ct_mul);
    fp12_println("Decrypted Message : ", &M_mul);

    fp12_t M3;fp12_init(&M3);
    fp12_mul(&M3, &M, &M2);
    fp12_println("M3 : ", &M3);


    // 時間測定
    struct timespec ts;
    struct tm tm;
    long double start_sec, start_nsec, end_sec, end_nsec, exe_time_buf0;
    int loop = 100;

    // efp12_t R;efp12_init(&R);
    // G_random
    // clock_gettime(CLOCK_REALTIME, &ts);
    // start_sec = ts.tv_sec;
    // start_nsec = ts.tv_nsec;
    // for (int i = 0; i < loop; i ++){
    //     G_random(&R, gen1, gen2);
    // }
    // clock_gettime(CLOCK_REALTIME, &ts);
    // end_sec = ts.tv_sec;
    // end_nsec = ts.tv_nsec;
    // exe_time_buf0 = ((end_sec - start_sec) + (end_nsec - start_nsec) / 1000000000.0)*1000.0/loop;
    // printf("Execution time for KeyGen!! : %Lf[ms]\n", exe_time_buf0);

    // KeyGen
    clock_gettime(CLOCK_REALTIME, &ts);
    start_sec = ts.tv_sec;
    start_nsec = ts.tv_nsec;
    for (int i = 0; i < loop; i ++){
        KeyGen(&pk, &sk, gen1, gen2);
    }
    clock_gettime(CLOCK_REALTIME, &ts);
    end_sec = ts.tv_sec;
    end_nsec = ts.tv_nsec;
    exe_time_buf0 = ((end_sec - start_sec) + (end_nsec - start_nsec) / 1000000000.0)*1000.0/loop;
    printf("Execution time for KeyGen!! : %Lf[ms]\n", exe_time_buf0);

    // // HomKeyGen
    // clock_gettime(CLOCK_REALTIME, &ts);
    // start_sec = ts.tv_sec;
    // start_nsec = ts.tv_nsec;
    // for (int i = 0; i < loop; i ++){
    //     HomKeyGen(&hk, pk, sk, omega);
    // }
    // clock_gettime(CLOCK_REALTIME, &ts);
    // end_sec = ts.tv_sec;
    // end_nsec = ts.tv_nsec;
    // exe_time_buf0 = ((end_sec - start_sec) + (end_nsec - start_nsec) / 1000000000.0)*1000.0/loop;
    // printf("Execution time for HomKeyGen!! : %Lf[ms]\n", exe_time_buf0);

    // // Enc
    // clock_gettime(CLOCK_REALTIME, &ts);
    // start_sec = ts.tv_sec;
    // start_nsec = ts.tv_nsec;
    // for (int i = 0; i < loop; i ++){
    //     Enc(&ct, pk, M, omega, pv);
    // }
    // clock_gettime(CLOCK_REALTIME, &ts);
    // end_sec = ts.tv_sec;
    // end_nsec = ts.tv_nsec;
    // exe_time_buf0 = ((end_sec - start_sec) + (end_nsec - start_nsec) / 1000000000.0)*1000.0/loop;
    // printf("Execution time for Enc!! : %Lf[ms]\n", exe_time_buf0);

    // // Dec
    // clock_gettime(CLOCK_REALTIME, &ts);
    // start_sec = ts.tv_sec;
    // start_nsec = ts.tv_nsec;
    // for (int i = 0; i < loop; i ++){
    //     Dec(&M, pk, sk, omega, ct);
    // }
    // clock_gettime(CLOCK_REALTIME, &ts);
    // end_sec = ts.tv_sec;
    // end_nsec = ts.tv_nsec;
    // exe_time_buf0 = ((end_sec - start_sec) + (end_nsec - start_nsec) / 1000000000.0)*1000.0/loop;
    // printf("Execution time for Dec!! : %Lf[ms]\n", exe_time_buf0);

    // // Test
    // clock_gettime(CLOCK_REALTIME, &ts);
    // start_sec = ts.tv_sec;
    // start_nsec = ts.tv_nsec;
    // for (int i = 0; i < loop; i ++){
    //     Test(pk, hk, ct);
    // }
    // clock_gettime(CLOCK_REALTIME, &ts);
    // end_sec = ts.tv_sec;
    // end_nsec = ts.tv_nsec;
    // exe_time_buf0 = ((end_sec - start_sec) + (end_nsec - start_nsec) / 1000000000.0)*1000.0/loop;
    // printf("Execution time for Test!! : %Lf[ms]\n", exe_time_buf0);

    // // Eval
    // clock_gettime(CLOCK_REALTIME, &ts);
    // start_sec = ts.tv_sec;
    // start_nsec = ts.tv_nsec;
    // for (int i = 0; i < loop; i ++){
    //     Eval(&ct_mul, pk, hk, ct, ct2, pv);
    // }
    // clock_gettime(CLOCK_REALTIME, &ts);
    // end_sec = ts.tv_sec;
    // end_nsec = ts.tv_nsec;
    // exe_time_buf0 = ((end_sec - start_sec) + (end_nsec - start_nsec) / 1000000000.0)*1000.0/loop;
    // printf("Execution time for Eval!! : %Lf[ms]\n", exe_time_buf0);
}
