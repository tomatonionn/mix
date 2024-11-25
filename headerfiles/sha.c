#include "../miller_header.h"

void Gamma(mpz_t gamma, const struct efp12 c1, const struct fp12 c2, const struct fp12 c3, const struct fp12 c4, const mpz_t p){

    // (c1, c2, c3, c4) → c_sum_fp12 → c_sum
    mpz_t c_sum;mpz_init(c_sum);
    struct fp12 c_sum_fp12;fp12_init(&c_sum_fp12);

    // c1.x +c1.y + c2 + c3 + c4 → c_sum_fp12
    fp12_add(&c_sum_fp12, c1.x, c1.y, p);
    fp12_add(&c_sum_fp12, c_sum_fp12, c2, p);
    fp12_add(&c_sum_fp12, c_sum_fp12, c3, p);
    fp12_add(&c_sum_fp12, c_sum_fp12, c4, p);

    // fp12_sum.x0.x0.x0.x0 + fp12_sum.x0.x0.x0.x0 + ... → mpz_t c_sum
    mpz_add(c_sum, c_sum_fp12.x0.x0.x0.x0, c_sum_fp12.x0.x0.x1.x0);
    mpz_add(c_sum, c_sum, c_sum_fp12.x0.x1.x0.x0);
    mpz_add(c_sum, c_sum, c_sum_fp12.x0.x1.x1.x0);
    mpz_add(c_sum, c_sum, c_sum_fp12.x0.x2.x0.x0);
    mpz_add(c_sum, c_sum, c_sum_fp12.x0.x2.x1.x0);
    mpz_add(c_sum, c_sum, c_sum_fp12.x1.x0.x0.x0);
    mpz_add(c_sum, c_sum, c_sum_fp12.x1.x0.x1.x0);
    mpz_add(c_sum, c_sum, c_sum_fp12.x1.x1.x0.x0);
    mpz_add(c_sum, c_sum, c_sum_fp12.x1.x1.x1.x0);
    mpz_add(c_sum, c_sum, c_sum_fp12.x1.x2.x0.x0);
    mpz_add(c_sum, c_sum, c_sum_fp12.x1.x2.x1.x0);


    // HASH(c_sum) → gamma
    // mpz_t型をバイト列に変換
    size_t c_sum_len = (mpz_sizeinbase(c_sum, 2) + 7) / 8; // バイト数を計算
    unsigned char *c_sum_bytes = (unsigned char *)malloc(c_sum_len);
    mpz_export(c_sum_bytes, NULL, 1, 1, 0, 0, c_sum);

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

    mpz_clear(c_sum);
}

void Function(mpz_t function, const struct fp12 c5){

    // c5.x0.x0.x0.x0 + c5.x0.x0.x0.x0 + ... → mpz_t c5_sum
    mpz_t c5_sum;mpz_init(c5_sum);
    mpz_add(c5_sum, c5.x0.x0.x0.x0, c5.x0.x0.x1.x0);
    mpz_add(c5_sum, c5_sum, c5.x0.x1.x0.x0);
    mpz_add(c5_sum, c5_sum, c5.x0.x1.x1.x0);
    mpz_add(c5_sum, c5_sum, c5.x0.x2.x0.x0);
    mpz_add(c5_sum, c5_sum, c5.x0.x2.x1.x0);
    mpz_add(c5_sum, c5_sum, c5.x1.x0.x0.x0);
    mpz_add(c5_sum, c5_sum, c5.x1.x0.x1.x0);
    mpz_add(c5_sum, c5_sum, c5.x1.x1.x0.x0);
    mpz_add(c5_sum, c5_sum, c5.x1.x1.x1.x0);
    mpz_add(c5_sum, c5_sum, c5.x1.x2.x0.x0);
    mpz_add(c5_sum, c5_sum, c5.x1.x2.x1.x0);


    // HASH(c5_sum) → function
    // mpz_t型をバイト列に変換
    size_t c5_sum_len = (mpz_sizeinbase(c5_sum, 2) + 7) / 8; // バイト数を計算
    unsigned char *c5_sum_bytes = (unsigned char *)malloc(c5_sum_len);
    mpz_export(c5_sum_bytes, NULL, 1, 1, 0, 0, c5_sum);

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

    mpz_clear(c5_sum);
}
