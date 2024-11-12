#include "../miller_header.h"

// f(x) = x^2 + 1
// f(α) = α^2 + 1 = 0
// α^2 = -1

// 初期化
void fp2_init(struct fp2 *X){
    fp_init(&X->x0);
    fp_init(&X->x1);
}

// 解放
void fp2_clear(struct fp2 *X){
    fp_clear(&X->x0);
    fp_clear(&X->x1);
}

// 表示
void fp2_printf(const struct fp2 X){
    printf("(");
    fp_printf(X.x0);
    printf(", ");
    fp_printf(X.x1);
    printf(")");
}

// 代入
void fp2_set(struct fp2 *S, const struct fp2 X){
    fp_set(&S->x0, X.x0);
    fp_set(&S->x1, X.x1);
}

// 生成
void fp2_random(struct fp2 *X, const mpz_t p, gmp_randstate_t state){
    fp_random(&X->x0, p, state);
    fp_random(&X->x1, p, state);
}

// 比較
int fp2_cmp(const struct fp2 X, const struct fp2 Y){
    if(fp_cmp(X.x0, Y.x0) == 0 && fp_cmp(X.x1, Y.x1) == 0){
        return 0;
    }
    else{
        return 1;
    }
}

// 負数
void fp2_neg(struct fp2 *S, struct fp2 X, const mpz_t p){
    fp_neg(&S->x0, X.x0, p);
    fp_neg(&S->x1, X.x1, p);
}

// 加算
void fp2_add(struct fp2 *S, const struct fp2 X, const struct fp2 Y, const mpz_t p){
    struct fp2 temp;fp2_init(&temp);

    fp_add(&temp.x0, X.x0, Y.x0, p);
    fp_add(&temp.x1, X.x1, Y.x1, p);
    fp2_set(S, temp);

    fp2_clear(&temp);
}

// 減算
void fp2_sub(struct fp2 *S, const struct fp2 X, const struct fp2 Y, const mpz_t p){
    struct fp2 temp;fp2_init(&temp);

    fp_sub(&temp.x0, X.x0, Y.x0, p);
    fp_sub(&temp.x1, X.x1, Y.x1, p);
    fp2_set(S, temp);

    fp2_clear(&temp);
}

// 乗算
void fp2_mul(struct fp2 *S, const struct fp2 X, const struct fp2 Y, const mpz_t p){
    struct fp T1, T2, T3, T4;
    fp_init(&T1);fp_init(&T2);fp_init(&T3);fp_init(&T4);
    struct fp2 temp;fp2_init(&temp);

    fp_mul(&T1, X.x0, Y.x0, p);
    fp_mul(&T2, X.x1, Y.x1, p);

    fp_add(&T3, X.x0, X.x1, p);
    fp_add(&T4, Y.x0, Y.x1, p);

    fp_sub(&temp.x0, T1, T2, p);
    fp_mul(&temp.x1, T3, T4, p);
    fp_sub(&temp.x1, temp.x1, T1, p);
    fp_sub(&temp.x1, temp.x1, T2, p);

    fp2_set(S, temp);

    fp_clear(&T1);fp_clear(&T2);fp_clear(&T3);fp_clear(&T4);
    fp2_clear(&temp);
}

// ２乗
void fp2_square(struct fp2 *S, const struct fp2 X, const mpz_t p){
    struct fp T1, T2, T3;
    fp_init(&T1);fp_init(&T2);fp_init(&T3);
    struct fp2 temp;fp2_init(&temp);

    fp_add(&T1, X.x0, X.x1, p);
    fp_sub(&T2, X.x0, X.x1, p);
    fp_mul(&T3, X.x0, X.x1, p);

    fp_mul(&temp.x0, T1, T2, p);
    fp_add(&temp.x1, T3, T3, p);

    fp2_set(S, temp);

    fp_clear(&T1);fp_clear(&T2);fp_clear(&T3);
    fp2_clear(&temp);
}

// 逆元
void fp2_inv(struct fp2 *S, const struct fp2 X, const mpz_t p){
    struct fp s, s0, s1, s_inv;
    fp_init(&s);fp_init(&s0);fp_init(&s1);fp_init(&s_inv);

    fp_mul(&s0, X.x0, X.x0, p);
    fp_mul(&s1, X.x1, X.x1, p);
    fp_add(&s, s0, s1, p);          // s = X0^2 + X1^2
    fp_inv(&s_inv, s, p);           // s = 1 / (X0^2 + X1^2)

    struct fp2 temp;fp2_init(&temp);
    fp_mul(&temp.x0, s_inv, X.x0, p);
    fp_mul(&temp.x1, s_inv, X.x1, p);
    mpz_neg(temp.x1.x0, temp.x1.x0);
    mpz_mod(temp.x1.x0, temp.x1.x0, p);
    fp2_set(S, temp);

    fp_clear(&s);fp_clear(&s0);fp_clear(&s1);fp_clear(&s_inv);
    fp2_clear(&temp);
}

// 冪乗算
void fp2_pow(struct fp2 *S, const struct fp2 X, const mpz_t s, const mpz_t p){
    struct fp2 temp;fp2_init(&temp);
    mpz_set_str(temp.x0.x0, "1", 10);
    char *scalar_binary = mpz_get_str (NULL, 2, s);
    size_t len = strlen(scalar_binary);

    for(size_t i = 0; i < len; i++){
        fp2_square(&temp, temp, p);
        if(scalar_binary[i] - '0' == 1){
            fp2_mul(&temp, temp, X, p);
        }
    }
    fp2_set(S, temp);

    free(scalar_binary); scalar_binary = NULL;
    fp2_clear(&temp);
}

// Frobenius写像
void fp2_Frobenius(struct fp2 *S, const struct fp2 X, const mpz_t p){
    struct fp2 temp;fp2_init(&temp);

    fp_set(&temp.x0, X.x0);
    fp_sub(&temp.x1, temp.x1, X.x1, p);

    fp2_set(S, temp);

    fp2_clear(&temp);
}

// 平方剰余判定
int fp2_legendre(const struct fp2 X, const mpz_t p){
    struct fp2 check;fp2_init(&check);
    mpz_t temp;mpz_init(temp);
    struct fp2 one;fp2_init(&one);mpz_set_str(one.x0.x0, "1", 10);
    mpz_t two;mpz_init_set_str(two, "2", 10);

    // X^(p^2 - 1)/2 = 1 or -1 
    mpz_mul(temp, p, p);
    mpz_sub(temp, temp, one.x0.x0);
    mpz_cdiv_q(temp, temp, two);
    fp2_pow(&check, X, temp, p);
    if(fp2_cmp(check, one) == 0){
        return 1;
    }
    else{
        return 0;
    }

    mpz_clears(temp, two, NULL);
    fp2_clear(&one);fp2_clear(&check);
}

void fp2_sqrt(struct fp2 *S, const struct fp2 X, const mpz_t p, gmp_randstate_t state){
    if(fp2_legendre(X, p) == 1){
        // Tonelli-Shanks
        mpz_t temp;
        mpz_init(temp);

        mpz_t one;
        mpz_init_set_str(one, "1", 10);

        mpz_t two;
        mpz_init_set_str(two, "2", 10);

        // STEP 1 : p^2-1を2で割れるまで計算
        // p^2-1 = Q*2^s
        mpz_t s;
        mpz_init_set_str(s, "0", 10);

        mpz_t p_minus_one;
        mpz_init(p_minus_one);
        mpz_mul(p_minus_one, p, p);
        mpz_sub(p_minus_one, p_minus_one, one);

        mpz_t surplus;
        mpz_init(surplus);

        mpz_t Q;
        mpz_init_set(Q, p_minus_one);

        while(true){
            mpz_cdiv_qr(temp, surplus, Q, two);
            if(mpz_sgn(surplus) != 0){
                break;
            }
            mpz_add(s, s, one);
            mpz_set (Q,temp);
        }

        gmp_printf("p^2 - 1 = Q * 2^s = %Zd * 2^%Zd\n", Q, s);

        // STEP 2
        struct fp2 z;
        mpz_inits(z.x0.x0, z.x1.x0, NULL);
        fp2_random(&z, p, state);
        while(fp2_legendre(z, p) == 1){
            fp2_random(&z, p, state);
        }

        gmp_printf("z = %Zdw + %Zd\n",z.x1.x0, z.x0.x0);

        // STEP 3 : 初期値
        // M
        mpz_t M_0;
        mpz_init_set(M_0, s);
        gmp_printf("M_0 = %Zd, ",M_0);

        // c
        struct fp2 c_0;
        mpz_inits(c_0.x0.x0, c_0.x1.x0, NULL);
        fp2_pow(&c_0, z, Q, p);
        gmp_printf("c_0 = %Zdw + %Zd, ",c_0.x1.x0, c_0.x0.x0);

        // t
        struct fp2 t_0;
        mpz_inits(t_0.x0.x0, t_0.x1.x0, NULL);
        fp2_pow(&t_0, X, Q, p);
        gmp_printf("t_0 = %Zdw + %Zd, ",t_0.x1.x0, t_0.x0.x0);

        // R
        mpz_t index;
        mpz_init(index);
        mpz_add(index, Q, one);
        mpz_cdiv_q (index, index, two);
        

        struct fp2 R_0;
        mpz_inits(R_0.x0.x0, R_0.x1.x0, NULL);
        fp2_pow(&R_0, X, index, p);
        gmp_printf("R_0 = %Zdw + %Zd\n",R_0.x1.x0, R_0.x0.x0);

        // その他
        struct fp2 floor;
        mpz_inits(floor.x0.x0, floor.x1.x0, NULL);

        mpz_t index2;
        mpz_init(index2);

        // // STEP 4 : ループ
        while(mpz_cmp(t_0.x0.x0, one) != 0 || mpz_sgn(t_0.x1.x0) != 0){
            // M_i+1
            mpz_t j;
            mpz_init(j);

            for(mpz_set(j,one); mpz_cmp(j,M_0) < 0; mpz_add(j,j,one)){
                mpz_t j2;
                mpz_init(j2);

                mpz_powm(j2, two, j, p_minus_one);

                struct fp2 ch;
                mpz_inits(ch.x0.x0, ch.x1.x0, NULL);
                fp2_pow(&ch, t_0, j2, p);

                if(mpz_cmp(ch.x0.x0, one) == 0 && mpz_sgn(ch.x1.x0) == 0){
                    mpz_set(index, M_0);
                    mpz_set(index2, M_0);
                    mpz_set(M_0, j);
                    break;

                    mpz_clears(j, j2, ch, NULL);
                }
            }

            gmp_printf("M_0 = %Zd, ",M_0);

            // c_i+1
            mpz_set(floor.x0.x0,c_0.x0.x0);
            mpz_set(floor.x1.x0,c_0.x1.x0);
            mpz_sub(index, index, M_0);
            mpz_powm(index, two, index, p);
            fp2_pow(&c_0, floor, index, p);

            gmp_printf("c_0 = %Zdw + %Zd, ",c_0.x1.x0, c_0.x0.x0);

            // t_i+1
            fp2_mul(&t_0, t_0, c_0, p);

            gmp_printf("t_0 = %Zdw + %Zd, ",t_0.x1.x0, t_0.x0.x0);

            // R_i+1
            struct fp2 index3;
            mpz_inits(index3.x0.x0, index3.x1.x0, NULL);

            mpz_sub(index2, index2, M_0);
            mpz_sub(index2, index2, one);
            mpz_powm(index2, two, index2, p);
            fp2_pow(&index3, floor, index2, p);
            fp2_mul(&R_0, R_0, index3, p);

            gmp_printf("R_0 = %Zdw + %Zd\n",R_0.x1.x0, R_0.x0.x0);

            fp2_clear(&index3);

        }
        mpz_set(S->x0.x0, R_0.x0.x0);
        mpz_set(S->x1.x0, R_0.x1.x0);

        // printf("rresidue\n");
        mpz_clears(temp, one, two, s, p_minus_one, surplus, Q, M_0, index, index2, NULL);
        fp2_clear(&z);
        fp2_clear(&c_0);
        fp2_clear(&t_0);
        fp2_clear(&R_0);
        fp2_clear(&floor);
        
    }
    else if(fp2_legendre(X, p) == 0){
        printf("nonresidue\n");
    }
}