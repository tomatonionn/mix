#include "../miller_header.h"

// f(x) = x^2 + 1
// f(α) = α^2 + 1 = 0
// g(x) = x^3 - (α + 1)
// g(β) = β^3 - (α + 1)

// 初期化
void fp6_init(struct fp6 *X){
    fp2_init(&X->x0);
    fp2_init(&X->x1);
    fp2_init(&X->x2);
}

// 解放
void fp6_clear(struct fp6 *X){
    fp2_clear(&X->x0);
    fp2_clear(&X->x1);
    fp2_clear(&X->x2);
}

// 表示
void fp6_printf(const struct fp6 X){
    printf("(");
    fp2_printf(X.x0);
    printf(", ");
    fp2_printf(X.x1);
    printf(", ");
    fp2_printf(X.x2);
    printf(")");
}

// 代入
void fp6_set(struct fp6 *S, struct fp6 temp){
    fp2_set(&S->x0, temp.x0);
    fp2_set(&S->x1, temp.x1);
    fp2_set(&S->x2, temp.x2);
}

// 生成
void fp6_random(struct fp6 *X, const mpz_t p, gmp_randstate_t state){
    fp2_random(&X->x0, p, state);
    fp2_random(&X->x1, p, state);
    fp2_random(&X->x2, p, state);
}

// 比較
int fp6_cmp(const struct fp6 X, const struct fp6 Y){
    if(fp2_cmp(X.x0, Y.x0) == 0 && fp2_cmp(X.x1, Y.x1) == 0 && fp2_cmp(X.x2, Y.x2) == 0){
        return 0;
    }
    else{
        return 1;
    }
}

// (α + 1)X
void a1_xi(struct fp2 *S, const struct fp2 X, const mpz_t p){
    struct fp2 temp;
    fp2_init(&temp);
    fp_sub(&temp.x0, X.x0, X.x1, p);
    fp_add(&temp.x1, X.x0, X.x1, p);
    fp2_set(S, temp);

    fp2_clear(&temp);
}

// 負数
void fp6_neg(struct fp6 *S, struct fp6 X, const mpz_t p){
    fp2_neg(&S->x0, X.x0, p);
    fp2_neg(&S->x1, X.x1, p);
    fp2_neg(&S->x2, X.x2, p);
}

// 加算
void fp6_add(struct fp6 *S, const struct fp6 X, const struct fp6 Y, const mpz_t p){
    struct fp6 temp;fp6_init(&temp);

    fp2_add(&temp.x0, X.x0, Y.x0, p);
    fp2_add(&temp.x1, X.x1, Y.x1, p);
    fp2_add(&temp.x2, X.x2, Y.x2, p);

    fp6_set(S, temp);

    fp6_clear(&temp);
}

// 減算
void fp6_sub(struct fp6 *S, const struct fp6 X, const struct fp6 Y, const mpz_t p){
    struct fp6 temp;fp6_init(&temp);

    fp2_sub(&temp.x0, X.x0, Y.x0, p);
    fp2_sub(&temp.x1, X.x1, Y.x1, p);
    fp2_sub(&temp.x2, X.x2, Y.x2, p);

    fp6_set(S, temp);

    fp6_clear(&temp);
}

// 乗算
void fp6_mul(struct fp6 *S, const struct fp6 X, const struct fp6 Y, const mpz_t p){
    struct fp2 T0, T1, T2, T3, T3_1, T3_2, T4, T4_1, T4_2, T5, T5_1, T5_2;
    fp2_init(&T0);fp2_init(&T1);fp2_init(&T2);
    fp2_init(&T3);fp2_init(&T3_1);fp2_init(&T3_2);
    fp2_init(&T4);fp2_init(&T4_1);fp2_init(&T4_2);
    fp2_init(&T5);fp2_init(&T5_1);fp2_init(&T5_2);

    // (X0, X1, X2)*(Y0, Y1, Y2)
    fp2_mul(&T0, X.x0, Y.x0, p);    // T0 = X0*Y0
    fp2_mul(&T1, X.x1, Y.x1, p);    // T1 = X1*Y1
    fp2_mul(&T2, X.x2, Y.x2, p);    // T2 = X2*Y2

    fp2_add(&T3_1, X.x1, X.x2, p);
    fp2_add(&T3_2, Y.x1, Y.x2, p);
    fp2_mul(&T3, T3_1, T3_2, p);    // T3 = (X1 + X2)*(Y1 + Y2)

    fp2_add(&T4_1, X.x0, X.x1, p);
    fp2_add(&T4_2, Y.x0, Y.x1, p);
    fp2_mul(&T4, T4_1, T4_2, p);    // T4 = (X0 + X1)*(Y0 + Y1)

    fp2_add(&T5_1, X.x0, X.x2, p);
    fp2_add(&T5_2, Y.x0, Y.x2, p);
    fp2_mul(&T5, T5_1, T5_2, p);    // T5 = (X0 + X2)*(Y0 + Y2)

    struct fp6 temp;fp6_init(&temp);

    fp2_sub(&temp.x0, T3, T1, p);
    fp2_sub(&temp.x0, temp.x0, T2, p);
    a1_xi(&temp.x0, temp.x0, p);
    fp2_add(&temp.x0, temp.x0, T0, p);  // S0 = T0 + (α + 1)*(T3 - T1 - T2)

    a1_xi(&temp.x1, T2, p);
    fp2_add(&temp.x1, temp.x1, T4, p);
    fp2_sub(&temp.x1, temp.x1, T0, p);
    fp2_sub(&temp.x1, temp.x1, T1, p);  // S1 = T4 - T0 - T1 + (α + 1)*T2

    fp2_add(&temp.x2, T1, T5, p);
    fp2_sub(&temp.x2, temp.x2, T0, p);
    fp2_sub(&temp.x2, temp.x2, T2, p);  // S2 = T5 - T0 - T2 + T1

    fp6_set(S, temp);

    fp2_clear(&T0);fp2_clear(&T1);fp2_clear(&T2);
    fp2_clear(&T3);fp2_clear(&T3_1);fp2_clear(&T3_2);
    fp2_clear(&T4);fp2_clear(&T4_1);fp2_clear(&T4_2);
    fp2_clear(&T5);fp2_clear(&T5_1);fp2_clear(&T5_2);
    fp6_clear(&temp);
}

// ２乗算
void fp6_square(struct fp6 *S, const struct fp6 X, const mpz_t p){
    struct fp2 T1, T2, T3, T4, T5, T6;
    fp2_init(&T1);fp2_init(&T2);fp2_init(&T3);
    fp2_init(&T4);fp2_init(&T5);fp2_init(&T6);

    // (X0, X1, X2)^2
    fp2_add(&T1, X.x1, X.x1, p);    // T1 = 2*x1
    fp2_square(&T2, X.x0, p);       // T2 = X0^2
    fp2_square(&T3, X.x2, p);       // T3 = X2^2
    fp2_mul(&T4, T1, X.x2, p);      // T4 = T1*X2
    fp2_mul(&T5, T1, X.x0, p);      // T5 = T1*X0
    fp2_add(&T6, X.x0, X.x1, p);    
    fp2_add(&T6, T6, X.x2, p);
    fp2_square(&T6, T6, p);         // T6 = (X0 + X1 + X2)^2

    struct fp6 temp;fp6_init(&temp);

    a1_xi(&temp.x0, T4, p);
    fp2_add(&temp.x0, temp.x0, T2, p);  // S0 = T2 + (α + 1)*T4
    a1_xi(&temp.x1, T3, p);
    fp2_add(&temp.x1, temp.x1, T5, p);  // S1 = T5 + (α + 1)*T3
    fp2_sub(&temp.x2, T6, T2, p);
    fp2_sub(&temp.x2, temp.x2, T3, p);
    fp2_sub(&temp.x2, temp.x2, T4, p);
    fp2_sub(&temp.x2, temp.x2, T5, p);  // S2 = T6 - (T2 + T3 + T4 + T5)

    fp6_set(S, temp);

    fp2_clear(&T1);fp2_clear(&T2);fp2_clear(&T3);
    fp2_clear(&T4);fp2_clear(&T5);fp2_clear(&T6);
    fp6_clear(&temp);
}

void make_epsilon(mpz_t epsilon, const mpz_t p){
    mpz_t one, two, three;
    mpz_init_set_str(one, "1", 10);
    mpz_init_set_str(two, "2", 10);
    mpz_init_set_str(three, "3", 10);

    mpz_sub(epsilon, p, one);
    mpz_cdiv_q (epsilon, epsilon, three);
    mpz_powm(epsilon, two, epsilon, p);

    mpz_clears(one, two, three, NULL);
}

// Frobenius写像
void fp6_Frobenius(struct fp6 *S, const struct fp6 X, const mpz_t p){
    struct fp6 tempS;fp6_init(&tempS);
    struct fp temp;fp_init(&temp);
    struct fp factor1;fp_init(&factor1);
    mpz_t two;mpz_init_set_str(two, "2", 10);

    // factor2 = 2^(p-1/6)
    mpz_sub_ui(factor1.x0, p, 1);
    mpz_cdiv_q_ui(factor1.x0, factor1.x0, 6);
    mpz_powm(factor1.x0, two, factor1.x0, p);

    // factor2 = 2^(p-1/3)
    struct fp factor2;fp_init(&factor2);
    mpz_sub_ui(factor2.x0, p, 1);
    mpz_cdiv_q_ui(factor2.x0, factor2.x0, 3);
    mpz_powm(factor2.x0, two, factor2.x0, p);

    mpz_set(tempS.x0.x0.x0, X.x0.x0.x0);
    fp_sub(&tempS.x0.x1, tempS.x0.x1, X.x0.x1, p);
    fp_mul(&tempS.x1.x0, X.x1.x1, factor1, p);fp_neg(&tempS.x1.x0, tempS.x1.x0, p);
    fp_mul(&tempS.x1.x1, X.x1.x0, factor1, p);fp_neg(&tempS.x1.x1, tempS.x1.x1, p);
    fp_mul(&tempS.x2.x0, X.x2.x0, factor2, p);fp_neg(&tempS.x2.x0, tempS.x2.x0, p);
    fp_mul(&tempS.x2.x1, X.x2.x1, factor2, p);

    fp6_set(S, tempS);

    mpz_clear(two);
    fp_clear(&temp);fp_clear(&factor1);fp_clear(&factor2);
    fp6_clear(&tempS);
}

// 逆元
void fp6_inv(struct fp6 *S, const struct fp6 X, const mpz_t p){
    struct fp6 XP, XP2, XP3, XP4, XP5, s, T;
    fp6_init(&XP);fp6_init(&XP2);fp6_init(&XP3);fp6_init(&XP4);fp6_init(&XP5);
    fp6_init(&s);fp6_init(&T);

    fp6_Frobenius(&XP, X, p);
    fp6_Frobenius(&XP2, XP, p);
    fp6_Frobenius(&XP3, XP2, p);
    fp6_Frobenius(&XP4, XP3, p);
    fp6_Frobenius(&XP5, XP4, p);

    fp6_mul(&T, XP, XP2, p);
    fp6_mul(&T, T, XP3, p);
    fp6_mul(&T, T, XP4, p);
    fp6_mul(&T, T, XP5, p);
    fp6_mul(&s, X, T, p);
    fp_inv(&s.x0.x0, s.x0.x0, p);

    struct fp6 temp;
    fp6_init(&temp);

    fp6_mul(&temp, s, T, p);
    fp6_set(S, temp);

    fp6_clear(&XP);fp6_clear(&XP2);fp6_clear(&XP3);fp6_clear(&XP4);fp6_clear(&XP5);
    fp6_clear(&s);fp6_clear(&T);fp6_clear(&temp);
}

void fp6_pow(struct fp6 *S, const struct fp6 X, const mpz_t s, const mpz_t p){
    struct fp6 temp;
    fp6_init(&temp);
    mpz_set_str(temp.x0.x0.x0, "1", 10);
    char *scalar_binary = mpz_get_str (NULL, 2, s);
    size_t len = strlen(scalar_binary);

    for(size_t i = 0; i < len; i++){
        fp6_square(&temp, temp, p);
        if(scalar_binary[i] - '0' == 1){
            fp6_mul(&temp, temp, X, p);
        }
    }
    fp6_set(S, temp);

    free(scalar_binary); scalar_binary = NULL;
    fp6_clear(&temp);
}

int fp6_legendre(const struct fp6 X, const mpz_t p){
    struct fp6 check;
    fp6_init(&check);
    mpz_t temp;
    mpz_init(temp);
    struct fp6 one;
    fp6_init(&one);
    mpz_set_str(one.x0.x0.x0, "1", 10);
    mpz_t two;
    mpz_init_set_str(two, "2", 10);
    mpz_mul(temp, p, p);
    mpz_mul(temp, temp, temp);
    mpz_mul(temp, temp, p);
    mpz_mul(temp, temp, p);
    mpz_sub(temp, temp, one.x0.x0.x0);
    mpz_cdiv_q(temp, temp, two);
    fp6_pow(&check, X, temp, p);
    if(fp6_cmp(check, one) == 0){
        return 1;
    }
    else{
        return 0;
    }

    mpz_clears(temp, one, two, NULL);
    fp6_clear(&check);
}

void fp6_sqrt(struct fp6 *S, const struct fp6 X, const mpz_t p, gmp_randstate_t state){
    if(fp6_legendre(X, p) == 1){
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

        // STEP 2
        struct fp6 z;fp6_init(&z);
        fp6_random(&z, p, state);
        while(fp6_legendre(z, p) == 1){
            fp6_random(&z, p, state);
        }

        // STEP 3 : 初期値
        // M
        mpz_t M_0;
        mpz_init_set(M_0, s);

        // c
        struct fp6 c_0;
        fp6_init(&c_0);
        fp6_pow(&c_0, z, Q, p);

        // t
        struct fp6 t_0;
        fp6_init(&t_0);
        fp6_pow(&t_0, X, Q, p);

        // R
        mpz_t index;
        mpz_init(index);
        mpz_add(index, Q, one);
        mpz_cdiv_q (index, index, two);
        

        struct fp6 R_0;fp6_init(&R_0);
        fp6_pow(&R_0, X, index, p);

        // その他
        struct fp6 floor;fp6_init(&floor);

        mpz_t index2;
        mpz_init(index2);

        // // STEP 4 : ループ
        while(mpz_cmp(t_0.x0.x0.x0, one) != 0 || mpz_sgn(t_0.x0.x1.x0) != 0 || mpz_sgn(t_0.x1.x0.x0) != 0 || mpz_sgn(t_0.x1.x1.x0) != 0 || mpz_sgn(t_0.x2.x0.x0) != 0 || mpz_sgn(t_0.x2.x1.x0) != 0){
            // M_i+1
            mpz_t j;
            mpz_init(j);

            for(mpz_set(j,one); mpz_cmp(j,M_0) < 0; mpz_add(j,j,one)){
                mpz_t j2;
                mpz_init(j2);

                mpz_powm(j2, two, j, p_minus_one);

                struct fp6 ch;
                fp6_init(&ch);
                fp6_pow(&ch, t_0, j2, p);

                if(mpz_cmp(ch.x0.x0.x0, one) == 0 && mpz_sgn(ch.x0.x1.x0) == 0 && mpz_sgn(ch.x1.x0.x0) == 0 && mpz_sgn(ch.x1.x1.x0) == 0 && mpz_sgn(ch.x2.x0.x0) == 0 && mpz_sgn(ch.x2.x1.x0) == 0){
                    mpz_set(index, M_0);
                    mpz_set(index2, M_0);
                    mpz_set(M_0, j);
                    printf("\nhello!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
                    break;

                    mpz_clears(j, j2, ch, NULL);
                }
            }

            gmp_printf("M_0 = %Zd, ",M_0);

            // c_i+1
            fp6_set(&floor,c_0);
            mpz_sub(index, index, M_0);
            mpz_powm(index, two, index, p);
            fp6_pow(&c_0, floor, index, p);

            gmp_printf("c_0 = {(%Zd , %Zd) , (%Zd , %Zd) , (%Zd , %Zd)}\n",
            c_0.x0.x0.x0, c_0.x0.x1.x0, c_0.x1.x0.x0, c_0.x1.x1.x0, c_0.x2.x0.x0, c_0.x2.x1.x0);

            // t_i+1
            fp6_mul(&t_0, t_0, c_0, p);

            gmp_printf("t_0 = {(%Zd , %Zd) , (%Zd , %Zd) , (%Zd , %Zd)}\n",
            t_0.x0.x0.x0, t_0.x0.x1.x0, t_0.x1.x0.x0, t_0.x1.x1.x0, t_0.x2.x0.x0, t_0.x2.x1.x0);

            // R_i+1
            struct fp6 index3;
            fp6_init(&index3);

            mpz_sub(index2, index2, M_0);
            mpz_sub(index2, index2, one);
            mpz_powm(index2, two, index2, p);
            fp6_pow(&index3, floor, index2, p);
            fp6_mul(&R_0, R_0, index3, p);

            gmp_printf("R_0 = {(%Zd , %Zd) , (%Zd , %Zd) , (%Zd , %Zd)}\n",
            R_0.x0.x0.x0, R_0.x0.x1.x0, R_0.x1.x0.x0, R_0.x1.x1.x0, R_0.x2.x0.x0, R_0.x2.x1.x0);

            fp6_clear(&index3);

        }
        fp6_set(S, R_0);

        // printf("rresidue\n");
        mpz_clears(temp, one, two, s, p_minus_one, surplus, Q, M_0, index, index2, NULL);
        fp6_clear(&z);fp6_clear(&c_0);fp6_clear(&t_0);fp6_clear(&R_0);fp6_clear(&floor);
        
    }
    else if(fp6_legendre(X, p) == 0){
        printf("nonresidue\n");
    }
}