#include "../miller_header.h"

// f(x) = x^2 + 1
// f(α) = α^2 + 1 = 0
// α^2 = -1
// g(x) = x^3 - (α + 1)
// g(β) = β^3 - (α + 1) = 0
// β^3 = (α + 1)
// h(x) = x^2 - β
// f(γ) = γ^2 - β = 0
// γ^2 = β

void fp12_init(struct fp12 *X){
    fp6_init(&X->x0);
    fp6_init(&X->x1);
}

void fp12_clear(struct fp12 *X){
    fp6_clear(&X->x0);
    fp6_clear(&X->x1);
}

void fp12_printf(const struct fp12 X){
    printf("(");
    fp6_printf(X.x0);
    printf(", ");
    fp6_printf(X.x1);
    printf(")");
}

void fp12_set(struct fp12 *S, struct fp12 X){
    fp6_set(&S->x0, X.x0);
    fp6_set(&S->x1, X.x1);
}

void fp12_random(struct fp12 *X, const mpz_t p, gmp_randstate_t state){
    fp6_random(&X->x0, p, state);
    fp6_random(&X->x1, p, state);
}

int fp12_cmp(const struct fp12 X, const struct fp12 Y){
    if(fp6_cmp(X.x0, Y.x0) == 0 && fp6_cmp(X.x1, Y.x1) == 0){
        return 0;
    }
    else{
        return 1;
    }
}

void mulb(struct fp6 *S, const struct fp6 X, const mpz_t p){
    struct fp6 temp;
    fp6_init(&temp);

    fp2_set(&temp.x1, X.x0);
    fp2_set(&temp.x2, X.x1);
    a1_xi(&temp.x0, X.x2, p);

    fp6_set(S, temp);

    fp6_clear(&temp);
}

// 負数
void fp12_neg(struct fp12 *S, struct fp12 X, const mpz_t p){
    fp6_neg(&S->x0, X.x0, p);
    fp6_neg(&S->x1, X.x1, p);
}

void fp12_add(struct fp12 *S, const struct fp12 X, const struct fp12 Y, const mpz_t p){
    struct fp12 temp;
    fp12_init(&temp);

    fp6_add(&temp.x0, X.x0, Y.x0, p);
    fp6_add(&temp.x1, X.x1, Y.x1, p);

    fp12_set(S, temp);

    fp12_clear(&temp);
}

void fp12_sub(struct fp12 *S, const struct fp12 X, const struct fp12 Y, const mpz_t p){
    struct fp12 temp;
    fp12_init(&temp);

    fp6_sub(&temp.x0, X.x0, Y.x0, p);
    fp6_sub(&temp.x1, X.x1, Y.x1, p);

    fp12_set(S, temp);

    fp12_clear(&temp);
}

void fp12_mul(struct fp12 *S, const struct fp12 X, const struct fp12 Y, const mpz_t p){
    struct fp6 T1, T2, T2b, T3, T4;
    fp6_init(&T1);fp6_init(&T2);fp6_init(&T2b);fp6_init(&T3);fp6_init(&T4);

    fp6_mul(&T1, X.x0, Y.x0, p);
    fp6_mul(&T2, X.x1, Y.x1, p);
    mulb(&T2b, T2, p);
    fp6_add(&T3, X.x0, X.x1, p);
    fp6_add(&T4, Y.x0, Y.x1, p);

    struct fp12 temp;
    fp12_init(&temp);

    fp6_add(&temp.x0, T1, T2b, p);
    fp6_mul(&temp.x1, T3, T4, p);
    fp6_sub(&temp.x1, temp.x1, T1, p);
    fp6_sub(&temp.x1, temp.x1, T2, p);

    fp12_set(S, temp);

    fp6_clear(&T1);fp6_clear(&T2);fp6_clear(&T2b);fp6_clear(&T3);fp6_clear(&T4);
    fp12_clear(&temp);
}

void fp12_square(struct fp12 *S, const struct fp12 X, const mpz_t p){
    struct fp12 temp;
    fp12_init(&temp);
    struct fp6 T1, T2;
    fp6_init(&T1);fp6_init(&T2);

    fp6_square(&T1, X.x0, p);
    fp6_square(&T2, X.x1, p);
    mulb(&T2, T2, p);
    fp6_add(&temp.x0, T1, T2, p);
    fp6_mul(&temp.x1, X.x0, X.x1, p);
    fp6_add(&temp.x1, temp.x1, temp.x1, p);

    fp12_set(S, temp);

    fp6_clear(&T1);fp6_clear(&T2);
    fp12_clear(&temp);
}

void fp12_Frobenius(struct fp12 *S, const struct fp12 X, const mpz_t p){
    // S = [{x0},{(y0),(y1),(y2)}]
    // x0, x1 = fp6_Frobenius(X)
    // y0, y1, y2 = x1 * (1 + α)^{(p - 1) / 6}
    struct fp6 x0, x1;fp6_init(&x0);fp6_init(&x1);
    struct fp2 y0, y1, y2;fp2_init(&y0);fp2_init(&y1);fp2_init(&y2);
    struct fp2 tmp;fp2_init(&tmp);mpz_set_ui(tmp.x0.x0, 1);mpz_set_ui(tmp.x1.x0, 1);   // tmp = (1 + α)
    mpz_t s;mpz_init(s);mpz_sub_ui(s, p, 1);mpz_cdiv_q_ui(s, s, 6);     // s = (p - 1) / 6
    fp2_pow(&tmp, tmp, s, p);       // tmp = (1 + α)^{(p - 1) / 6}

    fp6_Frobenius(&x0, X.x0, p);
    fp6_Frobenius(&x1, X.x1, p);
    fp2_set(&y0, x1.x0);
    fp2_set(&y1, x1.x1);
    fp2_set(&y2, x1.x2);

    fp2_mul(&y0, y0, tmp, p);
    fp2_mul(&y1, y1, tmp, p);
    fp2_mul(&y2, y2, tmp, p);

    fp6_set(&S->x0, x0);
    fp2_set(&S->x1.x0, y0);fp2_set(&S->x1.x1, y1);fp2_set(&S->x1.x2, y2);
}


void fp12_inv(struct fp12 *S, const struct fp12 X, const mpz_t p){
    struct fp12 X2;fp12_init(&X2);
    fp6_set(&X2.x0, X.x0);
    fp6_sub(&X2.x1, X2.x1, X.x1, p);

    struct fp6 Y,Yinv, temp1, temp2, b;
    fp6_init(&Y);fp6_init(&temp1);fp6_init(&temp2);fp6_init(&b);fp6_init(&Yinv);
    mpz_set_str(b.x1.x0.x0, "1", 10);
    fp6_square(&temp1, X.x0, p);
    fp6_square(&temp2, X.x1, p);
    fp6_mul(&temp2, temp2, b, p);
    fp6_sub(&Y, temp1, temp2, p);
    fp6_inv(&Yinv, Y, p);

    struct fp12 Y2, temp;fp12_init(&Y2);fp12_init(&temp);
    fp6_set(&Y2.x0, Yinv);
    fp12_mul(&temp, X2, Y2, p);

    fp12_set(S, temp);

    fp6_clear(&Y);fp6_clear(&Yinv);fp6_clear(&temp1);fp6_clear(&temp2);
    fp12_clear(&X2);fp12_clear(&Y2);fp12_clear(&temp);
}

void fp12_pow(struct fp12 *S, const struct fp12 X, const mpz_t s, const mpz_t p){
    struct fp12 temp;
    fp12_init(&temp);
    mpz_set_str(temp.x0.x0.x0.x0, "1", 10);
    char *scalar_binary = mpz_get_str(NULL, 2, s);
    size_t len = strlen(scalar_binary);

    for(size_t i = 0; i < len; i++){
        fp12_square(&temp, temp, p);
        if(scalar_binary[i] - '0' == 1){
            fp12_mul(&temp, temp, X, p);
        }
    }
    fp12_set(S, temp);

    free(scalar_binary); scalar_binary = NULL;
    fp12_clear(&temp);
}

int fp12_legendre(const struct fp12 X, const mpz_t p){
    struct fp12 check;
    fp12_init(&check);
    mpz_t temp;
    mpz_init(temp);
    struct fp12 one;
    fp12_init(&one);
    mpz_set_str(one.x0.x0.x0.x0, "1", 10);

    mpz_mul(temp, p, p);
    mpz_mul(temp, temp, temp);
    mpz_mul(temp, temp, temp);
    mpz_mul(temp, temp, p);
    mpz_mul(temp, temp, p);
    mpz_mul(temp, temp, p);
    mpz_mul(temp, temp, p);
    mpz_sub_ui(temp, temp, 1);
    mpz_cdiv_q_ui(temp, temp, 2);
    fp12_pow(&check, X, temp, p);
    if(fp12_cmp(check, one) == 0){
        return 1;
    }
    else{
        return 0;
    }

    mpz_clear(temp);
    fp12_clear(&check);
    fp12_clear(&one);
}

void fp12_sqrt(struct fp12 *S, const struct fp12 X, const mpz_t p, gmp_randstate_t state){
    if(fp12_legendre(X, p) == 1){
        // printf("rresidue!\n");
        
        mpz_t check;mpz_init(check);
        mpz_t three;mpz_init_set_str(three, "3", 10);
        mpz_t four;mpz_init_set_str(four, "4", 10);
        mpz_mul(check, p, p);
        mpz_mul(check, check, check);
        mpz_mul(check, check, check);
        mpz_mul(check, check, p);
        mpz_mul(check, check, p);
        mpz_mul(check, check, p);
        mpz_mul(check, check, p);
        mpz_mod(check, check, four);
        if(mpz_cmp(check, three) == 0){
            mpz_t index;
            mpz_init(index);
            mpz_t one;
            mpz_init_set_str(one, "1", 10);
            struct fp12 temp;
            fp12_init(&temp);

            mpz_add(index, p, one);
            mpz_cdiv_q(index, index, four);
            fp12_pow(&temp, X, index, p);
            fp12_set(S, temp);

            mpz_clears(index, one, NULL);
            fp12_clear(&temp);
        }
        else{
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
            mpz_mul(p_minus_one, p_minus_one, p_minus_one);
            mpz_mul(p_minus_one, p_minus_one, p_minus_one);
            mpz_mul(p_minus_one, p_minus_one, p);
            mpz_mul(p_minus_one, p_minus_one, p);
            mpz_mul(p_minus_one, p_minus_one, p);
            mpz_mul(p_minus_one, p_minus_one, p);
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
            struct fp12 z;
            fp12_init(&z);
            fp12_random(&z, p, state);
            while(fp12_legendre(z, p) == 1){
                fp12_random(&z, p, state);
            }

            // STEP 3 : 初期値
            // M
            mpz_t M_0;
            mpz_init_set(M_0, s);

            // c
            struct fp12 c_0;
            fp12_init(&c_0);
            fp12_pow(&c_0, z, Q, p);

            // t
            struct fp12 t_0;
            fp12_init(&t_0);
            fp12_pow(&t_0, X, Q, p);
           
            // R
            mpz_t index;
            mpz_init(index);
            mpz_add(index, Q, one);
            mpz_cdiv_q (index, index, two);


            struct fp12 R_0;
            fp12_init(&R_0);
            fp12_pow(&R_0, X, index, p);
            
            // その他
            struct fp12 floor;
            fp12_init(&floor);

            mpz_t index2;
            mpz_init(index2);

            // // STEP 4 : ループ
            while(mpz_cmp(t_0.x0.x0.x0.x0, one) != 0 || mpz_sgn(t_0.x0.x0.x1.x0) != 0 || mpz_sgn(t_0.x0.x1.x0.x0) != 0 || mpz_sgn(t_0.x0.x1.x1.x0) != 0 || mpz_sgn(t_0.x0.x2.x0.x0) != 0 || mpz_sgn(t_0.x0.x2.x1.x0) != 0
                    || mpz_sgn(t_0.x1.x0.x0.x0) != 0 || mpz_sgn(t_0.x1.x0.x1.x0) != 0 || mpz_sgn(t_0.x1.x1.x0.x0) != 0 || mpz_sgn(t_0.x1.x1.x1.x0) != 0 || mpz_sgn(t_0.x1.x2.x0.x0) != 0 || mpz_sgn(t_0.x1.x2.x1.x0) != 0){
                // M_i+1
                mpz_t j;mpz_init(j);

                for(mpz_set(j,one); mpz_cmp(j,M_0) < 0; mpz_add(j,j,one)){
                    mpz_t j2;
                    mpz_init(j2);

                    mpz_powm(j2, two, j, p_minus_one);

                    struct fp12 ch;
                    fp12_init(&ch);
                    fp12_pow(&ch, t_0, j2, p);

                    if(mpz_cmp(ch.x0.x0.x0.x0, one) == 0 && mpz_sgn(ch.x0.x0.x1.x0) == 0 && mpz_sgn(ch.x0.x1.x0.x0) == 0 && mpz_sgn(ch.x0.x1.x1.x0) == 0 && mpz_sgn(ch.x0.x2.x0.x0) == 0 && mpz_sgn(ch.x0.x2.x1.x0) == 0
                         && mpz_sgn(ch.x1.x0.x0.x0) == 0 && mpz_sgn(ch.x1.x0.x1.x0) == 0 && mpz_sgn(ch.x1.x1.x0.x0) == 0 && mpz_sgn(ch.x1.x1.x1.x0) == 0 && mpz_sgn(ch.x1.x2.x0.x0) == 0 && mpz_sgn(ch.x1.x2.x1.x0) == 0){
                        mpz_set(index, M_0);
                        mpz_set(index2, M_0);
                        mpz_set(M_0, j);

                        mpz_clear(j2);
                        fp12_clear(&ch);
                        break;
                    }
                }
                mpz_clear(j);
                

                // gmp_printf("M_0 = %Zd, ",M_0);

                // c_i+1
                fp12_set(&floor,c_0);
                mpz_sub(index, index, M_0);
                mpz_powm(index, two, index, p);
                fp12_pow(&c_0, floor, index, p);

                // t_i+1
                fp12_mul(&t_0, t_0, c_0, p);

                // R_i+1
                struct fp12 index3;
                fp12_init(&index3);

                mpz_sub(index2, index2, M_0);
                mpz_sub(index2, index2, one);
                mpz_powm(index2, two, index2, p);
                fp12_pow(&index3, floor, index2, p);
                fp12_mul(&R_0, R_0, index3, p);

                fp12_clear(&index3);
            }
            fp12_set(S, R_0);

            mpz_clears(temp, one, two, s, p_minus_one, surplus, Q, M_0, index, index2, NULL);
            fp12_clear(&z);
            fp12_clear(&c_0);
            fp12_clear(&t_0);
            fp12_clear(&R_0);
            fp12_clear(&floor);
            

        mpz_clears(check, three, four, NULL);
        }
    }
    else if(fp12_legendre(X, p) == 0){
        printf("nonresidue\n");
    }
}