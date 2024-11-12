#include "../miller_header.h"

void efp2_init(struct efp2 *X){
    fp2_init(&X->x);
    fp2_init(&X->y);
    X->inf = 0;
}

void efp2_clear(struct efp2 *X){
    fp2_clear(&X->x);
    fp2_clear(&X->y);
}

void efp2_set(struct efp2 *X, struct efp2 Y){
    fp2_set(&X->x, Y.x);
    fp2_set(&X->y, Y.y);
    X->inf = Y.inf;
}

// y^2 = x^3 + ax + b
void efp2_random(struct efp2 A, struct fp2 b, mpz_t p,gmp_randstate_t state){
    struct fp2 temp, temp1, temp2;
    fp2_init(&temp);
    fp2_init(&temp1);
    fp2_init(&temp2);

    while(true){
        fp2_random(&A.x, p, state);
        fp2_mul(&temp1, A.x, A.x, p);
        fp2_mul(&temp1, temp1, A.x, p);
        fp2_add(&temp, temp1, temp2, p);
        fp2_add(&temp, temp, b, p);

        // 平方剰余判定
        if(fp2_legendre(temp, p) == 1){
            fp2_sqrt(&A.y, temp, p, state);
 
            // y, -y のどちらを出力するか
            gmp_randstate_t state;
            unsigned long int seed; // 乱数のシード
            mpz_t rand;
            mpz_init(rand);
            seed = time(NULL); // 現在時刻をシードにする例
            gmp_randinit_default(state); // gmp_randinit_defaultにstateを渡します
            gmp_randseed_ui(state, seed); // 乱数のシードを設定します
            mpz_urandomb(rand, state, 1);
            if(mpz_sgn(rand) == 0){
                mpz_neg(A.y.x0.x0, A.y.x0.x0);
                mpz_mod(A.y.x0.x0, A.y.x0.x0, p);
                mpz_neg(A.y.x1.x0, A.y.x1.x0);
                mpz_mod(A.y.x1.x0, A.y.x1.x0, p);
            }
            break;
        }
    }
}

// ２倍算
void efp2_ecd(struct efp2 *R, struct efp2 P, mpz_t p){
    // 無限遠点処理
    if(P.inf == 1){
        return;
    }

    // 例外処理
    if(mpz_sgn(P.y.x0.x0) == 0 && mpz_sgn(P.y.x1.x0) == 0){
        gmp_printf("error : Yp = 0\n");
        return;
    }
 
    struct fp three;
    mpz_init_set_str(three.x0, "3", 10);
    struct fp2 lambda, lambda_numerator, lambda_denominator;     // 分子lambda_numerator, 分母lambda_denominator
    mpz_inits(lambda.x0.x0, lambda.x1.x0, NULL);
    mpz_inits(lambda_numerator.x0.x0, lambda_numerator.x1.x0, NULL);
    mpz_inits(lambda_denominator.x0.x0, lambda_denominator.x1.x0, NULL);
   
    fp2_mul(&lambda_numerator, P.x, P.x, p);
    fp_mul(&lambda_numerator.x0, lambda_numerator.x0, three, p);
    fp_mul(&lambda_numerator.x1, lambda_numerator.x1, three, p);
    fp2_add(&lambda_denominator, P.y, P.y, p);
    fp2_inv(&lambda_denominator, lambda_denominator, p);
    fp2_mul(&lambda, lambda_numerator, lambda_denominator, p);
 
    struct fp2 temp_Rx;
    mpz_inits(temp_Rx.x0.x0, temp_Rx.x1.x0, NULL);
    fp2_mul(&temp_Rx, lambda, lambda, p);
    fp2_sub(&temp_Rx, temp_Rx, P.x, p);
    fp2_sub(&temp_Rx, temp_Rx, P.x, p);
 
    struct fp2 temp_Ry;
    mpz_inits(temp_Ry.x0.x0, temp_Ry.x1.x0, NULL);
    fp2_sub(&temp_Ry, P.x, temp_Rx, p);
    fp2_mul(&temp_Ry, temp_Ry, lambda, p);
    fp2_sub(&temp_Ry, temp_Ry, P.y, p);
 
    mpz_set(R->x.x0.x0, temp_Rx.x0.x0);
    mpz_set(R->x.x1.x0, temp_Rx.x1.x0);
    mpz_set(R->y.x0.x0, temp_Ry.x0.x0);
    mpz_set(R->y.x1.x0, temp_Ry.x1.x0);
}
 
// 加算
void efp2_eca(struct efp2 *R, struct efp2 P, struct efp2 Q, mpz_t p){
    if(P.inf == 1){
        fp2_set(&R->x, Q.x);
        fp2_set(&R->y, Q.y);
        R->inf = 0;
        return;
    }

    else if(Q.inf == 1){
        fp2_set(&R->x, Q.x);
        fp2_set(&R->y, Q.y);
        R->inf = 0; 
        return;
    }

    else if(mpz_cmp(P.x.x0.x0, Q.x.x0.x0) == 0 && mpz_cmp(P.x.x1.x0, Q.x.x1.x0) == 0 && mpz_cmp(P.y.x0.x0, Q.y.x0.x0) == 0 && mpz_cmp(P.y.x1.x0, Q.y.x1.x0) == 0){
        efp2_ecd(R, P, p);
        return;
    }

    else if(mpz_cmp(P.x.x0.x0, Q.x.x0.x0) == 0 && mpz_cmp(P.x.x1.x0, Q.x.x1.x0) == 0){
        R->inf = 1;
        return;
    }
    
    struct fp2 lambda, lambda_numerator, lambda_denominator;
    fp2_init(&lambda);fp2_init(&lambda_numerator);fp2_init(&lambda_denominator);
 
    fp2_sub(&lambda_numerator, Q.y, P.y, p);
    fp2_sub(&lambda_denominator, Q.x, P.x, p);
    fp2_inv(&lambda_denominator, lambda_denominator, p);
    fp2_mul(&lambda, lambda_numerator, lambda_denominator, p);
 
    struct fp2 temp_Rx;fp2_init(&temp_Rx);
    fp2_mul(&temp_Rx, lambda, lambda, p);
    fp2_sub(&temp_Rx, temp_Rx, P.x, p);
    fp2_sub(&temp_Rx, temp_Rx, Q.x, p);
 
    struct fp2 temp_Ry;fp2_init(&temp_Ry);
    fp2_sub(&temp_Ry, P.x, temp_Rx, p);
    fp2_mul(&temp_Ry, temp_Ry, lambda, p);
    fp2_sub(&temp_Ry, temp_Ry, P.y, p);
 
    fp2_set(&R->x, temp_Rx);
    fp2_set(&R->y, temp_Ry);
    R->inf = 0;

    fp2_clear(&temp_Rx);fp2_clear(&temp_Ry);fp2_clear(&lambda);fp2_clear(&lambda_numerator);fp2_clear(&lambda_denominator);
}
 
// スカラー倍算
void efp2_scm(struct efp2 R, struct efp2 P, mpz_t s, mpz_t p){
    R.inf = 1;
    char *scalar_binary = mpz_get_str (NULL, 2, s);
    size_t len = strlen(scalar_binary);
 
    for(size_t i = 0; i < len; i++){
        efp2_ecd(&R, R, p);
        if(scalar_binary[i] - '0' == 1){
            efp2_eca(&R, R, P, p);
        }
    }
}

// int main(void){
//     printf("hello");
// }