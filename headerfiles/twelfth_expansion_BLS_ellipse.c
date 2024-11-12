#include "../miller_header.h"

void efp12_init(struct efp12 *X){
    fp12_init(&X->x);
    fp12_init(&X->y);
    X->inf = 0;
}

void efp12_clear(struct efp12 *X){
    fp12_clear(&X->x);
    fp12_clear(&X->y);
}

void efp12_printf(const struct efp12 X){
    printf("(");
    fp12_printf(X.x);
    printf(", ");
    fp12_printf(X.y);
    printf(")");
}

void efp12_set(struct efp12 *X, struct efp12 Y){
    fp12_set(&X->x, Y.x);
    fp12_set(&X->y, Y.y);
    X->inf = Y.inf;
}

int efp12_cmp(const struct efp12 X, const struct efp12 Y){
    if(fp12_cmp(X.x, Y.x) == 0 && fp12_cmp(X.y, Y.y) == 0){
        return 0;
    }
    else{
        return 1;
    }
}

// y^2 = x^3 + b
void efp12_random(struct efp12 *A, struct fp b, mpz_t p, gmp_randstate_t state){
    struct fp12 temp;fp12_init(&temp);
    struct efp12 tempA;fp12_init(&tempA.x);fp12_init(&tempA.y);tempA.inf = 0;

    while(true){
        fp12_random(&tempA.x, p, state);
        fp12_mul(&temp, tempA.x, tempA.x, p);
        fp12_mul(&temp, temp, tempA.x, p);
        fp_add(&temp.x0.x0.x0, temp.x0.x0.x0, b, p);
 
        // 平方剰余判定
        if(fp12_legendre(temp, p) == 1){
            fp12_sqrt(&tempA.y, temp, p, state);
            fp12_set(&A->x, tempA.x);fp12_set(&A->y, tempA.y);

            fp12_clear(&temp);fp12_clear(&tempA.x);fp12_clear(&tempA.y);
            break;
        }
    }

}

// ２倍算
void efp12_ecd(struct efp12 *R, struct efp12 P, mpz_t p){
    // 無限遠点処理
    if(P.inf == 1){
        return;
    }

    // 例外処理 Yp = 0
    struct fp12 zero;fp12_init(&zero); 
    if(fp12_cmp(P.y, zero) == 0){
        R->inf = 1;
        return;
    }
    fp12_clear(&zero);

    struct fp12 three;fp12_init(&three);mpz_set_str(three.x0.x0.x0.x0, "3", 10);
    struct fp12 lambda, lambda_numerator, lambda_denominator;     // 分子lambda_numerator, 分母lambda_denominator
    fp12_init(&lambda);fp12_init(&lambda_numerator);fp12_init(&lambda_denominator);

    fp12_mul(&lambda_numerator, P.x, P.x, p);
    fp12_mul(&lambda_numerator, lambda_numerator, three, p);
    fp12_add(&lambda_denominator, P.y, P.y, p);
    fp12_inv(&lambda_denominator, lambda_denominator, p);
    fp12_mul(&lambda, lambda_numerator, lambda_denominator, p);
 
    struct fp12 temp_Rx;fp12_init(&temp_Rx);
    fp12_mul(&temp_Rx, lambda, lambda, p);
    fp12_sub(&temp_Rx, temp_Rx, P.x, p);
    fp12_sub(&temp_Rx, temp_Rx, P.x, p);
 
    struct fp12 temp_Ry;fp12_init(&temp_Ry);
    fp12_sub(&temp_Ry, P.x, temp_Rx, p);
    fp12_mul(&temp_Ry, temp_Ry, lambda, p);
    fp12_sub(&temp_Ry, temp_Ry, P.y, p);
 
    fp12_set(&R->x, temp_Rx);
    fp12_set(&R->y, temp_Ry);
    R->inf = 0;

    fp12_clear(&three);
    fp12_clear(&temp_Rx);fp12_clear(&temp_Ry);
    fp12_clear(&lambda);fp12_clear(&lambda_numerator);fp12_clear(&lambda_denominator);
}

// 加算
void efp12_eca(struct efp12 *R, struct efp12 P, struct efp12 Q, mpz_t p){
    //　無限遠点処理
    if(P.inf == 1){
        fp12_set(&R->x, Q.x);
        fp12_set(&R->y, Q.y);
        R->inf = 0;
        return;
    }

    else if(Q.inf == 1){
        fp12_set(&R->x, P.x);
        fp12_set(&R->y, P.y);
        R->inf = 0;
        return;
    }

    // 同点処理 P = Q
    else if(efp12_cmp(P, Q) == 0){
        efp12_ecd(R, P, p);
        return;
    }

    // 例外処理 Xp = Xq
    else if(fp12_cmp(P.x, Q.x) == 0){
        R->inf = 1;
        return;
    }
    
    struct fp12 lambda, lambda_numerator, lambda_denominator;
    fp12_init(&lambda);fp12_init(&lambda_numerator);fp12_init(&lambda_denominator);

    fp12_sub(&lambda_numerator, Q.y, P.y, p);
    fp12_sub(&lambda_denominator, Q.x, P.x, p);
    fp12_inv(&lambda_denominator, lambda_denominator, p);
    fp12_mul(&lambda, lambda_numerator, lambda_denominator, p);
 
    struct fp12 temp_Rx;fp12_init(&temp_Rx);
    fp12_mul(&temp_Rx, lambda, lambda, p);
    fp12_sub(&temp_Rx, temp_Rx, P.x, p);
    fp12_sub(&temp_Rx, temp_Rx, Q.x, p);
 
    struct fp12 temp_Ry;fp12_init(&temp_Ry);
    fp12_sub(&temp_Ry, P.x, temp_Rx, p);
    fp12_mul(&temp_Ry, temp_Ry, lambda, p);
    fp12_sub(&temp_Ry, temp_Ry, P.y, p);
 
    fp12_set(&R->x, temp_Rx);
    fp12_set(&R->y, temp_Ry);
    R->inf = 0;

    fp12_clear(&temp_Rx);fp12_clear(&temp_Ry);fp12_clear(&lambda);fp12_clear(&lambda_numerator);fp12_clear(&lambda_denominator);
}
 

// スカラー倍算
void efp12_scm(struct efp12 *R, struct efp12 P, mpz_t s, mpz_t p){
    struct efp12 tempR;fp12_init(&tempR.x);fp12_init(&tempR.y);tempR.inf = 1;
    char *scalar_binary = mpz_get_str (NULL, 2, s);
    size_t len = strlen(scalar_binary);

    for(size_t i = 0; i < len; i++){
        efp12_ecd(&tempR, tempR, p);
        if(scalar_binary[i] - '0' == 1){
            efp12_eca(&tempR, tempR, P, p);
        }
    }

    fp12_set(&R->x, tempR.x);fp12_set(&R->y, tempR.y);R->inf = tempR.inf;

    free(scalar_binary); scalar_binary = NULL;
    fp12_clear(&tempR.x);fp12_clear(&tempR.y);
}