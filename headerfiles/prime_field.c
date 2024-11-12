#include "../miller_header.h"

extern int Mcounter;

// 初期化
void fp_init(struct fp *X){
    mpz_init(X->x0);
}

// 解放
void fp_clear(struct fp *X){
    mpz_clear(X->x0);
}

// 表示
void fp_printf(const struct fp X){
    gmp_printf("%Zd",X.x0);
}

// 代入
void fp_set(struct fp *S, const struct fp X){
    mpz_set(S->x0, X.x0);
}

// シード設定
void make_state(gmp_randstate_t state){
    gmp_randinit_mt (state);
    gmp_randseed_ui(state, time(NULL));
}

// 元の生成
void fp_random(struct fp *X, const mpz_t p, gmp_randstate_t state){
    mpz_t temp;mpz_init(temp);
    mpz_urandomm(temp, state, p);
    mpz_set(X->x0, temp);
    mpz_clear(temp);
}

// 比較
int fp_cmp(const struct fp X, const struct fp Y){
    if(mpz_cmp(X.x0, Y.x0) == 0){
        return 0;
    }
    else{
        return 1;
    }
}

// 負数
void fp_neg(struct fp *S, struct fp X, const mpz_t p){
    mpz_neg(S->x0, X.x0);
    mpz_mod(S->x0, S->x0, p);
}

// 和 S = X + Y
void fp_add(struct fp *S, const struct fp X, const struct fp Y, const mpz_t p){
    mpz_t temp;mpz_init(temp);
    mpz_add(temp, X.x0, Y.x0);
    mpz_mod(temp, temp, p);
    mpz_set(S->x0, temp);
    mpz_clear(temp);
}

// 差 S = X - Y
void fp_sub(struct fp *S, const struct fp X, const struct fp Y, const mpz_t p){
    mpz_t temp;mpz_init(temp);
    mpz_sub(temp, X.x0, Y.x0);
    mpz_mod(temp, temp, p);
    mpz_set(S->x0, temp);
    mpz_clear(temp);
}

// 積 S = X * Y
void fp_mul(struct fp *S, const struct fp X, const struct fp Y, const mpz_t p){
    mpz_t temp;mpz_init(temp);
    mpz_mul(temp, X.x0, Y.x0);
    mpz_mod(temp, temp, p);
    mpz_mod(S->x0, temp, p);
    mpz_clear(temp);
    Mcounter ++;
}

// 逆数 S = 1 / X
void fp_inv(struct fp *S, const struct fp X, const mpz_t p){
    mpz_t temp;
    mpz_init(temp);
    mpz_invert(temp, X.x0, p);
    mpz_set(S->x0, temp);
    mpz_clear(temp);
    Mcounter += 25;
}

// 冪乗 S = X ** s
void fp_pow(struct fp *S, const struct fp X, const mpz_t s, const mpz_t p){
    mpz_t temp;mpz_init(temp);
    mpz_powm(temp, X.x0, s, p);
    mpz_set(S->x0, temp);
    mpz_clear(temp);
}

// 平方剰余判定
int fp_legendre(const struct fp X, const mpz_t p){
    return mpz_legendre(X.x0, p);
}

// 1.2.7 平方根計算 b = √a
void fp_sqrt(struct fp *b, struct fp a, mpz_t p, gmp_randstate_t state){
    if(fp_legendre(a, p) == 1){
        mpz_t temp;
        mpz_init(temp);
        mpz_t zero;
        mpz_init_set_str(zero, "0", 10);
        mpz_t one;
        mpz_init_set_str(one, "1", 10);
        mpz_t two;
        mpz_init_set_str(two, "2", 10);
        mpz_t three;
        mpz_init_set_str(three, "3", 10);
        mpz_t four;
        mpz_init_set_str(four, "4", 10);
        mpz_mod(temp, p, four);
        // p が 4 で割って 3 余る素数の場合
        if(!mpz_cmp(temp,three)){
            mpz_t p_plus_one;
            mpz_init(p_plus_one);
            mpz_add(p_plus_one, p, one);

            mpz_t index;
            mpz_init(index);
            mpz_cdiv_q (index, p_plus_one, four);

            mpz_powm(temp, a.x0, index, p);
            mpz_set(b->x0, temp);
            
            mpz_clear(p_plus_one);
            mpz_clear(index);
        }

        // p が 4 で割って 1 余る素数の場合 Tonelli-Shanks
        else if(!mpz_cmp(temp,one)){
            // STEP 1 : p-1を2で割れるまで計算
            // p-1 = Q*2^s
            mpz_t s;
            mpz_init_set_str(s, "0", 10);

            mpz_t p_minus_one;
            mpz_init(p_minus_one);
            mpz_sub(p_minus_one, p, one);

            mpz_t surplus;
            mpz_init(surplus);

            mpz_t Q;
            mpz_init_set(Q, p_minus_one);

            while(true){
                mpz_cdiv_qr(temp, surplus, Q, two);
                if(mpz_cmp (surplus, zero) != 0){
                    break;
                }
                mpz_add(s, s, one);
                mpz_set (Q,temp);
            }

            gmp_printf("%Zd - 1 = %Zd * 2^%Zd\n",p,Q,s);

            // STEP 2
            struct fp z;
            mpz_init(z.x0);
            fp_random(&z, p, state);
            while(fp_legendre(z, p) != -1){
                fp_random(&z, p, state);
            }

            gmp_printf("z = %Zd\n",z);

            // STEP 3 : 初期値
            // M
            mpz_t M_0;
            mpz_init_set(M_0, s);

            // c
            mpz_t c_0;
            mpz_init(c_0);
            mpz_powm(c_0, z.x0, Q, p);

            // t
            mpz_t t_0;
            mpz_init(t_0);
            mpz_powm(t_0, a.x0, Q, p);

            // R
            mpz_t index;
            mpz_init(index);
            mpz_add(index, Q, one);
            mpz_cdiv_q (index, index, two);

            mpz_t R_0;
            mpz_init(R_0);
            mpz_powm(R_0, a.x0, index, p);

            // その他
            mpz_t floor;
            mpz_init(floor);

            mpz_t index2;
            mpz_init(index2);

            gmp_printf("M_0 = %Zd, ",M_0);
            gmp_printf("c_0 = %Zd, ",c_0);
            gmp_printf("t_0 = %Zd, ",t_0);
            gmp_printf("R_0 = %Zd\n",R_0);

            // // STEP 4 : ループ

            while(mpz_cmp(t_0, one) != 0){
                // M_i+1
                mpz_t j;
                mpz_init(j);

                for(mpz_set(j,one); mpz_cmp(j,M_0) < 0; mpz_add(j,j,one)){
                    mpz_t j2;
                    mpz_init(j2);
                    mpz_powm(j2,two,j,p);

                    mpz_t ch;
                    mpz_init(ch);
                    mpz_powm(ch, t_0, j2, p);

                    if(mpz_cmp(ch,one) == 0){
                        mpz_set(index, M_0);
                        mpz_set(index2, M_0);
                        mpz_set(M_0, j);
                        mpz_clear(j);
                        mpz_clear(j2);
                        mpz_clear(ch);
                        break;
                    }
                }                
                
                // c_i+1
                mpz_set(floor,c_0);
                mpz_sub(index, index, M_0);
                mpz_powm(index, two, index, p);
                mpz_powm(c_0, floor, index, p);

                // t_i+1
                mpz_mul(t_0, t_0, c_0);
                mpz_mod(t_0,t_0,p);

                // R_i+1
                mpz_sub(index2, index2, M_0);
                mpz_sub(index2, index2, one);
                mpz_powm(index2, two, index2, p);
                mpz_powm(index2, floor, index2, p);
                mpz_mul(R_0, R_0, index2);
                mpz_mod(R_0, R_0, p);
            }
            mpz_set(b->x0, R_0);
            mpz_clear(s);
            mpz_clear(p_minus_one);
            mpz_clear(surplus);
            mpz_clear(Q);
            fp_clear(&z);
            mpz_clear(M_0);
            mpz_clear(c_0);
            mpz_clear(t_0);
            mpz_clear(R_0);
            mpz_clear(index);
            mpz_clear(index2);
            mpz_clear(floor);

        }
        mpz_clear(temp);
        mpz_clear(zero);
        mpz_clear(one);
        mpz_clear(two);
        mpz_clear(three);
        mpz_clear(four);
    }
    else if(fp_legendre(a, p) == 0){
        mpz_set_str(b->x0, "0", 10);
    }
    else{
        printf("no sqrt\n");
        mpz_set_str(b->x0, "-1", 10);
    }
}