#include "../miller_header.h"

int Mcounter = 0;

void efp_random(struct efp *A, struct fp b, mpz_t p, gmp_randstate_t state){
    struct fp temp, temp1, temp2;
    struct efp tempA;
    mpz_inits(temp.x0, temp1.x0, temp2.x0, tempA.x.x0, tempA.y.x0, NULL);
    tempA.inf = 0;
 
    while(true){
        fp_random(&tempA.x, p, state);
        fp_mul(&temp1, tempA.x, tempA.x, p);
        fp_mul(&temp1, temp1, tempA.x, p);
        fp_add(&temp, temp1, b, p);
 
        // 平方剰余判定
        if(fp_legendre(temp, p) == 1){
            fp_sqrt(&tempA.y, temp, p, state);
            break;
        }
    }
    mpz_set(A->x.x0, tempA.x.x0);
    mpz_set(A->y.x0, tempA.y.x0);
    A->inf = tempA.inf;
}

void efp_ecd(struct efp *R, struct efp P, mpz_t p){
    // 無限遠点処理
    if(P.inf == 1){
        return;
    }

    else if(mpz_sgn(P.y.x0) == 0){
        R->inf = 1;
        return;
    }
    else{
        
        struct fp three;mpz_init_set_str(three.x0, "3", 10);
        struct fp lambda, lambda_numerator, lambda_denominator;     // 分子lambda_numerator, 分母lambda_denominator
        mpz_inits(lambda.x0, lambda_numerator.x0, lambda_denominator.x0, NULL);
    
        fp_mul(&lambda_numerator, P.x, P.x, p);
        fp_mul(&lambda_numerator, lambda_numerator, three, p);
        fp_add(&lambda_denominator, P.y, P.y, p);
        fp_inv(&lambda_denominator, lambda_denominator, p);
        fp_mul(&lambda, lambda_numerator, lambda_denominator, p);
    
        struct fp temp_Rx, temp_Ry;
        mpz_inits(temp_Rx.x0, temp_Ry.x0, NULL);

        fp_mul(&temp_Rx, lambda, lambda, p);
        fp_sub(&temp_Rx, temp_Rx, P.x, p);
        fp_sub(&temp_Rx, temp_Rx, P.x, p);

        fp_sub(&temp_Ry, P.x, temp_Rx, p);
        fp_mul(&temp_Ry, temp_Ry, lambda, p);
        fp_sub(&temp_Ry, temp_Ry, P.y, p);
    
        mpz_set(R->x.x0, temp_Rx.x0);
        mpz_set(R->y.x0, temp_Ry.x0);
        R->inf = 0;
    }
}
 
void efp_eca(struct efp *R, struct efp P, struct efp Q, mpz_t p){
    if(P.inf == 1){
        mpz_set(R->x.x0, Q.x.x0);
        mpz_set(R->y.x0, Q.y.x0);
        R->inf = 0;
        return;
    }

    else if(Q.inf == 1){
        mpz_set(R->x.x0, P.x.x0);
        mpz_set(R->y.x0, P.y.x0);
        R->inf = 0;
        return;
    }

    else if(mpz_cmp(P.x.x0,Q.x.x0) == 0 && mpz_cmp(P.y.x0,Q.y.x0) == 0){
        efp_ecd(R, P, p);
        return;
    }

    else if(mpz_cmp(P.x.x0,Q.x.x0) == 0){
        R->inf = 1;
        return;
    }
    
    struct fp lambda, lambda_numerator, lambda_denominator;     // 分子 : lambda_numerator, 分母 : lambda_denominator
    mpz_inits(lambda.x0, lambda_numerator.x0, lambda_denominator.x0, NULL);
 
    fp_sub(&lambda_numerator, Q.y, P.y, p);
    fp_sub(&lambda_denominator, Q.x, P.x, p);
    fp_inv(&lambda_denominator, lambda_denominator, p);
    fp_mul(&lambda, lambda_numerator, lambda_denominator, p);
 
    struct fp temp_Rx, temp_Ry;
    mpz_inits(temp_Rx.x0, temp_Ry.x0, NULL);
    fp_mul(&temp_Rx, lambda, lambda, p);
    fp_sub(&temp_Rx, temp_Rx, P.x, p);
    fp_sub(&temp_Rx, temp_Rx, Q.x, p);
 
    fp_sub(&temp_Ry, P.x, temp_Rx, p);
    fp_mul(&temp_Ry, temp_Ry, lambda, p);
    fp_sub(&temp_Ry, temp_Ry, P.y, p);
 
    mpz_set(R->x.x0, temp_Rx.x0);
    mpz_set(R->y.x0, temp_Ry.x0);
    R->inf = 0;
    
}
 
void efp_scm(struct efp *R, struct efp P, mpz_t s, mpz_t p){
    struct efp tempR;mpz_inits(tempR.x.x0, tempR.y.x0, NULL);tempR.inf = 1;
    char *scalar_binary = mpz_get_str (NULL, 2, s);
    size_t len = strlen(scalar_binary);
 
    for(size_t i = 0; i < len; i++){
        efp_ecd(&tempR, tempR, p);
        if(scalar_binary[i] - '0' == 1){
            efp_eca(&tempR, tempR, P, p);
        }
    }

    mpz_set(R->x.x0, tempR.x.x0);
    mpz_set(R->y.x0, tempR.y.x0);
    R->inf = tempR.inf;
}

void gen_z(mpz_t z) {
    mpz_t tmp; mpz_init(tmp);

    // set z <- 2^77+2^11-2^9-2^6
    mpz_ui_pow_ui(tmp, 2, 77);
    mpz_add(z, z, tmp);
    mpz_ui_pow_ui(tmp, 2, 11);
    mpz_add(z, z, tmp);
    mpz_ui_pow_ui(tmp, 2, 9);
    mpz_sub(z, z, tmp);
    mpz_ui_pow_ui(tmp, 2, 6);
    mpz_sub(z, z, tmp);

    mpz_clear(tmp);
}

void print_params(mpz_t z, mpz_t prime, mpz_t r, mpz_t t) {
    gmp_printf("z = %Zd\n", z);
    gmp_printf("p = %Zd\n", prime);
    gmp_printf("r = %Zd\n", r);
    gmp_printf("t = %Zd\n", t);
}

void gen_params(mpz_t z, mpz_t prime, mpz_t r, mpz_t t) {
    mpz_t tmp0; mpz_init(tmp0);
    mpz_t tmp1; mpz_init(tmp1);
    mpz_t tmp2; mpz_init(tmp2);
    // int bit_num = 70;

    gen_z(z);

    // set prime and r
    mpz_sub_ui(tmp0, z, 1); // set tmp0 <- z-1
    mpz_pow_ui(tmp0, tmp0, 2); // set tmp0 <- tmp0^2 <- (z-1)^2
    mpz_pow_ui(tmp1, z, 4); // set tmp1 <- z^4
    mpz_pow_ui(tmp2, z, 2); // set tmp2 <- z^2
    mpz_sub(tmp1, tmp1, tmp2); // set tmp1 <- tmp1-tmp2 <- z^4-z^2
    mpz_add_ui(r, tmp1, 1); // set r <- tmp1+1 <- z^4-z^2+1
    mpz_mul(tmp0, tmp0, r); // set tmp0 <- tmp0*r <- (z-1)^2*(z^4-z^2+1)
    mpz_tdiv_q_ui(tmp0, tmp0, 3); // set tmp0 <- tmp0/3 <- (z-1)^2*(z^4-z^2+1)/3
    mpz_add(prime, tmp0, z); // set prime <- tmp0+z <- (z-1)^2*(z^4-z^2+1)/3+z
    // set t
    mpz_add_ui(t, z, 1); // set t <- z+1
    mpz_clears(tmp0, tmp1, tmp2, NULL);
}

void rank_number(mpz_t E, int n, mpz_t t, mpz_t p){
    mpz_t tn, tn_1, tn_2;
    mpz_inits(tn, tn_1, tn_2, NULL);

    // 初期値設定
    mpz_set(tn_2, t);   // t1 = t
    mpz_mul(tn_1, t, t);mpz_sub(tn_1, tn_1, p);mpz_sub(tn_1, tn_1, p);  // t2 = t1^2 -2p

    // 例外処理
    if(n == 1 || n == 2){
        mpz_pow_ui(E, p, n);
        mpz_add_ui(E, E, 1);
        if(n == 1){
            mpz_sub(E, E, tn_2);
        }
        else if(n == 2){
            mpz_sub(E, E, tn_1);
        }
        return;
    }

    mpz_t temp1, temp2;
    mpz_inits(temp1, temp2, NULL);

    // tn = tn-1 * t1 - p * tn-2
    for(int i = n - 2; i > 0; i--){
        mpz_mul(temp1, tn_1, t);
        mpz_mul(temp2, p, tn_2);
        mpz_sub(tn, temp1, temp2);
        mpz_set(tn_2, tn_1);
        mpz_set(tn_1, tn);
    }
    // gmp_printf("t : %Zd\n", tn);
    
    mpz_pow_ui(E, p, n);
    mpz_add_ui(E, E, 1);
    mpz_sub(E, E, tn);

    mpz_clears(tn, tn_1, tn_2, temp1, temp2, NULL);
}

void l_TP(struct fp12 *S, struct efp12 P, struct efp12 Q, struct efp12 T, mpz_t p){
    struct fp12 temp, temp1, temp2, temp3, temp4;
    fp12_init(&temp);fp12_init(&temp1);fp12_init(&temp2);fp12_init(&temp3);fp12_init(&temp4);

    fp12_sub(&temp1, Q.y, P.y, p);  //yQ-yP
    fp12_sub(&temp2, P.y, T.y, p);  //yP-yT
    fp12_sub(&temp3, P.x, T.x, p);  //xP-xT
    fp12_inv(&temp3, temp3, p); 
    fp12_sub(&temp4, Q.x, P.x, p);  //Qx-Px

    fp12_mul(&temp2, temp2, temp3, p);  //(yP-yT)/(xP-xT)
    fp12_mul(&temp2, temp2, temp4, p);  //
    fp12_sub(&temp, temp1, temp2, p);

    fp12_set(S, temp);

    fp12_clear(&temp);fp12_clear(&temp1);fp12_clear(&temp2);fp12_clear(&temp3);fp12_clear(&temp4);
}

void l_TQ_twist(struct fp12 *S, struct efp12 P, struct efp12 Q, struct efp12 T, mpz_t p){
    // twist set
    struct fp2 Qx_twist, Qy_twist;fp2_init(&Qx_twist);fp2_init(&Qy_twist);
    struct fp2 Tx_twist, Ty_twist;fp2_init(&Tx_twist);fp2_init(&Ty_twist);
    fp2_set(&Qx_twist, Q.x.x0.x1);fp2_set(&Qy_twist, Q.y.x1.x1);
    fp2_set(&Tx_twist, T.x.x0.x1);fp2_set(&Ty_twist, T.y.x1.x1);

    struct fp2 tmp;fp2_init(&tmp);  // tmp = (Yq - Yt) / (Xq - Xt)
    struct fp2 tmp1;fp2_init(&tmp1);
    fp2_sub(&tmp, Qy_twist, Ty_twist, p);     // tmp = Yq - Yt
    fp2_sub(&tmp1, Qx_twist, Tx_twist, p);    // tmp1 = Xq - Xt
    fp2_inv(&tmp1, tmp1, p);
    fp2_mul(&tmp, tmp, tmp1, p);

    struct fp2 S0, S1, S2, tmp_S;
    fp2_init(&S0);fp2_init(&S1);fp2_init(&S2);fp2_init(&tmp_S);

    fp2_set(&S0, P.y.x0.x0);

    fp2_mul(&S1, tmp, Qx_twist, p);
    fp2_sub(&S1, S1, Qy_twist, p);

    fp2_mul(&tmp_S, tmp, P.x.x0.x0, p);
    fp2_sub(&S2, S2, tmp_S, p);

    fp2_set(&S->x0.x0, S0);
    fp2_set(&S->x1.x0, S2);
    fp2_set(&S->x1.x1, S1);

    fp2_clear(&Qx_twist);fp2_clear(&Qy_twist);
    fp2_clear(&Tx_twist);fp2_clear(&Ty_twist);
    fp2_clear(&tmp);fp2_clear(&tmp1);
    fp2_clear(&S0);fp2_clear(&S1);fp2_clear(&S2);fp2_clear(&tmp_S);

}

void l_TT(struct fp12 *S, struct efp12 Q, struct efp12 T, mpz_t p){
    struct fp12 temp, temp1, temp2, temp3, temp4;
    fp12_init(&temp);fp12_init(&temp1);fp12_init(&temp2);fp12_init(&temp3);fp12_init(&temp4);

    fp12_sub(&temp1, Q.y, T.y, p);
    fp12_square(&temp2, T.x, p);
    fp12_add(&temp3, temp2, temp2, p);
    fp12_add(&temp2, temp3, temp2, p);//3Xt^2
    fp12_add(&temp3, T.y, T.y, p);//2yt
    fp12_inv(&temp3, temp3, p);
    fp12_sub(&temp4, Q.x, T.x, p);

    fp12_mul(&temp2, temp2, temp3, p);// 3Xt^2/2yt
    fp12_mul(&temp2, temp2, temp4, p);
    fp12_sub(&temp, temp1, temp2, p);
    
    fp12_set(S, temp);

    fp12_clear(&temp);fp12_clear(&temp1);fp12_clear(&temp2);fp12_clear(&temp3);fp12_clear(&temp4);
}

void l_TT_twist(struct fp12 *S, struct efp12 P, struct efp12 T, mpz_t p){
    // twist set
    struct fp2 Tx_twist, Ty_twist;fp2_init(&Tx_twist);fp2_init(&Ty_twist);
    fp2_set(&Tx_twist, T.x.x0.x1);fp2_set(&Ty_twist, T.y.x1.x1);

    // l_TT
    struct fp2 tmp;fp2_init(&tmp);  // tmp = 3Xt^2 / 2Yt
    struct fp2 tmp1;fp2_init(&tmp1);
    fp2_square(&tmp1, Tx_twist, p);     //tmp1 = Xt^2
    fp2_add(&tmp, tmp1, tmp1, p);
    fp2_add(&tmp, tmp, tmp1, p);
    fp2_add(&tmp1, Ty_twist, Ty_twist, p);   // tmp1 = 2Yt
    fp2_inv(&tmp1, tmp1, p);
    fp2_mul(&tmp, tmp, tmp1, p);

    struct fp2 S0, S1, S2, tmp_S;
    fp2_init(&S0);fp2_init(&S1);fp2_init(&S2);fp2_init(&tmp_S);

    fp2_set(&S0, P.y.x0.x0);

    fp2_mul(&S1, tmp, Tx_twist, p);
    fp2_sub(&S1, S1, Ty_twist, p);

    fp2_mul(&tmp_S, tmp, P.x.x0.x0, p);
    fp2_sub(&S2, S2, tmp_S, p);

    fp2_set(&S->x0.x0, S0);
    fp2_set(&S->x1.x0, S2);
    fp2_set(&S->x1.x1, S1);

    fp2_clear(&Tx_twist);fp2_clear(&Ty_twist);
    fp2_clear(&tmp);fp2_clear(&tmp1);
    fp2_clear(&S0);fp2_clear(&S1);fp2_clear(&S2);fp2_clear(&tmp_S);
}

void generate1(struct efp12 *P, struct fp b, mpz_t z, mpz_t r, mpz_t p, gmp_randstate_t state){
    struct efp tempP;fp_init(&tempP.x);fp_init(&tempP.y);tempP.inf = 0;

    efp_random(&tempP, b, p, state);
    // gmp_printf("tempP : (%Zd, %Zd)\n", tempP.x.x0, tempP.y.x0);
    struct fp c1, c2, c3;mpz_inits(c1.x0, c2.x0, c3.x0, NULL);
    fp_mul(&c1, tempP.x, tempP.x, p);
    fp_mul(&c1, c1, tempP.x, p);
    fp_add(&c1, c1, b, p);
    // gmp_printf("x : %Zd\n", c1.x0);

    fp_mul(&c2, tempP.y, tempP.y, p);
    // gmp_printf("y : %Zd\n", c2.x0);

    mpz_t temp_E;mpz_init(temp_E);
    mpz_sub(temp_E, p, z);mpz_cdiv_q (temp_E, temp_E, r);
    // gmp_printf("tempE : %Zd\n", temp_E);
    efp_scm(&tempP, tempP, temp_E, p);//
    mpz_set(P->x.x0.x0.x0.x0, tempP.x.x0);mpz_set(P->y.x0.x0.x0.x0, tempP.y.x0);

    mpz_clear(temp_E);
    fp_clear(&tempP.x);fp_clear(&tempP.y);
}

void generate2(struct efp12 *Q, struct fp b, mpz_t E, mpz_t r, mpz_t p, gmp_randstate_t state){
    struct efp12 P;
    fp12_init(&P.x);fp12_init(&P.y);P.inf = 0;
    efp12_random(&P, b, p, state);

    // B定義
    struct efp12 B;fp12_init(&B.x);fp12_init(&B.y);B.inf = 0;
    mpz_t temp_r;mpz_init(temp_r);
    mpz_t temp_E;mpz_init(temp_E);
    mpz_mul(temp_r, r, r);
    mpz_cdiv_q (temp_E, E, temp_r);
    efp12_scm(&B, P, temp_E, p);
    
    // R定義
    struct efp12 R;fp12_init(&R.x);fp12_init(&R.y);R.inf = 0;
    // fp12_pow(&R.x, B.x, p, p);
    // fp12_pow(&R.y, B.y, p, p);
    fp12_Frobenius(&R.x, B.x, p);
    fp12_Frobenius(&R.y, B.y, p);
    
    // S定義
    struct efp12 S;fp12_init(&S.x);fp12_init(&S.y);S.inf = 0;
    fp12_set(&S.x, B.x);
    fp12_set(&S.y, B.y);
    fp12_neg(&S.y, S.y, p);
    
    // Q計算
    struct efp12 tempQ;fp12_init(&tempQ.x);fp12_init(&tempQ.y);tempQ.inf = 0;
    efp12_eca(&tempQ, R, S, p);
    fp12_set(&Q->x, tempQ.x);fp12_set(&Q->y, tempQ.y);

    mpz_clears(temp_r, temp_E, NULL);
    fp12_clear(&P.x);fp12_clear(&P.y);
    fp12_clear(&B.x);fp12_clear(&B.y);
    fp12_clear(&R.x);fp12_clear(&R.y);
    fp12_clear(&S.x);fp12_clear(&S.y);
    fp12_clear(&tempQ.x);fp12_clear(&tempQ.y);

}

void efp12_eca_twist(struct efp12 *R, struct efp12 P, struct efp12 Q, mpz_t p){
    struct efp2 P_twist;efp2_init(&P_twist);
    fp2_set(&P_twist.x, P.x.x0.x1);fp2_set(&P_twist.y, P.y.x1.x1);

    struct efp2 Q_twist;efp2_init(&Q_twist);
    fp2_set(&Q_twist.x, Q.x.x0.x1);fp2_set(&Q_twist.y, Q.y.x1.x1);

    struct efp2 R_tmp;efp2_init(&R_tmp);

    efp2_eca(&R_tmp, P_twist, Q_twist, p);

    fp2_set(&R->x.x0.x1, R_tmp.x);fp2_set(&R->y.x1.x1, R_tmp.y);   

    efp2_clear(&P_twist);efp2_clear(&Q_twist);efp2_clear(&R_tmp); 
}

void efp12_ecd_twist(struct efp12 *R, struct efp12 P, mpz_t p){
    struct efp2 P_twist;efp2_init(&P_twist);
    fp2_set(&P_twist.x, P.x.x0.x1);fp2_set(&P_twist.y, P.y.x1.x1);

    struct efp2 R_tmp;efp2_init(&R_tmp);
    efp2_ecd(&R_tmp, P_twist, p);

    fp2_set(&R->x.x0.x1, R_tmp.x);fp2_set(&R->y.x1.x1, R_tmp.y);

    efp2_clear(&P_twist);efp2_clear(&R_tmp);
}

void optimal_ate_miller(struct fp12 *ans, mpz_t z, struct efp12 P, struct efp12 Q, mpz_t p){
	// 初期化
	struct fp12 tempf;fp12_init(&tempf);mpz_set_str(tempf.x0.x0.x0.x0, "1", 10);
	struct efp12 T;fp12_init(&T.x);fp12_init(&T.y);T.inf = 0;
	fp12_set(&T.x, Q.x);fp12_set(&T.y, Q.y);

    // if(mpz_sgn(z) == -1){
    //     mpz_neg(z, z);
    // }
    // gmp_printf("|z| : %Zd\n", z);
	char *z_binary = mpz_get_str(NULL, 2, z);
	size_t len = strlen(z_binary);
    // printf("len : %ld\n", len);
    // printf("z -> %s\n", z_binary);

    struct fp12 temp1;
    fp12_init(&temp1);
	for(size_t i = 0; i < (len-1); i ++){
        // printf("%d", z_binary[i+1] - '0');
		fp12_square(&tempf, tempf, p);
       	// l_TT(&temp1, P, T, p);
        l_TT_twist(&temp1, P, T, p);
        
		fp12_mul(&tempf, tempf, temp1, p);
		// efp12_ecd(&T, T, p);
		efp12_ecd_twist(&T, T, p);
		if(z_binary[i+1] - '0' == 1){
			// l_TP(&temp1, Q, P, T, p);
            l_TQ_twist(&temp1, P, Q, T, p);
            // gmp_printf("\nlTQ : [{(%Zd,%Zd),(%Zd,%Zd),(%Zd,%Zd)},{(%Zd,%Zd),(%Zd,%Zd),(%Zd,%Zd)}]\n",
            // temp1.x0.x0.x0.x0, temp1.x0.x0.x1.x0, temp1.x0.x1.x0.x0, temp1.x0.x1.x1.x0, temp1.x0.x2.x0.x0, temp1.x0.x2.x1.x0, temp1.x1.x0.x0.x0, temp1.x1.x0.x1.x0, temp1.x1.x1.x0.x0, temp1.x1.x1.x1.x0, temp1.x1.x2.x0.x0, temp1.x1.x2.x1.x0);
			fp12_mul(&tempf, tempf, temp1, p);
			// efp12_eca(&T, T, Q, p);
            efp12_eca_twist(&T, T, Q, p);
		}
	}
    // printf("\n");
	fp12_set(ans, tempf);

    fp12_clear(&tempf);fp12_clear(&temp1);fp12_clear(&T.x);fp12_clear(&T.y);
    free(z_binary); z_binary = NULL;
}

void final_exponentiation(struct fp12 *S, struct fp12 f, mpz_t r, mpz_t p){
	mpz_t temp;mpz_init(temp);
	// mpz_t temp1;mpz_init(temp1);
	struct fp12 ans;fp12_init(&ans);
	mpz_pow_ui(temp, p, 12);
	mpz_sub_ui(temp, temp, 1);
    // mpz_cdiv_r(temp1, temp, r);
	mpz_cdiv_q(temp, temp, r);
    // gmp_printf("p^12 - 1/r : %Zd\n", temp);
	fp12_pow(&ans, f, temp, p);

    fp12_set(S, ans);

    mpz_clear(temp);fp12_clear(&ans);
}

void easy_final_exponentiation(struct fp12 *S, struct fp12 f, mpz_t r, mpz_t p){
    //// easy part
    struct fp12 F, tmp_f, tmp_f1, tmp_f2, tmp_f3;
    fp12_init(&F);fp12_init(&tmp_f);fp12_init(&tmp_f1);fp12_init(&tmp_f2);fp12_init(&tmp_f3);
    mpz_t s;mpz_init_set_ui(s, 6);

    // tmp_f1 = (f^p)^6
    fp12_Frobenius(&tmp_f1, f, p);
    fp12_pow(&tmp_f1, tmp_f1, s, p);

    // tmp_f2 = 1/f
    fp12_inv(&tmp_f2, f, p);

    // F = f^(p^6 - 1)
    fp12_mul(&F, tmp_f1, tmp_f2, p);

    // tmp_f3 = F^(p^2)
    fp12_Frobenius(&tmp_f3, F, p);
    fp12_mul(&tmp_f3, tmp_f3, tmp_f3, p);

    // tmp_f = F^(p^2 + 1)
    fp12_mul(&tmp_f, tmp_f3, F, p);


    //// hard part
    // tmp_p = (p^4 - p^2 + 1) / r
    mpz_t tmp_p, tmp_p1, tmp_p2;mpz_inits(tmp_p, tmp_p1, tmp_p2, NULL);
    mpz_pow_ui(tmp_p1, p, 4);
    mpz_mul(tmp_p2, p, p);
    mpz_sub(tmp_p, tmp_p1, tmp_p2);
    mpz_add_ui(tmp_p, tmp_p, 1);
    mpz_cdiv_q(tmp_p, tmp_p, r);

    fp12_pow(&tmp_f, tmp_f, tmp_p, p);

    fp12_set(S, tmp_f);

    mpz_clears(s, tmp_p, tmp_p1, tmp_p2, NULL);
    fp12_clear(&F);fp12_clear(&tmp_f);fp12_clear(&tmp_f1);fp12_clear(&tmp_f2);fp12_clear(&tmp_f3);
}

void hard_final_exponentiation(struct fp12 *S, struct fp12 f, mpz_t r, mpz_t p, mpz_t z){
    //// easy part
    struct fp12 F, tmp_f, tmp_f1, tmp_f2, tmp_f3;
    fp12_init(&F);fp12_init(&tmp_f);fp12_init(&tmp_f1);fp12_init(&tmp_f2);fp12_init(&tmp_f3);
    mpz_t s;mpz_init_set_ui(s, 6);

    // tmp_f1 = (f^p)^6
    fp12_Frobenius(&tmp_f1, f, p);
    fp12_pow(&tmp_f1, tmp_f1, s, p);

    // tmp_f2 = 1/f
    fp12_inv(&tmp_f2, f, p);

    // F = f^(p^6 - 1)
    fp12_mul(&F, tmp_f1, tmp_f2, p);

    // tmp_f3 = F^(p^2)
    fp12_Frobenius(&tmp_f3, F, p);
    fp12_mul(&tmp_f3, tmp_f3, tmp_f3, p);

    // tmp_f = F^(p^2 + 1)
    fp12_mul(&tmp_f, tmp_f3, F, p);

    //// hard part
    // hard part = t3p^3 + t2p^2 + t1p + t0
    mpz_t t0, t1, t2, t3;mpz_inits(t0, t1, t2, t3, NULL);

    // t3 = m = (χ - 1)^2 / r
    mpz_sub_ui(t3, z, 1);
    mpz_mul(t3, t3, t3);
    mpz_cdiv_q(t3, t3, r);

    // t2 = χ *　t3
    mpz_mul(t2, z, t3);

    // t1 = χ * t2 - t3
    mpz_mul(t1, z, t2);
    mpz_sub(t1, t1, t3);

    // t0 = χ * t1 + 1
    mpz_mul(t0, z, t1);
    mpz_add_ui(t0, t0, 1);

    mpz_t tmp_p, tmp_p1, tmp_p2, tmp_p3;
    mpz_inits(tmp_p, tmp_p1, tmp_p2, tmp_p3, NULL);

    mpz_pow_ui(tmp_p1, p, 3);
    mpz_mul(tmp_p1, tmp_p1, t3);

    mpz_mul(tmp_p2, p, p);
    mpz_mul(tmp_p2, tmp_p2, t2);

    mpz_mul(tmp_p3, t1, p);

    mpz_add(tmp_p, tmp_p1, tmp_p2);
    mpz_add(tmp_p, tmp_p, tmp_p3);
    mpz_add(tmp_p, tmp_p, t0);

    fp12_pow(&tmp_f, tmp_f, tmp_p, p);

    fp12_set(S, tmp_f);

    mpz_clears(s, tmp_p, tmp_p1, tmp_p2, NULL);
    fp12_clear(&F);fp12_clear(&tmp_f);fp12_clear(&tmp_f1);fp12_clear(&tmp_f2);fp12_clear(&tmp_f3);
}

void bilinearity(struct efp12 P, struct efp12 Q, mpz_t a, mpz_t b, mpz_t z, mpz_t r, mpz_t p){
    struct efp12 tempP, tempQ;
    fp12_init(&tempP.x);fp12_init(&tempP.y);tempP.inf = 0; 
    fp12_init(&tempQ.x);fp12_init(&tempQ.y);tempQ.inf = 0;
    struct fp12 f1, f2, f3;
    fp12_init(&f1);fp12_init(&f2);fp12_init(&f3);

    // １回目
    printf("\nf1start\n");
    efp12_scm(&tempP, P, b, p);
    efp12_scm(&tempQ, Q, a, p);
    optimal_ate_miller(&f1, z, tempP, tempQ, p);
    final_exponentiation(&f1, f1, r, p);
    gmp_printf("pair1 : ");
    fp12_printf(f1);

    // ２回目
    printf("\nf2start\n");
    efp12_scm(&tempP, P, a, p);
    efp12_scm(&tempQ, Q, b, p);
    optimal_ate_miller(&f2, z, tempP, tempQ, p);
    final_exponentiation(&f2, f2, r, p);
    gmp_printf("pair2 : ");
    fp12_printf(f2);

    // ３回目
    printf("\nf3start\n");
    mpz_t s;mpz_init(s);mpz_mul(s, a, b);
    gmp_printf("ab : %Zd\n", s);
    optimal_ate_miller(&f3, z, P, Q, p);
    final_exponentiation(&f3, f3, r, p);
    fp12_pow(&f3, f3, s, p);
    gmp_printf("pair3 : ");
    fp12_printf(f3);

    if(fp12_cmp(f1, f2) == 0 && fp12_cmp(f2, f3) == 0){
        printf("\nSuccess!!\n");
        return;
    }
    else{
        printf("\nMiss...\n");
        return;
    }

    fp12_clear(&tempP.x);fp12_clear(&tempP.y);fp12_clear(&tempQ.x);fp12_clear(&tempQ.y);
    fp12_clear(&f1);fp12_clear(&f2);fp12_clear(&f3);
    mpz_clear(s);
}