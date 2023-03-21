#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <gmp.h>

/*FLT: The number of loops defined by bitlen(p) with p224, p256, p384.*/
#define P_224_LOOPS_FLT 224 
#define P_256_LOOPS_FLT 256
#define P_384_LOOPS_FLT 384

/*FCTMI: The number of loops defined by (49d+57)//17 when d = 224, 256, or 384.*/
#define P_224_LOOPS_FCTMI 634//649 
#define P_256_LOOPS_FCTMI 724//741
#define P_384_LOOPS_FCTMI 1086//1110

/*hdFCTMI: The number of loops when d = 224, 256, or 384.*/
#define P_224_LOOPS_hdFCTMI 517 
#define P_256_LOOPS_hdFCTMI 590
#define P_384_LOOPS_hdFCTMI 885

/*Precomputed numbers of FCTMI.*/
#define P_224_PRE_FCIMI "2008672555230201739572940614232879648470511825043254648242560"//"26644009792206268056205042062915547495932996820470046226332801564673"
#define P_256_PRE_FCIMI "115789218083876703796405708104922954266844366797366597945313458537426472636415"//"90462569673685862684289949614566734656231367827980882643939160353509669339137"
#define P_384_PRE_FCIMI "39401102631331629224042747884736582634501236456156400192182914717222615670170491983949444240325277863502310224573143"//"39124746653188674882124692525769249793887199085144063008056343542510001139863950072793525947872276405658285025530879"

/*Precomputed numbers of hdFCTMI.*/
#define P_224_PRE_hdFCIMI "3369993333197670545056871185872896170421996906849403135584645414912"
#define P_256_PRE_hdFCIMI "49471717823667788566512707997982133781718357895554504641135267859606011904"
#define P_384_PRE_hdFCIMI "36092402410143488406830006651847393538951028178043571563308995486226118841864038653366404110975726264988796028928"

/*ECTMI: The number of loops defined by 2*bitlen(p) with p224, p256, p384.*/
#define P_224_LOOPS_ECTMI 448 
#define P_256_LOOPS_ECTMI 512
#define P_384_LOOPS_ECTMI 768

/*Precomputed numbers of ECTMI.*/
#define P_224_PRE_ECIMI "18831305204698540654176465351577770740545525012632715657217"
#define P_256_PRE_ECIMI "161759680002903838760694582340907787320905632487091026015073739997175"
#define P_384_PRE_ECIMI "39401997719623594004609818748773390158923402225834345130666249913583939917682844212927918355281735000658345862363656"

/*Definition of the field.*/
#define P_224 "26959946667150639794667015087019630673557916260026308143510066298881"
#define P_256 "115792089210356248762697446949407573530086143415290314195533631308867097853951"
#define P_384 "39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319"

/*FLT: q = a^-1 mod p.*/
void FLT(mpz_t q, mpz_t a, mpz_t p, int LOOPS){
	char k[LOOPS];
	mpz_t p_minus_two;
	mpz_init(p_minus_two);
	mpz_sub_ui(p_minus_two, p, 2);
	mpz_get_str(k, 2, p_minus_two);

	mpz_set(q, a);
	
	for(int i = 1; i < LOOPS; ++i){
		mpz_mul(q, q, q);
		mpz_mod(q, q, p);
		if(k[i] == '1'){
			mpz_mul(q, q, a);
			mpz_mod(q, q, p);
		}
	}
	return;
}

void FCTMI(mpz_t q, mpz_t a, mpz_t p, mpz_t pre_comp, int LOOPS){
	mpz_t u, v;
	mpz_init_set(u, a); mpz_init_set(v, p);

	mpz_t r; 
	mpz_set_ui(q, 0); mpz_init_set_ui(r, 1);

	int delta = 1; 
	int s, z; 
	int temp1; mpz_t temp2, temp3; mpz_init(temp2); mpz_init(temp3);

	for (int i = 0; i < LOOPS; ++i){
		z = mpz_tstbit(u, 0); //z = LSB(u)
		s = (delta >> 31) + 1; //s = signbit(-delta)

		/*delta = 1 + (1 - 2sz) delta*/
		s = s & z; // s = sz
		temp1 = (s << 1); // temp1 = 2sz
		temp1 = 1 - temp1; // temp1 = 1-2sz
		delta *= temp1; // delta = (1-2sz) delta
		delta += 1; // delta = delta + 1

		/*part 1 u = (u + (1-2sz)z v)/2*/
		z *= temp1; // z = (1-2sz)z
		mpz_set(temp2, v); temp2->_mp_size *= z;
		//mpz_mul_si(temp2, v, z); // temp2 = (1-2sz)zv
		
		/*v = v xor sz(v xor u)*/
		mpz_xor(temp3, v, u); // temp3 = v xor u
		temp3->_mp_size *= s;
		//mpz_mul_ui(temp3, temp3, s); // temp3 = sz(v xor u)
		mpz_xor(v, v, temp3); // v = v xor sz(v xor u)

		/*part 2 u = (u + (1-2sz)z v)/2*/
		mpz_add(u, u, temp2); // u = u + (1-2sz)zv
		mpz_tdiv_q_2exp(u, u, 1); // u = (u + (1-2sz)zv)/2

		/*part 1 r = (1-2sz) z q + r*/
		mpz_set(temp2, q); temp2->_mp_size *= z;
		//mpz_mul_si(temp2, q, z); // temp2 = (1-2sz)z q

		/*q = 2(q xor sz(q xor r)) q = 2(sz r+(sz XOR 1)q)*/
		mpz_xor(temp3, q, r); // temp3 = q xor r
		temp3->_mp_size *= s;
		//mpz_mul_ui(temp3, temp3, s); // temp3 = sz(q xor r)
		mpz_xor(q, q, temp3); // q = q xor sz(q xor r)
		mpz_mul_2exp(q, q, 1); // q = 2(q xor sz(q xor r)) 
		/*
		temp1 = s ^ 1; // temp1 = sz XOR 1
		mpz_mul_si(q, q, temp1); // q = q (sz XOR 1)
		mpz_mul_si(temp3, r, s); // temp3 = r sz
		mpz_add(q, q, temp3); // q = q (sz XOR 1) + r sz
		mpz_mul_2exp(q, q, 1); // q = 2(q (sz XOR 1) + r sz)
		*/

		/*part 2 r = (1-2sz) z q + r*/
		mpz_add(r, r, temp2);		
	}
	q->_mp_size *= mpz_sgn(v);
	//mpz_mul_si(q, q, mpz_sgn(v)); 
	mpz_mul(q, q, pre_comp); mpz_mod(q, q, p); 
	mpz_clear(u); mpz_clear(v); mpz_clear(r); mpz_clear(temp2); mpz_clear(temp3);
	return;
}

void hdFCTMI(mpz_t q, mpz_t a, mpz_t p, mpz_t pre_comp, int LOOPS){
	mpz_t u, v;
	mpz_init_set(u, a); mpz_init_set(v, p);

	mpz_t r; 
	mpz_set_ui(q, 0); mpz_init_set_ui(r, 1);

	int delta = 1; 
	int s, z; 
	int temp1; mpz_t temp2, temp3; mpz_init(temp2); mpz_init(temp3);

	for (int i = 0; i < LOOPS; ++i){
		z = mpz_tstbit(u, 0); //z = LSB(u)
		s = (delta >> 31) + 1; //s = signbit(-delta)

		/*delta = 2 + (1 - 2sz) delta*/
		s = s & z; // s = sz
		temp1 = (s << 1); // temp1 = 2sz
		temp1 = 1 - temp1; // temp1 = 1-2sz
		delta *= temp1; // delta = (1-2sz) delta
		delta += 2; // delta = delta + 2

		/*part 1 u = (u + (1-2sz)z v)/2*/
		z *= temp1; // z = (1-2sz)z
		mpz_set(temp2, v); temp2->_mp_size *= z;
		//mpz_mul_si(temp2, v, z); // temp2 = (1-2sz)zv
		
		/*v = v xor sz(v xor u)*/
		mpz_xor(temp3, v, u); // temp3 = v xor u
		temp3->_mp_size *= s;
		//mpz_mul_ui(temp3, temp3, s); // temp3 = sz(v xor u)
		mpz_xor(v, v, temp3); // v = v xor sz(v xor u)

		/*part 2 u = (u + (1-2sz)z v)/2*/
		mpz_add(u, u, temp2); // u = u + (1-2sz)zv
		mpz_tdiv_q_2exp(u, u, 1); // u = (u + (1-2sz)zv)/2

		/*part 1 r = (1-2sz) z q + r*/
		mpz_set(temp2, q); temp2->_mp_size *= z;
		//mpz_mul_si(temp2, q, z); // temp2 = (1-2sz)z q

		/*q = 2(q xor sz(q xor r)) q = 2(sz r+(sz XOR 1)q)*/
		mpz_xor(temp3, q, r); // temp3 = q xor r
		temp3->_mp_size *= s;
		//mpz_mul_ui(temp3, temp3, s); // temp3 = sz(q xor r)
		mpz_xor(q, q, temp3); // q = q xor sz(q xor r)
		mpz_mul_2exp(q, q, 1); // q = 2(q xor sz(q xor r)) 
		/*
		temp1 = s ^ 1; // temp1 = sz XOR 1
		mpz_mul_si(q, q, temp1); // q = q (sz XOR 1)
		mpz_mul_si(temp3, r, s); // temp3 = r sz
		mpz_add(q, q, temp3); // q = q (sz XOR 1) + r sz
		mpz_mul_2exp(q, q, 1); // q = 2(q (sz XOR 1) + r sz)
		*/

		/*part 2 r = (1-2sz) z q + r*/
		mpz_add(r, r, temp2);		
	}
	q->_mp_size*=mpz_sgn(v);
	//mpz_mul_si(q, q, mpz_sgn(v)); 
	mpz_mul(q, q, pre_comp); mpz_mod(q, q, p); 
	mpz_clear(u); mpz_clear(v); mpz_clear(r); mpz_clear(temp2); mpz_clear(temp3);
	return;
}

/*Our ECTMI*/
/*
void ECTMI_IF(mpz_t q, mpz_t a, mpz_t p, mpz_t pre_comp, int LOOPS){
	mpz_t u, v;
	mpz_init_set(u, a); mpz_init_set(v, p);

	mpz_t r; 
	mpz_set_str(q, "0", 10); mpz_init_set_str(r, "1", 10);
	
	int u_lsb, v_lsb;
	mpz_t temp1, temp2; mpz_init(temp1); mpz_init(temp2);
	
	for (int i = 0; i < LOOPS; ++i){
		u_lsb = mpz_tstbit(u, 0);
		v_lsb = mpz_tstbit(v, 0);
		mpz_sub(temp1, v, u); //sub temp1 = v - u 
		if(u_lsb == 1 && v_lsb == 1){
			mpz_tdiv_q_2exp(temp1, temp1, 1); //shift
			mpz_sub(temp2, q, r); //sub
			mpz_mul_2exp(v, r, 1); //shift
			if(mpz_cmp(u, temp1) > 0){
				mpz_set(q, v); mpz_set(r, temp2); mpz_set(temp2, u); mpz_set(v, temp2); mpz_set(temp2, temp1); mpz_set(u, temp2);// 6assign
			}
			else{
				mpz_set(q, temp2); mpz_set(r, v); mpz_set(temp2, u); mpz_set(u, temp1); mpz_set(v, u); mpz_set(u, temp2);// 6assign
			}
		}
		else if(u_lsb == 1 && v_lsb == 0){
			mpz_tdiv_q_2exp(u, v, 1); //shift
			mpz_sub(temp2, q, r); //sub
			mpz_mul_2exp(temp2, temp2, 1); //shift
			if(mpz_cmp(u, temp1) > 0){
				mpz_set(v, u); mpz_set(u, temp1); mpz_set(temp1, q); mpz_set(q, temp2); mpz_set(r, q); mpz_set(q, temp1); // 6assign
			}
			else{
				mpz_set(r, q); mpz_set(q, temp2); mpz_set(temp2, u); mpz_set(u, temp1); mpz_set(v, u); mpz_set(u, temp2);// 6assign
			}
		}
		else{
			mpz_tdiv_q_2exp(v, u, 1); //shift
			mpz_sub(temp2, q, r); //sub
			mpz_mul_2exp(temp2, temp2, 1); //shift
			if(mpz_cmp(v, temp1) > 0){
				mpz_set(q, r); mpz_set(r, temp2); mpz_set(temp2, v); mpz_set(v, temp1); mpz_set(u, v); mpz_set(v, temp2);// 6assign
			}
			else{
				mpz_set(u, v); mpz_set(v, temp1); mpz_set(temp1, r); mpz_set(r, temp2); mpz_set(q, r); mpz_set(r, temp1); // 6assign
			}
		}
	}
	mpz_mul(q, q, pre_comp); mpz_mod(q, q, p); 
	mpz_clear(u); mpz_clear(v); mpz_clear(r); mpz_clear(temp1); mpz_clear(temp2);
	return;
}
*/

void ECTMI2(mpz_t q, mpz_t a, mpz_t p, mpz_t pre_comp, int LOOPS){
	mpz_t u, v, r;
	mpz_init_set(u, a); mpz_init_set(v, p); mpz_set_ui(q, 0); mpz_init_set_ui(r, 1);
	
	mpz_t temp1; mpz_init(temp1);
	mpz_t temp2; mpz_init(temp2);
	mpz_t temp3; mpz_init(temp3);
	bool u_lsb, v_lsb, t2;
	int t1;
	
	mpz_t *sort1[2]; sort1[0] = &temp1; sort1[1] = &temp2;
	mpz_t *sort2[2]; sort2[0] = &v; sort2[1] = &temp3;
	for (int i = 0; i < LOOPS; ++i){
		u_lsb = mpz_tstbit(u, 0);
		v_lsb = mpz_tstbit(v, 0);
		
		mpz_set(temp2, v); temp2->_mp_size *= u_lsb;
		//mpz_mul_ui(temp2, v, u_lsb);
		t1 = 2-(u_lsb<<1)-v_lsb;
		mpz_set(temp1, u); temp1->_mp_size *= t1;
		//mpz_mul_si(temp1, u, t1);
		mpz_add(temp2, temp2, temp1);// add
		mpz_tdiv_q_2exp(temp2, temp2, 1);// shift
		
		mpz_set(temp1, q); temp1->_mp_size *= u_lsb;
		//mpz_mul_ui(temp1, q, u_lsb);
		mpz_set(temp3, r); temp3->_mp_size *= t1;
		//mpz_mul_si(temp3, r, t1);
		mpz_add(temp3, temp1, temp3);//add
		
		t2 = v_lsb ^ u_lsb;
		v->_mp_size *= t2;
		//mpz_mul_ui(v, v, t2);
		t1 = ((v_lsb & u_lsb)<<1)-1;
		u->_mp_size *= t1;
		//mpz_mul_si(u, u, t1);
		mpz_add(temp1, v, u);//add
		
		q->_mp_size *= t2;
		//mpz_mul_ui(q, q, t2);
		r->_mp_size *= t1;
		//mpz_mul_si(r, r, t1);
		mpz_add(v, q, r);//add
		mpz_mul_2exp(v, v, 1);//shift

		u_lsb = mpz_cmp(temp2, temp1) + 1;
		t2 = !u_lsb;
		mpz_set(q, *sort2[u_lsb]);
		mpz_set(r, *sort2[t2]);
		mpz_set(v, *sort1[u_lsb]);
		mpz_set(u, *sort1[t2]);
		/*
		mpz_xor(u, temp1, temp2); //xor
		mpz_mul_ui(u, u, u_lsb); 
		mpz_xor(v, u, temp1); //xor
		mpz_xor(u, u, temp2); //xor
		*/
		/*
		mpz_mul_ui(v, temp2, b1);
		mpz_mul_ui(u, temp1, !b1);
		mpz_xor(v, v, u); //xor
		mpz_mul_ui(u, temp2, !b1);
		mpz_mul_ui(temp1, temp1, b1);
		mpz_xor(u, u, temp1); //xor
		*/
	}
	mpz_mul(q, q, pre_comp); mpz_mod(q, q, p); 
	mpz_clear(temp1); mpz_clear(temp2); mpz_clear(temp3); 
	mpz_clear(u); mpz_clear(v); mpz_clear(r);
	return;
}

void ECTMI3(mpz_t *q, mpz_t a, mpz_t p, mpz_t pre_comp, int LOOPS){
	mpz_t u, v, r, temp1, temp2;
	mpz_init_set(u, a); mpz_init_set(v, p); mpz_set_ui(*q, 0); mpz_init_set_ui(r, 1);
	mpz_init(temp1); mpz_init(temp2);

	bool u_lsb, v_lsb, bf, bs;
	mpz_t temp3; mpz_init(temp3);
	
	mpz_t *A[3][2]; mpz_t *B[3][2];

	A[0][0] = A[1][1] = A[2][1] = &temp1;
	A[1][0] = &v;
	A[2][0] = A[0][1] = &u;
	
	B[0][0] = B[1][1] = B[2][1] = &temp2;
	B[1][0] = q;
	B[2][0] = B[0][1] = &r;
	
	uint8_t C[2][2];
	C[0][0] = C[1][1] = 0;
	C[0][1] = 2; C[1][0] = 1;

	for (int i = 0; i < LOOPS; ++i){
		mpz_sub(temp1, v, u); //sub temp1 = v - u
		mpz_sub(temp2, *q, r); //sub temp2 = q - r

		u_lsb = mpz_tstbit(u, 0);
		v_lsb = mpz_tstbit(v, 0);
		
		C[0][0] = C[u_lsb][v_lsb]; 

		mpz_tdiv_q_2exp(*A[C[0][0]][0], *A[C[0][0]][0], 1); //shift

		mpz_mul_2exp(*B[C[0][0]][1], *B[C[0][0]][1], 1); //shift

		bs = mpz_cmp(*A[C[0][0]][0], *A[C[0][0]][1]) + 1; //32bit_add
		bf = !bs;

		mpz_set(temp3, *A[C[0][0]][bf]); // new v
		mpz_set(u, *A[C[0][0]][bs]); 
		mpz_set(v, temp3);

		mpz_set(temp3, *B[C[0][0]][bf]); //new q
		mpz_set(r, *B[C[0][0]][bs]);
		mpz_set(*q, temp3);
	}
	mpz_mul(*q, *q, pre_comp); mpz_mod(*q, *q, p); 
	mpz_clear(temp1); mpz_clear(temp2); mpz_clear(temp3); 
	mpz_clear(u); mpz_clear(v); mpz_clear(r);
	return;
}

static __inline__ uint64_t GetCC(){
	unsigned int a, d;
	asm volatile("mfence;rdtsc" : "=a" (a), "=d" (d)::"memory");
	return ((uint64_t)a) | (((uint64_t)d)<<32);
}
