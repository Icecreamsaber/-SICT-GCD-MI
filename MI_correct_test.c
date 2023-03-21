#include <time.h>
#include "Modular_Inversion.h"

int main(){
	gmp_randstate_t state;
	gmp_randinit_mt(state);
	gmp_randseed_ui(state, time(NULL));

	mpz_t p, pre_fctmi, pre_hdfctmi, pre_ectmi, result, result1;
	mpz_init_set_str(p, P_256, 10);
	mpz_init_set_str(pre_fctmi, P_256_PRE_FCIMI, 10);
	mpz_init_set_str(pre_hdfctmi, P_256_PRE_hdFCIMI, 10);
	mpz_init_set_str(pre_ectmi, P_256_PRE_ECIMI, 10);
	mpz_init(result);
	mpz_init(result1);
	
	mpz_t test_num; mpz_init(test_num);
	for(int i = 0; i < 10000; ++i){
		mpz_urandomm(test_num, state, p);

		mpz_invert(result, test_num, p);
		//gmp_printf("mpz_invert: %Zd\n", result);

		FLT(result1, test_num, p, P_256_LOOPS_FLT);
		//gmp_printf("FLT: %Zd\n", result1);
		if(mpz_cmp(result, result1) != 0){
			gmp_printf("FLT: %Zd!=%Zd\n", result1, result);
			return 0;
		}
		
		FCTMI(result1, test_num, p, pre_fctmi, P_256_LOOPS_FCTMI);
		//gmp_printf("FCTMI: %Zd\n", result1);
		if(mpz_cmp(result, result1) != 0){
			gmp_printf("FCTMI: %Zd!=%Zd\n", result1, result);
			return 0;
		}
		
		hdFCTMI(result1, test_num, p, pre_hdfctmi, P_256_LOOPS_hdFCTMI);
		//gmp_printf("hdFCTMI: %Zd\n", result1);
		if(mpz_cmp(result, result1) != 0){
			gmp_printf("hdFCTMI: %Zd!=%Zd\n", result1, result);
			return 0;
		}

		ECTMI2(result1, test_num, p, pre_ectmi, P_256_LOOPS_ECTMI);
		//gmp_printf("ECTMI2: %Zd\n", result1);
		if(mpz_cmp(result, result1) != 0){
			gmp_printf("ECTMI2: %Zd!=%Zd\n", result1, result);
			return 0;
		}

		ECTMI3(&result1, test_num, p, pre_ectmi, P_256_LOOPS_ECTMI);
		//gmp_printf("ECTMI3: %Zd\n", result1);
		if(mpz_cmp(result, result1) != 0){
			gmp_printf("ECTMI3: %Zd!=%Zd\n", result1, result);
			return 0;
		}
	}
	printf("PASS!\n");	
	return 0;
}
