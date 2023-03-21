#include <time.h>
#include "GCD.h"

int main(){
	gmp_randstate_t state;
	gmp_randinit_mt(state);
	gmp_randseed_ui(state, time(NULL));
	
	mpz_t x1, y1, result, result1;
	
	mpz_init(x1);
	mpz_init(y1);
	mpz_init(result);
	mpz_init(result1);
	
	for(int i = 0; i < 10000; ++i){
		mpz_urandomb(x1, state, 256);
		while(mpz_tstbit(x1, 0) == 0){
			mpz_urandomb(x1, state, 256);
		}
		mpz_urandomm(y1, state, x1);
		
		FCTMI(result, y1, x1, P_256_LOOPS_FCTMI);
		//gmp_printf("FCTMI: %Zd\n", result);
		
		hdFCTMI(result1, y1, x1, P_256_LOOPS_hdFCTMI);
		//gmp_printf("hdFCTMI: %Zd\n", result);
		
		if(mpz_cmp(result, result1) != 0){
			gmp_printf("hdFCTMI: x1=%Zd, x2=%Zd, %Zd\n", x1, y1, result1);
		}
		
		ECTMI(result1, y1, x1, P_256_LOOPS_ECTMI);
		//gmp_printf("ECTMI: %Zd\n", result);
		if(mpz_cmp(result, result1) != 0){
			gmp_printf("ECTMI: x1=%Zd, x2=%Zd, %Zd\n", x1, y1, result1);
		}
		
		ECTMI2(result1, y1, x1, P_256_LOOPS_ECTMI);
		//gmp_printf("ECTsMI: %Zd\n", result);
		if(mpz_cmp(result, result1) != 0){
			gmp_printf("ECTMI2: x1=%Zd, x2=%Zd, %Zd\n", x1, y1, result1);
		}
		
		ECTMI3(&result1, y1, x1, P_256_LOOPS_ECTMI);
		//gmp_printf("ECTMI: %Zd\n", result);
		if(mpz_cmp(result, result1) != 0){
			gmp_printf("ECTMI3: x1=%Zd, x2=%Zd, %Zd\n", x1, y1, result1);
		}
	}
	printf("PASS!\n");
	return 0;
}
