#include <time.h>
#include <sys/time.h>
#include "Modular_Inversion.h"

int main(){
	gmp_randstate_t state;
	gmp_randinit_mt(state);
	gmp_randseed_ui(state, time(NULL));

	//struct timeval start, end;
	//long sec, usec;
	uint64_t CC1, CC2;
	
	mpz_t p, pre_fctmi, pre_hdfctmi, pre_ectmi, result;
	mpz_init_set_str(p, P_384, 10);
	mpz_init_set_str(pre_fctmi, P_384_PRE_FCIMI, 10);
	mpz_init_set_str(pre_hdfctmi, P_384_PRE_hdFCIMI, 10);
	mpz_init_set_str(pre_ectmi, P_384_PRE_ECIMI, 10);
	mpz_init(result);

	int N = 100000;
	mpz_t test_num[N];
	for(int i = 0; i < N; ++i){
		mpz_init(test_num[i]);
		mpz_urandomm(test_num[i], state, p);
	}

	//gettimeofday(&start, NULL);
	CC1 = GetCC();
	for(int i = 0; i < N; ++i){
		FLT(result, test_num[i], p, P_384_LOOPS_FLT);
	}
	//gettimeofday(&end, NULL);
    	CC2 = GetCC();
    	printf("FLT: %ld\n", CC2-CC1);
    	//sec = end.tv_sec-start.tv_sec;
    	//usec = end.tv_usec-start.tv_usec;
	//printf("FLT: %ld\n", (sec*1000000+usec));
	
	//gettimeofday(&start, NULL);
	CC1 = GetCC();
	for(int i = 0; i < N; ++i){
		FCTMI(result, test_num[i], p, pre_fctmi, P_384_LOOPS_FCTMI);
	}
	//gettimeofday(&end, NULL);
    	//sec = end.tv_sec-start.tv_sec;
    	//usec = end.tv_usec-start.tv_usec;
	//printf("FCTMI: %ld\n", (sec*1000000+usec));
	CC2 = GetCC();
    	printf("FCTMI: %ld\n", CC2-CC1);
    	
    	//gettimeofday(&start, NULL);
	CC1 = GetCC();
	for(int i = 0; i < N; ++i){
		hdFCTMI(result, test_num[i], p, pre_hdfctmi, P_384_LOOPS_hdFCTMI);
	}
	//gettimeofday(&end, NULL);
    	//sec = end.tv_sec-start.tv_sec;
    	//usec = end.tv_usec-start.tv_usec;
	//printf("FCTMI: %ld\n", (sec*1000000+usec));
	CC2 = GetCC();
    	printf("hdFCTMI: %ld\n", CC2-CC1);

	//gettimeofday(&start, NULL);
	CC1 = GetCC();
	for(int i = 0; i < N; ++i){
		ECTMI2(result, test_num[i], p, pre_ectmi, P_384_LOOPS_ECTMI);
	}
	//gettimeofday(&end, NULL);
    	//sec = end.tv_sec-start.tv_sec;
    	//usec = end.tv_usec-start.tv_usec;
	//printf("ECTMI_IF: %ld\n", (sec*1000000+usec));
	CC2 = GetCC();
    	printf("ECTMI2: %ld\n", CC2-CC1);
	
	//gettimeofday(&start, NULL);
	CC1 = GetCC();
	for(int i = 0; i < N; ++i){
		ECTMI3(&result, test_num[i], p, pre_ectmi, P_384_LOOPS_ECTMI);
	}
	//gettimeofday(&end, NULL);
    	//sec = end.tv_sec-start.tv_sec;
    	//usec = end.tv_usec-start.tv_usec;
	//printf("ECTMI: %ld\n", (sec*1000000+usec));
	CC2 = GetCC();
    	printf("ECTMI3: %ld\n", CC2-CC1);

	return 0;
}
