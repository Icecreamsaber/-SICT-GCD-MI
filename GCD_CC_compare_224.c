#include <time.h>
#include <sys/time.h>
#include "GCD.h"

int main(){
	gmp_randstate_t state;
	gmp_randinit_mt(state);
	gmp_randseed_ui(state, time(NULL));

	//struct timeval start, end;
	//long sec, usec;
	uint64_t CC1, CC2;
	
	mpz_t p, result;
	mpz_init_set_str(p, P_224, 10); mpz_sub_ui(p, p, 1);
	mpz_init(result);

	int N = 100000;
	mpz_t test_num[N];
	for(int i = 0; i < N; ++i){
		mpz_init(test_num[i]);
		mpz_urandomm(test_num[i], state, p);
		while(mpz_tstbit(test_num[i], 0) == 0){
			mpz_urandomm(test_num[i], state, p);
		}
	}
	
	//gettimeofday(&start, NULL);
	CC1 = GetCC();
	for(int i = 0; i < N; ++i){
		FCTMI(result, p, test_num[i], P_224_LOOPS_FCTMI);
	}
	//gettimeofday(&end, NULL);
    	CC2 = GetCC();
    	printf("FGCD: %ld\n", CC2-CC1);
    	//sec = end.tv_sec-start.tv_sec;
    	//usec = end.tv_usec-start.tv_usec;
	//printf("FGCD: %ld\n", (sec*1000000+usec));
	
	//gettimeofday(&start, NULL);
	CC1 = GetCC();
	for(int i = 0; i < N; ++i){
		hdFCTMI(result, p, test_num[i], P_224_LOOPS_hdFCTMI);
	}
	CC2 = GetCC();
    	printf("hdFGCD: %ld\n", CC2-CC1);
	//gettimeofday(&end, NULL);
    	//sec = end.tv_sec-start.tv_sec;
    	//usec = end.tv_usec-start.tv_usec;
	//printf("hdFGCD: %ld\n", (sec*1000000+usec));

	//gettimeofday(&start, NULL);
	CC1 = GetCC();
	for(int i = 0; i < N; ++i){
		ECTMI(result, test_num[i], p, P_224_LOOPS_ECTMI);
	}
	CC2 = GetCC();
    	printf("EGCD: %ld\n", CC2-CC1);
	//gettimeofday(&end, NULL);
    	//sec = end.tv_sec-start.tv_sec;
    	//usec = end.tv_usec-start.tv_usec;
	//printf("EGCD: %ld\n", (sec*1000000+usec));
	
	//gettimeofday(&start, NULL);
	CC1 = GetCC();
	for(int i = 0; i < N; ++i){
		ECTMI2(result, test_num[i], p, P_224_LOOPS_ECTMI);
	}
	CC2 = GetCC();
    	printf("EGCD2: %ld\n", CC2-CC1);
	//gettimeofday(&end, NULL);
    	//sec = end.tv_sec-start.tv_sec;
    	//usec = end.tv_usec-start.tv_usec;
	//printf("EGCD2: %ld\n", (sec*1000000+usec));
	
	//gettimeofday(&start, NULL);
	CC1 = GetCC();
	for(int i = 0; i < N; ++i){
		ECTMI3(&result, test_num[i], p, P_224_LOOPS_ECTMI);
	}
	CC2 = GetCC();
    	printf("EGCD3: %ld\n", CC2-CC1);
	//gettimeofday(&end, NULL);
    	//sec = end.tv_sec-start.tv_sec;
    	//usec = end.tv_usec-start.tv_usec;
	//printf("EGCD3: %ld\n", (sec*1000000+usec));

	return 0;
}
