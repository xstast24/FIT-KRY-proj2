/**
 * Kry projekt 2 - Implementace RSA
 * Filip Stastny (xstast24)
 */
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <time.h>
#include <gmp.h>

// TODO ... -Wall -Wextra -pedantic TODO pridat?, odebrat spousteni z make
using namespace std;

class Vars {
   public:
    int bits; // requested size of public modulus in bits
    mpz_t p; // p prime number
    mpz_t q; // q prime number
    mpz_t n; // n = p * q
    gmp_randstate_t rngState;
};

void generateKey(Vars vars);
void getRandomOddNumber(mpz_t* primeNumber, unsigned long bits, Vars vars);
bool isOddNumberPrimeNumber(mpz_t *number);
void exitWithMsg(string message);

void exitWithMsg(string message) {
    cerr << message;
    exit(1);
}

void generateKey(Vars vars) {
    mpz_init(vars.p);
    mpz_init(vars.q);
    mpz_init(vars.n);

    // init random number generator
    srand(time(NULL));
    unsigned long seed = rand();
    gmp_randinit_mt(vars.rngState);
    gmp_randseed_ui(vars.rngState, seed);

    // generate p & q
    if (vars.bits % 2 == 0) {
        // even length of key
        unsigned long bitsCount = vars.bits/2;
        
        getRandomOddNumber(&vars.p, bitsCount, vars);
        getRandomOddNumber(&vars.q, bitsCount, vars);
        gmp_printf ("P: %Zd\n", vars.p);
        gmp_printf ("Q: %Zd\n", vars.q);
        cout << endl << isOddNumberPrimeNumber(&vars.p) << endl;
        cout << endl << isOddNumberPrimeNumber(&vars.q) << endl;
    } else {
        // odd length of key
    }

}

void getRandomPrimeNumber(mpz_t *primeNumber, unsigned long bitsLength, Vars vars) {
    
}

/* Generate random odd number.*/
void getRandomOddNumber(mpz_t* oddNumber, unsigned long bitsLength, Vars vars) {
    mpz_urandomb(*oddNumber, vars.rngState, bitsLength);
    // set last bit to 1 = odd number
    mpz_t lastBit;
    mpz_init(lastBit);
    mpz_set_ui(lastBit, 1);
    mpz_ior(*oddNumber, *oddNumber, lastBit);
}

/*Check if given ODD number is a prime number.*/
bool isOddNumberPrimeNumber(mpz_t *number) {
    // check modulo for incrementing number until increment gets bigger than half of the number
    // after that, it is sure thaat number is prime number
    mpz_t i;
    mpz_init(i);
    mpz_set_ui(i, 3);
    mpz_t halfNumber;
    mpz_init(halfNumber);
    mpz_cdiv_q_ui(halfNumber, *number, 2);
    // add 1 to be sure it works if 1 gets rounded when int division odd number
    mpz_add_ui(halfNumber, halfNumber, 1);

    // cycle until the divider gets bigger than the number/2
    // increment by 2 > testing only odd dividers (number is odd originally, so not dividable by 2)
    mpz_t modulo;
    mpz_init(modulo);
    for(i; mpz_cmp(halfNumber, i) > 0; mpz_add_ui(i, i, 2)){
        mpz_cdiv_r(modulo, *number, i);
        if(mpz_cmp_ui(modulo, 0) == 0) {
            return false;
        } 
    }
    return true;
}

int main(int argc, char *argv[])
{
    Vars vars; //object holding commonly used variables

    if (argc < 3) {
        exitWithMsg("Invalid arguments count.");
    } else {
        // generating keys
        if (argc == 3 && argv[1] == string("-g")) {
            vars.bits = stoi(argv[2]);
            generateKey(vars);
        } else {
            exitWithMsg("Invalid arguments combination.");
        }
    }
    
    return 0;
}
