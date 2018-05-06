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
#include <ctime>
#include <gmp.h>

using namespace std;

#define MILLER_RABIN_STEPS 20

gmp_randstate_t rngState; // mpi rng (MT) state

class Vars {
   public:
    int bits; // requested size of public modulus in bits
    mpz_t p; // p prime number
    mpz_t q; // q prime number
    mpz_t n; // n = p * q
    mpz_t phi; // phi(n) = (p - 1) * (q - 1)
    mpz_t e; // e > 1  &&  e < phi(n)  &&  gcd(e, phi = 1
    mpz_t d;
    // tmp helper variables
    mpz_t tmpa; // helper vars for multiplicative inverse counting
    mpz_t tmpb; // helper vars for multiplicative inverse counting
    mpz_t tmps0; // helper vars for multiplicative inverse counting
    mpz_t tmps1; // helper vars for multiplicative inverse counting
};

Vars generateKey(Vars vars);
Vars generatePandQ(Vars vars);
void countMinPossibleOddPrimeValue(mpz_t* minValue, mpz_t* msbBits);
void getRandomOddNumber(mpz_t* primeNumber, unsigned long bits, Vars vars);
bool isPrimeNumberMillerRabin(mpz_t *number); // doesn't work properly, ignores many primes for unknown reasons, seems to work for more than 20 bit, but takes too long
bool isOddNumberPrimeNumber(mpz_t *number); // checking primes by dividing all numbers step by step, too slow, works up to +-60 bits = 10 seconds, 50 bits or less = instant
Vars findE(Vars vars);
Vars obtainMultiplicativeInverse(Vars vars);
void exitWithMsg(string message);


int main(int argc, char *argv[])
{
    Vars vars; //object holding commonly used variables

    if (argc < 3) {
        exitWithMsg("Invalid arguments count.");
    } else {
        // generating keys
        if (argc == 3 && argv[1] == string("-g")) {
            vars.bits = stoi(argv[2]);
            vars = generateKey(vars);
            // output results to stdout
            gmp_printf ("%Zd ", vars.p);
            gmp_printf ("%Zd ", vars.q);
            gmp_printf ("%Zd ", vars.n);
            gmp_printf ("%Zd ", vars.e);
            gmp_printf ("%Zd", vars.d);
        } else {
            exitWithMsg("Invalid arguments combination.");
        }
    }
    
    return 0;
}

/* Generate RSA parameters - keys */
Vars generateKey(Vars vars) {
    // init random number generator
    srand(clock()); // init srand with proc tick count
    unsigned long seed = rand();
    gmp_randinit_mt(rngState);
    gmp_randseed_ui(rngState, seed);

    // generate p and q values
    mpz_init(vars.p);
    mpz_init(vars.q);
    vars = generatePandQ(vars);

    // count n value
    mpz_init(vars.n);
    mpz_mul(vars.n, vars.p, vars.q);

    // count phi(n) = (p - 1) * (q - 1) ...phi is always even
    mpz_t pMinusOne;
    mpz_init(pMinusOne);
    mpz_sub_ui(pMinusOne, vars.p, 1);
    mpz_t qMinusOne;
    mpz_init(qMinusOne);
    mpz_sub_ui(qMinusOne, vars.q, 1);
    mpz_init(vars.phi);
    mpz_mul(vars.phi, pMinusOne, qMinusOne);

    // count e
    mpz_init(vars.e);
    vars = findE(vars);

    // count d
    mpz_init(vars.d);
    mpz_init(vars.tmpa);
    mpz_init(vars.tmpb);
    mpz_init(vars.tmps0);
    mpz_init(vars.tmps1);
    mpz_set(vars.tmpa, vars.e);
    mpz_set(vars.tmpb, vars.phi);
    mpz_set_ui(vars.tmps0, 1);
    mpz_set_ui(vars.tmps1, 0);
    vars = obtainMultiplicativeInverse(vars);
    
    return vars;
}

/* Generate random prime number with needed bit length. */
Vars generatePandQ(Vars vars) {
    // generate p & q
    if (vars.bits % 2 == 0) {
        // even length of key
        unsigned long bitsCount = vars.bits/2;
        // prepare bit shifting, so the result N will have needed length
        unsigned long twoBits = 3UL << bitsCount-2;
        mpz_t twoBitsShifted;
        mpz_init(twoBitsShifted);
        mpz_set_ui(twoBitsShifted, twoBits);
        // count minimal possible value, to check this when generating prime number
        mpz_t minValue;
        mpz_init(minValue);
        countMinPossibleOddPrimeValue(&minValue, &twoBitsShifted);

        // generate Q with two first bits set to 1, so final N will have first bit 1
        getRandomOddNumber(&vars.p, bitsCount, vars);
        mpz_ior(vars.p, vars.p, twoBitsShifted); // shift "1 1" to first two bits
        while (!isOddNumberPrimeNumber(&vars.p)) {
            // try another odd number by substracting 2, if reached minValue, generate new number
            if (mpz_cmp(vars.p, minValue) >= 0) { // positive value if op1 > op2
                mpz_sub_ui(vars.p, vars.p, 2); // substract 2 and try again
            } else {
                // cant substract anymore, min value reached, generate a new number
                getRandomOddNumber(&vars.p, bitsCount, vars);
                countMinPossibleOddPrimeValue(&minValue, &twoBitsShifted);
            }
        }

        // generate Q with two first bits set to 1, so final N will have first bit 1
        getRandomOddNumber(&vars.q, bitsCount, vars);
        mpz_ior(vars.q, vars.q, twoBitsShifted);
        while (!isOddNumberPrimeNumber(&vars.q)) {
            // try another odd number by substracting 2, if reached minValue, generate new number
            if (mpz_cmp(vars.q, minValue) >= 0) { // positive value if op1 > op2
                mpz_sub_ui(vars.q, vars.q, 2); // substract 2 and try again
            } else {
                // cant substract anymore, min value reached, generate a new number
                getRandomOddNumber(&vars.q, bitsCount, vars);
                countMinPossibleOddPrimeValue(&minValue, &twoBitsShifted);
            }
        }
        mpz_clear(twoBitsShifted); mpz_clear(minValue);
    } else {
        // odd length of key
        unsigned long bitsCount1 = vars.bits/2 + 1; // one number is 1 bit longer
        unsigned long bitsCount2 = vars.bits/2;

        mpz_t lastBit; mpz_init(lastBit);
        mpz_set_ui(lastBit, 1);
        // generate Q with two first bits set to 1, so final N will have first bit 1
        mpz_urandomb(vars.p, rngState, bitsCount1);
        mpz_setbit(vars.p, bitsCount1 - 1);
        mpz_setbit(vars.p, bitsCount1 - 2); // shift "1 1" to first two bits
        mpz_ior(vars.p, vars.p, lastBit);
        while (!isOddNumberPrimeNumber(&vars.p)) {
            mpz_urandomb(vars.p, rngState, bitsCount1);
            mpz_setbit(vars.p, bitsCount1 - 1);
            mpz_setbit(vars.p, bitsCount1 - 2); // shift "1 1" to first two bits
            mpz_ior(vars.p, vars.p, lastBit);
        }
        
        // generate Q with two first bits set to 1, so final N will have first bit 1
        mpz_urandomb(vars.q, rngState, bitsCount1);
        mpz_setbit(vars.q, bitsCount2 - 1);
        mpz_setbit(vars.q, bitsCount2 - 2); // shift "1 1" to first two bits
        mpz_ior(vars.q, vars.q, lastBit);
        while (!isOddNumberPrimeNumber(&vars.q)) {
            mpz_urandomb(vars.q, rngState, bitsCount2);
            mpz_setbit(vars.q, bitsCount2 - 1);
            mpz_setbit(vars.q, bitsCount2 - 2);
            mpz_ior(vars.q, vars.q, lastBit);
        }
        mpz_clear(lastBit);
    }
    return vars;
}

/* Milelr-rabin test, algorithm based on pseudocode from wikipedia */
// TODO doesnt work properly, ignores many primes
bool isPrimeNumberMillerRabin(mpz_t *number) {

    mpz_t n; mpz_t d; mpz_t x; mpz_t a; mpz_t r; mpz_t s; mpz_t tmp; mpz_t mod; mpz_t tmp2; mpz_t nMinusOne;
    mpz_init(n);
    mpz_init(d);
    mpz_init(x);
    mpz_init(a);
    mpz_init(r);
    mpz_init(s);
    mpz_init(tmp);
    mpz_init(mod);
    mpz_init(tmp2);
    mpz_init(nMinusOne);
    bool continueLoop = false;

    mpz_set(n, *number);
    mpz_set(tmp, *number);
    mpz_sub_ui(tmp, tmp, 1); // n - 1
    mpz_set(nMinusOne, tmp);
    mpz_cdiv_r_ui (mod, tmp, 2);
    mpz_set_ui(s, 0);
    while (mpz_cmp_ui(mod, 0) == 0) {
        mpz_div_ui(tmp, tmp, 2);
        mpz_cdiv_r_ui(mod, tmp, 2);
        mpz_add_ui(s, s, 1); // s ... (2 ^ s) * d
    }
    
    // (2^s) * d ... d is modulo from previous cycle
    mpz_abs(d, mod);
    mpz_sub_ui(tmp2, nMinusOne, 3); // n - 4, used as random generator limit
    for (unsigned long i=0; i < MILLER_RABIN_STEPS; i++) {
        mpz_urandomm(a, rngState, tmp2); // gen 0 to n-4
        mpz_add_ui(a, a, 2); // make it 2 to n-2
        mpz_powm(x, a, d, n);
        if ((mpz_cmp_ui(x, 1) == 0) || (mpz_cmp(x, nMinusOne) == 0)) {
            continue;
        }
        
        for (unsigned long j = 1; mpz_cmp_ui(s, j+1) > 0; j++) {
            mpz_powm_ui(x, x, 2, n);
            if (mpz_cmp_ui(x, 0) == 0) {
                mpz_clear(n);mpz_clear(d);mpz_clear(x); mpz_clear(a);mpz_clear(r);mpz_clear(s);mpz_clear(tmp);mpz_clear(mod);mpz_clear(tmp2);mpz_clear(nMinusOne);
                return false;
            }
            if (mpz_cmp(x, nMinusOne) == 0) {
                continueLoop = true;
                break;
            }
        }
        if (continueLoop == true) {
            continueLoop = false;
            continue;
        }
        mpz_clear(n);mpz_clear(d);mpz_clear(x); mpz_clear(a);mpz_clear(r);mpz_clear(s);mpz_clear(tmp);mpz_clear(mod);mpz_clear(tmp2);mpz_clear(nMinusOne);
        return false;
    }
    mpz_clear(n);mpz_clear(d);mpz_clear(x); mpz_clear(a);mpz_clear(r);mpz_clear(s);mpz_clear(tmp);mpz_clear(mod);mpz_clear(tmp2);mpz_clear(nMinusOne);
    return true;
}

/* Count minimal possible value that a needed prime number P/Q can have for given bit length.
   Does OR for given msbBits and the last lsbBit (odd value), then returns the result. */
void countMinPossibleOddPrimeValue(mpz_t* minValue, mpz_t* msbBits) {
    // OR last bit
    mpz_t lastBit;
    mpz_init(lastBit);
    mpz_set_ui(lastBit, 1);
    mpz_ior(*minValue, *minValue, lastBit);
    // OR msbBits
    mpz_ior(*minValue, *minValue, *msbBits);
    mpz_clear(lastBit);
}

/* Generate random odd number.*/
void getRandomOddNumber(mpz_t* oddNumber, unsigned long bitsLength, Vars vars) {
    mpz_urandomb(*oddNumber, rngState, bitsLength);
    // set last bit to 1 = odd number
    mpz_t lastBit;
    mpz_init(lastBit);
    mpz_set_ui(lastBit, 1);
    mpz_ior(*oddNumber, *oddNumber, lastBit);
    mpz_clear(lastBit);
}

/*Check if given ODD number is a prime number.*/
bool isOddNumberPrimeNumber(mpz_t *number) {
    if (mpz_cmp_ui(*number, 1) == 0) {
        return false;
    }
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
            mpz_clear(i);mpz_clear(halfNumber);mpz_clear(modulo);
            return false;
        } 
    }
    mpz_clear(i);mpz_clear(halfNumber);mpz_clear(modulo);
    return true;
}

/* Find GCD using Eucl. algorithm */
void countGcd(mpz_t* gcd, mpz_t* x, mpz_t* y) {
    mpz_t mod; mpz_t tmp1; mpz_t tmp2;
    mpz_init(mod); mpz_init(tmp1); mpz_init(tmp2);
    mpz_set(tmp1, *x);
    mpz_set(tmp2, *y);

    mpz_cdiv_r(mod, tmp1, tmp2);
    while (mpz_cmp_ui(mod, 0) != 0) {
        mpz_set(tmp1, tmp2);
        mpz_set(tmp2, mod);
        mpz_cdiv_r(mod, tmp1, tmp2);
    }

    mpz_set(*gcd, tmp2);

    mpz_clear(mod); mpz_clear(tmp1); mpz_clear(tmp2);
}


/* Find e using Eucl. algorithm */
Vars findE(Vars vars) {
    // gen random e <= phi
    mpz_urandomm(vars.e, rngState, vars.phi); // gen 0 to phi-1
    mpz_add_ui(vars.e, vars.e, 1); // make it 1 to phi

    // GCD is 1, the E is okay
    mpz_t gcd;
    mpz_init(gcd);
    countGcd(&gcd, &vars.e, &vars.phi);
    if (mpz_cmp_ui(gcd, 1) == 0) {
        return vars;
    } else {
        mpz_t tmp;
        mpz_init(tmp);
        mpz_set(tmp, vars.e);
        // decrement e until get to 1 and try to GCD
        while (mpz_cmp_ui(vars.e, 1) > 0) {
            mpz_sub_ui(vars.e, vars.e, 1);
            countGcd(&gcd, &vars.e, &vars.phi);
            if (mpz_cmp_ui(gcd, 1) == 0) {
                mpz_clear(tmp); mpz_clear(gcd);
                return vars;
            }
        }
        // has not found e going down to 1, try going up to phi
        mpz_set(vars.e, tmp);
        while (mpz_cmp(vars.e, vars.phi) < 0) {
            mpz_add_ui(vars.e, vars.e, 1);
            countGcd(&gcd, &vars.e, &vars.phi);
            if (mpz_cmp_ui(gcd, 1) == 0) {
                mpz_clear(tmp); mpz_clear(gcd);
                return vars;
            }
        }
    }
    
    // should not get here, but if yes, return 3 (given by a project definition)
    mpz_set_ui(vars.e, 3);
    mpz_clear(gcd);
    return vars;
}

// multiplicative inverse, algorithm based on publicly available code from rosettacode.org
//Vars obtainMultiplicativeInverse(int a, int b, int s0 = 1, int s1 = 0) {
Vars obtainMultiplicativeInverse(Vars vars) {
        if (mpz_cmp_ui(vars.tmpb, 0) == 0) {
            mpz_set(vars.d, vars.tmps0);
            return vars;
        } else {
            // obtainMultiplicativeInverse(tmpa, tmpb, s0, s1);
            // turn row above into row below:
            // obtainMultiplicativeInverse(tmpb, tmpa%tmpb, s1, s0 - s1*(par1/par2))
            mpz_t tmp1; mpz_t tmp2; mpz_t tmp3; mpz_t tmp4;
            mpz_init(tmp1); mpz_init(tmp2); mpz_init(tmp3); mpz_init(tmp4);
            mpz_set(tmp1, vars.tmpa);
            mpz_set(tmp2, vars.tmpb);
            mpz_set(tmp3, vars.tmps0);
            mpz_set(tmp4, vars.tmps1);
            // par1
            mpz_set(vars.tmpa, tmp2);

            // par2
            mpz_cdiv_r (vars.tmpb, tmp1, tmp2);

            // s0
            mpz_set(vars.tmps0, tmp4);

            // s1
            mpz_t tmp;
            mpz_init(tmp);
            mpz_div(tmp, tmp1, tmp2);
            mpz_mul(tmp, tmp4, tmp);
            mpz_sub(vars.tmps1, tmp3, tmp);

            mpz_clear(tmp1); mpz_clear(tmp2); mpz_clear(tmp3); mpz_clear(tmp4); mpz_clear(tmp);
            return obtainMultiplicativeInverse(vars);
        }
}

/* exit program with message to std::err and result code 1 */
void exitWithMsg(string message) {
    cerr << message;
    exit(1);
}
