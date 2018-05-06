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
#include <bitset> // TODO

// TODO ... -Wall -Wextra -pedantic TODO pridat?, odebrat spousteni z make
using namespace std;

class Vars {
   public:
    int bits; // requested size of public modulus in bits
    mpz_t p; // p prime number
    mpz_t q; // q prime number
    mpz_t n; // n = p * q
    mpz_t phi; // phi(n) = (p - 1) * (q - 1)
    mpz_t e; // e > 1  &&  e < phi(n)  &&  gcd(e, phi = 1
    mpz_t d;
    gmp_randstate_t rngState; // mpi rng (MT) state
    mpz_t tmpa; // helper vars for multiplicative inverse counting
    mpz_t tmpb; // helper vars for multiplicative inverse counting
    mpz_t tmps0; // helper vars for multiplicative inverse counting
    mpz_t tmps1; // helper vars for multiplicative inverse counting
};

Vars generateKey(Vars vars);
Vars generatePandQ(Vars vars);
void countMinPossibleOddPrimeValue(mpz_t* minValue, mpz_t* msbBits);
void getRandomOddNumber(mpz_t* primeNumber, unsigned long bits, Vars vars);
bool isOddNumberPrimeNumber(mpz_t *number);
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
    gmp_randinit_mt(vars.rngState);
    gmp_randseed_ui(vars.rngState, seed);

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

    // count e TODO
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
    } else {
        // odd length of key
        unsigned long bitsCount1 = vars.bits/2 + 1; // one number is 1 bit longer
        unsigned long bitsCount2 = vars.bits/2;
        cout << bitsCount1 << " " << bitsCount2 << endl; // TODO
        // prepare bit shifting, so the result N will have needed length
        unsigned long twoBits = 3UL << bitsCount1-2; // one number is 1 bit longer
        mpz_t twoBitsShifted1;
        mpz_init(twoBitsShifted1);
        mpz_set_ui(twoBitsShifted1, twoBits);

        twoBits = 3UL << bitsCount2-2; // one number is 1 bit lonshorter
        mpz_t twoBitsShifted2;
        mpz_init(twoBitsShifted2);
        mpz_set_ui(twoBitsShifted2, twoBits);

        // generate Q with two first bits set to 1, so final N will have first bit 1
        getRandomOddNumber(&vars.p, bitsCount1, vars);
        mpz_ior(vars.p, vars.p, twoBitsShifted1); // shift "1 1" to first two bits
        while (!isOddNumberPrimeNumber(&vars.p)) {
            getRandomOddNumber(&vars.p, bitsCount1, vars);
            mpz_ior(vars.p, vars.p, twoBitsShifted1); // shift "1 1" to first two bits
        }

        // generate Q with two first bits set to 1, so final N will have first bit 1
        getRandomOddNumber(&vars.q, bitsCount2, vars);
        mpz_ior(vars.q, vars.q, twoBitsShifted2);
        while (!isOddNumberPrimeNumber(&vars.q)) {
            getRandomOddNumber(&vars.q, bitsCount2, vars);
            mpz_ior(vars.q, vars.q, twoBitsShifted2);
        }
    }
    return vars;
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
            return false;
        } 
    }
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
}


/* Find e using Eucl. algorithm */
Vars findE(Vars vars) {
    // gen random e <= phi
    mpz_urandomm(vars.e, vars.rngState, vars.phi); // gen 0 to phi-1
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
                return vars;
            }
        }
        // has not found e going down to 1, try going up to phi
        mpz_set(vars.e, tmp);
        while (mpz_cmp(vars.e, vars.phi) < 0) {
            mpz_add_ui(vars.e, vars.e, 1);
            countGcd(&gcd, &vars.e, &vars.phi);
            if (mpz_cmp_ui(gcd, 1) == 0) {
                return vars;
            }
        }
    }
    // should not get here, but if yes, return 3 (given by a project definition)
    mpz_set_ui(vars.e, 3);
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

            return obtainMultiplicativeInverse(vars);
        }
}

/* exit program with message to std::err and result code 1 */
void exitWithMsg(string message) {
    cerr << message;
    exit(1);
}
