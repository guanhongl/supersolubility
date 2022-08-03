#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_factor.h"

typedef struct {
    slong count; 
    slong MIN;
    slong MAX;
} myarg_t;

// max is 2**63 - 1 ~ 10**18
static slong pow10[19] = {
    1, 10, 100, 1000, 10000, 
    100000, 1000000, 10000000, 100000000, 1000000000, 
    10000000000, 100000000000, 1000000000000, 10000000000000, 100000000000000, 
    1000000000000000, 10000000000000000, 100000000000000000, 1000000000000000000 
};

/**
 * returns 10**n
 * https://stackoverflow.com/a/18581693
 */
slong quick_pow10(int n)
{
    return pow10[n]; 
}

/**
 * Sets f to the greatest common divisor of g and h. The result is always positive, even if one of g and h is negative
 * http://flintlib.org/doc/fmpz.html#c.fmpz_gcd
 */
void gcd(fmpz_t f, slong g, slong h)
{
    fmpz_t a, b;
    fmpz_init_set_si(a, g);
    fmpz_init_set_si(b, h);

    fmpz_gcd(f, a, b);

    fmpz_clear(a);
    fmpz_clear(b);
}

/**
 * returns 1 if n is an ss number, else return 0
 */
int is_ss(slong n)
{
    int pass = 1;
    // fmpz_factor_t consists of two fmpz vectors representing bases and exponents
    fmpz_factor_t factors;
    // Initialises an fmpz_factor_t structure
    fmpz_factor_init(factors);
    // takes a machine integer n as input
    fmpz_factor_si(factors, n);
    // the number of factors
    slong limit = factors->num;

    slong i, k, j; // the indices
    slong p_i, a_i, p_k, a_k, p_j;
    slong psi;
    // condition 1: for all i, 
    // the distinct prime factors of gcd(n, psi) and gcd(n, p_i - 1) are equal
    for (i = 0; i < limit; i++) {
        p_i = fmpz_get_si(factors->p + i);
        a_i = fmpz_get_si(factors->exp + i);
        // psi(p**k) = (p**k - 1)(p**(k-1) - 1)...(p - 1)
        psi = 1;
        for (k = a_i; k > 0; k--) {
            fmpz_t product;
            fmpz_init(product);
            fmpz_pow_ui(product, (factors->p + i), k);

            psi *= fmpz_get_si(product) - 1;

            fmpz_clear(product);
        }
        // get gcd(n, psi)
        fmpz_t gcd_1;
        fmpz_init(gcd_1);
        gcd(gcd_1, n, psi);
        // get gcd(n, p_i - 1)
        fmpz_t gcd_2;
        fmpz_init(gcd_2);
        gcd(gcd_2, n, p_i - 1);
        // get the distinct prime factors of gcd_1, gcd_2
        fmpz_factor_t fac_1, fac_2;
        fmpz_factor_init(fac_1);
        fmpz_factor_init(fac_2);
        fmpz_factor(fac_1, gcd_1);
        fmpz_factor(fac_2, gcd_2);
        // condition 1
        if (fac_1->num == fac_2->num) {
            // for any i, if the factors do not equal, then fail
            for (j = 0; j < fac_1-> num; j++) {
                if (fmpz_get_si(fac_1->p + j) != fmpz_get_si(fac_2->p + j)) {
                    pass = 0;
                    i = limit; // break
                    j = fac_1->num;
                }
            }
        }
        else {
            pass = 0;
            i = limit; // break
        }

        fmpz_clear(gcd_1);
        fmpz_clear(gcd_2);
        fmpz_factor_clear(fac_1);
        fmpz_factor_clear(fac_2);
    }
    // condition 2: if there exists p_i <= a_k (i != k)
    if (pass == 1) {
        // try all p_i <= a_k
        for (i = 0; i < limit; i++) {
            p_i = fmpz_get_si(factors->p + i);
            a_i = fmpz_get_si(factors->exp + i);
            for (k = 0; k < limit; k++) {
                p_k = fmpz_get_si(factors->p + k);
                a_k = fmpz_get_si(factors->exp + k);
                // condition 2
                if (i != k && p_i <= a_k) {
                    // condition 2a: there does not exist a prime p_j
                    // s.t. p_i divides p_j - 1 and p_j divides p_k - 1
                    for (j = 0; j < limit; j++) {
                        p_j = fmpz_get_si(factors->p + j);
                        if 
                        (
                            ((p_j - 1) % p_i == 0) &&
                            ((p_k - 1) % p_j == 0)
                        )
                        {
                            pass = 0;
                            j = limit; // break
                        }
                    }
                    // condition 2b: a_i <= 2
                    // if a_i == 2, then p_i**2 divides p_k - 1
                    if (pass == 1) {
                        if
                        (
                            (a_i == 1) ||
                            ( a_i == 2 && ((p_k - 1) % (p_i * p_i) == 0) )
                        )
                        {
                            // do nothing
                        }
                        else {
                            pass = 0;
                        }
                    }
                    
                    // break loop
                    if (pass == 0) {
                        k = limit;
                        i = limit;
                    }
                }
            }
        }
    }
    
    // Clears an fmpz_factor_t structure
    fmpz_factor_clear(factors);

    return pass;
}

/**
 * the function invoked by a thread, which takes arg, a pointer, and returns a pointer
 */
void * thread(void * arg) {	
    myarg_t * myarg = (myarg_t *) arg; // define myarg, a myarg_t pointer, and point it to the value at arg

    for (slong n = myarg->MIN; n <= myarg->MAX; n++) {
        if (is_ss(n) == 1) {
            myarg->count++;
        }
    }

    return NULL; // return NULL; terminate the thread
}

// slong (signed long) max is 9223372036854775807 or 2**63 - 1
/**
 * cmd line args: $./ss EXP NUM_THREADS or $./ss MIN MAX NUM_THREADS
 * e.g. $./ss 8 8 (10**8 max, 8 threads)
 * e.g. $./ss 2 1000 8 (2 min, 1000 max, 8 threads)
 */
int main(int argc, char* argv[]) 
{
    // the default is 10**3
    int EXP = 3; 
    slong MIN = 2;
    slong MAX = 1000;
    int NUM_THREADS = 1;

    if (argc == 2 || argc > 4) {
        printf("[ERROR] incorrect number of command line arguments.\n");

        return 1;
    }
    // get the EXP from the cmd line args
    if (argc == 3) {
        EXP = strtol(argv[1], NULL, 10); // ignore leftover
        MAX = quick_pow10(EXP);
        printf("EXP %d\n", EXP);
        flint_printf("MAX %wd\n", MAX);

        NUM_THREADS = strtol(argv[2], NULL, 10);
    }
    // get the MIN, MAX
    if (argc == 4) {
        MIN = strtol(argv[1], NULL, 10);
        MAX = strtol(argv[2], NULL, 10);
        flint_printf("MIN %wd\n", MIN);
        flint_printf("MAX %wd\n", MAX);
        EXP = 1;

        NUM_THREADS = strtol(argv[3], NULL, 10);
    }

    flint_set_num_threads(NUM_THREADS);
    printf("num_threads %d\n", flint_get_num_threads());

    // error check for max
    if (argc == 3 && EXP > 18) {
        printf("[ERROR] MAX cannot exceed 10**18.\n");

        return 1;
    }
    if (argc == 4 && MIN > MAX) {
        printf("[ERROR] MIN cannot be greater than MAX.\n");

        return 1;
    }
    fflush(stdout);

    FILE* fp = fopen("output.txt", "w");

    pthread_t threads[NUM_THREADS];
    myarg_t myargs[NUM_THREADS];

    struct timespec start, end;
    double cpu_time = 0.0;
    int e, t; // the indices
    slong n, count, max, step;
    n = MIN - 1; 
    count = 0; // the total

    flint_fprintf(fp, "MIN %wd, MAX %wd\n", MIN, MAX);
    fprintf(fp, "N\t\t\t\tcount\t\t\t\ttime (s)\n");
    // for each exponent
    for (e = 1; e <= EXP; e++) {
        max = (argc == 4) ? MAX : quick_pow10(e);
        step = (max - n) / NUM_THREADS;

        clock_gettime(CLOCK_MONOTONIC, &start);

        // for each thread
        for (t = 0; t < NUM_THREADS; t++) {
            myargs[t].count = 0;
            myargs[t].MIN = n + 1;
            myargs[t].MAX = (t == NUM_THREADS - 1) ?
                max : n + step;

            pthread_create(&threads[t], NULL, thread, &myargs[t]);

            n += step;
        }
        
        for (t = 0; t < NUM_THREADS; t++) {	
            pthread_join(threads[t], NULL); // wait for the specified thread to terminate
            count += myargs[t].count; // sum
        }

        clock_gettime(CLOCK_MONOTONIC, &end);
        cpu_time += end.tv_sec - start.tv_sec;
        cpu_time += (end.tv_nsec - start.tv_nsec) / 1000000000.0;

        fflush(fp);
        if (argc == 4) {
            flint_fprintf(fp, "%wd\t\t\t\t%wd\t\t\t\t%f\n", MAX, count, cpu_time);
        }
        else {
            flint_fprintf(fp, "10**%d\t\t\t\t%wd\t\t\t\t%f\n", e, count, cpu_time);
        }

        n = max;
    }

    flint_printf("count %wd\n", count);
    printf("cpu_time %f\n", cpu_time);

    fclose(fp);
    return 0;
}