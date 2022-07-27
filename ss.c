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
 * returns 2**n
 */
// slong quick_pow2(int n)
// {
//     // max is 2**62
//     // slong pow2 = 1 << n;
//     slong pow2 = 0x4000000000000000; // 0100 ... 0000
//     pow2 = pow2 >> (62 - n);

//     return pow2;
// }

/**
 * Sets f to the greatest common divisor of g and h. The result is always positive, even if one of g and h is negative
 * http://flintlib.org/doc/fmpz.html#c.fmpz_gcd
 */
void get_gcd(fmpz_t f, slong g, slong h)
{
    fmpz_t a, b;
    fmpz_init_set_si(a, g);
    fmpz_init_set_si(b, h);

    fmpz_gcd(f, a, b);

    fmpz_clear(a);
    fmpz_clear(b);
}

/**
 * if any product in products is divisble by n, return n
 * else return 1
 */
// slong list_gcd(int *products, slong n, slong size)
// {
//     for (slong i = 0; i < size; i++) {
//         if (products[i] % n == 0) {
//             return n;
//         }
//     }

//     return 1;
// }

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

    // try all combinations of p_i and p_j
    slong i, j, e, k; // the indices
    slong e_i, e_j, p_i, p_j, p_k;
    slong list_gcd;
    for (i = 0; i < limit; i++) {
        for (j = 0; j < limit; j++) {
            if (j != i) {
                e_i = fmpz_get_si(factors->exp + i); // the ith exponent
                e_j = fmpz_get_si(factors->exp + j); // the jth exponent
                p_i = fmpz_get_si(factors->p + i); // the ith factor
                p_j = fmpz_get_si(factors->p + j); // the jth factor
                // allocate mem for products
                // e.g. products for 7**2 = [6, 48]
                // int *products = malloc(e_j * sizeof(slong));
                list_gcd = 1;

                for (e = 1; e <= e_j; e++) {
                    // product = the jth factor ** e
                    fmpz_t product;
                    fmpz_init(product);
                    fmpz_pow_ui(product, (factors->p + j), e);
                    // append the product - 1
                    // products[e - 1] = fmpz_get_si(product) - 1;
                    if ((fmpz_get_si(product) - 1) % p_i == 0) {
                        list_gcd = p_i;
                        e = e_j + 1; // break
                    }

                    fmpz_clear(product);
                }

                fmpz_t gcd;
                fmpz_init(gcd);
                get_gcd(gcd, p_i, p_j - 1);
                // CONDITION ONE
                // if (list_gcd(products, p_i, e_j) == fmpz_get_si(gcd)) {
                if (list_gcd == fmpz_get_si(gcd)) {
                    // CONDITION TWO
                    if (p_i <= e_j) {
                        // CONDITION THREE
                        if (e_i == 1 || e_i == 2) {
                            // product = the ith factor ** the ith exp
                            fmpz_t product;
                            fmpz_init(product);
                            fmpz_pow_ui(product, (factors->p + i), e_i);
                            // CONDITION FOUR
                            if ((p_j - 1) % fmpz_get_si(product) == 0) {
                                // CONDITION FIVE
                                for (k = i + 1; k < j; k++) {
                                    p_k = fmpz_get_si(factors->p + k); // the (i + 1)th factor

                                    if 
                                    (
                                        ((p_k - 1) % p_i == 0) &&
                                        ((p_j - 1) % p_k == 0)
                                    )
                                    {
                                        pass = 0;
                                        k = j; // break
                                    }
                                }
                            }
                            // CONDITION FOUR
                            else {
                                pass = 0;
                            }

                            fmpz_clear(product);
                        } 
                        // CONDITION THREE
                        else {
                            pass = 0;
                        }
                    }
                    // CONDITION TWO
                    // continue if CON ONE is true and CON TWO is false
                } 
                // CONDITION ONE
                else {
                    pass = 0;
                }

                fmpz_clear(gcd);
                // free(products);

                // break loops if fail
                if (pass == 0) {
                    i = limit;
                    j = limit;
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
        EXP = atoi(argv[1]);
        MAX = quick_pow10(EXP);
        printf("EXP %d\n", EXP);
        flint_printf("MAX %wd\n", MAX);

        NUM_THREADS = atoi(argv[2]);
    }
    // get the MIN, MAX
    if (argc == 4) {
        MIN = atol(argv[1]);
        MAX = atol(argv[2]);
        flint_printf("MIN %wd\n", MIN);
        flint_printf("MAX %wd\n", MAX);
        EXP = 1;

        NUM_THREADS = atoi(argv[3]);
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

            // flint_printf("t %d, MIN %wd, MAX %wd\n", t, myargs[t].MIN, myargs[t].MAX);

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

        // flint_fprintf(fp, "%wd\t\t\t\t%wd\t\t\t\t%f\n", quick_pow10(e), count, cpu_time);
        if (argc == 4) {
            flint_fprintf(fp, "%wd\t\t\t\t%wd\t\t\t\t%f\n", MAX, count, cpu_time);
        }
        else {
            flint_fprintf(fp, "10**%d\t\t\t\t%wd\t\t\t\t%f\n", e, count, cpu_time);
        }
        fflush(fp);

        n = max;
    }

    flint_printf("count %wd\n", count);
    printf("cpu_time %f\n", cpu_time);

    fclose(fp);
    return 0;
}
