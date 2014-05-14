#include "m_rand.h"

double mode_rand() {
    static int old_iseed, new_iseed, alfa;

    time_t now;
    struct tm *l_time;
    time(&now);
    l_time = localtime(&now);

    int ry, rm, rd, rh, rn, rs;
    ry = l_time->tm_year;
    rm = l_time->tm_mon;
    rd = l_time->tm_mday;
    rh = l_time->tm_hour;
    rn = l_time->tm_min;
    rs = l_time->tm_sec;

    int iseed = ry+70*(rm+12*(rd+31*(rh+23*(rn+59*rs))));
    if ((iseed%2) == 0) {
        iseed-=1;
    }

    if (old_iseed == iseed) {
        if (!alfa) {
            new_iseed = (16807*iseed+2147483647)%2147483647;
            alfa = 1;
        } else {
            new_iseed = (16807*new_iseed+2147483647)%2147483647;
        }
        iseed = new_iseed;
    } else {
        old_iseed = iseed;
        alfa = 0;
    }

    int ia = 16807;
    int ic = 2147483647;
    int iq = 127773;
    int ir = 2836;

    int ih = iseed/iq;
    int il = iseed%iq;
    int it = ia*il-ir*ih;

    if (it > 0) {
        iseed = it;
    } else {
        iseed = ic+it;
    }
    double rr = iseed/((ic/1.0));

    return rr;
}
