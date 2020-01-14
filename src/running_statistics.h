#ifndef _RUNNING_STATISTICS_H
#define _RUNNING_STATISTICS_H

typedef struct runningstat {
    long m_n;
    double m_oldM;
    double m_newM;
    double m_oldS;
    double m_newS;
    double m_min;
    double m_max;
} runningstat;

runningstat init_runningstat() {
    runningstat r;
    r.m_n = 0;
    r.m_oldM = 0.0;
    r.m_newM = 0.0;
    r.m_oldS = 0.0;
    r.m_newS = 0.0;
    r.m_min = 0.0;
    r.m_max = 0.0;

    return r;
}

void clear(runningstat *r) {
    r->m_n = 0;
    r->m_oldM = r->m_newM = r->m_oldS = r->m_newS = r->m_min = r->m_max = 0.0;
}

void push(runningstat *r, double x) {
    r->m_n += 1;

    // See Knuth TAOCP vol 2, 3rd edition, page 232
    if (r->m_n == 1) {
        r->m_oldM = r->m_newM = x;
        r->m_oldS = 0.0;
        r->m_min = r->m_max = x;
    } else {
        r->m_newM = r->m_oldM + (x - r->m_oldM) / r->m_n;
        r->m_newS = r->m_oldS + (x - r->m_oldM) * (x - r->m_newM);

        // Set up for next iteration
        r->m_oldM = r->m_newM;
        r->m_oldS = r->m_newS;

        // Min and max
        if (x < r->m_min) {
            r->m_min = x;
        } else if (x > r->m_max) {
            r->m_max = x;
        }
    }
}

long num(runningstat *r) {
    return r->m_n;
}

double mean(runningstat *r) {
    return (r->m_n > 0) ? r->m_newM : 0.0;
}

double variance(runningstat *r) {
    return (r->m_n > 1) ? (r->m_newS / (r->m_n - 1)) : 0.0;
}

double least(runningstat *r) {
    return r->m_min;
}

double greatest(runningstat *r) {
    return r->m_max;
}

#endif
