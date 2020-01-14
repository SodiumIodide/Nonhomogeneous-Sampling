#ifndef _CONFIG_H
#include <config.h>
#endif

#ifdef S2
    #ifndef _S2_H
        #include "s2.h"
    #endif
    #define SOLVE() do {\
        s2();\
    } while (0)
#endif

#ifdef PUREABS
    #ifndef _PURE_ABS_H
        #include "pure_abs.h"
    #endif
    #define SOLVE() do {\
        pure_abs();\
    } while (0)
#endif

#ifdef GEOMETRYTIME
    #ifndef _GEOMETRY_TIME_H
        #include "geometry_time.h"
    #endif
    #define SOLVE() do {\
        geometry_time();\
    } while (0)
#endif

#ifdef ANALYTICAL
    #ifndef _ANALYTICAL_H
        #include "analytical.h"
    #endif
    #define SOLVE() do {\
        analytical();\
    } while (0)
#endif

#ifdef VOLFRAC
    #ifndef _VOLUME_FRACTION_H
        #include "volume_fraction.h"
    #endif
    #define SOLVE() do {\
        volume_fraction();\
    } while (0)
#endif

int main(int argc, char *argv[]) {
    (void)argc;
    (void)argv;

    SOLVE();

    return 0;
}
