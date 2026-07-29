#include <floattetwild/get_mem.h>

extern "C" size_t getPeakRSS();

size_t floatTetWild::get_peak_mem(){
    return getPeakRSS()/(1024*1024);
}
