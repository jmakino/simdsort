#include <algorithm>
#include <functional>
extern "C"{
    void sort_int64_c(int64_t * a, int length)
    {
	std::sort(a, a+length);
    }
}
