#ifndef GLOBAL_H
#define GLOBAL_H

#include <string>
#include "random_generator.h"

namespace Global {

inline std::string project_directory;
inline Random_generator prng;

inline constexpr double pi {3.14159265358979323846};
}


#endif // GLOBAL_H
