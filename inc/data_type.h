#ifndef MATRICE_DATA_TYPE_H
#define MATRICE_DATA_TYPE_H

#ifdef __cplusplus
   #include <cstdint>
#else
   #include <stdint.h>
#endif

#ifdef __cplusplus
extern "C"{
#endif
// Data type
typedef enum data_type : uint32_t { knone = 0, kint = 1, kfloat = 2, kdouble = 3, kint32 = 4, kfloat32 = 5, kfloat64 = 6, kint64=7} data_type;

#ifdef __cplusplus
}
#endif

#endif
