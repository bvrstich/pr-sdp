#ifdef PQG

#define __G_CON

#endif

#ifdef PQGT1

#define __G_CON
#define __T1_CON

#endif

#ifdef PQGT2

#define __G_CON
#define __T2_CON

#endif

#ifdef PQGT

#define __G_CON
#define __T1_CON
#define __T2_CON

#endif

#include "lapack.h"
#include "Matrix.h"
#include "Vector.h"
#include "TPM.h"
#include "SPM.h"
#include "PHM.h"
#include "DPM.h"
#include "PPHM.h"

#include "SUP.h"
#include "EIG.h"
