//nog enkele definities:
#ifdef PQ

#define __Q_CON

#endif

#ifdef PQG

#define __Q_CON
#define __G_CON

#endif

#ifdef PQGT1

#define __Q_CON
#define __G_CON
#define __T1_CON

#endif

#include "lapack.h"
#include "Matrix.h"
#include "TPM.h"
#include "SPM.h"
#include "PHM.h"
#include "DPM.h"

#include "SUP.h"
#include "EIG.h"
