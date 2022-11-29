#ifndef EXPORT_WANNIER90_INCLUDED
#define EXPORT_WANNIER90_INCLUDED

//#include "StdFace_vals.h"

#if defined(_UHF)

void ExportGeometry(struct StdIntList *StdI);
void ExportInteraction(struct StdIntList *StdI);

#endif

#endif /* EXPORT_WANNIER90_INCLUDED */
