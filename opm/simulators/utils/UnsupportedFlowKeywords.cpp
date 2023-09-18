/*
  Copyright 2021 Equinor.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/simulators/utils/UnsupportedFlowKeywords.hpp>

namespace Opm::FlowKeywordValidation
{

const KeywordValidation::UnsupportedKeywords& unsupportedKeywords()
{
    static const KeywordValidation::UnsupportedKeywords unsupported_keywords = {
        {"ACTION", {true, std::nullopt}},
        {"ACTIONG", {true, std::nullopt}},
        {"ACTIONR", {true, std::nullopt}},
        {"ACTIONS", {true, std::nullopt}},
        {"ACTIONW", {true, std::nullopt}},
        {"ACTPARAM", {true, std::nullopt}},
        {"ADSALNOD", {true, std::nullopt}},
        {"ADDZCORN", {true, std::nullopt}},
        {"ADSORP", {true, std::nullopt}},
        {"AITS", {false, std::nullopt}},
        {"AITSOFF", {false, std::nullopt}},
        {"ALKADS", {true, std::nullopt}},
        {"ALKALINE", {true, std::nullopt}},
        {"ALKROCK", {true, std::nullopt}},
        {"API", {true, std::nullopt}},
        {"ALPOLADS", {true, std::nullopt}},
        {"ALSURFAD", {true, std::nullopt}},
        {"ALSURFST", {true, std::nullopt}},
        {"AMALGAM", {true, std::nullopt}},
        {"APIGROUP", {true, std::nullopt}},
        {"APILIM", {true, std::nullopt}},
        {"APIVD", {true, std::nullopt}},
        {"AQANCONL", {true, std::nullopt}},
        {"AQANNC", {true, std::nullopt}},
        {"AQANTRC", {true, std::nullopt}},
        {"AQUALIST", {true, std::nullopt}},
        {"AQUCHGAS", {true, std::nullopt}},
        {"AQUCHWAT", {true, std::nullopt}},
        {"AQUCWFAC", {true, std::nullopt}},
        {"AQUFET", {true, std::string{"Use the AQUFETP keyword instead"}}},
        {"AQUNNC", {true, std::nullopt}},
        {"AUTOCOAR", {true, std::nullopt}},
        {"AUTOREF", {true, std::nullopt}},
        {"BIGMODEL", {false, std::nullopt}},
        {"BDENSITY", {true, std::nullopt}},
        {"BGGI", {true, std::nullopt}},
        {"BOGI", {true, std::nullopt}},
        {"BOUNDARY", {true, std::nullopt}},
        {"BPARA", {true, std::nullopt}},
        {"BPIDIMS", {true, std::nullopt}},
        {"BTOBALFA", {true, std::nullopt}},
        {"BTOBALFV", {true, std::nullopt}},
        {"CALTRAC", {true, std::nullopt}},
        {"CARFIN", {true, std::nullopt}},
        {"CART", {true, std::nullopt}},
        {"CBMOPTS", {true, std::nullopt}},
        {"CECON", {true, std::nullopt}},
        {"CECONT", {true, std::nullopt}},
        {"COAL", {true, std::nullopt}},
        {"COALADS", {true, std::nullopt}},
        {"COALNUM", {true, std::nullopt}},
        {"COALPP", {true, std::nullopt}},
        {"COARSEN", {true, std::nullopt}},
        {"COLLAPSE", {true, std::nullopt}},
        {"COLUMNS", {true, std::nullopt}},
        {"CBMOPTS", {true, std::nullopt}},
        {"COMPDATX", {true, std::nullopt}},
        {"COMPDATL", {true, std::nullopt}},
        {"COMPDATM", {true, std::nullopt}},
        {"COMPDATL", {true, std::nullopt}},
        {"COMPIMB", {true, std::nullopt}},
        {"COMPFLSH", {true, std::nullopt}},
        {"COMPINJK", {true, std::nullopt}},
        {"COMPLMPL", {true, std::nullopt}},
        {"COMPOFF", {true, std::nullopt}},
        {"COMPRIV", {true, std::nullopt}},
        {"COMPRP", {true, std::nullopt}},
        {"COMPRPL", {true, std::nullopt}},
        {"COMPSEGL", {true, std::nullopt}},
        {"COMPVE", {true, std::nullopt}},
        {"COMPVEL", {true, std::nullopt}},
        {"CPIFACT", {true, std::nullopt}},
        {"CPIFACTL", {true, std::nullopt}},
        {"CSKIN", {true, std::nullopt}},
        {"CONNECTION", {true, std::nullopt}},
        {"CONNECTION_PROBE", {true, std::nullopt}},
        {"COORDSYS", {true, std::nullopt}},
        {"COPYBOX", {true, std::nullopt}},
        {"CRITPERM", {true, std::nullopt}},
        {"DATUMR", {true, std::nullopt}},
        {"DATUMRX", {true, std::nullopt}},
        {"DCQDEFN", {true, std::nullopt}},
        {"DEBUG", {false, std::nullopt}},
        {"DELAYACT", {true, std::nullopt}},
        {"DEPTHTAB", {true, std::nullopt}},
        {"DIAGDISP", {true, std::nullopt}},
        {"DIFF", {true, std::nullopt}},
        {"DIFFCOAL", {true, std::nullopt}},
        {"DIFFDP", {true, std::nullopt}},
        {"DIFFMMF", {true, std::nullopt}},
        {"DIFFMR", {true, std::nullopt}},
        {"DIFFMR-", {true, std::nullopt}},
        {"DIFFMTHT", {true, std::nullopt}},
        {"DIFFMTH-", {true, std::nullopt}},
        {"DIFFMX", {true, std::nullopt}},
        {"DIFFMX-", {true, std::nullopt}},
        {"DIFFMY", {true, std::nullopt}},
        {"DIFFMY-", {true, std::nullopt}},
        {"DIFFMZ", {true, std::nullopt}},
        {"DIFFMZ-", {true, std::nullopt}},
        {"DIFFR", {true, std::nullopt}},
        {"DIFFTHT", {true, std::nullopt}},
        {"DIFFX", {true, std::nullopt}},
        {"DIFFY", {true, std::nullopt}},
        {"DIFFZ", {true, std::nullopt}},
        {"DIMPES", {true, std::nullopt}},
        {"DIMPLICT", {true, std::nullopt}},
        {"DISPDIMS", {true, std::nullopt}},
        {"DISPERSE", {true, std::nullopt}},
        {"DOMAINS", {true, std::nullopt}},
        {"DPGRID", {true, std::nullopt}},
        {"DPKRMOD", {true, std::nullopt}},
        {"DPNUM", {true, std::nullopt}},
        {"DR", {true, std::string{"Use the DRV keyword instead"}}},
        {"DRILPRI", {true, std::nullopt}},
        {"DSPDEINT", {true, std::nullopt}},
        {"DTHETA", {true, std::string{"Use the DTHETAV keyword instead"}}},
        {"DUALPERM", {true, std::nullopt}},
        {"DUALPORO", {true, std::nullopt}},
        {"DUMPCUPL", {true, std::nullopt}},
        {"DUMPFLUX", {true, std::nullopt}},
        {"DYNAMICR", {true, std::nullopt}},
        {"DYNRDIMS", {true, std::nullopt}},
        {"DYNAMICR", {true, std::nullopt}},
        {"DZMATRIX", {true, std::nullopt}},
        {"DZMTRX", {true, std::nullopt}},
        {"DZMTRXV", {true, std::nullopt}},
        {"DZNET", {true, std::nullopt}},
        {"ECHO", {false, std::nullopt}},
        {"ECLMC", {true, std::nullopt}},
        {"EHYSTRR", {true, std::nullopt}},
        {"ENDDYN", {true, std::nullopt}},
        {"ENDFIN", {true, std::nullopt}},
        {"ENDNUM", {true, std::nullopt}},
        {"ENDSKIP", {true, std::nullopt}},
        {"ENKRVD", {true, std::nullopt}},
        {"ENKSRVD", {true, std::nullopt}},
        {"ENPCVD", {true, std::nullopt}},
        {"ENPTVD", {true, std::nullopt}},
        {"ENSPCVD", {true, std::nullopt}},
        {"EPSDBGS", {false, std::nullopt}},
        {"EPSDEBUG", {false, std::nullopt}},
        {"EQLZCORN", {true, std::nullopt}},
        {"EQUALREG", {true, std::nullopt}},
        {"ESSNODE", {true, std::nullopt}},
        {"EXCAVATE", {true, std::nullopt}},
        {"EXCEL", {false, std::nullopt}},
        {"EXTFIN", {true, std::nullopt}},
        {"EXTHOST", {true, std::nullopt}},
        {"EXTRAPMS", {false, std::nullopt}},
        {"EXTREPGL", {true, std::nullopt}},
        {"FBHPDEF", {true, std::nullopt}},
        {"FHERCHBL", {true, std::nullopt}},
        {"FRICTION", {true, std::nullopt}},
        {"FIPOWG", {false, std::string{"Report not available"}}},
        {"FIPSEP", {true, std::nullopt}},
        {"FLUXREG", {true, std::nullopt}},
        {"FLUXTYPE", {true, std::nullopt}},
        {"FMTHMD", {true, std::nullopt}},
        {"FOAMDCYO", {true, std::nullopt}},
        {"FOAMDCYW", {true, std::nullopt}},
        {"FOAMFCN", {true, std::nullopt}},
        {"FOAMFRM", {true, std::nullopt}},
        {"FOAMFSO", {true, std::nullopt}},
        {"FOAMFST", {true, std::nullopt}},
        {"FOAMFSW", {true, std::nullopt}},
        {"FOAMMOBP", {true, std::nullopt}},
        {"FOAMMOBS", {true, std::nullopt}},
        {"FORMFEED", {true, std::nullopt}},
        {"FULLIMP", {false, std::nullopt}},
        {"GEGONT", {true, std::nullopt}},
        {"GETDATA", {true, std::nullopt}},
        {"GASBEGIN", {true, std::nullopt}},
        {"GASCONC", {true, std::nullopt}},
        {"GASEND", {true, std::nullopt}},
        {"GASFCOMP", {true, std::nullopt}},
        {"GASFDECR", {true, std::nullopt}},
        {"GASFDELC", {true, std::nullopt}},
        {"GASFIELD", {true, std::nullopt}},
        {"GASFTARG", {true, std::nullopt}},
        {"GASMONTH", {true, std::nullopt}},
        {"GASPERIO", {true, std::nullopt}},
        {"GASSATC", {true, std::nullopt}},
        {"GASYEAR", {true, std::nullopt}},
        {"GCALECON", {true, std::nullopt}},
        {"GCONCAL", {true, std::nullopt}},
        {"GCONENG", {true, std::nullopt}},
        {"GCONPRI", {true, std::nullopt}},
        {"GCONTOL", {true, std::nullopt}},
        {"GCUTBACK", {true, std::nullopt}},
        {"GCUTBACT", {true, std::nullopt}},
        {"GCVD", {true, std::nullopt}},
        {"GDCQ", {true, std::nullopt}},
        {"GDCQECON", {true, std::nullopt}},
        {"GDIMS", {true, std::nullopt}},
        {"GDORIENT", {false, std::nullopt}},
        {"GDRILPOT", {true, std::nullopt}},
        {"GECONT", {true, std::nullopt}},
        {"GETGLOB", {true, std::nullopt}},
        {"GI", {true, std::nullopt}},
        {"GIALL", {true, std::nullopt}},
        {"GIMODEL", {true, std::nullopt}},
        {"GINODE", {true, std::nullopt}},
        {"GLIFTLIM", {true, std::nullopt}},
        {"GNETDP", {true, std::nullopt}},
        {"GNETINJE", {true, std::nullopt}},
        {"GNETPUMP", {true, std::nullopt}},
        {"GRADGRUP", {true, std::nullopt}},
        {"GRADRESV", {true, std::nullopt}},
        {"GRADRFT", {true, std::nullopt}},
        {"GRADWELL", {true, std::nullopt}},
        {"GRAVCONS", {true, std::nullopt}},
        {"GRAVDR", {true, std::nullopt}},
        {"GRAVDRB", {true, std::nullopt}},
        {"GRAVDRM", {true, std::nullopt}},
        {"GRDREACH", {true, std::nullopt}},
        {"GRUPMAST", {true, std::nullopt}},
        {"GRUPRIG", {true, std::nullopt}},
        {"GRUPSLAV", {true, std::nullopt}},
        {"GRUPTARG", {true, std::nullopt}},
        {"GSATINJE", {true, std::nullopt}},
        {"GSEPCOND", {true, std::nullopt}},
        {"GSSCPTST", {true, std::nullopt}},
        {"GSWINGF", {true, std::nullopt}},
        {"GTADD", {true, std::nullopt}},
        {"GTMULT", {true, std::nullopt}},
        {"GUIDECAL", {true, std::nullopt}},
        {"GSATPROD", {true, std::nullopt}},
        {"GUPFREQ", {true, std::nullopt}},
        {"GWRTWCV", {true, std::nullopt}},
        {"HALFTRAN", {true, std::nullopt}},
        {"HAxxxxxx", {true, std::nullopt}},
        {"HBNUM", {true, std::nullopt}},
        {"HDISP", {true, std::nullopt}},
        {"HMAQUCT", {true, std::nullopt}},
        {"HMAQUFET", {true, std::nullopt}},
        {"HMAQUNUM", {true, std::nullopt}},
        {"HMDIMS", {true, std::nullopt}},
        {"HMFAULTS", {true, std::nullopt}},
        {"HMMLAQUN", {true, std::nullopt}},
        {"HMMLCTAQ", {true, std::nullopt}},
        {"HMMLFTAQ", {true, std::nullopt}},
        {"HMMLTWCN", {true, std::nullopt}},
        {"HMMULTxx", {true, std::nullopt}},
        {"HMMULTFT", {true, std::nullopt}},
        {"HMMULTSG", {true, std::nullopt}},
        {"HMPROPS", {true, std::nullopt}},
        {"HMROCK", {true, std::nullopt}},
        {"HMROCKT", {true, std::nullopt}},
        {"HMRREF", {true, std::nullopt}},
        {"HMWELCON", {true, std::nullopt}},
        {"HMWPIMLT", {true, std::nullopt}},
        {"HMxxxxxx", {true, std::nullopt}},
        {"HRFIN", {true, std::nullopt}},
        {"HWKRO", {true, std::nullopt}},
        {"HWKRORG", {true, std::nullopt}},
        {"HWKRORW", {true, std::nullopt}},
        {"HWKRW", {true, std::nullopt}},
        {"HWKRWR", {true, std::nullopt}},
        {"HWPCW", {true, std::nullopt}},
        {"HWSNUM", {true, std::nullopt}},
        {"HWSOGCR", {true, std::nullopt}},
        {"HWSOWCR", {true, std::nullopt}},
        {"HWSWCR", {true, std::nullopt}},
        {"HWSWL", {true, std::nullopt}},
        {"HWSWLPC", {true, std::nullopt}},
        {"HWSWU", {true, std::nullopt}},
        {"HXFIN", {true, std::nullopt}},
        {"HYDRHEAD", {true, std::nullopt}},
        {"HYFIN", {true, std::nullopt}},
        {"HYMOBGDR", {true, std::nullopt}},
        {"HYST", {true, std::nullopt}},
        {"HYSTCHCK", {true, std::nullopt}},
        {"HZFIN", {true, std::nullopt}},
        {"IHOST", {true, std::nullopt}},
        {"IMBNUMMF", {true, std::nullopt}},
        {"IMKRVD", {true, std::nullopt}},
        {"IMPCVD", {true, std::nullopt}},
        {"IMPES", {true, std::nullopt}},
        {"IMPLICIT", {true, std::nullopt}},
        {"IMPTVD", {true, std::nullopt}},
        {"IMSPCVD", {true, std::nullopt}},
        {"INSPEC", {true, std::nullopt}},
        {"INTPC", {true, std::nullopt}},
        {"IONROCK", {true, std::nullopt}},
        {"IONXROCK", {true, std::nullopt}},
        {"IONXSURF", {true, std::nullopt}},
        {"ISOLNUM", {true, std::nullopt}},
        {"JFUNCR", {true, std::nullopt}},
        {"KRNUMMF", {true, std::nullopt}},
        {"KRNUMR-", {true, std::nullopt}},
        {"KRNUMT-", {true, std::nullopt}},
        {"KRNUMX-", {true, std::nullopt}},
        {"KRNUMY-", {true, std::nullopt}},
        {"KRNUMZ-", {true, std::nullopt}},
        {"LANGMPL", {true, std::nullopt}},
        {"LANGMUIR", {true, std::nullopt}},
        {"LANGSOLV", {true, std::nullopt}},
        {"LCUNIT", {true, std::nullopt}},
        {"LGR", {true, std::nullopt}},
        {"LGRCOPY", {true, std::nullopt}},
        {"LGRFREE", {true, std::nullopt}},
        {"LGRLOCK", {true, std::nullopt}},
        {"LGROFF", {true, std::nullopt}},
        {"LGRON", {true, std::nullopt}},
        {"LICENSES", {true, std::nullopt}},
        {"LINCOM", {true, std::nullopt}},
        {"LINKPERM", {true, std::nullopt}},
        {"LKRO", {true, std::nullopt}},
        {"LKRORG", {true, std::nullopt}},
        {"LKRORW", {true, std::nullopt}},
        {"LKRW", {true, std::nullopt}},
        {"LKRWR", {true, std::nullopt}},
        {"LOAD", {true, std::nullopt}},
        {"LOWSALT", {true, std::nullopt}},
        {"LPCW", {true, std::nullopt}},
        {"LSALTFNC", {true, std::nullopt}},
        {"LSLTWNUM", {true, std::nullopt}},
        {"LSNUM", {true, std::nullopt}},
        {"LSOGCR", {true, std::nullopt}},
        {"LSOWCR", {true, std::nullopt}},
        {"LSWCR", {true, std::nullopt}},
        {"LSWL", {true, std::nullopt}},
        {"LSWLPC", {true, std::nullopt}},
        {"LSWU", {true, std::nullopt}},
        {"LTOSIGMA", {true, std::nullopt}},
        {"LWKRO", {true, std::nullopt}},
        {"LWKRORG", {true, std::nullopt}},
        {"LWKRORW", {true, std::nullopt}},
        {"LWKRW", {true, std::nullopt}},
        {"LWKRWR", {true, std::nullopt}},
        {"LWPCW", {true, std::nullopt}},
        {"LWSLTNUM", {true, std::nullopt}},
        {"LWSNUM", {true, std::nullopt}},
        {"LWSOGCR", {true, std::nullopt}},
        {"LWSOWCR", {true, std::nullopt}},
        {"LWSWCR", {true, std::nullopt}},
        {"LWSWL", {true, std::nullopt}},
        {"LWSWLPC", {true, std::nullopt}},
        {"LWSWU", {true, std::nullopt}},
        {"LX", {true, std::nullopt}},
        {"LXFIN", {true, std::nullopt}},
        {"LY", {true, std::nullopt}},
        {"LYFIN", {true, std::nullopt}},
        {"LZ", {true, std::nullopt}},
        {"LZFIN", {true, std::nullopt}},
        {"MASSFLOW", {true, std::nullopt}},
        {"MATCORR", {true, std::nullopt}},
        {"MEMORY", {false, std::nullopt}},
        {"MESSAGE", {false, std::nullopt}},
        {"MESSOPTS", {false, std::nullopt}},
        {"MESSSRVC", {true, std::nullopt}},
        {"MINNNCT", {true, std::nullopt}},
        {"MLANG", {true, std::nullopt}},
        {"MLANGSLV", {true, std::nullopt}},
        {"MONITOR", {false, std::nullopt}},
        {"MPFANUM", {true, std::nullopt}},
        {"MPFNNC", {true, std::nullopt}},
        {"MSGFILE", {false, std::nullopt}},
        {"MULSGGD", {true, std::nullopt}},
        {"MULSGGDV", {true, std::nullopt}},
        {"MULTOUTS", {true, std::nullopt}},
        {"MULTREAL", {true, std::nullopt}},
        {"MULTREGD", {true, std::nullopt}},
        {"MULTREGH", {true, std::nullopt}},
        {"MULTSIG", {true, std::nullopt}},
        {"MULTSIGV", {true, std::nullopt}},
        {"MULT_XYZ", {true, std::nullopt}},
        {"NARROW", {true, std::nullopt}},
        {"NCONSUMP", {true, std::nullopt}},
        {"NEFAC", {true, std::nullopt}},
        {"NETCOMPA", {true, std::nullopt}},
        {"NEXT", {false, std::nullopt}},
        {"NEXTSTPL", {true, std::nullopt}},
        {"NINENUM", {true, std::nullopt}},
        {"NINEPOIN", {true, std::nullopt}},
        {"NMATOPTS", {true, std::nullopt}},
        {"NMATRIX", {true, std::nullopt}},
        {"NODPPM", {true, std::nullopt}},
        {"NOECHO", {false, std::nullopt}},
        {"NOHMD", {true, std::nullopt}},
        {"NOHMO", {true, std::nullopt}},
        {"NOHYST", {true, std::nullopt}},
        {"NOWARNEP", {false, std::nullopt}},
        {"NRSOUT", {true, std::nullopt}},
        {"NNEWTF", {true, std::nullopt}},
        {"NOCASC", {true, std::nullopt}},
        {"NOGGF", {true, std::nullopt}},
        {"NOINSPEC", {false, std::nullopt}},
        {"NOMONITO", {true, std::nullopt}},
        {"NONNC", {true, std::nullopt}},
        {"NORSSPEC", {false, std::nullopt}},
        {"NOWARN", {false, std::nullopt}},
        {"NSTACK", {false, std::nullopt}},
        {"NWATREM", {true, std::nullopt}},
        {"NXFIN", {true, std::nullopt}},
        {"NYFIN", {true, std::nullopt}},
        {"NZFIN", {true, std::nullopt}},
        {"OFM", {true, std::nullopt}},
        {"OILAPI", {true, std::nullopt}},
        {"OLDTRAN", {true, std::nullopt}},
        {"OLDTRANR", {true, std::nullopt}},
        {"OPTIONS", {true, std::nullopt}},
        {"OUTRAD", {true, std::string{"Use the DRV keyword instead"}}},
        {"OUTSOL", {false, std::nullopt}},
        {"PARAOPTS", {false, std::nullopt}},
        {"PCG32D", {true, std::nullopt}},
        {"PCW32D", {true, std::nullopt}},
        {"PERMJFUN", {true, std::nullopt}},
        {"PETOPTS", {true, std::nullopt}},
        {"PLYESAL", {true, std::nullopt}},
        {"PLYKRRF", {true, std::nullopt}},
        {"PLYOPTS", {true, std::nullopt}},
        {"PLYRMDEN", {true, std::nullopt}},
        {"PLYROCKM", {true, std::nullopt}},
        {"PLYTRRF", {true, std::nullopt}},
        {"PLYTRRFA", {true, std::nullopt}},
        {"PLYVISCS", {true, std::nullopt}},
        {"PLYVISCT", {true, std::nullopt}},
        {"PLYVSCST", {true, std::nullopt}},
        {"PVZG", {true, std::nullopt}},
        {"PMAX", {true, std::nullopt}},
        {"PRIORITY", {true, std::nullopt}},
        {"PSTEADY", {true, std::nullopt}},
        {"PSWRG", {true, std::nullopt}},
        {"PSWRO", {true, std::nullopt}},
        {"PVCO", {true, std::nullopt}},
        {"PVZG", {true, std::nullopt}},
        {"QDRILL", {true, std::nullopt}},
        {"QDRILL", {true, std::nullopt}},
        {"QHRATING", {true, std::nullopt}},
        {"QMOBIL", {true, std::nullopt}},
        {"PARALLEL", {false, std::nullopt}},
        {"PARTTRAC", {true, std::nullopt}},
        {"PBUB", {true, std::nullopt}},
        {"PCG", {true, std::nullopt}},
        {"PCW", {true, std::nullopt}},
        {"PDEW", {true, std::nullopt}},
        {"PEBI", {true, std::nullopt}},
        {"PECOEFS", {true, std::nullopt}},
        {"PEDIMS", {true, std::nullopt}},
        {"PEGTABX", {true, std::nullopt}},
        {"PEKTABX", {true, std::nullopt}},
        {"PENUM", {true, std::nullopt}},
        {"PERMAVE", {true, std::nullopt}},
        {"PERMXY", {true, std::nullopt}},
        {"PERMYZ", {true, std::nullopt}},
        {"PERMZX", {true, std::nullopt}},
        {"PETGRID", {true, std::nullopt}},
        {"PICOND", {true, std::nullopt}},
        {"PIMULTAB", {true, std::nullopt}},
        {"PINCHNUM", {true, std::nullopt}},
        {"PINCHOUT", {true, std::nullopt}},
        {"PINCHREG", {true, std::nullopt}},
        {"PINCHXY", {true, std::nullopt}},
        {"PLYADSS", {true, std::nullopt}},
        {"PLYATEMP", {true, std::nullopt}},
        {"PLYCAMAX", {true, std::nullopt}},
        {"PLYDHFLF", {true, std::nullopt}},
        {"PRORDER", {true, std::nullopt}},
        {"PRVD", {true, std::nullopt}},
        {"PVTGWO", {true, std::nullopt}},
        {"RAINFALL", {true, std::nullopt}},
        {"RBEDCONT", {true, std::nullopt}},
        {"RADFIN", {true, std::nullopt}},
        {"RADFIN4", {true, std::nullopt}},
        {"RCMASTS", {true, std::nullopt}},
        {"REACACT", {true, std::nullopt}},
        {"REACHES", {true, std::nullopt}},
        {"READDATA", {true, std::nullopt}},
        {"RESIDNUM", {true, std::nullopt}},
        {"RESVNUM", {true, std::nullopt}},
        {"RIVDEBUG", {true, std::nullopt}},
        {"RIVRXSEC", {true, std::nullopt}},
        {"RIVERSYS", {true, std::nullopt}},
        {"RIVRDIMS", {true, std::nullopt}},
        {"RIVRPROP", {true, std::nullopt}},
        {"RIVRXSE", {true, std::nullopt}},
        {"RIVSALT", {true, std::nullopt}},
        {"RIVTRACE", {true, std::nullopt}},
        {"ROCKFRAC", {true, std::nullopt}},
        {"ROCKPAMA", {true, std::nullopt}},
        {"ROCKTABH", {true, std::nullopt}},
        {"ROCKTABW", {true, std::nullopt}},
        {"ROCKTHSG", {true, std::nullopt}},
        {"ROCKTSIG", {true, std::nullopt}},
        {"ROCKV", {true, std::nullopt}},
        {"RPTCPL", {false, std::nullopt}},
        {"RPTGRIDL", {false, std::nullopt}},
        {"RPTHM", {false, std::nullopt}},
        {"RPTHMG", {false, std::nullopt}},
        {"RPTHMD", {false, std::nullopt}},
        {"RPTHMW", {false, std::nullopt}},
        {"RPTINIT", {false, std::nullopt}},
        {"RPTISOL", {false, std::nullopt}},
        {"RPTPROPS", {false, std::nullopt}},
        {"RPTREGS", {false, std::nullopt}},
        {"RPTSOL", {false, std::nullopt}},
        {"RSGI", {true, std::nullopt}},
        {"RSSPE", {true, std::nullopt}},
        {"RSSSPEC", {false, std::nullopt}},
        {"RVCONST", {true, std::nullopt}},
        {"RVCONSTT", {true, std::nullopt}},
        {"RVGI", {true, std::nullopt}},
        {"REFINE", {true, std::nullopt}},
        {"RADFIN4", {true, std::nullopt}},
        {"RHO", {true, std::nullopt}},
        {"RKTRMDIR", {true, std::nullopt}},
        {"RPTGRID", {false, std::nullopt}},
        {"RPTPROS", {true, std::nullopt}},
        {"PRTRST", {true, std::nullopt}},
        {"RPTRUNSP", {false, std::nullopt}},
        {"RPTSMRY", {false, std::nullopt}},
        {"RSCONST", {true, std::nullopt}},
        {"RSCONSTT", {true, std::nullopt}},
        {"RSSPEC", {true, std::nullopt}},
        {"SAMG", {true, std::nullopt}},
        {"SAVE", {false, std::nullopt}},
        {"SKIP", {true, std::nullopt}},
        {"SKIP100", {true, std::nullopt}},
        {"SKIP300", {true, std::nullopt}},
        {"SALTNODE", {true, std::nullopt}},
        {"SALTREST", {true, std::nullopt}},
        {"SCALELIM", {true, std::nullopt}},
        {"SCDATAB", {true, std::nullopt}},
        {"SCDETAB", {true, std::nullopt}},
        {"SCDPTAB", {true, std::nullopt}},
        {"SCDPTRAC", {true, std::nullopt}},
        {"SCDPDIMS", {true, std::nullopt}},
        {"SCVD", {true, std::nullopt}},
        {"SEPVALS", {true, std::nullopt}},
        {"SFOAM", {true, std::nullopt}},
        {"SGF32D", {true, std::nullopt}},
        {"SHRATE", {true, std::string{"See the PLYSHEAR keyword instead"}}},
        {"SIGMA", {true, std::nullopt}},
        {"SIGMAGD", {true, std::nullopt}},
        {"SIGMAGDV", {true, std::nullopt}},
        {"SIGMATH", {true, std::nullopt}},
        {"SIGMAV", {true, std::nullopt}},
        {"SIMULATE", {true, std::nullopt}},
        {"SKRO", {true, std::nullopt}},
        {"SKRORG", {true, std::nullopt}},
        {"SKRORW", {true, std::nullopt}},
        {"SKRW", {true, std::nullopt}},
        {"SKRWR", {true, std::nullopt}},
        {"SLAVES", {true, std::nullopt}},
        {"SMULTX", {true, std::nullopt}},
        {"SMULTY", {true, std::nullopt}},
        {"SMULTZ", {true, std::nullopt}},
        {"SOCRS", {true, std::nullopt}},
        {"SOF32D", {true, std::nullopt}},
        {"SOLVCONC", {true, std::nullopt}},
        {"SOLVDIMS", {true, std::nullopt}},
        {"SOLVDIRS", {false, std::nullopt}},
        {"SOLVFRAC", {true, std::nullopt}},
        {"SOLVNUM", {true, std::nullopt}},
        {"SOLWNUM", {true, std::nullopt}},
        {"SOMGAS", {true, std::nullopt}},
        {"SOMWAT", {true, std::nullopt}},
        {"SSGCR", {true, std::nullopt}},
        {"SSGL", {true, std::nullopt}},
        {"SSOGCR", {true, std::nullopt}},
        {"SSOWCR", {true, std::nullopt}},
        {"SSWCR", {true, std::nullopt}},
        {"SSWL", {true, std::nullopt}},
        {"SSWU", {true, std::nullopt}},
        {"STOG", {true, std::nullopt}},
        {"STOW", {true, std::nullopt}},
        {"STWG", {true, std::nullopt}},
        {"SURF", {true, std::nullopt}},
        {"SURFACT", {true, std::nullopt}},
        {"SURFACTW", {true, std::nullopt}},
        {"SURFADDW", {true, std::nullopt}},
        {"SURFADS", {true, std::nullopt}},
        {"SURFCAPD", {true, std::nullopt}},
        {"SURFESAL", {true, std::nullopt}},
        {"SURFNUM", {true, std::nullopt}},
        {"SURFOPTS", {true, std::nullopt}},
        {"SURFROCK", {true, std::nullopt}},
        {"SURFST", {true, std::nullopt}},
        {"SURFSTES", {true, std::nullopt}},
        {"SURFVISC", {true, std::nullopt}},
        {"SURFWNUM", {true, std::nullopt}},
        {"SWF32D", {true, std::nullopt}},
        {"SWINGFAC", {true, std::nullopt}},
        {"TEMPNODE", {true, std::nullopt}},
        {"TEMPTVD", {true, std::nullopt}},
        {"TIGHTEN", {true, std::nullopt}},
        {"TIGHTENP", {true, std::nullopt}},
        {"TIME", {true, std::nullopt}},
        {"TNUM", {true, std::nullopt}},
        {"TPAMEPS", {true, std::nullopt}},
        {"TPAMEPSS", {true, std::nullopt}},
        {"TRACERKM", {true, std::nullopt}},
        {"TRACERKP", {true, std::nullopt}},
        {"TRACITVD", {true, std::nullopt}},
        {"TRACTVD", {true, std::nullopt}},
        {"TRACITVD", {true, std::nullopt}},
        {"TRADS", {true, std::nullopt}},
        {"TRANGL", {true, std::nullopt}},
        {"TRANR", {true, std::nullopt}},
        {"TRANTHT", {true, std::nullopt}},
        {"TRDCY", {true, std::nullopt}},
        {"TRDIF", {true, std::nullopt}},
        {"TRDIS", {true, std::nullopt}},
        {"TRKPF", {true, std::nullopt}},
        {"TRNHD", {true, std::nullopt}},
        {"TRPLPORO", {true, std::nullopt}},
        {"TRROCK", {true, std::nullopt}},
        {"TUNINGDP", {false, std::nullopt}},
        {"TUNINGH", {false, std::nullopt}},
        {"TUNINGL", {false, std::nullopt}},
        {"TUNINGS", {false, std::nullopt}},
        {"TZONE", {true, std::nullopt}},
        {"UDT", {true, std::nullopt}},
        {"UDTDIMS", {true, std::nullopt}},
        {"UNCODHMD", {true, std::nullopt}},
        {"UNIFOUTS", {false, std::nullopt}},
        {"UNIFSAVE", {false, std::nullopt}},
        {"USECUPL", {true, std::nullopt}},
        {"USEFLUX", {true, std::nullopt}},
        {"USENOFLO", {true, std::nullopt}},
        {"VDFLOW", {true, std::nullopt}},
        {"VDFLOWR", {true, std::nullopt}},
        {"VE", {true, std::nullopt}},
        {"VEDEBUG", {true, std::nullopt}},
        {"VEFIN", {true, std::nullopt}},
        {"VEFRAC", {true, std::nullopt}},
        {"VEFRACP", {true, std::nullopt}},
        {"VEFRACPV", {true, std::nullopt}},
        {"VEFRACV", {true, std::nullopt}},
        {"VFPCHK", {true, std::nullopt}},
        {"VFPTABL", {true, std::nullopt}},
        {"VISAGE", {true, std::nullopt}},
        {"VISCD", {true, std::nullopt}},
        {"VISDATES", {true, std::nullopt}},
        {"VISOPTS", {true, std::nullopt}},
        {"WAITBAL", {true, std::nullopt}},
        {"WALKALIN", {true, std::nullopt}},
        {"WALQCALC", {true, std::nullopt}},
        {"WAPI", {true, std::nullopt}},
        {"WARN", {false, std::nullopt}},
        {"WBHGLR", {true, std::nullopt}},
        {"WBOREVOL", {true, std::nullopt}},
        {"WCALCVAL", {true, std::nullopt}},
        {"WCONINJ", {true, std::nullopt}},
        {"WCONINJP", {true, std::nullopt}},
        {"WCUTBACK", {true, std::nullopt}},
        {"WCUTBACT", {true, std::nullopt}},
        {"WCYCLE", {true, std::nullopt}},
        {"WDFACCOR", {true, std::nullopt}},
        {"WDFAC", {true, std::nullopt}},
        {"WDRILTIM", {true, std::nullopt}},
        {"WDRILPRI", {true, std::nullopt}},
        {"WDRILRES", {true, std::nullopt}},
        {"WECONINJ", {true, std::nullopt}},
        {"WECONT", {true, std::nullopt}},
        {"WELCNTL", {true, std::nullopt}},
        {"WELDEBUG", {true, std::nullopt}},
        {"WELDRAW", {true, std::nullopt}},
        {"WELEVNT", {true, std::nullopt}},
        {"WELMOVEL", {true, std::nullopt}},
        {"WELOPENL", {true, std::nullopt}},
        {"WELPRI", {true, std::nullopt}},
        {"WELSOMIN", {true, std::nullopt}},
        {"WELSPECL", {true, std::nullopt}},
        {"WFRICSEG", {true, std::nullopt}},
        {"WFRICSGL", {true, std::nullopt}},
        {"WFRICTN", {true, std::nullopt}},
        {"WFRICTNL", {true, std::nullopt}},
        {"WGASPROD", {true, std::nullopt}},
        {"WGORPEN", {true, std::nullopt}},
        {"WH2NUM", {true, std::nullopt}},
        {"WH3NUM", {true, std::nullopt}},
        {"WHEDREFD", {true, std::nullopt}},
        {"WHTEMP", {true, std::nullopt}},
        {"WLIMTOL", {true, std::nullopt}},
        {"WLIFT", {true, std::nullopt}},
        {"WLISTARG", {true, std::nullopt}},
        {"WLISTNAM", {true, std::nullopt}},
        {"WLISTOPT", {true, std::nullopt}},
        {"WNETCTRL", {true, std::nullopt}},
        {"WNETDP", {true, std::nullopt}},
        {"WORKLIM", {true, std::nullopt}},
        {"WORKTHP", {true, std::nullopt}},
        {"WPIMULTL", {true, std::nullopt}},
        {"WPITAB", {true, std::nullopt}},
        {"WPLUG", {true, std::nullopt}},
        {"WPOLYRED", {true, std::nullopt}},
        {"WPOTCALC", {false, std::nullopt}},
        {"WREGROUP", {true, std::nullopt}},
        {"WSCCLEAN", {true, std::nullopt}},
        {"WSCCLENL", {true, std::nullopt}},
        {"WSCTAB", {true, std::nullopt}},
        {"WSEGDFIN", {true, std::nullopt}},
        {"WSEGDFMD", {true, std::nullopt}},
        {"WSEGDFPA", {true, std::nullopt}},
        {"WSEGEXSS", {true, std::nullopt}},
        {"WSEGFLIM", {true, std::nullopt}},
        {"WSEGFMOD", {true, std::nullopt}},
        {"WSEGINIT", {true, std::nullopt}},
        {"WSEGITER", {false, std::nullopt}},
        {"WSEGLABY", {true, std::nullopt}},
        {"WSEGLINK", {true, std::nullopt}},
        {"WSEGMULT", {true, std::nullopt}},
        {"WSEGPROP", {true, std::nullopt}},
        {"WSEGPULL", {true, std::nullopt}},
        {"WSEGSEP", {true, std::nullopt}},
        {"WSEGSOLV", {true, std::nullopt}},
        {"WSEGTABL", {true, std::nullopt}},
        {"WSURFACT", {true, std::nullopt}},
        {"WTADD", {true, std::nullopt}},
        {"WTEMPQ", {true, std::nullopt}},
        {"WTHPMAX", {true, std::nullopt}},
        {"WTMULT", {true, std::nullopt}},
        {"ZIPPY2", {false, std::nullopt}},
        {"ZIPP2OFF", {false, std::nullopt}},
    };

    return unsupported_keywords;
}

}
