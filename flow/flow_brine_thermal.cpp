/*
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
#include "config.h"

#include <opm/material/common/ResetLocale.hpp>
#include <opm/simulators/flow/Main.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

namespace Opm {
namespace Properties {
namespace TTag {
struct EclFlowBrineProblem {
    using InheritsFrom = std::tuple<EclFlowProblem>;
};
}
template<class TypeTag>
struct EnableBrine<TypeTag, TTag::EclFlowBrineProblem> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::EclFlowBrineProblem> {
    static constexpr bool value = true;
};
}}

int main(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::EclFlowBrineProblem;
    auto mainObject = Opm::Main(argc, argv);
    return mainObject.runStatic<TypeTag>();
}