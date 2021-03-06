/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2018,2019,2020, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements test of some pulling routines
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_pulling
 */
#include "gmxpre.h"

#include "gromacs/pulling/pull.h"

#include <cmath>

#include <algorithm>

#include <gtest/gtest.h>

#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull_internal.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "../../mdtypes/md_enums.h"

namespace gmx
{
namespace
{

using gmx::test::defaultRealTolerance;

class PullTest : public ::testing::Test
{
protected:
    PullTest() {}

    static void test(PbcType pbcType, matrix box)
    {
        t_pbc pbc;

        // PBC stuff
        set_pbc(&pbc, pbcType, box);

        GMX_ASSERT(pbc.ndim_ePBC >= 1 && pbc.ndim_ePBC <= DIM,
                   "Tests only support PBC along at least x and at most x, y, and z");

        real boxSizeZSquared;
        if (pbc.ndim_ePBC > ZZ)
        {
            boxSizeZSquared = gmx::square(box[ZZ][ZZ]);
        }
        else
        {
            boxSizeZSquared = GMX_REAL_MAX;
        }

        {
            // Distance pulling in all 3 dimensions
            t_pull_coord params;
            params.eGeom   = epullgDIST;
            params.dim[XX] = 1;
            params.dim[YY] = 1;
            params.dim[ZZ] = 1;
            pull_coord_work_t pcrd(params);
            clear_dvec(pcrd.spatialData.vec);

            real minBoxSize2 = GMX_REAL_MAX;
            for (int d = 0; d < pbc.ndim_ePBC; d++)
            {
                minBoxSize2 = std::min(minBoxSize2, norm2(box[d]));
            }
            EXPECT_REAL_EQ_TOL(0.25 * minBoxSize2, max_pull_distance2(&pcrd, &pbc),
                               defaultRealTolerance());
        }

        {
            // Distance pulling along Z
            t_pull_coord params;
            params.eGeom   = epullgDIST;
            params.dim[XX] = 0;
            params.dim[YY] = 0;
            params.dim[ZZ] = 1;
            pull_coord_work_t pcrd(params);
            clear_dvec(pcrd.spatialData.vec);
            EXPECT_REAL_EQ_TOL(0.25 * boxSizeZSquared, max_pull_distance2(&pcrd, &pbc),
                               defaultRealTolerance());
        }

        {
            // Directional pulling along Z
            t_pull_coord params;
            params.eGeom   = epullgDIR;
            params.dim[XX] = 1;
            params.dim[YY] = 1;
            params.dim[ZZ] = 1;
            pull_coord_work_t pcrd(params);
            clear_dvec(pcrd.spatialData.vec);
            pcrd.spatialData.vec[ZZ] = 1;
            EXPECT_REAL_EQ_TOL(0.25 * boxSizeZSquared, max_pull_distance2(&pcrd, &pbc),
                               defaultRealTolerance());
        }

        {
            // Directional pulling along X
            t_pull_coord params;
            params.eGeom   = epullgDIR;
            params.dim[XX] = 1;
            params.dim[YY] = 1;
            params.dim[ZZ] = 1;
            pull_coord_work_t pcrd(params);
            clear_dvec(pcrd.spatialData.vec);
            pcrd.spatialData.vec[XX] = 1;

            real minDist2 = square(box[XX][XX]);
            for (int d = XX + 1; d < DIM; d++)
            {
                minDist2 -= square(box[d][XX]);
            }
            EXPECT_REAL_EQ_TOL(0.25 * minDist2, max_pull_distance2(&pcrd, &pbc), defaultRealTolerance());
        }

        {
            // transformation pull coordinate test
            // Create standard pull coordinate
            t_pull_coord params;
            params.eGeom = epullgDIST;
            pull_coord_work_t x1_pcrd(params);
            // Create transformation pull coordinate
            params.eGeom      = epullgTRANSFORMATION;
            const char * expression = "x1^2 + 3";
            params.expression = expression;
            pull_coord_work_t transformation_pcrd(params);
            transformation_pcrd.expressionParser.initialize(1);

            pull_t pull;
            pull.coord = std::vector<pull_coord_work_t>({x1_pcrd, transformation_pcrd});
            for (double v = 0; v < 10; v++)
            {
                // transformation pull coord value
                x1_pcrd.spatialData.value = v;
                get_pull_coord_value(&pull, 1, &pbc);
                EXPECT_DOUBLE_EQ(v * v + 3, transformation_pcrd.spatialData.value)
                        << "transformation coordinate value does not match the expected expression.";

                // force and derivative
                double transformation_force     = v + 0.5;
                transformation_pcrd.scalarForce = transformation_force;
                double x1_force = compute_force_from_transformation_coord(&pull, 1, 0);
                EXPECT_DOUBLE_EQ(2 * v * transformation_force, x1_force)
                        << "Force distributed from transformation coordinate "
                           "should be the derivative times the transformation force.";
            }
        }
    }
};

TEST_F(PullTest, MaxPullDistanceXyzScrewBox)
{
    matrix box = { { 10, 0, 0 }, { 0, 10, 0 }, { 0, 0, 10 } };

    test(PbcType::Screw, box);
}

TEST_F(PullTest, MaxPullDistanceXyzCubicBox)
{
    matrix box = { { 10, 0, 0 }, { 0, 10, 0 }, { 0, 0, 10 } };

    test(PbcType::Xyz, box);
}

TEST_F(PullTest, MaxPullDistanceXyzTricBox)
{
    matrix box = { { 10, 0, 0 }, { 3, 10, 0 }, { 3, 4, 10 } };

    test(PbcType::Xyz, box);
}

TEST_F(PullTest, MaxPullDistanceXyzLongBox)
{
    matrix box = { { 10, 0, 0 }, { 0, 10, 0 }, { 0, 0, 30 } };

    test(PbcType::Xyz, box);
}

TEST_F(PullTest, MaxPullDistanceXySkewedBox)
{
    matrix box = { { 10, 0, 0 }, { 5, 8, 0 }, { 0, 0, 0 } };

    test(PbcType::XY, box);
}

} // namespace

} // namespace gmx
