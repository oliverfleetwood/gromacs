/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
/*! \file
 * \brief
 * Contains code for calling GROMACS routines from an external program
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \inpublicapi
 * \ingroup module_qmmm
 */
#include <stdlib.h>
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/mtop_util.h"
#include "gromacs/legacyheaders/sim_util.h"
#include "gromacs/legacyheaders/vsite.h"
#include "gromacs/legacyheaders/constr.h"
#include "gromacs/legacyheaders/nrnb.h"
#include "gromacs/legacyheaders/vec.h"
#include "gromacs/legacyheaders/mdatoms.h"
#include "gromacs/legacyheaders/force.h"
#include "gromacs/legacyheaders/oenv.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/mdlib/minimize.h"
#include "gromacs/mmslave.h"
#include "gromacs/mmslave/mmslave.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/timing/wallcycle.h"

/* Note: the C-interface is all the way down in this file */

namespace gmx
{

GromacsInABox::GromacsInABox(FILE             *fplog,
                             const t_commrec  *cr,
                             const gmx_mtop_t *mtop,
                             const t_inputrec *ir,
                             matrix            box)
{
    int koko = 0;
    // Initiate everything
    printf("MOMO %d\n", koko++);
    bFirst_ = TRUE;

    printf("MOMO %d\n", koko++);

    s_min_ = init_em_state();
    copy_mat(box, s_min_->s.box);
    printf("MOMO %d\n", koko++);

    top_ = gmx_mtop_generate_local_top(mtop, ir);;

    printf("MOMO %d\n", koko++);
    //! Flops
    init_nrnb(&nrnb_);
    printf("MOMO %d\n", koko++);

    //! CPU Accounting
    wcycle_ = wallcycle_init(fplog, 1, (t_commrec *)cr, cr->nnodes, 0);
    printf("MOMO %d\n", koko++);

    //! Energetics
    gstat_ = global_stat_init((t_inputrec *)ir);

    printf("MOMO %d\n", koko++);
    //! Virtual sites
    vsite_ = NULL; // init_vsite((gmx_mtop_t *)mtop, (t_commrec *)cr, FALSE);
    printf("MOMO %d\n", koko++);

    //! Constraints
    //gmx_edsam_t ed = NULL;
    //init_state(&state_, mtop->natoms, 1, 0, 0, 0);
    constr_ = NULL; //init_constraints(fplog, (gmx_mtop_t *)mtop, (t_inputrec *)ir, ed, &state_, (t_commrec *)cr);
    printf("MOMO %d\n", koko++);

    //! FC Data
    fcd_ = (t_fcdata *)calloc(1, sizeof(*fcd_));
    printf("MOMO %d\n", koko++);

    //! Molecular graph
    graph_ = mk_graph(fplog, &(top_->idef), 0,
                      mtop->natoms, FALSE, FALSE);
    printf("MOMO %d\n", koko++);

    //! MD Atoms
    mdatoms_ = init_mdatoms(fplog, (gmx_mtop_t *)mtop, FALSE);
    atoms2md((gmx_mtop_t *)mtop, (t_inputrec *)ir,
             0, NULL,
             0, mtop->natoms,
             mdatoms_);
    printf("MOMO %d\n", koko++);

    //! Force record
    output_env_init_default(&oenv_);
    fr_ = mk_forcerec();
    init_forcerec(fplog, oenv_, fr_, fcd_, ir, mtop, cr, box,
                  NULL, NULL, NULL, NULL, NULL,
                  FALSE, 0.0);
    fr_->qr->QMMMscheme = eQMMMschemeslave;
    printf("MOMO %d\n", koko++);

    //! Energy data
    enerd_ = (gmx_enerdata_t *)calloc(1, sizeof(*enerd_));

    //! Check lambda stuff
    int n_lambda = 0;
    init_enerdata(ir->opts.ngener, n_lambda, enerd_);
    printf("MOMO %d\n", koko++);

    //! Accounting
    count_ = 0;

}

GromacsInABox::~GromacsInABox()
{
    // Delete everything
}

MMSlave::MMSlave(const t_commrec *cr)
{
    x_    = NULL;
    v_    = NULL;
    f_    = NULL;
    giab_ = NULL;
    cr_   = cr;
}

bool MMSlave::readTpr(const char *tpr)
{
    t_tpxheader tpx;
    int         version, generation, natoms;
    int         koko = 0;

    printf("KOKO %d\n", koko++);
    read_tpxheader(tpr, &tpx, FALSE, &version, &generation);
    natoms = tpx.natoms;
    x_     = (rvec *)calloc(natoms, sizeof(rvec));
    v_     = (rvec *)calloc(natoms, sizeof(rvec));
    f_     = (rvec *)calloc(natoms, sizeof(rvec));

    printf("KOKO %d\n", koko++);
    (void) read_tpx(tpr, &inputrec_, box_, &natoms, x_,
                    (tpx.bV ? v_ : NULL),
                    (tpx.bF ? f_ : NULL), &mtop_);

    gmx_mtop_atomloop_all_t aloop = gmx_mtop_atomloop_all_init(&mtop_);
    int                     at_global;
    t_atom                 *atom;

    groupSize_.resize(1+inputrec_.opts.ngQM);
    while (gmx_mtop_atomloop_all_next(aloop, &at_global, &atom))
    {
        groupSize_[ggrpnr(&(mtop_.groups), egcQMMM, at_global)]++;
    }
    GMX_RELEASE_ASSERT((natoms == nAtoms()),
                       "Total number of atoms not consistent with group indices");
    printf("KOKO %d\n", koko++);

    giab_ = new GromacsInABox(stdout, cr_, &mtop_, &inputrec_, box_);
    printf("KOKO %d\n", koko++);

    return true;
}

int MMSlave::nAtoms()
{
    int n = 0;
    for (unsigned int i = 0; (i < groupSize_.size()); i++)
    {
        n += groupSize_[i];
    }
    return n;
}

static bool copyIt(int         natoms_dst,
                   rvec       *x_dst,
                   int         natoms_src,
                   const rvec *x_src)
{
    if ((natoms_dst < natoms_src) || (NULL == x_dst))
    {
        return false;
    }
    for (int i = 0; (i < natoms_src); i++)
    {
        copy_rvec(x_src[i], x_dst[i]);
    }
    return true;
}

bool MMSlave::copyX(int natoms, rvec *x)
{
    return copyIt(natoms, x, groupSize_[0] + groupSize_[1], x_);
}

bool MMSlave::copyV(int natoms, rvec *v)
{
    return copyIt(natoms, v, groupSize_[0] + groupSize_[1], v_);
}

bool MMSlave::copyF(int natoms, rvec *f)
{
    return copyIt(natoms, f, groupSize_[0] + groupSize_[1], f_);
}

void MMSlave::cleanUp()
{
    if (NULL != x_)
    {
        free(x_);
    }
    if (NULL != v_)
    {
        free(v_);
    }
    if (NULL != f_)
    {
        free(f_);
    }
}

bool MMSlave::setAtomQ(atom_id id, double q)
{
    t_atom               *atom;

    gmx_mtop_atomlookup_t alook = gmx_mtop_atomlookup_init(&mtop_);

    gmx_mtop_atomnr_to_atom(alook, id, &atom);
    atom->q  = q;
    atom->qB = q;

    gmx_mtop_atomlookup_destroy(alook);

    return true;
}

bool MMSlave::getAtomQ(atom_id id, double *q)
{
    t_atom               *atom;

    gmx_mtop_atomlookup_t alook = gmx_mtop_atomlookup_init(&mtop_);

    gmx_mtop_atomnr_to_atom(alook, id, &atom);
    *q = atom->q;

    gmx_mtop_atomlookup_destroy(alook);

    return true;
}

bool MMSlave::getAtomNumber(atom_id id, int *atomNumber)
{
    t_atom               *atom;

    gmx_mtop_atomlookup_t alook = gmx_mtop_atomlookup_init(&mtop_);

    gmx_mtop_atomnr_to_atom(alook, id, &atom);
    *atomNumber = atom->atomnumber;
    if (NOTSET == *atomNumber)
    {
        GMX_THROW(InvalidInputError("No information about atom numbers in the tpr file"));
    }
    gmx_mtop_atomlookup_destroy(alook);

    return true;
}

bool MMSlave::getGroupID(atom_id id, int *groupID)
{
    *groupID = ggrpnr(&(mtop_.groups), egcQMMM, id);

    return true;
}

bool MMSlave::calcEnergy(FILE       *fplog,
                         const rvec *x,
                         rvec       *f,
                         double     *energy)
{
    int hoho = 0;
    clear_mat(giab_->vir_);
    clear_mat(giab_->pres_);
    clear_rvec(giab_->mu_tot_);

    // Copy the coordinates.
    if (!copyIt(nAtoms(), x_, nAtoms(), x))
    {
        GMX_THROW(InternalError("Copying coordinates"));
    }
    printf("HOHOHO %d\n", hoho++);
    // Make sure the coordinates are in the state too!
    giab_->s_min_->s.x = (rvec *)x;
    giab_->s_min_->f = (rvec *)f;
    evaluate_energy(fplog,
                    (t_commrec *)cr_,
                    &mtop_,
                    giab_->s_min_,
                    giab_->top_,
                    &inputrec_,
                    &giab_->nrnb_,
                    giab_->wcycle_,
                    giab_->gstat_,
                    giab_->vsite_,
                    giab_->constr_,
                    giab_->fcd_,
                    giab_->graph_,
                    giab_->mdatoms_,
                    giab_->fr_,
                    giab_->mu_tot_,
                    giab_->enerd_,
                    giab_->vir_,
                    giab_->pres_,
                    giab_->count_,
                    giab_->bFirst_);
    giab_->bFirst_ = FALSE;
    printf("HOHOHO %d\n", hoho++);

    // Copy the forces. Check!
    copyIt(nAtoms(), f, nAtoms(), giab_->s_min_->f);
    printf("HOHOHO %d\n", hoho++);

    // Copy the enery
    *energy = giab_->enerd_->term[F_EPOT];
    printf("HOHOHO %d\n", hoho++);

    return true;
}

}

//! Abstract type for the mmslave code
typedef struct gmx_mmslave {
    //! Embedded C++ class
    gmx::MMSlave *mms;
} gmx_mmslave;

/* Routines for C interface to the MMSlave class */
gmx_mmslave_t mmslave_init(const t_commrec *cr)
{
    gmx_mmslave *gms;

    gms      = (gmx_mmslave *) calloc(1, sizeof(gmx_mmslave));
    gms->mms = new gmx::MMSlave(cr);

    return gms;
}

void mmslave_done(gmx_mmslave_t gms)
{
    delete gms->mms;
    free(gms);
}

int mmslave_read_tpr(const char   *tpr,
                     gmx_mmslave_t gms)
{
    printf("Koko\n");
    if (gms->mms->readTpr(tpr))
    {
        return 1;
    }
    return 0;
}

int mmslave_natoms(gmx_mmslave_t gms)
{
    return gms->mms->nAtoms();
}

int mmslave_ngroups(gmx_mmslave_t gms)
{
    return gms->mms->nGroups();
}

int mmslave_copyX(gmx_mmslave_t gms, int natoms, rvec *x)
{
    if (gms->mms->copyX(natoms, x))
    {
        return 1;
    }
    return 0;
}

int mmslave_copyV(gmx_mmslave_t gms, int natoms, rvec *v)
{
    if (gms->mms->copyX(natoms, v))
    {
        return 1;
    }
    return 0;
}

int mmslave_copyF(gmx_mmslave_t gms, int natoms, rvec *f)
{
    if (gms->mms->copyX(natoms, f))
    {
        return 1;
    }
    return 0;
}

void mmslave_clean(gmx_mmslave_t gms)
{
    gms->mms->cleanUp();
    delete gms->mms;
    free(gms);
}

int mmslave_set_q(gmx_mmslave_t gms,
                  atom_id       id,
                  double        q)
{
    if (gms->mms->setAtomQ(id, q))
    {
        return 1;
    }
    return 0;
}

int mmslave_get_q(gmx_mmslave_t gms,
                  atom_id       id,
                  double       *q)
{
    if (gms->mms->getAtomQ(id, q))
    {
        return 1;
    }
    return 0;
}

int mmslave_get_atomnumber(gmx_mmslave_t gms,
                           atom_id       id)
{
    int atomNumber;

    if (gms->mms->getAtomNumber(id, &atomNumber))
    {
        return atomNumber;
    }
    return NOTSET;
}

int mmslave_get_group_id(gmx_mmslave_t gms,
                         atom_id       id)
{
    int groupID;

    if (gms->mms->getGroupID(id, &groupID))
    {
        return groupID;
    }
    return NOTSET;
}

int mmslave_calc_energy(gmx_mmslave_t gms,
                        FILE         *fplog,
                        const rvec   *x,
                        rvec         *f,
                        double       *energy)
{
    if (gms->mms->calcEnergy(fplog, x, f, energy))
    {
        return 1;
    }
    return 0;
}
